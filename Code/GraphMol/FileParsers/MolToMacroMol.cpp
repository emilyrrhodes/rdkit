//
//  Copyright (C) 2025 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/Conversions.h>
#include <GraphMol/MacroMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace RDKit {
namespace {

class BondConnectionDef {
 public:
  BondConnectionDef(unsigned int atomIdx1Init, unsigned int atomIdx2Init)
      : atomIdx1(atomIdx1Init), atomIdx2(atomIdx2Init) {}

  unsigned int atomIdx1;
  unsigned int atomIdx2;

 public:
  bool operator<(const BondConnectionDef &other) const {
    if (atomIdx1 != other.atomIdx1) {
      return atomIdx1 < other.atomIdx1;
    }
    return atomIdx2 < other.atomIdx2;
  }
};

bool isTemplateMatchAHit(
    const MatchVectType &match, const ROMol &mol, const ROMol *templateMol,
    const RWMol *queryMol, std::map<unsigned int, unsigned int> &atomMap,
    const std::vector<SubstanceGroup::AttachPoint> &supAttachPoints,
    std::map<BondConnectionDef, std::string> &bondConnectionMap,
    std::vector<unsigned int> &atomsInMatch) {
  // get a new set of match points to be used for the new template ref atoms
  atomsInMatch.clear();
  // newAttachOrds.clear();
  std::map<unsigned int, unsigned int> molToQueryMap;

  // check to see if any of the atoms or bonds are already used
  for (auto pair : match) {
    if (atomMap.contains(pair.second)) {
      // atom already used in another match - skip this match
      return false;
    }

    molToQueryMap[pair.second] = pair.first;
    atomsInMatch.push_back(pair.second);
  }
  // for each atom in the match, check its bonds.  If it has a bond to
  // another atom in the hit,
  //  make sure that bond is the query.  If not, skip this match.

  // If the atom is an attachment point for the template, see if there is
  // a bond to an atom NOT in the match. If so, record the attachment
  // order for later.

  // if the atom is not an attachment point,  it has an external bond skip this
  // match.

  std::vector<bool> attachPointUsed(supAttachPoints.size(), false);
  for (auto molToQueryMapItem : molToQueryMap) {
    auto atom = mol.getAtomWithIdx(molToQueryMapItem.first);
    auto queryAtom = queryMol->getAtomWithIdx(molToQueryMapItem.second);
    unsigned int templateAtomIdx =
        queryAtom->getProp<unsigned int>("origAtomId");
    for (const auto bond : mol.atomBonds(atom)) {
      unsigned int nbrAtomIdx = bond->getOtherAtomIdx(atom->getIdx());
      if (molToQueryMap.contains(nbrAtomIdx)) {
        // make sure the bond is also in the query mol - if not skip this
        // match

        auto queryBond = queryMol->getBondBetweenAtoms(
            molToQueryMap[bond->getBeginAtomIdx()],
            molToQueryMap[bond->getEndAtomIdx()]);

        if (!queryBond) {
          return false;
        }

        // bond to another atom is in the query - ok
        continue;
      } else {
        // bond to an atom NOT in the main SUP - could be a leaving group

        bool sapFound = false;
        for (unsigned int attachPointIndex = 0;
             attachPointIndex < supAttachPoints.size(); ++attachPointIndex) {
          auto attachPoint = supAttachPoints[attachPointIndex];
          if (attachPointUsed[attachPointIndex]) {
            continue;
          }

          if (attachPoint.aIdx == templateAtomIdx) {
            // see if the attachment point is an H or OH, and the nbr atom is
            // also an H or OH if so, this leaving group can be considered part
            // of the template

            sapFound = true;
            attachPointUsed[attachPointIndex] = true;

            auto nbrAtom = mol.getAtomWithIdx(nbrAtomIdx);
            const Atom *templateNbrAtom = nullptr;
            if (attachPoint.lvIdx >= 0) {
              templateNbrAtom = templateMol->getAtomWithIdx(attachPoint.lvIdx);
            }
            // if the template nbr atom is a single atom and so is the matched
            // atom, and they have the same atomic num and numHs, treat this as
            // part of the template
            if (attachPoint.lvIdx >= 0 && nbrAtom->getDegree() == 1 &&
                templateNbrAtom->getDegree() == 1 &&
                templateNbrAtom->getAtomicNum() == nbrAtom->getAtomicNum() &&
                templateNbrAtom->getTotalNumHs() == nbrAtom->getTotalNumHs()) {
              // ok - treat this as part of the template

              atomsInMatch.push_back(nbrAtomIdx);
            } else {
              // save the attach point relative to the full atom molecule being
              // converted

              bondConnectionMap[BondConnectionDef(atom->getIdx(), nbrAtomIdx)] =
                  attachPoint.id;
            }
            break;
          }
        }

        if (!sapFound) {
          // external bond that is not allowed - skip this match
          return false;
        }
      }
    }
  }

  return true;
}

}  // namespace

std::unique_ptr<RDKit::MacroMol> MolToMacroMol(
    const ROMol &mol, const RDKit::MacroMolTemplateLib &templates,
    MolToMacroMolParams molToMacroMolParams) {
  auto macro_mol = std::make_unique<MacroMol>();

  Conformer *conf = nullptr;
  if (mol.getNumConformers() > 0 && templates.doesLibHaveCoords()) {
    conf = new Conformer();
    conf->set3D(false);

    macro_mol->addConformer(conf, true);
  }

  std::map<unsigned int, unsigned int> atomMap;
  std::map<BondConnectionDef, std::string> bondConnectionMap;
  for (const auto &templateEntry : templates) {
    auto templateMol = templateEntry.get();
    templateMol->updatePropertyCache(false);
    std::vector<std::string> templateNames;

    std::string templateAtomClass;
    std::string templateNameToUse;
    templateMol->getPropIfPresent<std::string>(common_properties::molAtomClass,
                                               templateAtomClass);
    templateMol->getPropIfPresent<std::vector<std::string>>(
        common_properties::templateNames, templateNames);
    switch (molToMacroMolParams.macroUseTemplateName) {
      case MacroMolUseTemplateName::UseFirstName:
        templateNameToUse = templateNames[0];
        break;
      case MacroMolUseTemplateName::UseSecondName:
        templateNameToUse = templateNames.back();
        break;
    }

    // get the sgroup that is the base for this template

    const RDKit::SubstanceGroup *sgroup = templateMol->getMainSgroup();

    if (sgroup == nullptr) {
      continue;  // skip this template
    }
    auto supAttachPoints = sgroup->getAttachPoints();

    std::unique_ptr<RWMol> queryMol(new RWMol(*templateMol));

    queryMol->beginBatchEdit();
    for (auto atom : queryMol->atoms()) {
      atom->setProp("origAtomId", atom->getIdx());

      if (std::find(sgroup->getAtoms().begin(), sgroup->getAtoms().end(),
                    atom->getIdx()) == sgroup->getAtoms().end()) {
        queryMol->removeAtom(atom->getIdx());
      }
    }
    queryMol->commitBatchEdit();

    // find all occurrences of the template in the molecule

    MolOps::sanitizeMol(*queryMol);

    queryMol->updatePropertyCache(false);

    std::vector<MatchVectType> matchVect;

    SubstructMatchParameters params;
    params.recursionPossible = false;
    params.useChirality = true;
    params.useQueryQueryMatches = false;
    params.maxMatches = 0;  // find all matches

    matchVect = SubstructMatch(mol, *queryMol, params);

    if (!matchVect.size()) {
      continue;
    }

    bool templateCopied = false;
    for (const auto &match : matchVect) {
      std::vector<unsigned int> atomsInMatch;

      // add this match to the MacroMol if it is a valid hit

      if (!isTemplateMatchAHit(match, mol, templateMol, queryMol.get(), atomMap,
                               supAttachPoints, bondConnectionMap,
                               atomsInMatch)) {
        continue;
      }

      if (!templateCopied) {
        // add the template to the MacroMol
        std::unique_ptr<MacroMolTemplate> tempTemplate(
            new MacroMolTemplate(*templateMol));
        macro_mol->addTemplate(tempTemplate);
        templateCopied = true;
      }

      auto newAtomIdx = macro_mol->getNumAtoms();
      for (const auto atomMatch : atomsInMatch) {
        atomMap[atomMatch] = newAtomIdx;
      }

      // create the macro atom reference

      macro_mol->addAtom(new Atom(0), true, true);
      auto macro_mol_atom = macro_mol->getAtomWithIdx(newAtomIdx);

      macro_mol_atom->setAtomicNum(0);
      macro_mol_atom->setProp(common_properties::dummyLabel, templateNameToUse);
      macro_mol_atom->setProp(common_properties::molAtomClass,
                              templateAtomClass);

      if (conf) {
        conf->resize(macro_mol->getNumAtoms());

        RDGeom::Point3D pos;
        for (const auto atomMatch : atomsInMatch) {
          RDGeom::Point3D atomPos = mol.getConformer().getAtomPos(atomMatch);
          pos += atomPos;
        }
        pos /= static_cast<double>(atomsInMatch.size());

        conf->setAtomPos(macro_mol_atom->getIdx(), pos);
      }
    }
  }

  // add all atoms that are not in a template

  for (const auto atom : mol.atoms()) {
    if (!atomMap.contains(atom->getIdx())) {
      auto newAtomIdx = macro_mol->getNumAtoms();
      macro_mol->addAtom(new Atom(*atom), true, true);
      atomMap[atom->getIdx()] = newAtomIdx;
      if (conf) {
        conf->resize(macro_mol->getNumAtoms());
        conf->setAtomPos(newAtomIdx,
                         mol.getConformer().getAtomPos(atom->getIdx()));
      }
    }
  }

  // bonds:  if the atoms of the original bonds are mapped to the same
  // atom, they are in the template def

  for (auto bond : mol.bonds()) {
    unsigned int begAtomIdx = bond->getBeginAtomIdx();
    unsigned int endAtomIdx = bond->getEndAtomIdx();

    if (atomMap.contains(begAtomIdx) && atomMap.contains(endAtomIdx) &&
        atomMap[begAtomIdx] == atomMap[endAtomIdx]) {
      continue;
    }

    // bond'S atoms are in different new atoms - add the bond

    macro_mol->addBond(atomMap[begAtomIdx], atomMap[endAtomIdx],
                       bond->getBondType());
    auto newBond = macro_mol->getBondWithIdx(macro_mol->getNumBonds() - 1);

    if (bondConnectionMap.contains(BondConnectionDef(begAtomIdx, endAtomIdx))) {
      newBond->setProp(
          common_properties::_MolFileBondAttachPt1,
          bondConnectionMap[BondConnectionDef(begAtomIdx, endAtomIdx)]);
    }
    if (bondConnectionMap.contains(BondConnectionDef(endAtomIdx, begAtomIdx))) {
      newBond->setProp(
          common_properties::_MolFileBondAttachPt2,
          bondConnectionMap[BondConnectionDef(endAtomIdx, begAtomIdx)]);
    }
  }

  return macro_mol;
}

}  // end of namespace RDKit
