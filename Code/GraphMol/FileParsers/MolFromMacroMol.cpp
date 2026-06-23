//
//  Copyright (C) 2025 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/FileParseException.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/Conversions.h>
#include <GraphMol/MacroMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include "MolSGroupParsing.h"

#include <algorithm>
#include <climits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace RDKit::SGroupParsing;

namespace RDKit {
namespace {

// maps (main atom#, template atom#) to new atom#.  A template atom# of UINT_MAX
// marks an atom that is not part of a template.
using OriginAtomMap = std::map<std::pair<unsigned int, unsigned int>, unsigned int>;

// maps (main atom#, attach label) to template atom#
using AttachMap = std::map<std::pair<unsigned int, std::string>, unsigned int>;

unsigned int getNewAtomForBond(const Atom *atom, const std::string &lbl,
                               const OriginAtomMap &originAtomMap,
                               const AttachMap &attachMap) {
  unsigned int atomIdx = atom->getIdx();
  if (lbl == "") {  // not a template atom
    return originAtomMap.at(std::pair(atomIdx, UINT_MAX));
  }

  // if here , it is a template atom

  auto attachMapIt = attachMap.find(std::pair(atomIdx, lbl));
  if (attachMapIt == attachMap.end()) {
    throw FileParseException("Attachment ord not found");
  }
  return originAtomMap.at(std::pair(atomIdx, attachMapIt->second));
}

void copySgroupIntoResult(
    RWMol &mol, const MacroMol &macroMol, const unsigned int atomIdx,
    const RDKit::SubstanceGroup &sgroup, std::string sgroupName,
    std::vector<std::unique_ptr<SubstanceGroup>> &newSgroups,
    RDKit::Conformer *newConf, const RDGeom::Point3D &coordOffset,
    OriginAtomMap &originAtomMap) {
  // add a superatom sgroup to mark the atoms from this macro atom
  // expansion. These new superatom sgroups are not put into the output mol
  // yet, because the bonds have not be added to the mol nor to the sgroup.
  // They are saved in an array to be added later

  const std::string typ = "SUP";
  newSgroups.emplace_back(new SubstanceGroup((ROMol *)&mol, typ));
  auto newSgroup = newSgroups.back().get();
  newSgroup->setProp("LABEL", sgroupName);

  // copy the atoms of the sgroup into the new molecule

  if (newConf) {
    newConf->resize(newConf->getNumAtoms() +
                    macroMol.getTemplate(atomIdx)->getNumAtoms());
  }

  for (auto templateAtomIdx : sgroup.getAtoms()) {
    auto templateAtom =
        macroMol.getTemplate(atomIdx)->getAtomWithIdx(templateAtomIdx);
    auto newAtom = new Atom(*templateAtom);

    mol.addAtom(newAtom, true, true);
    newSgroup->addAtomWithIdx(newAtom->getIdx());

    originAtomMap[std::pair(atomIdx, templateAtomIdx)] = newAtom->getIdx();

    if (newConf) {
      newConf->setAtomPos(newAtom->getIdx(),
                          coordOffset + macroMol.getTemplate(atomIdx)
                                            ->getConformer()
                                            .getAtomPos(templateAtomIdx));
    }
  }
}

void processBondInMainMol(const Bond *bond, RWMol &mol,
                          const OriginAtomMap &originAtomMap,
                          const AttachMap &attachMap) {
  std::string lbl1 = "", lbl2 = "";
  bond->getPropIfPresent(common_properties::_MolFileBondAttachPt1, lbl1);
  bond->getPropIfPresent(common_properties::_MolFileBondAttachPt2, lbl2);

  unsigned int newBeginAtom =
      getNewAtomForBond(bond->getBeginAtom(), lbl1, originAtomMap, attachMap);
  if (newBeginAtom == UINT_MAX) {
    throw FileParseException("Error getting new atom for bond");
  }

  unsigned int newEndAtom =
      getNewAtomForBond(bond->getEndAtom(), lbl2, originAtomMap, attachMap);
  if (newEndAtom == UINT_MAX) {
    throw FileParseException("Error getting new atom for bond");
  }

  mol.addBond(newBeginAtom, newEndAtom, bond->getBondType());
}

}  // namespace

std::unique_ptr<RDKit::RWMol> MolFromMacroMol(
    const MacroMol *macroMol,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams,
    const RDKit::MolFromMacroMolParams &molFromMacroMolParams) {
  auto mol = std::make_unique<RWMol>();

  // maps main atom# and template atom# to new atom#
  OriginAtomMap originAtomMap;

  // maps main atom# and attach label to template atom#
  AttachMap attachMap;

  // first get some information from the templates to be used when
  // creating the coords for the new atoms. this is a dirty approach
  // that simply expands the orginal macro atom coords to be big
  // enough to hold any expanded macro atom. No attempt is made to
  // make this look nice, or to avoid overlaps.
  std::map<const MacroMolTemplate *, RDGeom::Point3D> templateCentroids;
  double maxSize = 0.0;

  const Conformer *conf = nullptr;
  std::unique_ptr<Conformer> newConf(nullptr);
  if (macroMol->getNumConformers() != 0) {
    conf = &macroMol->getConformer(0);
    newConf.reset(new Conformer(macroMol->getNumAtoms()));
    newConf->set3D(conf->is3D());

    // loop over all main mol atoms and make a list of the templates used
    // external libs could have thousands of tempaltes, most not used by the
    // current macro mol

    std::set<const MacroMolTemplate *> templatesInUse;
    for (unsigned int atomIdx = 0; atomIdx != macroMol->getNumAtoms();
         ++atomIdx) {
      auto templatePtr = macroMol->getTemplate(atomIdx);
      if (templatePtr != nullptr && !templatesInUse.contains(templatePtr)) {
        templatesInUse.insert(templatePtr);
      }
    }

    for (auto templateMol : templatesInUse) {
      RDGeom::Point3D sumOfCoords;
      const RDKit::Conformer *templateConf = nullptr;
      auto confCount = templateMol->getNumConformers();
      if (confCount == 0) {
        conf = nullptr;
        break;
      }
      templateConf = &templateMol->getConformer(0);
      RDGeom::Point3D maxCoord = templateConf->getAtomPos(0);
      RDGeom::Point3D minCoord = maxCoord;
      for (unsigned int atomIdx = 0; atomIdx < templateMol->getNumAtoms();
           ++atomIdx) {
        auto atomCoord = templateConf->getAtomPos(atomIdx);
        sumOfCoords += atomCoord;

        if (atomCoord.x > maxCoord.x) {
          maxCoord.x = atomCoord.x;
        }
        if (atomCoord.y > maxCoord.y) {
          maxCoord.y = atomCoord.y;
        }
        if (atomCoord.z > maxCoord.z) {
          maxCoord.z = atomCoord.z;
        }
        if (atomCoord.x < minCoord.x) {
          minCoord.x = atomCoord.x;
        }
        if (atomCoord.y < minCoord.y) {
          minCoord.y = atomCoord.y;
        }
        if (atomCoord.z < minCoord.z) {
          minCoord.z = atomCoord.z;
        }
      }
      templateCentroids[templateMol] =
          (sumOfCoords / templateMol->getNumAtoms());
      if (maxCoord.x - minCoord.x > maxSize) {
        maxSize = maxCoord.x - minCoord.x;
      }
      if (maxCoord.y - minCoord.y > maxSize) {
        maxSize = maxCoord.y - minCoord.y;
      }
      if (maxCoord.z - minCoord.z > maxSize) {
        maxSize = maxCoord.z - minCoord.z;
      }
    }
  }

  // loop over the bonds in the main mol.  For each bond, record the ATTACHMAP
  // entries.  We do not know the mapped template atom yet, so add a
  // placeholder for that.  We will fill in the real template atom after we
  // have processed the main atoms and know the mapping of main atoms to
  // template atoms.

  for (auto &bond : macroMol->bonds()) {
    std::string lbl;
    if (bond->getPropIfPresent(common_properties::_MolFileBondAttachPt1, lbl)) {
      auto key = std::pair(bond->getBeginAtomIdx(), lbl);
      if (attachMap.find(key) != attachMap.end()) {
        throw FileParseException("Duplicate attachment point label");
      }
      attachMap[key] = UINT_MAX;
    }
    if (bond->getPropIfPresent(common_properties::_MolFileBondAttachPt2, lbl)) {
      auto key = std::pair(bond->getEndAtomIdx(), lbl);
      if (attachMap.find(key) != attachMap.end()) {
        throw FileParseException("Duplicate attachment point label");
      }
      attachMap[key] = UINT_MAX;  // placeholder for now, will fill in real
                                  // template atom idx later
    }
  }

  // for each atom in the main mol, expand it to full atom form

  std::vector<StereoGroup> newStereoGroups;
  std::vector<Atom *> absoluteAtoms;
  std::vector<Bond *> absoluteBonds;

  std::vector<std::unique_ptr<SubstanceGroup>> newSgroups;

  unsigned int atomCount = macroMol->getNumAtoms();
  for (unsigned int atomIdx = 0; atomIdx < atomCount; ++atomIdx) {
    auto atom = macroMol->getAtomWithIdx(atomIdx);
    std::string dummyLabel = "";
    std::string atomClass = "";

    if (!atom->getPropIfPresent(common_properties::dummyLabel, dummyLabel) ||
        !atom->getPropIfPresent(common_properties::molAtomClass, atomClass) ||
        dummyLabel == "" || atomClass == "") {
      // NOT a template atom - just copy it
      auto newAtom = new Atom(*atom);
      mol->addAtom(newAtom, true, true);

      originAtomMap[std::pair(atomIdx, UINT_MAX)] = newAtom->getIdx();
      if (conf != nullptr) {
        newConf->setAtomPos(newAtom->getIdx(),
                            conf->getAtomPos(atomIdx) * maxSize);
      }
    } else {  // it is a macro atom - expand it

      unsigned int seqId = 0;
      std::string seqName = "";
      atom->getPropIfPresent(common_properties::molAtomSeqId, seqId);
      atom->getPropIfPresent(common_properties::molAtomSeqName, seqName);

      auto templateMol = macroMol->getTemplate(atomIdx);
      std::vector<std::string> templateNames;
      std::string templateNameToUse;

      templateMol->getProp<std::vector<std::string>>(
          common_properties::templateNames, templateNames);
      switch (molFromMacroMolParams.macroTemplateNames) {
        case MacroMolTemplateNames::UseFirstName:
          templateNameToUse = templateNames[0];
          break;
        case MacroMolTemplateNames::UseSecondName:
          templateNameToUse = templateNames.back();
          break;
        case MacroMolTemplateNames::AsEntered:
          templateNameToUse = dummyLabel;
          break;
        case MacroMolTemplateNames::All:
          templateNameToUse = "";
          for (const auto &nm : templateNames) {
            if (templateNameToUse != "") {
              templateNameToUse += "+";
            }
            templateNameToUse += nm;
          }
          break;
      }

      // first find the sgroup that is the base for this atom's
      // template

      const SubstanceGroup *sgroup = templateMol->getMainSgroup();

      // add the atoms from the main template to the new molecule

      std::string sgroupName = atomClass + "_";
      if (seqId != 0) {
        sgroupName += std::to_string(seqId) + "_";
      } else {
        sgroupName += "na_";
      }
      if (seqName != "") {
        sgroupName += "_" + seqName;
      }

      sgroupName += templateNameToUse;

      RDGeom::Point3D coordOffset;
      if (conf) {
        coordOffset = (conf->getAtomPos(atomIdx) * maxSize) -
                      templateCentroids[templateMol];
      }

      copySgroupIntoResult(*mol, *macroMol, atomIdx, *sgroup, sgroupName,
                           newSgroups, newConf.get(), coordOffset,
                           originAtomMap);

      // if we are including atoms from leaving groups, go through the
      // attachment points of the main sgroup. If the attach point is
      // not found in the attachMap, then find the sgroup for that
      // attach point and add its atoms the molecule

      for (auto attachPoint : sgroup->getAttachPoints()) {
        auto key = std::pair(atomIdx, attachPoint.id);
        if (attachMap.find(key) != attachMap.end()) {
          // fill in the tempate atom id for this attachPoint

          attachMap[key] = attachPoint.aIdx;
        } else if (molFromMacroMolParams.includeLeavingGroups) {
          // this attach point was not found, so the leaving group is
          // included in the output molecule (if there is one).

          for (auto lgSgroup : getSubstanceGroups(*templateMol)) {
            std::string lgSup;
            std::string lgSgroupAtomClass;
            if (lgSgroup.getPropIfPresent<std::string>("TYPE", lgSup) &&
                lgSup == "SUP" &&
                lgSgroup.getPropIfPresent<std::string>("CLASS",
                                                       lgSgroupAtomClass) &&
                lgSgroupAtomClass == "LGRP") {
              auto lgSgroupAtoms = lgSgroup.getAtoms();
              if (std::find(lgSgroupAtoms.begin(), lgSgroupAtoms.end(),
                            attachPoint.lvIdx) != lgSgroupAtoms.end()) {
                std::string sgroupName = dummyLabel;
                if (seqId != 0) {
                  sgroupName += "_" + std::to_string(seqId);
                }
                if (seqName != "") {
                  sgroupName += "_" + seqName;
                }
                sgroupName += "_" + attachPoint.id;
                copySgroupIntoResult(*mol, *macroMol, atomIdx, lgSgroup,
                                     sgroupName, newSgroups, newConf.get(),
                                     coordOffset, originAtomMap);

                break;
              }
            }
          }
        }
      }

      // copy the bonds of the template into the new molecule
      // if the bonds are between atoms in the new molecule
      // Bonds to atoms in leaving groups that "left" are NOT copied

      for (auto bond : templateMol->bonds()) {
        if (originAtomMap.find(std::pair(atomIdx, bond->getBeginAtomIdx())) ==
                originAtomMap.end() ||
            originAtomMap.find(std::pair(atomIdx, bond->getEndAtomIdx())) ==
                originAtomMap.end()) {
          continue;  // bond not in the new molecule
        }
        auto newBeginAtomIdx =
            originAtomMap[std::pair(atomIdx, bond->getBeginAtomIdx())];
        auto newEndAtomIdx =
            originAtomMap[std::pair(atomIdx, bond->getEndAtomIdx())];

        auto newBond = new Bond(bond->getBondType());
        newBond->setBeginAtomIdx(newBeginAtomIdx);
        newBond->setEndAtomIdx(newEndAtomIdx);
        newBond->updateProps(*bond, false);
        mol->addBond(newBond, true);
      }

      // take care of stereo groups in the template
      // abs groups are added to the list of abs atoms and bonds, so
      // that we can add ONE abs group later

      for (auto &sg : templateMol->getStereoGroups()) {
        std::vector<Atom *> newGroupAtoms;
        std::vector<Bond *> newGroupBonds;

        for (auto sgAtom : sg.getAtoms()) {
          auto originAtom = std::pair(atomIdx, sgAtom->getIdx());

          auto newAtomPtr = originAtomMap.find(originAtom);
          if (newAtomPtr != originAtomMap.end()) {
            newGroupAtoms.push_back(mol->getAtomWithIdx(newAtomPtr->second));
          }
        }

        for (auto sgBond : sg.getBonds()) {
          auto originBeginAtom =
              std::pair(atomIdx, sgBond->getBeginAtomIdx());
          auto originEndAtom = std::pair(atomIdx, sgBond->getEndAtomIdx());

          auto newBeginAtomPtr = originAtomMap.find(originBeginAtom);
          auto newEndAtomPtr = originAtomMap.find(originEndAtom);
          if (newBeginAtomPtr != originAtomMap.end() &&
              newEndAtomPtr != originAtomMap.end()) {
            auto newBond = mol->getBondBetweenAtoms(newBeginAtomPtr->second,
                                                    newEndAtomPtr->second);
            if (newBond != nullptr) {
              newGroupBonds.push_back(newBond);
            }
          }
        }

        if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
          absoluteAtoms.insert(absoluteAtoms.end(), newGroupAtoms.begin(),
                               newGroupAtoms.end());
          absoluteBonds.insert(absoluteBonds.end(), newGroupBonds.begin(),
                               newGroupBonds.end());
        } else {
          // make a new group

          newStereoGroups.emplace_back(sg.getGroupType(), newGroupAtoms,
                                       newGroupBonds);
        }
      }
    }
  }

  if (mol->getNumAtoms() == 0) {
    return std::move(mol);
  }

  if (conf != nullptr) {
    newConf->resize(mol->getNumAtoms());
    mol->addConformer(newConf.release());
  }

  // now deal with the bonds from the original mol.

  for (auto bond : macroMol->bonds()) {
    processBondInMainMol(bond, *mol, originAtomMap, attachMap);
  }

  // copy any attrs from the main mol

  for (auto &prop : macroMol->getPropList(false, false)) {
    std::string propVal;
    if (macroMol->getPropIfPresent(prop, propVal)) {
      mol->setProp(prop, propVal);
    }
  }

  // copy the sgroups from the main mol for atoms not in a template

  if (molFromMacroMolParams.outputSgroups) {
    for (auto &sg : getSubstanceGroups(*macroMol)) {
      if (sg.getIsValid()) {
        auto &atoms = sg.getAtoms();
        std::vector<unsigned int> newAtoms;
        for (auto atom : atoms) {
          auto originAtom = std::pair(atom, UINT_MAX);
          auto newAtomPtr = originAtomMap.find(originAtom);
          if (newAtomPtr != originAtomMap.end()) {
            newAtoms.push_back(newAtomPtr->second);
          } else {
            // some atoms were in templates and others were not - cannot
            // add this sgroup
            newAtoms.clear();
            break;
          }
        }
        if (newAtoms.size() > 0) {
          const std::string type = "SUP";
          newSgroups.emplace_back(new SubstanceGroup(mol.get(), type));
          auto newSg = newSgroups.back().get();

          newSg->updateProps(sg, false);
          newSg->setAtoms(newAtoms);
        }
      }
    }

    // now that we have all substance groups from the template and from
    // the non-template atoms, and we have all the bonds, find the
    // Xbonds for each substance group and add them

    for (auto bond : mol->bonds()) {
      for (auto &sg : newSgroups) {
        // if one atom of the bond is found and the other is not in the
        // sgroup, this is a Xbond
        auto sgAtoms = sg->getAtoms();
        if ((std::find(sgAtoms.begin(), sgAtoms.end(),
                       bond->getBeginAtomIdx()) == sgAtoms.end()) !=
            (std::find(sgAtoms.begin(), sgAtoms.end(), bond->getEndAtomIdx()) ==
             sgAtoms.end())) {
          sg->addBondWithIdx(bond->getIdx());
        }
      }
    }

    if (newSgroups.size() > 0) {
      for (auto &sg : newSgroups) {
        addSubstanceGroup(*mol, *sg.get());
      }
    }
  }
  newSgroups.clear();  // just tidy cleanup

  // take care of stereo groups in the main mol - for atoms that are
  // NOT template refs

  for (auto &sg : macroMol->getStereoGroups()) {
    std::vector<Atom *> newGroupAtoms;
    std::vector<Bond *> newGroupBonds;

    for (auto sgAtom : sg.getAtoms()) {
      auto originAtom = std::pair(sgAtom->getIdx(), UINT_MAX);
      auto newAtomPtr = originAtomMap.find(originAtom);
      if (newAtomPtr != originAtomMap.end()) {
        newGroupAtoms.push_back(mol->getAtomWithIdx(newAtomPtr->second));
      }
    }

    for (auto sgBond : sg.getBonds()) {
      auto originBeginAtom = std::pair(sgBond->getBeginAtomIdx(), UINT_MAX);
      auto originEndAtom = std::pair(sgBond->getEndAtomIdx(), UINT_MAX);
      auto newBeginAtomPtr = originAtomMap.find(originBeginAtom);
      auto newEndAtomPtr = originAtomMap.find(originEndAtom);

      if (newBeginAtomPtr != originAtomMap.end() &&
          newEndAtomPtr != originAtomMap.end()) {
        auto newBond = mol->getBondBetweenAtoms(newBeginAtomPtr->second,
                                                newEndAtomPtr->second);
        if (newBond != nullptr) {
          newGroupBonds.push_back(newBond);
        }
      }

      if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
        absoluteAtoms.insert(absoluteAtoms.end(), newGroupAtoms.begin(),
                             newGroupAtoms.end());
        absoluteBonds.insert(absoluteBonds.end(), newGroupBonds.begin(),
                             newGroupBonds.end());
      } else {
        // make a new group

        newStereoGroups.emplace_back(sg.getGroupType(), newGroupAtoms,
                                     newGroupBonds);
      }
    }
  }

  // make an absolute group that contains any absolute atoms or bonds
  // from either the main mol or the template instantiations

  if (!absoluteAtoms.empty() || !absoluteBonds.empty()) {
    newStereoGroups.emplace_back(StereoGroupType::STEREO_ABSOLUTE, absoluteAtoms,
                                 absoluteBonds);
  }

  if (newStereoGroups.size() > 0) {
    mol->setStereoGroups(newStereoGroups);
  }

  // finishMolProcessing reads getConformer() unconditionally.  When the
  // source MacroMol had no coordinates the result has no conformer, so give
  // it a default 2D conformer (positions are meaningless, but chirality is
  // preserved from the template atom tags rather than from coordinates).
  if (mol->getNumConformers() == 0) {
    auto *defaultConf = new Conformer(mol->getNumAtoms());
    defaultConf->set3D(false);
    mol->addConformer(defaultConf, true);
  }

  bool chiralityPossible = false;
  RDKit::FileParserUtils::finishMolProcessing(mol.get(), chiralityPossible,
                                              molFileParserParams);

  unsigned int failedOp = 0;
  MolOps::sanitizeMol(*(mol.get()), failedOp, MolOps::SANITIZE_KEKULIZE);
  MolOps::sanitizeMol(*(mol.get()), failedOp, MolOps::SANITIZE_SETAROMATICITY);

  return std::move(mol);
}

std::unique_ptr<RDKit::RWMol> MolFromMacroMol(const MacroMol &macroMol) {
  return MolFromMacroMol(&macroMol,
                         RDKit::v2::FileParsers::MolFileParserParams(),
                         RDKit::MolFromMacroMolParams());
}

}  // end of namespace RDKit
