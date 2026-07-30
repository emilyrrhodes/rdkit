//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/FileParsers/MolFromMacroMol.h>

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/FileParseException.h>

#include <map>
#include <memory>
#include <set>
#include <string>

namespace RDKit {
namespace {

struct MacroAtomDetails {
  MonomerClass monomerClass = MonomerClass::Other;
  std::string label;
};

struct ExpandedMacroAtom {
  std::map<unsigned int, unsigned int> templateToResultAtom;
  std::map<std::string, unsigned int> attachPointToTemplateAtom;
};

struct BuildState {
  RWMol &result;
  std::map<unsigned int, unsigned int> atomToResultAtom;
  std::map<unsigned int, ExpandedMacroAtom> macroAtoms;
};

bool getMacroAtomDetails(const Atom &atom, MacroAtomDetails &info) {
  const auto *macroInfo = atom.getMacroAtomInfo();
  if (!macroInfo) {
    return false;
  }
  info.monomerClass = macroInfo->getMonomerClass();
  info.label = macroInfo->getSymbol();
  return !info.label.empty();
}

std::string getAttachPointId(int attachPoint) {
  if (attachPoint < 0) {
    return "";
  }
  return std::to_string(attachPoint);
}

const MacroMolTemplate *findTemplate(
    const MacroMolTemplateLibrary &templates, const MacroAtomDetails &info) {
  const auto *templ = templates.getBySymbol(info.monomerClass, info.label);
  if (!templ) {
    templ = templates.getByName(info.monomerClass, info.label);
  }
  return templ;
}

void addBondCopy(RWMol &result, const Bond &oldBond,
                 unsigned int beginAtomIdx, unsigned int endAtomIdx,
                 Bond::BondType bondType) {
  auto *newBond = new Bond(bondType);
  newBond->setBeginAtomIdx(beginAtomIdx);
  newBond->setEndAtomIdx(endAtomIdx);
  newBond->setBondDir(oldBond.getBondDir());
  newBond->updateProps(oldBond, false);
  result.addBond(newBond, true);
}

void addBondCopy(RWMol &result, const Bond &oldBond,
                 unsigned int beginAtomIdx, unsigned int endAtomIdx) {
  addBondCopy(result, oldBond, beginAtomIdx, endAtomIdx,
              oldBond.getBondType());
}

void copyRegularAtom(BuildState &state, const Atom &atom) {
  auto *newAtom = new Atom(atom);
  state.atomToResultAtom[atom.getIdx()] =
      state.result.addAtom(newAtom, true, true);
}

std::set<unsigned int> findLeavingAtomsToSkip(const MacroMol &macroMol,
                                              const Atom &macroAtom,
                                              const MacroMolTemplate &templ) {
  std::set<unsigned int> atomsToSkip;
  const auto *mainSgroup = templ.getMainSgroup();

  for (const auto *macroBond : macroMol.bonds()) {
    const auto *macroBondInfo = macroBond->getMacroBondInfo();
    if (!macroBondInfo) {
      continue;
    }

    for (const auto &bondInfo : macroBondInfo->getBonds()) {
      std::string attachPointId;
      if (macroBond->getBeginAtomIdx() == macroAtom.getIdx()) {
        attachPointId = getAttachPointId(bondInfo.beginAttachPt);
      } else if (macroBond->getEndAtomIdx() == macroAtom.getIdx()) {
        attachPointId = getAttachPointId(bondInfo.endAttachPt);
      }

      if (attachPointId.empty()) {
        continue;
      }

      for (const auto &attachPoint : mainSgroup->getAttachPoints()) {
        if (attachPoint.id == attachPointId && attachPoint.lvIdx >= 0) {
          atomsToSkip.insert(static_cast<unsigned int>(attachPoint.lvIdx));
        }
      }
    }
  }

  return atomsToSkip;
}

void rememberTemplateAttachPoints(ExpandedMacroAtom &expanded,
                                  const MacroMolTemplate &templ) {
  const auto *mainSgroup = templ.getMainSgroup();
  for (const auto &attachPoint : mainSgroup->getAttachPoints()) {
    expanded.attachPointToTemplateAtom[attachPoint.id] = attachPoint.aIdx;
  }
}

void copyTemplateAtoms(BuildState &state, ExpandedMacroAtom &expanded,
                       const MacroMolTemplate &templ,
                       const std::set<unsigned int> &atomsToSkip) {
  for (const auto *oldAtom : templ.atoms()) {
    if (atomsToSkip.find(oldAtom->getIdx()) != atomsToSkip.end()) {
      continue;
    }

    auto *newAtom = new Atom(*oldAtom);
    expanded.templateToResultAtom[oldAtom->getIdx()] =
        state.result.addAtom(newAtom, true, true);
  }
}

void copyTemplateBonds(BuildState &state, const ExpandedMacroAtom &expanded,
                       const MacroMolTemplate &templ) {
  for (const auto *oldBond : templ.bonds()) {
    const auto beginAtom =
        expanded.templateToResultAtom.find(oldBond->getBeginAtomIdx());
    const auto endAtom =
        expanded.templateToResultAtom.find(oldBond->getEndAtomIdx());
    if (beginAtom == expanded.templateToResultAtom.end() ||
        endAtom == expanded.templateToResultAtom.end()) {
      continue;
    }
    addBondCopy(state.result, *oldBond, beginAtom->second, endAtom->second);
  }
}

void copyMacroAtom(BuildState &state, const MacroMol &macroMol,
                   const MacroMolTemplateLibrary &templates, const Atom &atom,
                   const MacroAtomDetails &info) {
  const auto *templ = findTemplate(templates, info);
  if (!templ || !templ->getMainSgroup()) {
    throw FileParseException("No template found for macro atom " +
                             std::to_string(atom.getIdx()));
  }

  ExpandedMacroAtom expanded;
  const auto atomsToSkip = findLeavingAtomsToSkip(macroMol, atom, *templ);

  rememberTemplateAttachPoints(expanded, *templ);
  copyTemplateAtoms(state, expanded, *templ, atomsToSkip);
  copyTemplateBonds(state, expanded, *templ);

  state.macroAtoms[atom.getIdx()] = expanded;
}

void copyAtoms(BuildState &state, const MacroMol &macroMol,
               const MacroMolTemplateLibrary &templates) {
  for (const auto *atom : macroMol.atoms()) {
    MacroAtomDetails info;
    if (getMacroAtomDetails(*atom, info)) {
      copyMacroAtom(state, macroMol, templates, *atom, info);
    } else {
      copyRegularAtom(state, *atom);
    }
  }
}

unsigned int getResultAtomForBond(const BuildState &state, const Atom &atom,
                                  const std::string &attachPointId) {
  MacroAtomDetails info;
  if (!getMacroAtomDetails(atom, info)) {
    return state.atomToResultAtom.at(atom.getIdx());
  }

  const auto &expanded = state.macroAtoms.at(atom.getIdx());
  const auto attachPoint =
      expanded.attachPointToTemplateAtom.find(attachPointId);
  if (attachPoint == expanded.attachPointToTemplateAtom.end()) {
    throw FileParseException("Macro atom bond is missing an attachment point");
  }

  return expanded.templateToResultAtom.at(attachPoint->second);
}

void copyMacroMolBonds(BuildState &state, const MacroMol &macroMol) {
  for (const auto *oldBond : macroMol.bonds()) {
    const auto *macroBondInfo = oldBond->getMacroBondInfo();
    if (!macroBondInfo) {
      const auto beginAtomIdx =
          getResultAtomForBond(state, *oldBond->getBeginAtom(), "");
      const auto endAtomIdx =
          getResultAtomForBond(state, *oldBond->getEndAtom(), "");
      addBondCopy(state.result, *oldBond, beginAtomIdx, endAtomIdx);
      continue;
    }

    for (const auto &bondInfo : macroBondInfo->getBonds()) {
      const auto beginAtomIdx = getResultAtomForBond(
          state, *oldBond->getBeginAtom(),
          getAttachPointId(bondInfo.beginAttachPt));
      const auto endAtomIdx = getResultAtomForBond(
          state, *oldBond->getEndAtom(), getAttachPointId(bondInfo.endAttachPt));
      addBondCopy(state.result, *oldBond, beginAtomIdx, endAtomIdx,
                  static_cast<Bond::BondType>(bondInfo.bondType));
    }
  }
}

}  // namespace

std::unique_ptr<RWMol> MolFromMacroMol(
    const MacroMol &macroMol, const MacroMolTemplateLibrary &templates) {
  auto result = std::make_unique<RWMol>();
  BuildState state{*result, {}, {}};

  copyAtoms(state, macroMol, templates);
  copyMacroMolBonds(state, macroMol);

  result->updatePropertyCache(false);
  return result;
}

}  // namespace RDKit
