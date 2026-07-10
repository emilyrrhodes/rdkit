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
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/FileParseException.h>

#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <vector>

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

struct BuiltinLeavingGroup {
  std::vector<unsigned int> atomIdxs;
  unsigned int attachAtomIdx = 0;
  unsigned int leavingAtomIdx = 0;
  int attachPoint = 0;
};

struct BuiltinTemplateDef {
  const char *symbol;
  const char *templateName;
  MonomerClass monomerClass;
  const char *smiles;
  std::vector<unsigned int> mainAtomIdxs;
  std::vector<BuiltinLeavingGroup> leavingGroups;
};

const std::array<BuiltinTemplateDef, 22> builtinTemplates{{
    {"A", "ALA", MonomerClass::AminoAcid,
     "C[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CB :1.pdbName. CA "
     ":2.pdbName. N  :3.pdbName. H  :4.pdbName. C  :5.pdbName. O  "
     ":6.pdbName. OXT|",
     {0, 1, 2, 4, 5},
     {{{3}, 2, 3, 1}, {{6}, 4, 6, 2}}},
    {"R", "ARG", MonomerClass::AminoAcid,
     "N=C(N)NCCC[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. "
     "NH2:1.pdbName. CZ :2.pdbName. NH1:3.pdbName. NE :4.pdbName. CD "
     ":5.pdbName. CG :6.pdbName. CB :7.pdbName. CA :8.pdbName. N  "
     ":9.pdbName. H  :10.pdbName. C  :11.pdbName. O  :12.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11},
     {{{9}, 8, 9, 1}, {{12}, 10, 12, 2}}},
    {"N", "ASN", MonomerClass::AminoAcid,
     "NC(=O)C[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. ND2:1.pdbName. CG "
     ":2.pdbName. OD1:3.pdbName. CB :4.pdbName. CA :5.pdbName. N  :6.pdbName. "
     "H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 7, 8},
     {{{6}, 5, 6, 1}, {{9}, 7, 9, 2}}},
    {"D", "ASP", MonomerClass::AminoAcid,
     "O=C([C@H](CC(=O)[OH])N[H])[OH] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
     "OD1:6.pdbName. OD2:7.pdbName. N  :8.pdbName. H  :9.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 7},
     {{{8}, 7, 8, 1}, {{9}, 1, 9, 2}, {{6}, 4, 6, 3}}},
    {"C", "CYS", MonomerClass::AminoAcid,
     "O=C([C@H](CS[H])N[H])[OH] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. SG "
     ":5.pdbName. HG :6.pdbName. N  :7.pdbName. H  :8.pdbName. OXT|",
     {0, 1, 2, 3, 4, 6},
     {{{7}, 6, 7, 1}, {{8}, 1, 8, 2}, {{5}, 4, 5, 3}}},
    {"Q", "GLN", MonomerClass::AminoAcid,
     "NC(=O)CC[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. NE2:1.pdbName. CD "
     ":2.pdbName. OE1:3.pdbName. CG :4.pdbName. CB :5.pdbName. CA :6.pdbName. "
     "N  :7.pdbName. H  :8.pdbName. C  :9.pdbName. O  :10.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 6, 8, 9},
     {{{7}, 6, 7, 1}, {{10}, 8, 10, 2}}},
    {"E", "GLU", MonomerClass::AminoAcid,
     "O=C([C@H](CCC(=O)[OH])N[H])[OH] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
     ":5.pdbName. CD :6.pdbName. OE1:7.pdbName. OE2:8.pdbName. N  "
     ":9.pdbName. H  :10.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 6, 8},
     {{{9}, 8, 9, 1}, {{10}, 1, 10, 2}, {{7}, 5, 7, 3}}},
    {"G", "GLY", MonomerClass::AminoAcid,
     "O=C(CN[H])[OH] |atomProp:0.pdbName. O  :1.pdbName. C  "
     ":2.pdbName. CA :3.pdbName. N  :4.pdbName. H  :5.pdbName. OXT|",
     {0, 1, 2, 3},
     {{{4}, 3, 4, 1}, {{5}, 1, 5, 2}}},
    {"H", "HIS", MonomerClass::AminoAcid,
     "O=C([C@H](Cc1cnc[nH]1)N[H])[OH] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
     ":5.pdbName. CD2:6.pdbName. NE2:7.pdbName. CE1:8.pdbName. "
     "ND1:9.pdbName. N  :10.pdbName. H  :11.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
     {{{10}, 9, 10, 1}, {{11}, 1, 11, 2}}},
    {"I", "ILE", MonomerClass::AminoAcid,
     "CC[C@H](C)[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CD1:1.pdbName. "
     "CG1:2.pdbName. CB :3.pdbName. CG2:4.pdbName. CA :5.pdbName. N  "
     ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 7, 8},
     {{{6}, 5, 6, 1}, {{9}, 7, 9, 2}}},
    {"L", "LEU", MonomerClass::AminoAcid,
     "CC(C)C[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CD1:1.pdbName. "
     "CG :2.pdbName. CD2:3.pdbName. CB :4.pdbName. CA :5.pdbName. N  "
     ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 7, 8},
     {{{6}, 5, 6, 1}, {{9}, 7, 9, 2}}},
    {"K", "LYS", MonomerClass::AminoAcid,
     "O=C([C@H](CCCCN[H])N[H])[OH] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
     ":5.pdbName. CD :6.pdbName. CE :7.pdbName. NZ :8.pdbName. "
     "HZ1:9.pdbName. N  :10.pdbName. H  :11.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 6, 7, 9},
     {{{10}, 9, 10, 1}, {{11}, 1, 11, 2}, {{8}, 7, 8, 3}}},
    {"M", "MET", MonomerClass::AminoAcid,
     "CSCC[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CE :1.pdbName. "
     "SD :2.pdbName. CG :3.pdbName. CB :4.pdbName. CA :5.pdbName. N  "
     ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 7, 8},
     {{{6}, 5, 6, 1}, {{9}, 7, 9, 2}}},
    {"F", "PHE", MonomerClass::AminoAcid,
     "O=C([C@H](Cc1ccccc1)N[H])[OH] |atomProp:0.pdbName. O  :1.pdbName. C "
     " :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
     "CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. CE2:9.pdbName. "
     "CD2:10.pdbName. N  :11.pdbName. H  :12.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
     {{{11}, 10, 11, 1}, {{12}, 1, 12, 2}}},
    {"P", "PRO", MonomerClass::AminoAcid,
     "O=C([C@@H]1CCCN1[H])[OH] |atomProp:0.pdbName. O  :1.pdbName. C "
     " :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD "
     ":6.pdbName. N  :7.pdbName. H  :8.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 6},
     {{{7}, 6, 7, 1}, {{8}, 1, 8, 2}}},
    {"S", "SER", MonomerClass::AminoAcid,
     "O=C([C@H](CO)N[H])[OH] |atomProp:0.pdbName. O  :1.pdbName. C  "
     ":2.pdbName. CA :3.pdbName. CB :4.pdbName. OG :5.pdbName. N  "
     ":6.pdbName. H  :7.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5},
     {{{6}, 5, 6, 1}, {{7}, 1, 7, 2}}},
    {"T", "THR", MonomerClass::AminoAcid,
     "C[C@@H](O)[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. "
     "CG2:1.pdbName. CB :2.pdbName. OG1:3.pdbName. CA :4.pdbName. N  "
     ":5.pdbName. H  :6.pdbName. C  :7.pdbName. O  :8.pdbName. OXT|",
     {0, 1, 2, 3, 4, 6, 7},
     {{{5}, 4, 5, 1}, {{8}, 6, 8, 2}}},
    {"W", "TRP", MonomerClass::AminoAcid,
     "O=C([C@H](Cc1c[nH]c2ccccc12)N[H])[OH] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
     ":5.pdbName. CD1:6.pdbName. NE1:7.pdbName. CE2:8.pdbName. "
     "CZ2:9.pdbName. CH2:10.pdbName. CZ3:11.pdbName. CE3:12.pdbName. "
     "CD2:13.pdbName. N  :14.pdbName. H  :15.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13},
     {{{14}, 13, 14, 1}, {{15}, 1, 15, 2}}},
    {"Y", "TYR", MonomerClass::AminoAcid,
     "O=C([C@H](Cc1ccc(O)cc1)N[H])[OH] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
     "CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. OH :9.pdbName. "
     "CE2:10.pdbName. CD2:11.pdbName. N  :12.pdbName. H  :13.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
     {{{12}, 11, 12, 1}, {{13}, 1, 13, 2}}},
    {"V", "VAL", MonomerClass::AminoAcid,
     "CC(C)[C@H](N[H])C(=O)[OH] |atomProp:0.pdbName. CG1:1.pdbName. "
     "CB :2.pdbName. CG2:3.pdbName. CA :4.pdbName. N  :5.pdbName. H  "
     ":6.pdbName. C  :7.pdbName. O  :8.pdbName. OXT|",
     {0, 1, 2, 3, 4, 6, 7},
     {{{5}, 4, 5, 1}, {{8}, 6, 8, 2}}},
    {"U", "SEC", MonomerClass::AminoAcid,
     "O=C([C@H](C[SeH])N[H])[OH] |atomProp:0.pdbName. O  :1.pdbName. "
     "C  :2.pdbName. CA :3.pdbName. CB :4.pdbName.SE  :5.pdbName. N  "
     ":6.pdbName. H  :7.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5},
     {{{6}, 5, 6, 1}, {{7}, 1, 7, 2}}},
    {"O", "PYL", MonomerClass::AminoAcid,
     "C[C@@H]1CC=N[C@H]1C(=O)NCCCC[C@H](N[H])C(=O)[OH] "
     "|atomProp:0.pdbName. CB2:1.pdbName. CG2:2.pdbName. CD2:3.pdbName. "
     "CE2:4.pdbName. N2 :5.pdbName. CA2:6.pdbName. C2 :7.pdbName. O2 "
     ":8.pdbName. NZ :9.pdbName. CE :10.pdbName. CD :11.pdbName. CG "
     ":12.pdbName. CB :13.pdbName. CA :14.pdbName. N  :15.pdbName. H  "
     ":16.pdbName. C  :17.pdbName. O  :18.pdbName. OXT|",
     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17},
     {{{15}, 14, 15, 1}, {{18}, 16, 18, 2}}},
}};

std::unique_ptr<RWMol> parseBuiltinTemplateMol(const char *smiles) {
  SmilesParserParams params;
  params.removeHs = false;
  auto mol = std::unique_ptr<RWMol>(SmilesToMol(smiles, params));
  if (!mol) {
    throw FileParseException("Could not parse built-in MacroMol template");
  }
  return mol;
}

void setTemplateGroups(MacroMolTemplate &templ,
                       const BuiltinTemplateDef &builtinTemplate) {
  templ.setMainGroup(builtinTemplate.mainAtomIdxs);
  for (const auto &leavingGroup : builtinTemplate.leavingGroups) {
    templ.addLeavingGroup(leavingGroup.atomIdxs, leavingGroup.attachAtomIdx,
                          leavingGroup.leavingAtomIdx,
                          leavingGroup.attachPoint);
  }
}

std::unique_ptr<MacroMolTemplate> makeBuiltinTemplate(
    const BuiltinTemplateDef &builtinTemplate) {
  auto parsedMol = parseBuiltinTemplateMol(builtinTemplate.smiles);

  auto templ = std::make_unique<MacroMolTemplate>(
      *parsedMol, builtinTemplate.monomerClass, builtinTemplate.templateName,
      builtinTemplate.symbol, builtinTemplate.smiles);
  setTemplateGroups(*templ, builtinTemplate);
  return templ;
}

const std::array<std::unique_ptr<const MacroMolTemplate>,
                 builtinTemplates.size()> &
getBuiltinTemplates() {
  static const auto templates = [] {
    std::array<std::unique_ptr<const MacroMolTemplate>,
               builtinTemplates.size()>
        result{};
    for (std::size_t i = 0; i < builtinTemplates.size(); ++i) {
      result[i] = makeBuiltinTemplate(builtinTemplates[i]);
    }
    return result;
  }();
  return templates;
}

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

MacroMolTemplateLibrary &getGlobalMacroMolTemplateLibrary() {
  static MacroMolTemplateLibrary templates;
  static std::once_flag initialized;
  std::call_once(initialized, []() { addBuiltinMacroMolTemplates(templates); });
  return templates;
}

void addBuiltinMacroMolTemplates(MacroMolTemplateLibrary &templates) {
  for (const auto &builtinTemplate : getBuiltinTemplates()) {
    templates.addTemplate(
        std::make_unique<MacroMolTemplate>(*builtinTemplate));
  }
}

std::unique_ptr<RWMol> MolFromMacroMol(const MacroMol &macroMol) {
  return MolFromMacroMol(macroMol, getGlobalMacroMolTemplateLibrary());
}

}  // namespace RDKit
