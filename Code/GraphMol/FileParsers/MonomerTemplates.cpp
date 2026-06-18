//
//  Copyright (C) 2026 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MonomerTemplates.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <GraphMol/RWMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/MacroMolTemplate.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace RDKit {

namespace {
// Built-in amino-acid definitions (symbol -> map-numbered SMILES).  Each
// map-numbered atom is a single-atom leaving group; the map number is the
// connection label.  The CXSMILES extensions carry PDB atom names.
const std::unordered_map<std::string, std::string> builtin_monomer_data({
    {"A", "C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CB :1.pdbName. CA "
          ":2.pdbName. N  :3.pdbName. H  :4.pdbName. C  :5.pdbName. O  "
          ":6.pdbName. OXT|"},
    {"C", "O=C([C@H](CS[H:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. SG "
          ":5.pdbName. HG :6.pdbName. N  :7.pdbName. H  :8.pdbName. OXT|"},
    {"D",
     "O=C([C@H](CC(=O)[OH:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
     "OD1:6.pdbName. OD2:7.pdbName. N  :8.pdbName. H  :9.pdbName. OXT|"},
    {"E", "O=C([C@H](CCC(=O)[OH:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
          ":5.pdbName. CD :6.pdbName. OE1:7.pdbName. OE2:8.pdbName. N  "
          ":9.pdbName. H  :10.pdbName. OXT|"},
    {"F",
     "O=C([C@H](Cc1ccccc1)N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C "
     " :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
     "CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. CE2:9.pdbName. "
     "CD2:10.pdbName. N  :11.pdbName. H  :12.pdbName. OXT|"},
    {"G", "O=C(CN[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C  "
          ":2.pdbName. CA :3.pdbName. N  :4.pdbName. H  :5.pdbName. OXT|"},
    {"H", "O=C([C@H](Cc1cnc[nH]1)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
          ":5.pdbName. CD2:6.pdbName. NE2:7.pdbName. CE1:8.pdbName. "
          "ND1:9.pdbName. N  :10.pdbName. H  :11.pdbName. OXT|"},
    {"I",
     "CC[C@H](C)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CD1:1.pdbName. "
     "CG1:2.pdbName. CB :3.pdbName. CG2:4.pdbName. CA :5.pdbName. N  "
     ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|"},
    {"K", "O=C([C@H](CCCCN[H:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
          ":5.pdbName. CD :6.pdbName. CE :7.pdbName. NZ :8.pdbName. "
          "HZ1:9.pdbName. N  :10.pdbName. H  :11.pdbName. OXT|"},
    {"L", "CC(C)C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CD1:1.pdbName. "
          "CG :2.pdbName. CD2:3.pdbName. CB :4.pdbName. CA :5.pdbName. N  "
          ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|"},
    {"M", "CSCC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CE :1.pdbName. "
          "SD :2.pdbName. CG :3.pdbName. CB :4.pdbName. CA :5.pdbName. N  "
          ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|"},
    {"N",
     "NC(=O)C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. ND2:1.pdbName. CG "
     ":2.pdbName. OD1:3.pdbName. CB :4.pdbName. CA :5.pdbName. N  :6.pdbName. "
     "H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|"},
    {"O", "C[C@@H]1CC=N[C@H]1C(=O)NCCCC[C@H](N[H:1])C(=O)[OH:2] "
          "|atomProp:0.pdbName. CB2:1.pdbName. CG2:2.pdbName. CD2:3.pdbName. "
          "CE2:4.pdbName. N2 :5.pdbName. CA2:6.pdbName. C2 :7.pdbName. O2 "
          ":8.pdbName. NZ :9.pdbName. CE :10.pdbName. CD :11.pdbName. CG "
          ":12.pdbName. CB :13.pdbName. CA :14.pdbName. N  :15.pdbName. H  "
          ":16.pdbName. C  :17.pdbName. O  :18.pdbName. OXT|"},
    {"P", "O=C([C@@H]1CCCN1[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C "
          " :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD "
          ":6.pdbName. N  :7.pdbName. H  :8.pdbName. OXT|"},
    {"Q",
     "NC(=O)CC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. NE2:1.pdbName. CD "
     ":2.pdbName. OE1:3.pdbName. CG :4.pdbName. CB :5.pdbName. CA :6.pdbName. "
     "N  :7.pdbName. H  :8.pdbName. C  :9.pdbName. O  :10.pdbName. OXT|"},
    {"R", "N=C(N)NCCC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. "
          "NH2:1.pdbName. CZ :2.pdbName. NH1:3.pdbName. NE :4.pdbName. CD "
          ":5.pdbName. CG :6.pdbName. CB :7.pdbName. CA :8.pdbName. N  "
          ":9.pdbName. H  :10.pdbName. C  :11.pdbName. O  :12.pdbName. OXT|"},
    {"S", "O=C([C@H](CO)N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C  "
          ":2.pdbName. CA :3.pdbName. CB :4.pdbName. OG :5.pdbName. N  "
          ":6.pdbName. H  :7.pdbName. OXT|"},
    {"T", "C[C@@H](O)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. "
          "CG2:1.pdbName. CB :2.pdbName. OG1:3.pdbName. CA :4.pdbName. N  "
          ":5.pdbName. H  :6.pdbName. C  :7.pdbName. O  :8.pdbName. OXT|"},
    {"U", "O=C([C@H](C[SeH])N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. "
          "C  :2.pdbName. CA :3.pdbName. CB :4.pdbName.SE  :5.pdbName. N  "
          ":6.pdbName. H  :7.pdbName. OXT|"},
    {"V", "CC(C)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CG1:1.pdbName. "
          "CB :2.pdbName. CG2:3.pdbName. CA :4.pdbName. N  :5.pdbName. H  "
          ":6.pdbName. C  :7.pdbName. O  :8.pdbName. OXT|"},
    {"W", "O=C([C@H](Cc1c[nH]c2ccccc12)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
          ":5.pdbName. CD1:6.pdbName. NE1:7.pdbName. CE2:8.pdbName. "
          "CZ2:9.pdbName. CH2:10.pdbName. CZ3:11.pdbName. CE3:12.pdbName. "
          "CD2:13.pdbName. N  :14.pdbName. H  :15.pdbName. OXT|"},
    {"Y",
     "O=C([C@H](Cc1ccc(O)cc1)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
     "CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. OH :9.pdbName. "
     "CE2:10.pdbName. CD2:11.pdbName. N  :12.pdbName. H  :13.pdbName. OXT|"}});
}  // anonymous namespace

std::unique_ptr<MacroMolTemplate> buildAminoAcidTemplate(
    const std::string &symbol, const std::string &smiles) {
  RWMol *m = SmilesToMol(smiles, 0, /*sanitize=*/false);
  if (m == nullptr) {
    throw ValueErrorException("could not parse builtin amino-acid SMILES for " +
                              symbol);
  }

  // collect the leaving atoms (those with a non-zero atom map number); each is
  // a single-atom leaving group connected to a single core "attach" atom.
  struct LeavingInfo {
    std::string id;
    unsigned int lvIdx;
    unsigned int aIdx;
  };
  std::vector<LeavingInfo> leavings;
  boost::dynamic_bitset<> isLeaving(m->getNumAtoms());

  for (auto atom : m->atoms()) {
    int mapNum = atom->getAtomMapNum();
    if (mapNum == 0) {
      continue;
    }
    unsigned int lvIdx = atom->getIdx();
    // the leaving atom has exactly one neighbor: the attach atom
    auto neighbors = m->atomNeighbors(atom);
    CHECK_INVARIANT(neighbors.begin() != neighbors.end(),
                    "leaving atom has no neighbors");
    unsigned int aIdx = (*neighbors.begin())->getIdx();

    leavings.push_back({std::to_string(mapNum), lvIdx, aIdx});
    isLeaving.set(lvIdx);
    // clear so the map number doesn't leak into output SMILES
    atom->setAtomMapNum(0);
  }

  std::vector<unsigned int> coreAtoms;
  coreAtoms.reserve(m->getNumAtoms() - leavings.size());
  for (auto atom : m->atoms()) {
    if (!isLeaving[atom->getIdx()]) {
      coreAtoms.push_back(atom->getIdx());
    }
  }

  // main SUP sgroup: core atoms + an attachment point per leaving atom
  {
    SubstanceGroup sg(m, "SUP");
    sg.setProp<std::string>("CLASS", "AminoAcid");
    sg.setAtoms(coreAtoms);
    for (const auto &lv : leavings) {
      sg.addAttachPoint(lv.aIdx, static_cast<int>(lv.lvIdx), lv.id);
    }
    addSubstanceGroup(*m, sg);
  }

  // one LGRP sgroup per leaving atom
  for (const auto &lv : leavings) {
    SubstanceGroup sg(m, "SUP");
    sg.setProp<std::string>("CLASS", "LGRP");
    sg.setAtoms(std::vector<unsigned int>{lv.lvIdx});
    addSubstanceGroup(*m, sg);
  }

  std::unique_ptr<RWMol> up(m);
  return std::make_unique<MacroMolTemplate>(
      up, std::string("AminoAcid"), symbol,
      std::vector<std::pair<std::string, std::string>>{});
}

void loadBuiltinMacroMolTemplates(MacroMolTemplateLib &lib) {
  for (const auto &[symbol, smiles] : builtin_monomer_data) {
    auto tmpl = buildAminoAcidTemplate(symbol, smiles);
    lib.addTemplate(tmpl);
  }
}

namespace {
// Register the builtin loader so the global library is populated lazily the
// first time MacroMolTemplateLib::getGlobalLibrary() is called.  GraphMol
// cannot call SmilesToMol itself, so the FileParsers layer supplies the loader.
struct BuiltinTemplateRegistrar {
  BuiltinTemplateRegistrar() {
    MacroMolTemplateLib::setGlobalLibraryLoader(&loadBuiltinMacroMolTemplates);
  }
};
const BuiltinTemplateRegistrar g_builtinTemplateRegistrar;
}  // anonymous namespace

}  // namespace RDKit
