//
// Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/Conformer.h>
#include <GraphMol/FileParsers/MolToMacroMol.h>
#include <GraphMol/MacroMolTemplate.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <catch2/catch_all.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace RDKit;

namespace {

struct LeavingGroupDef {
  std::vector<unsigned int> atoms;
  unsigned int attachAtomIdx;
  unsigned int leavingAtomIdx;
  int attachPoint;
};

std::shared_ptr<MacroMolEntry> makeMacroMolEntry(
    const std::string &templateName, const std::string &symbol,
    const std::string &smiles, const std::vector<unsigned int> &mainAtoms,
    const std::vector<LeavingGroupDef> &leavingGroups = {}) {
  auto parsed = std::unique_ptr<RWMol>(SmilesToMol(smiles));
  if (!parsed) {
    throw std::runtime_error("could not parse test template smiles");
  }

  auto macroTemplate = std::make_shared<MacroMolTemplate>(*parsed);
  macroTemplate->setMainGroup(mainAtoms, "AA");
  for (const auto &leavingGroup : leavingGroups) {
    macroTemplate->addLeavingGroup(leavingGroup.atoms,
                                   leavingGroup.attachAtomIdx,
                                   leavingGroup.leavingAtomIdx,
                                   leavingGroup.attachPoint);
  }

  auto entry = std::make_shared<MacroMolEntry>();
  entry->monomerClass = "AA";
  entry->templateName = templateName;
  entry->symbol = symbol;
  entry->molTemplate = macroTemplate;
  return entry;
}

std::unique_ptr<RWMol> smilesToMolForMacroMolTest(const std::string &smiles) {
  auto mol = std::unique_ptr<RWMol>(SmilesToMol(smiles));
  if (!mol) {
    throw std::runtime_error("could not parse test molecule smiles");
  }
  return mol;
}

}  // namespace

TEST_CASE("MolToMacroMol minimal converter", "[MolToMacroMol]") {
  SECTION("larger template is tried before a smaller overlapping template") {
    MacroMolTemplateLibrary templates;
    auto small =
        makeMacroMolEntry("SMALL", "S", "CCO", {0, 1}, {{{2}, 1, 2, 1}});
    auto large = makeMacroMolEntry("LARGE", "L", "CCC", {0, 1, 2});

    templates.addEntry(small);
    templates.addEntry(large);

    REQUIRE(templates.entries().size() == 2);
    CHECK(templates.entries()[0] == large);
    CHECK(templates.entries()[1] == small);

    auto mol = smilesToMolForMacroMolTest("CCC");
    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    CHECK(macroMol->getNumAtoms() == 1);
    CHECK(macroMol->getNumBonds() == 0);
    const auto *macroInfo = macroMol->getAtomWithIdx(0)->getMacroAtomInfo();
    REQUIRE(macroInfo);
    CHECK(macroInfo->getSymbol() == "L");
    CHECK(macroInfo->getMonomerClass() == "AA");
  }

  SECTION("equal-size templates preserve insertion order") {
    MacroMolTemplateLibrary templates;
    auto first = makeMacroMolEntry("FIRST", "F", "CC", {0, 1});
    auto second = makeMacroMolEntry("SECOND", "S", "CC", {0, 1});

    templates.addEntry(first);
    templates.addEntry(second);

    REQUIRE(templates.entries().size() == 2);
    CHECK(templates.entries()[0] == first);
    CHECK(templates.entries()[1] == second);

    auto mol = smilesToMolForMacroMolTest("CC");
    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    REQUIRE(macroMol->getNumAtoms() == 1);
    const auto *macroInfo = macroMol->getAtomWithIdx(0)->getMacroAtomInfo();
    REQUIRE(macroInfo);
    CHECK(macroInfo->getSymbol() == "F");
    CHECK(macroInfo->getMonomerClass() == "AA");
  }

  SECTION("single hit creates a macro atom and copies unmatched atoms") {
    MacroMolTemplateLibrary templates;
    templates.addEntry(
        makeMacroMolEntry("ETHYL", "Et", "CCO", {0, 1}, {{{2}, 1, 2, 1}}));

    auto mol = smilesToMolForMacroMolTest("CCC");
    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    CHECK(macroMol->getNumAtoms() == 2);
    CHECK(macroMol->getNumBonds() == 1);

    const auto *macroAtom = macroMol->getAtomWithIdx(0);
    const auto *macroInfo = macroAtom->getMacroAtomInfo();
    REQUIRE(macroInfo);
    CHECK(macroInfo->getSymbol() == "Et");
    CHECK(macroInfo->getMonomerClass() == "AA");

    const auto *copiedAtom = macroMol->getAtomWithIdx(1);
    CHECK(copiedAtom->getMacroAtomInfo() == nullptr);
    CHECK(copiedAtom->getAtomicNum() == 6);
  }

  SECTION("external bonds between template hits become macro bonds") {
    MacroMolTemplateLibrary templates;
    templates.addEntry(makeMacroMolEntry("CO_UNIT", "X", "NCON", {1, 2},
                                         {{{0}, 1, 0, 1}, {{3}, 2, 3, 2}}));

    auto mol = smilesToMolForMacroMolTest("COCO");
    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    CHECK(macroMol->getNumAtoms() == 2);
    REQUIRE(macroMol->getNumBonds() == 1);

    const auto *bond = macroMol->getBondWithIdx(0);
    CHECK(bond->getProp<int>(common_properties::_MacroMolBeginAttachPt) == 2);
    CHECK(bond->getProp<int>(common_properties::_MacroMolEndAttachPt) == 1);
  }

  SECTION("entries without usable templates are skipped") {
    MacroMolTemplateLibrary templates;

    auto nullTemplate = std::make_shared<MacroMolEntry>();
    nullTemplate->monomerClass = "AA";
    nullTemplate->templateName = "NULL";
    nullTemplate->symbol = "N";
    templates.addEntry(nullTemplate);

    auto noMainSgroup = std::make_shared<MacroMolEntry>();
    noMainSgroup->monomerClass = "AA";
    noMainSgroup->templateName = "NO_MAIN";
    noMainSgroup->symbol = "M";
    noMainSgroup->molTemplate = std::make_shared<MacroMolTemplate>();
    templates.addEntry(noMainSgroup);

    auto mol = smilesToMolForMacroMolTest("C");
    mol->addConformer(new Conformer(mol->getNumAtoms()));

    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    CHECK(macroMol->getNumAtoms() == 1);
    CHECK(macroMol->getAtomWithIdx(0)->getAtomicNum() == 6);
    CHECK(macroMol->getNumConformers() == 0);
  }
}
