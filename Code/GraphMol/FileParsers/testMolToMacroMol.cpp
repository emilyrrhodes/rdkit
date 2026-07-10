//
// Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/MacroMol.h>
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

  auto macroTemplate = std::make_shared<MacroMolTemplate>(*parsed);
  macroTemplate->setMainGroup(mainAtoms, "OTHER");
  for (const auto &leavingGroup : leavingGroups) {
    macroTemplate->addLeavingGroup(
        leavingGroup.atoms, leavingGroup.attachAtomIdx,
        leavingGroup.leavingAtomIdx, leavingGroup.attachPoint);
  }

  auto entry = std::make_shared<MacroMolEntry>();
  entry->monomerClass = MonomerClass::OTHER;
  entry->templateName = templateName;
  entry->symbol = symbol;
  entry->molTemplate = macroTemplate;
  return entry;
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

    auto mol = SmilesToMol("CCC");
    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    CHECK(macroMol->getNumAtoms() == 1);
    CHECK(macroMol->getNumBonds() == 0);
    const auto *macroInfo = macroMol->getAtomWithIdx(0)->getMacroAtomInfo();
    REQUIRE(macroInfo);
    CHECK(macroInfo->getSymbol() == "L");
    CHECK(macroInfo->getMonomerClass() == MonomerClass::OTHER);
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

    auto mol = SmilesToMol("CC");
    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    REQUIRE(macroMol->getNumAtoms() == 1);
    const auto *macroInfo = macroMol->getAtomWithIdx(0)->getMacroAtomInfo();
    REQUIRE(macroInfo);
    CHECK(macroInfo->getSymbol() == "F");
    CHECK(macroInfo->getMonomerClass() == MonomerClass::OTHER);
  }

  SECTION("single hit creates a macro atom and copies unmatched atoms") {
    MacroMolTemplateLibrary templates;
    templates.addEntry(
        makeMacroMolEntry("ETHYL", "Et", "CCO", {0, 1}, {{{2}, 1, 2, 1}}));

    auto mol = SmilesToMol("CCC");
    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    CHECK(macroMol->getNumAtoms() == 2);
    CHECK(macroMol->getNumBonds() == 1);

    const auto *macroAtom = macroMol->getAtomWithIdx(0);
    const auto *macroInfo = macroAtom->getMacroAtomInfo();
    REQUIRE(macroInfo);
    CHECK(macroInfo->getSymbol() == "Et");
    CHECK(macroInfo->getMonomerClass() == MonomerClass::OTHER);

    const auto *copiedAtom = macroMol->getAtomWithIdx(1);
    CHECK(copiedAtom->getMacroAtomInfo() == nullptr);
    CHECK(copiedAtom->getAtomicNum() == 6);

    const auto *bondInfo = macroMol->getBondWithIdx(0)->getMacroBondInfo();
    REQUIRE(bondInfo);
    REQUIRE(bondInfo->getNumBonds() == 1);
    CHECK(bondInfo->getBond(0).beginAttachPt == 1);
    CHECK(bondInfo->getBond(0).endAttachPt == -1);
    CHECK(bondInfo->getBond(0).bondType ==
          static_cast<unsigned int>(Bond::BondType::SINGLE));
  }

  SECTION(
      "symmetric template orientations are kept for attachment validation") {
    MacroMolTemplateLibrary templates;
    templates.addEntry(
        makeMacroMolEntry("ETHYL", "Et", "CCO", {0, 1}, {{{2}, 1, 2, 1}}));

    auto mol = SmilesToMol("NCC");
    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    REQUIRE(macroMol->getNumAtoms() == 2);
    REQUIRE(macroMol->getNumBonds() == 1);

    CHECK(macroMol->getAtomWithIdx(0)->getAtomicNum() == 7);
    const auto *macroInfo = macroMol->getAtomWithIdx(1)->getMacroAtomInfo();
    REQUIRE(macroInfo);
    CHECK(macroInfo->getSymbol() == "Et");

    const auto *bondInfo = macroMol->getBondWithIdx(0)->getMacroBondInfo();
    REQUIRE(bondInfo);
    REQUIRE(bondInfo->getNumBonds() == 1);
    const auto macroBond = bondInfo->getBond(0);
    CHECK(((macroBond.beginAttachPt == -1 && macroBond.endAttachPt == 1) ||
           (macroBond.beginAttachPt == 1 && macroBond.endAttachPt == -1)));
  }

  SECTION("unannotated external bonds reject a template hit") {
    MacroMolTemplateLibrary templates;
    templates.addEntry(makeMacroMolEntry("ETHYL", "Et", "CC", {0, 1}));

    auto mol = SmilesToMol("CCC");
    std::unique_ptr<MacroMol> macroMol;
    CHECK_NOTHROW(macroMol = MolToMacroMol(*mol, templates));

    REQUIRE(macroMol);
    CHECK(macroMol->getNumAtoms() == 3);
    CHECK(macroMol->getNumBonds() == 2);
    for (const auto *atom : macroMol->atoms()) {
      CHECK(atom->getMacroAtomInfo() == nullptr);
    }
    for (const auto *bond : macroMol->bonds()) {
      CHECK(bond->getMacroBondInfo() == nullptr);
    }
  }

  SECTION("external bonds between template hits become macro bonds") {
    MacroMolTemplateLibrary templates;
    templates.addEntry(makeMacroMolEntry("CO_UNIT", "X", "NCON", {1, 2},
                                         {{{0}, 1, 0, 1}, {{3}, 2, 3, 2}}));

    auto mol = SmilesToMol("COCO");
    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    CHECK(macroMol->getNumAtoms() == 2);
    REQUIRE(macroMol->getNumBonds() == 1);

    const auto *bond = macroMol->getBondWithIdx(0);
    const auto *bondInfo = bond->getMacroBondInfo();
    REQUIRE(bondInfo);
    REQUIRE(bondInfo->getNumBonds() == 1);
    CHECK(bondInfo->getBond(0).beginAttachPt == 2);
    CHECK(bondInfo->getBond(0).endAttachPt == 1);
    CHECK(bondInfo->getBond(0).bondType ==
          static_cast<unsigned int>(Bond::BondType::SINGLE));
  }

  SECTION("entries without usable templates are skipped") {
    MacroMolTemplateLibrary templates;

    auto nullTemplate = std::make_shared<MacroMolEntry>();
    nullTemplate->monomerClass = MonomerClass::AA;
    nullTemplate->templateName = "NULL";
    nullTemplate->symbol = "N";
    templates.addEntry(nullTemplate);

    auto noMainSgroup = std::make_shared<MacroMolEntry>();
    noMainSgroup->monomerClass = MonomerClass::AA;
    noMainSgroup->templateName = "NO_MAIN";
    noMainSgroup->symbol = "M";
    noMainSgroup->molTemplate = std::make_shared<MacroMolTemplate>();
    templates.addEntry(noMainSgroup);

    auto mol = SmilesToMol("C");
    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    CHECK(macroMol->getNumAtoms() == 1);
    // The unmatched atom is carried through unchanged.
    CHECK(macroMol->getAtomWithIdx(0)->getAtomicNum() == 6);
  }

  SECTION("re-adding the same entry does not duplicate it") {
    MacroMolTemplateLibrary templates;
    auto entry = makeMacroMolEntry("DUP", "D", "CC", {0, 1});

    templates.addEntry(entry);
    templates.addEntry(entry);

    CHECK(templates.entries().size() == 1);
  }

  SECTION("an attachment atom may have at most one external bond") {
    // A single-atom main group with one attachment point. When it matches an
    // atom that carries two bonds leaving the monomer, the hit is rejected
    // rather than reusing the single attach point for both crossing bonds.
    MacroMolTemplateLibrary templates;
    templates.addEntry(
        makeMacroMolEntry("METHYL", "Me", "CO", {0}, {{{1}, 0, 1, 1}}));

    auto mol = SmilesToMol("CCC");
    std::unique_ptr<MacroMol> macroMol;
    CHECK_NOTHROW(macroMol = MolToMacroMol(*mol, templates));

    REQUIRE(macroMol);
    REQUIRE(macroMol->getNumAtoms() == 3);
    // The two terminal carbons (one external bond each) become macro atoms; the
    // central carbon (two external bonds, one attach point) is copied plainly.
    unsigned int macroAtomCount = 0;
    const Atom *plainAtom = nullptr;
    for (const auto *atom : macroMol->atoms()) {
      if (atom->getMacroAtomInfo() != nullptr) {
        ++macroAtomCount;
      } else {
        plainAtom = atom;
      }
    }
    CHECK(macroAtomCount == 2);
    REQUIRE(plainAtom != nullptr);
    CHECK(plainAtom->getAtomicNum() == 6);
    CHECK(plainAtom->getDegree() == 2);
  }

  SECTION("regular bonds keep their metadata when nothing matches") {
    MacroMolTemplateLibrary templates;  // empty: everything copied as-is

    auto mol = std::unique_ptr<RWMol>(SmilesToMol("C/C=C/C"));
    REQUIRE(mol);
    // Tag the double bond with a custom property to confirm it survives.
    auto *doubleBond = mol->getBondBetweenAtoms(1, 2);
    REQUIRE(doubleBond);
    REQUIRE(doubleBond->getStereo() != Bond::BondStereo::STEREONONE);
    doubleBond->setProp<int>("customTag", 7);

    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    CHECK(macroMol->getNumAtoms() == 4);
    CHECK(macroMol->getNumBonds() == 3);
    for (const auto *bond : macroMol->bonds()) {
      CHECK(bond->getMacroBondInfo() == nullptr);
    }

    const auto *copiedDouble = macroMol->getBondBetweenAtoms(1, 2);
    REQUIRE(copiedDouble);
    CHECK(copiedDouble->getBondType() == Bond::BondType::DOUBLE);
    CHECK(copiedDouble->getStereo() == doubleBond->getStereo());
    REQUIRE(copiedDouble->getStereoAtoms().size() == 2);
    int copiedTag = 0;
    REQUIRE(copiedDouble->getPropIfPresent<int>("customTag", copiedTag));
    CHECK(copiedTag == 7);
  }

  SECTION("copied regular atoms preserve query state") {
    MacroMolTemplateLibrary templates;  // empty

    auto mol = std::unique_ptr<RWMol>(SmartsToMol("[#6]-[#8]"));
    REQUIRE(mol);
    REQUIRE(mol->getAtomWithIdx(0)->hasQuery());
    REQUIRE(mol->getAtomWithIdx(1)->hasQuery());

    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    REQUIRE(macroMol->getNumAtoms() == 2);
    CHECK(macroMol->getAtomWithIdx(0)->hasQuery());
    CHECK(macroMol->getAtomWithIdx(1)->hasQuery());
  }

  SECTION("copied regular bonds preserve query state") {
    MacroMolTemplateLibrary templates;  // empty

    auto mol = std::unique_ptr<RWMol>(SmartsToMol("CC"));
    REQUIRE(mol);
    REQUIRE(mol->getBondWithIdx(0)->hasQuery());

    auto macroMol = MolToMacroMol(*mol, templates);

    REQUIRE(macroMol);
    REQUIRE(macroMol->getNumBonds() == 1);
    CHECK(macroMol->getBondWithIdx(0)->hasQuery());
  }

  SECTION("addLeavingGroup rejects an attachment atom outside the main group") {
    auto parsed = std::unique_ptr<RWMol>(SmilesToMol("CCO"));
    MacroMolTemplate templ(*parsed);
    templ.setMainGroup({0, 1}, "OTHER");
    // Atom 2 is not part of the main group {0, 1}.
    CHECK_THROWS(templ.addLeavingGroup({2}, 2, 2, 1));
  }
}
