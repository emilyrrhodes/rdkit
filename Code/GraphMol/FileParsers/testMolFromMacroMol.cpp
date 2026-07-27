//
// Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/Atom.h>
#include <GraphMol/FileParsers/MolFromMacroMol.h>
#include <GraphMol/MacroMol.h>
#include <GraphMol/MacroMolTemplate.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <catch2/catch_all.hpp>

#include <memory>

using namespace RDKit;

namespace {

constexpr bool doIsomericSmiles = true;
constexpr bool doKekule = false;
constexpr int rootedAtAtom = 0;
constexpr bool canonical = false;

std::shared_ptr<MacroMolTemplate> makeAlanineTemplate() {
  SmilesParserParams params;
  params.removeHs = false;
  auto alanine = std::unique_ptr<RWMol>(
      SmilesToMol("C[C@H](N[H:1])C(=O)[OH:2]", params));
  auto alanineTemplate = std::make_shared<MacroMolTemplate>(*alanine);
  alanineTemplate->setMainGroup({0, 1, 2, 4, 5}, MonomerClass::AminoAcid);
  alanineTemplate->addLeavingGroup({3}, 2, 3, 1);
  alanineTemplate->addLeavingGroup({6}, 4, 6, 2);
  return alanineTemplate;
}

std::shared_ptr<MacroMolTemplate> makeGlycineTemplate() {
  SmilesParserParams params;
  params.removeHs = false;
  auto glycine =
      std::unique_ptr<RWMol>(SmilesToMol("O=C(CN[H:1])[OH:2]", params));
  auto glycineTemplate = std::make_shared<MacroMolTemplate>(*glycine);
  glycineTemplate->setMainGroup({0, 1, 2, 3}, MonomerClass::AminoAcid);
  glycineTemplate->addLeavingGroup({4}, 3, 4, 1);
  glycineTemplate->addLeavingGroup({5}, 1, 5, 2);
  return glycineTemplate;
}

void addMacroMolTemplateEntry(
    MacroMolTemplateLibrary &templates, const char *templateName,
    const char *symbol, const std::shared_ptr<MacroMolTemplate> &macroTemplate) {
  auto entry = std::make_shared<MacroMolEntry>();
  entry->monomerClass = MonomerClass::AminoAcid;
  entry->templateName = templateName;
  entry->symbol = symbol;
  entry->molTemplate = macroTemplate;
  templates.addEntry(entry);
}

}  // namespace

TEST_CASE("MolFromMacroMol converts a single amino acid macro atom",
          "[MolFromMacroMol]") {
  MacroMolTemplateLibrary templates;
  addMacroMolTemplateEntry(templates, "ALA", "A", makeAlanineTemplate());

  MacroMol macroMol;
  macroMol.addMacroAtom("A", MonomerClass::AminoAcid);

  auto mol = MolFromMacroMol(macroMol, templates);

  REQUIRE(mol);
  CHECK(mol->getNumAtoms() == 7);
  CHECK(mol->getNumBonds() == 6);
  CHECK(MolToSmiles(*mol, doIsomericSmiles, doKekule, rootedAtAtom,
                    canonical) ==
        "C[C@H](N[H:1])C(=O)[OH:2]");
}

TEST_CASE("MolFromMacroMol converts connected amino acid macro atoms",
          "[MolFromMacroMol]") {
  MacroMolTemplateLibrary templates;
  addMacroMolTemplateEntry(templates, "ALA", "A", makeAlanineTemplate());
  addMacroMolTemplateEntry(templates, "GLY", "G", makeGlycineTemplate());

  MacroMol macroMol;
  const auto alanineMacroAtom = macroMol.addMacroAtom("A", MonomerClass::AminoAcid);
  const auto glycineMacroAtom = macroMol.addMacroAtom("G", MonomerClass::AminoAcid);
  macroMol.addMacroBond(alanineMacroAtom, glycineMacroAtom, 2, 1);

  auto mol = MolFromMacroMol(macroMol, templates);

  REQUIRE(mol);
  CHECK(mol->getNumAtoms() == 11);
  CHECK(mol->getNumBonds() == 10);
  CHECK(MolToSmiles(*mol, doIsomericSmiles, doKekule, rootedAtAtom,
                    canonical) ==
        "C[C@H](N[H:1])C(=O)NCC(=O)[OH:2]");
}

TEST_CASE("MolFromMacroMol converts mixed macro and atomistic atoms",
          "[MolFromMacroMol]") {
  MacroMolTemplateLibrary templates;
  addMacroMolTemplateEntry(templates, "ALA", "A", makeAlanineTemplate());

  MacroMol macroMol;
  const auto macroAtom = macroMol.addMacroAtom("A", MonomerClass::AminoAcid);
  const auto atomisticAtom = macroMol.addAtom(new Atom(6), false, true);
  macroMol.addMacroAtomToAtomBond(macroAtom, atomisticAtom, 1);

  auto mol = MolFromMacroMol(macroMol, templates);

  REQUIRE(mol);
  CHECK(mol->getNumAtoms() == 7);
  CHECK(mol->getNumBonds() == 6);
  CHECK(MolToSmiles(*mol, doIsomericSmiles, doKekule, rootedAtAtom,
                    canonical) ==
        "C[C@H](NC)C(=O)[OH:2]");
}
