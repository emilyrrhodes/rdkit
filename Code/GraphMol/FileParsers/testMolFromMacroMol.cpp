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

#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace RDKit;

namespace {

constexpr bool doIsomericSmiles = true;
constexpr bool doKekule = false;
constexpr int rootedAtAtom = 0;
constexpr bool canonical = false;

struct ExpectedBuiltinTemplate {
  const char *symbol;
  const char *templateName;
  MonomerClass monomerClass;
};

const std::array<ExpectedBuiltinTemplate, 22> expectedBuiltinTemplates{{
    {"A", "ALA", MonomerClass::AminoAcid},
    {"R", "ARG", MonomerClass::AminoAcid},
    {"N", "ASN", MonomerClass::AminoAcid},
    {"D", "ASP", MonomerClass::AminoAcid},
    {"C", "CYS", MonomerClass::AminoAcid},
    {"Q", "GLN", MonomerClass::AminoAcid},
    {"E", "GLU", MonomerClass::AminoAcid},
    {"G", "GLY", MonomerClass::AminoAcid},
    {"H", "HIS", MonomerClass::AminoAcid},
    {"I", "ILE", MonomerClass::AminoAcid},
    {"L", "LEU", MonomerClass::AminoAcid},
    {"K", "LYS", MonomerClass::AminoAcid},
    {"M", "MET", MonomerClass::AminoAcid},
    {"F", "PHE", MonomerClass::AminoAcid},
    {"P", "PRO", MonomerClass::AminoAcid},
    {"S", "SER", MonomerClass::AminoAcid},
    {"T", "THR", MonomerClass::AminoAcid},
    {"W", "TRP", MonomerClass::AminoAcid},
    {"Y", "TYR", MonomerClass::AminoAcid},
    {"V", "VAL", MonomerClass::AminoAcid},
    {"U", "SEC", MonomerClass::AminoAcid},
    {"O", "PYL", MonomerClass::AminoAcid},
}};

std::unique_ptr<MacroMolTemplate> makeAlanineTemplate() {
  SmilesParserParams params;
  params.removeHs = false;
  const std::string smiles = "C[C@H](N[H:1])C(=O)[OH:2]";
  auto alanine = std::unique_ptr<RWMol>(SmilesToMol(smiles, params));
  auto alanineTemplate = std::make_unique<MacroMolTemplate>(
      *alanine, MonomerClass::AminoAcid, "ALA", "A", smiles);
  alanineTemplate->setMainGroup({0, 1, 2, 4, 5});
  alanineTemplate->addLeavingGroup({3}, 2, 3, 1);
  alanineTemplate->addLeavingGroup({6}, 4, 6, 2);
  return alanineTemplate;
}

std::unique_ptr<MacroMolTemplate> makeGlycineTemplate() {
  SmilesParserParams params;
  params.removeHs = false;
  const std::string smiles = "O=C(CN[H:1])[OH:2]";
  auto glycine = std::unique_ptr<RWMol>(SmilesToMol(smiles, params));
  auto glycineTemplate = std::make_unique<MacroMolTemplate>(
      *glycine, MonomerClass::AminoAcid, "GLY", "G", smiles);
  glycineTemplate->setMainGroup({0, 1, 2, 3});
  glycineTemplate->addLeavingGroup({4}, 3, 4, 1);
  glycineTemplate->addLeavingGroup({5}, 1, 5, 2);
  return glycineTemplate;
}

}  // namespace

TEST_CASE("MolFromMacroMol converts a single amino acid macro atom",
          "[MolFromMacroMol]") {
  MacroMolTemplateLibrary templates;
  templates.addTemplate(makeAlanineTemplate());

  MacroMol macroMol;
  macroMol.addMacroAtom("A", MonomerClass::AminoAcid);

  auto mol = MolFromMacroMol(macroMol, templates);

  REQUIRE(mol);
  CHECK(mol->getNumAtoms() == 7);
  CHECK(mol->getNumBonds() == 6);
  CHECK(MolToSmiles(*mol, doIsomericSmiles, doKekule, rootedAtAtom,
                    canonical) == "C[C@H](N[H:1])C(=O)[OH:2]");
}

TEST_CASE("MolFromMacroMol converts connected amino acid macro atoms",
          "[MolFromMacroMol]") {
  MacroMolTemplateLibrary templates;
  templates.addTemplate(makeAlanineTemplate());
  templates.addTemplate(makeGlycineTemplate());

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
  templates.addTemplate(makeAlanineTemplate());

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

TEST_CASE("Global MacroMol template library contains built-in amino acids",
          "[MolFromMacroMol]") {
  const auto &templates = getGlobalMacroMolTemplateLibrary();

  for (const auto &[symbol, templateName, monomerClass] :
       expectedBuiltinTemplates) {
    const auto *bySymbol = templates.getBySymbol(monomerClass, symbol);
    INFO("symbol: " << symbol);
    REQUIRE(bySymbol);
    CHECK(bySymbol->getSymbol() == symbol);
    CHECK(bySymbol->getTemplateName() == templateName);
    CHECK(bySymbol->getMonomerClass() == monomerClass);

    const auto *byTemplateName =
        templates.getByTemplateName(monomerClass, templateName);
    INFO("templateName: " << templateName);
    REQUIRE(byTemplateName);
    CHECK(byTemplateName == bySymbol);
  }
}

TEST_CASE("addBuiltinMacroMolTemplates populates another library",
          "[MolFromMacroMol]") {
  MacroMolTemplateLibrary templates;
  addBuiltinMacroMolTemplates(templates);

  const auto *alanine =
      templates.getBySymbol(MonomerClass::AminoAcid, "A");
  REQUIRE(alanine);
  CHECK(alanine->getTemplateName() == "ALA");

  const auto *glycine =
      templates.getByTemplateName(MonomerClass::AminoAcid, "GLY");
  REQUIRE(glycine);
  CHECK(glycine->getSymbol() == "G");

  MacroMolTemplateLibrary secondTemplates;
  addBuiltinMacroMolTemplates(secondTemplates);
  const auto *secondAlanine =
      secondTemplates.getBySymbol(MonomerClass::AminoAcid, "A");
  REQUIRE(secondAlanine);
  CHECK(secondAlanine != alanine);
}

TEST_CASE("MolFromMacroMol uses the global amino-acid template library",
          "[MolFromMacroMol]") {
  MacroMol macroMol;
  macroMol.addMacroAtom("A", MonomerClass::AminoAcid);

  auto mol = MolFromMacroMol(macroMol);

  REQUIRE(mol);
  CHECK(mol->getNumAtoms() == 7);
  CHECK(mol->getNumBonds() == 6);
  CHECK(MolToSmiles(*mol, doIsomericSmiles, doKekule, rootedAtAtom,
                    canonical) == "C[C@H](N[H])C(=O)O");
}

TEST_CASE("MolFromMacroMol expands peptides with global built-ins",
          "[MolFromMacroMol]") {
  MacroMol macroMol;
  const auto alanineMacroAtom = macroMol.addMacroAtom("A", MonomerClass::AminoAcid);
  const auto glycineMacroAtom = macroMol.addMacroAtom("G", MonomerClass::AminoAcid);
  macroMol.addMacroBond(alanineMacroAtom, glycineMacroAtom, 2, 1);

  auto mol = MolFromMacroMol(macroMol);

  REQUIRE(mol);
  CHECK(mol->getNumAtoms() == 11);
  CHECK(mol->getNumBonds() == 10);
  CHECK(MolToSmiles(*mol, doIsomericSmiles, doKekule, rootedAtAtom,
                    canonical) == "C[C@H](N[H])C(=O)NCC(=O)O");
}

TEST_CASE("Built-in MacroMol templates use direct group definitions",
          "[MolFromMacroMol]") {
  const auto &templates = getGlobalMacroMolTemplateLibrary();
  const auto *alanineTemplate =
      templates.getByTemplateName(MonomerClass::AminoAcid, "ALA");
  REQUIRE(alanineTemplate);

  for (const auto *atom : alanineTemplate->atoms()) {
    CHECK(atom->getAtomMapNum() == 0);
  }

  const auto *mainSgroup = alanineTemplate->getMainSgroup();
  REQUIRE(mainSgroup != nullptr);
  CHECK(mainSgroup->getProp<std::string>("TYPE") == "SUP");
  CHECK(mainSgroup->getProp<std::string>("CLASS") == "AminoAcid");
  CHECK(mainSgroup->getAtoms() == std::vector<unsigned int>({0, 1, 2, 4, 5}));

  auto leavingGroups = alanineTemplate->getLeavingGroups();
  REQUIRE(leavingGroups.size() == 2);
  CHECK(leavingGroups[0]->getProp<std::string>("TYPE") == "SUP");
  CHECK(leavingGroups[0]->getProp<std::string>("CLASS") == "LGRP");
  CHECK(leavingGroups[0]->getAtoms() == std::vector<unsigned int>({3}));
  CHECK(leavingGroups[1]->getProp<std::string>("TYPE") == "SUP");
  CHECK(leavingGroups[1]->getProp<std::string>("CLASS") == "LGRP");
  CHECK(leavingGroups[1]->getAtoms() == std::vector<unsigned int>({6}));

  const auto &attachPoints = mainSgroup->getAttachPoints();
  REQUIRE(attachPoints.size() == 2);
  CHECK(attachPoints[0].aIdx == 2);
  CHECK(attachPoints[0].lvIdx == 3);
  CHECK(attachPoints[0].id == "1");
  CHECK(attachPoints[1].aIdx == 4);
  CHECK(attachPoints[1].lvIdx == 6);
  CHECK(attachPoints[1].id == "2");
}

TEST_CASE("Built-in MacroMol templates include side-chain leaving groups",
          "[MolFromMacroMol]") {
  const auto &templates = getGlobalMacroMolTemplateLibrary();
  const auto *cysteineTemplate =
      templates.getByTemplateName(MonomerClass::AminoAcid, "CYS");
  REQUIRE(cysteineTemplate);

  auto leavingGroups = cysteineTemplate->getLeavingGroups();
  REQUIRE(leavingGroups.size() == 3);
  CHECK(leavingGroups[2]->getAtoms() == std::vector<unsigned int>({5}));

  const auto *mainSgroup = cysteineTemplate->getMainSgroup();
  REQUIRE(mainSgroup != nullptr);
  const auto &attachPoints = mainSgroup->getAttachPoints();
  REQUIRE(attachPoints.size() == 3);
  CHECK(attachPoints[2].aIdx == 4);
  CHECK(attachPoints[2].lvIdx == 5);
  CHECK(attachPoints[2].id == "3");
}
