//
// Copyright (C) 2026 Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/Conversions.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MacroMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <catch2/catch_all.hpp>

using namespace RDKit;

TEST_CASE("testBuildMacroMol") {
  // Build a simple MacroMol with three macro atoms and two bonds, and check
  // that the MacroMol has the expected number of atoms and bonds, and that the
  // "sequence" of template names in the macro atoms is as expected.
  auto macro_mol = std::make_unique<MacroMol>();
  auto macro_atom_1 = macro_mol->addMacroAtom(MonomerClass::AA, "A");
  auto macro_atom_2 = macro_mol->addMacroAtom(MonomerClass::AA, "C");
  auto macro_atom_3 = macro_mol->addMacroAtom(MonomerClass::AA, "D");
  macro_mol->addMacroBond(macro_atom_1, macro_atom_2, 2, 1);
  macro_mol->addMacroBond(macro_atom_2, macro_atom_3, 2, 1);
  CHECK(macro_mol->getNumAtoms() == 3);
  CHECK(macro_mol->getNumBonds() == 2);
  std::string sequence;
  for (const auto &atom : macro_mol->atoms()) {
    std::string templateName =
        atom->getProp<std::string>(common_properties::dummyLabel);
    sequence += templateName;
    CHECK(atom->getProp<std::string>(common_properties::molAtomClass) == "AA");
  }
  CHECK(sequence == "ACD");
}

TEST_CASE("testSubstructureSearchWithMacroMols") {
  // Load in a PDB file as a MacroMol, build a simple MacroMol that is a subset
  // of the first MacroMol, and check that a substructure search finds the
  // expected match.
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/1DNG.pdb";
  std::unique_ptr<RWMol> mol(PDBFileToMol(fname));
  std::unique_ptr<RDKit::MacroMol> pdb_macro_mol =
      RDKit::MolToMacroMol(*mol, MacroMolTemplateLib::getGlobalLibrary());

  for (const auto &atom : pdb_macro_mol->atoms()) {
    if (!atom->hasProp(common_properties::dummyLabel) ||
        !atom->hasProp(common_properties::molAtomClass)) {
      std::cerr << "Atom " << atom->getIdx()
                << " is missing required properties for MacroMol: "
                << common_properties::dummyLabel << " or "
                << common_properties::molAtomClass << std::endl;
      continue;
    }
    std::string templateName =
        atom->getProp<std::string>(common_properties::dummyLabel);
    std::string className =
        atom->getProp<std::string>(common_properties::molAtomClass);
    std::cerr << "Atom " << atom->getIdx() << " template name: " << templateName
              << " class name: " << className << std::endl;
  }

  auto simple_macro_mol = std::make_unique<MacroMol>();
  auto macro_atom_1 = simple_macro_mol->addMacroAtom(MonomerClass::AA, "A");
  auto macro_atom_2 = simple_macro_mol->addMacroAtom(MonomerClass::AA, "Y");
  auto macro_atom_3 = simple_macro_mol->addMacroAtom(MonomerClass::AA, "E");
  simple_macro_mol->addMacroBond(macro_atom_1, macro_atom_2, 2, 1);
  simple_macro_mol->addMacroBond(macro_atom_2, macro_atom_3, 2, 1);

  RDKit::GeneralizedSubstruct::ExtendedQueryMol query(
      std::make_unique<RWMol>(*simple_macro_mol));
  auto match =
      RDKit::GeneralizedSubstruct::hasSubstructMatch(*pdb_macro_mol, query);
  CHECK(match);
}

TEST_CASE("testMacroMolSubstructureMatching") {
  // Build MacroMols with the same number of atoms and bonds, but different
  // attachment points, template names, and monomer classes, and check that
  // substructure searches behave as expected.
  auto macro_mol_1 = std::make_unique<MacroMol>();
  auto mm_1_atom_1 = macro_mol_1->addMacroAtom(MonomerClass::AA, "A");
  auto mm_1_atom_2 = macro_mol_1->addMacroAtom(MonomerClass::AA, "C");
  auto mm_1_atom_3 = macro_mol_1->addMacroAtom(MonomerClass::AA, "D");
  macro_mol_1->addMacroBond(mm_1_atom_1, mm_1_atom_2, 2, 1);
  macro_mol_1->addMacroBond(mm_1_atom_2, mm_1_atom_3, 2, 1);

  auto macro_mol_2 = std::make_unique<MacroMol>();
  auto mm_2_atom_1 = macro_mol_2->addMacroAtom(MonomerClass::AA, "A");
  auto mm_2_atom_2 = macro_mol_2->addMacroAtom(MonomerClass::AA, "C");
  macro_mol_2->addMacroBond(
      mm_2_atom_1, mm_2_atom_2, 2,
      1);  // NOTE: SAME attachment point location as macro_mol_1

  auto macro_mol_3 = std::make_unique<MacroMol>();
  auto mm_3_atom_1 = macro_mol_3->addMacroAtom(MonomerClass::AA, "A");
  auto mm_3_atom_2 = macro_mol_3->addMacroAtom(MonomerClass::AA, "C");
  macro_mol_3->addMacroBond(
      mm_3_atom_1, mm_3_atom_2, 2,
      3);  // NOTE: DIFFERENT attachment point location as macro_mol_1

  auto macro_mol_4 = std::make_unique<MacroMol>();
  auto mm_4_atom_1 = macro_mol_4->addMacroAtom(MonomerClass::AA, "G");
  // NOTE: DIFFERENT template name compared to macro_mol_1
  auto mm_4_atom_2 = macro_mol_4->addMacroAtom(MonomerClass::AA, "C");
  macro_mol_4->addMacroBond(mm_4_atom_1, mm_4_atom_2, 2, 1);

  auto macro_mol_5 = std::make_unique<MacroMol>();
  auto mm_5_atom_1 = macro_mol_5->addMacroAtom(MonomerClass::NA, "A");
  // NOTE: DIFFERENT monomer class compared to macro_mol_1
  auto mm_5_atom_2 = macro_mol_5->addMacroAtom(MonomerClass::AA, "C");
  macro_mol_5->addMacroBond(mm_5_atom_1, mm_5_atom_2, 2, 1);

  RDKit::GeneralizedSubstruct::ExtendedQueryMol query_2(
      std::make_unique<RWMol>(*macro_mol_2));
  auto match_1_2 =
      RDKit::GeneralizedSubstruct::hasSubstructMatch(*macro_mol_1, query_2);
  CHECK(match_1_2);
  RDKit::GeneralizedSubstruct::ExtendedQueryMol query_3(
      std::make_unique<RWMol>(*macro_mol_3));
  auto match_1_3 =
      RDKit::GeneralizedSubstruct::hasSubstructMatch(*macro_mol_1, query_3);
  CHECK(
      !match_1_3);  // This should fail because a change in the attachment point
                    // location should be considered a change chemical structure
  RDKit::GeneralizedSubstruct::ExtendedQueryMol query_4(
      std::make_unique<RWMol>(*macro_mol_4));
  auto match_1_4 =
      RDKit::GeneralizedSubstruct::hasSubstructMatch(*macro_mol_1, query_4);
  CHECK(!match_1_4);  // This should fail because a change in the template name
                      // should be considered a change chemical structure
  RDKit::GeneralizedSubstruct::ExtendedQueryMol query_5(
      std::make_unique<RWMol>(*macro_mol_5));
  auto match_1_5 =
      RDKit::GeneralizedSubstruct::hasSubstructMatch(*macro_mol_1, query_5);
  CHECK(!match_1_5);  // This should fail because a change in the monomer class
                      // should be considered a change chemical structure
}

TEST_CASE("testMacroMolToAtomisticMol") {
  // Build a simple MacroMol with three macro atoms and two bonds, convert it to
  // an atomistic mol, and check that the resulting atomistic mol has the
  // expected number of atoms and bonds.
  auto macro_mol = std::make_unique<MacroMol>();
  auto macro_atom_1 = macro_mol->addMacroAtom(MonomerClass::AA, "A");
  auto macro_atom_2 = macro_mol->addMacroAtom(MonomerClass::AA, "C");
  auto macro_atom_3 = macro_mol->addMacroAtom(MonomerClass::AA, "D");
  macro_mol->addMacroBond(macro_atom_1, macro_atom_2, 2, 1);
  macro_mol->addMacroBond(macro_atom_2, macro_atom_3, 2, 1);

  // outputSgroups=false so the leaving-group hydrogens at unused attachment
  // points (the N-terminal amine H and the cysteine thiol H) are collapsed to
  // implicit Hs by default H removal rather than being protected by the
  // per-residue superatom sgroups.
  MolFromMacroMolParams molFromMacroMolParams;
  molFromMacroMolParams.outputSgroups = false;
  auto atomistic_mol = RDKit::MolFromMacroMol(
      macro_mol.get(), RDKit::v2::FileParsers::MolFileParserParams(),
      molFromMacroMolParams);
  CHECK(atomistic_mol->getNumAtoms() == 20);
  CHECK(atomistic_mol->getNumBonds() == 19);
  CHECK(atomistic_mol->getNumConformers() == 0);
  auto smiles = RDKit::MolToSmiles(*atomistic_mol);
  // compare canonical-to-canonical so the check is independent of the input
  // SMILES atom ordering
  std::unique_ptr<ROMol> expected(RDKit::SmilesToMol(
      "N[C@@H](C)C(=O)N[C@@H](CS)C(=O)N[C@@H](CC(=O)O)C(=O)O"));
  CHECK(smiles == RDKit::MolToSmiles(*expected));
}

TEST_CASE("testMacroMolCanonicalize") {
  // Build two MacroMols with the same template atoms and bonds, but different
  // orders of the macro atoms, canonicalize them, and check that the canonical
  // SMILES for the two MacroMols are the same after canonicalization.
  auto macro_mol_1 = std::make_unique<MacroMol>();
  auto mm_1_atom_1 = macro_mol_1->addMacroAtom(MonomerClass::AA, "A");
  auto mm_1_atom_2 = macro_mol_1->addMacroAtom(MonomerClass::AA, "C");
  auto mm_1_atom_3 = macro_mol_1->addMacroAtom(MonomerClass::AA, "D");
  macro_mol_1->addMacroBond(mm_1_atom_1, mm_1_atom_2, 2, 1);
  macro_mol_1->addMacroBond(mm_1_atom_2, mm_1_atom_3, 2, 1);

  auto macro_mol_2 = std::make_unique<MacroMol>();
  auto mm_2_atom_1 = macro_mol_2->addMacroAtom(MonomerClass::AA, "C");
  auto mm_2_atom_2 = macro_mol_2->addMacroAtom(MonomerClass::AA, "D");
  auto mm_2_atom_3 = macro_mol_2->addMacroAtom(MonomerClass::AA, "A");
  macro_mol_2->addMacroBond(mm_2_atom_3, mm_2_atom_1, 2, 1);
  macro_mol_2->addMacroBond(mm_2_atom_1, mm_2_atom_2, 2, 1);

  auto can_smiles_1 = RDKit::MolToSmiles(*macro_mol_1);
  auto can_smiles_2 = RDKit::MolToSmiles(*macro_mol_2);
  CHECK(can_smiles_1 == can_smiles_2);
}
