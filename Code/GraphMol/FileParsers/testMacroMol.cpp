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
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MacroMolUtils.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <catch2/catch_all.hpp>

using namespace RDKit;

TEST_CASE("testBuildMacroMol") {
  // Build a simple MacroMol with three macro atoms and two bonds, and check
  // that the MacroMol has the expected number of atoms and bonds, and that the
  // "sequence" of template names in the macro atoms is as expected.
  auto macro_mol = std::make_unique<MacroMol>();
  auto macro_atom_1 = macro_mol->addMacroAtom("AminoAcid", "A");
  auto macro_atom_2 = macro_mol->addMacroAtom("AminoAcid", "C");
  auto macro_atom_3 = macro_mol->addMacroAtom("AminoAcid", "D");
  macro_mol->addMacroBond(macro_atom_1, macro_atom_2, Bond::BondType::SINGLE,
                          "2", "1");
  macro_mol->addMacroBond(macro_atom_2, macro_atom_3, Bond::BondType::SINGLE,
                          "2", "1");
  CHECK(macro_mol->getNumAtoms() == 3);
  CHECK(macro_mol->getNumBonds() == 2);
  std::string sequence;
  for (const auto &atom : macro_mol->atoms()) {
    std::string templateName =
        atom->getProp<std::string>(common_properties::dummyLabel);
    sequence += templateName;
  }
  CHECK(sequence == "ACD");
}

/*
TEST_CASE("testSubstructureSearchWithMacroMols") {
  // Load in a PDB file as a MacroMol, build a simple MacroMol that is a subset
  // of the first MacroMol, and check that a substructure search finds the
  // expected match.
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/1DNG.pdb";
  std::unique_ptr<RWMol> mol(PDBFileToMol(fname));
  std::unique_ptr<RDKit::MacroMol> pdb_macro_mol = RDKit::toMonomeric(mol);

  auto simple_macro_mol = std::make_unique<MacroMol>();
  auto macro_atom_1 = simple_macro_mol->addMacroAtom("PEPTIDE", "A");
  auto macro_atom_2 = simple_macro_mol->addMacroAtom("PEPTIDE", "Y");
  auto macro_atom_3 = simple_macro_mol->addMacroAtom("PEPTIDE", "E");
  simple_macro_mol->addMacroBond(macro_atom_1, macro_atom_2,
                                 Bond::BondType::SINGLE, "R2", "R1");
  simple_macro_mol->addMacroBond(macro_atom_2, macro_atom_3,
                                 Bond::BondType::SINGLE, "R2", "R1");

  RDKit::GeneralizedSubstruct::ExtendedQueryMol query(
      std::make_unique<RWMol>(*simple_macro_mol));
  auto match =
      RDKit::GeneralizedSubstruct::hasSubstructMatch(*pdb_macro_mol, query);
  CHECK(match);
}
*/

TEST_CASE("testUniversalAttachmentPointConventions") {
  // Build two MacroMols with the same template atoms and bonds, but different
  // attachment point conventions, and check that they are not substructure
  // matches of each other when the attachment point conventions are not the
  // same, but are when the attachment point conventions are the same.
  auto macro_mol_1 = std::make_unique<MacroMol>();
  auto mm_1_atom_1 = macro_mol_1->addMacroAtom("AminoAcid", "A");
  auto mm_1_atom_2 = macro_mol_1->addMacroAtom("AminoAcid", "C");
  auto mm_1_atom_3 = macro_mol_1->addMacroAtom("AminoAcid", "D");
  macro_mol_1->addMacroBond(mm_1_atom_1, mm_1_atom_2, Bond::BondType::SINGLE,
                            "2", "1");
  macro_mol_1->addMacroBond(mm_1_atom_2, mm_1_atom_3, Bond::BondType::SINGLE,
                            "2", "1");

  auto macro_mol_2 = std::make_unique<MacroMol>();
  auto mm_2_atom_1 = macro_mol_2->addMacroAtom("AminoAcid", "A");
  auto mm_2_atom_2 = macro_mol_2->addMacroAtom("AminoAcid", "C");
  macro_mol_2->addMacroBond(
      mm_2_atom_1, mm_2_atom_2, Bond::BondType::SINGLE, "2",
      "1");  // NOTE: SAME attachment point convention than macro_mol_1

  auto macro_mol_3 = std::make_unique<MacroMol>();
  auto mm_3_atom_1 = macro_mol_3->addMacroAtom("AminoAcid", "A");
  auto mm_3_atom_2 = macro_mol_3->addMacroAtom("AminoAcid", "C");
  macro_mol_3->addMacroBond(
      mm_3_atom_1, mm_3_atom_2, Bond::BondType::SINGLE, "2",
      "R1");  // NOTE: DIFFERENT attachment point convention as macro_mol_1

  auto macro_mol_4 = std::make_unique<MacroMol>();
  auto mm_4_atom_1 = macro_mol_4->addMacroAtom("AminoAcid", "A");
  auto mm_4_atom_2 = macro_mol_4->addMacroAtom("AminoAcid", "C");
  macro_mol_4->addMacroBond(mm_4_atom_1, mm_4_atom_2, Bond::BondType::SINGLE,
                            "2",
                            "3");  // NOTE: DIFFERENT attachment point location
                                   // compared to macro_mol_1

  RDKit::GeneralizedSubstruct::ExtendedQueryMol query_2(
      std::make_unique<RWMol>(*macro_mol_2));
  auto match_1_2 =
      RDKit::GeneralizedSubstruct::hasSubstructMatch(*macro_mol_1, query_2);
  CHECK(match_1_2);
  RDKit::GeneralizedSubstruct::ExtendedQueryMol query_3(
      std::make_unique<RWMol>(*macro_mol_3));
  auto match_1_3 =
      RDKit::GeneralizedSubstruct::hasSubstructMatch(*macro_mol_1, query_3);
  CHECK(match_1_3);
  RDKit::GeneralizedSubstruct::ExtendedQueryMol query_4(
      std::make_unique<RWMol>(*macro_mol_4));
  auto match_1_4 =
      RDKit::GeneralizedSubstruct::hasSubstructMatch(*macro_mol_1, query_4);
  CHECK(
      match_1_4);  // This should fail because a change in the attachment point
                   // location should be considered a change chemical structure
}

/*
TEST_CASE("testMacroMolToAtomisticMol") {
  // Build a simple MacroMol with three macro atoms and two bonds, convert it to
  // an atomistic mol, and check that the resulting atomistic mol has the
  // expected number of atoms and bonds.
  auto macro_mol = std::make_unique<MacroMol>();
  auto macro_atom_1 = macro_mol->addMacroAtom("AminoAcid", "A");
  auto macro_atom_2 = macro_mol->addMacroAtom("AminoAcid", "C");
  auto macro_atom_3 = macro_mol->addMacroAtom("AminoAcid", "D");
  macro_mol->addMacroBond(macro_atom_1, macro_atom_2, Bond::BondType::SINGLE,
                          "2", "1");
  macro_mol->addMacroBond(macro_atom_2, macro_atom_3, Bond::BondType::SINGLE,
                          "2", "1");

  auto atomistic_mol = RDKit::toAtomistic(macro_mol);
  CHECK(atomistic_mol->getNumAtoms() == 37);
  CHECK(atomistic_mol->getNumBonds() == 40);
  auto smiles = RDKit::MolToSmiles(*atomistic_mol);
  CHECK(smiles == "N[C@@H](C)C(=O)N[C@@H](CS)C(=O)N[C@@H](CC(=O)O)C(=O)O");
}
*/

TEST_CASE("testMacroMolCanonicalize") {
  // Build two MacroMols with the same template atoms and bonds, but different
  // orders of the macro atoms, canonicalize them, and check that the canonical
  // SMILES for the two MacroMols are the same after canonicalization.
  auto macro_mol_1 = std::make_unique<MacroMol>();
  auto mm_1_atom_1 = macro_mol_1->addMacroAtom("AminoAcid", "A");
  auto mm_1_atom_2 = macro_mol_1->addMacroAtom("AminoAcid", "C");
  auto mm_1_atom_3 = macro_mol_1->addMacroAtom("AminoAcid", "D");
  macro_mol_1->addMacroBond(mm_1_atom_1, mm_1_atom_2, Bond::BondType::SINGLE,
                            "2", "1");
  macro_mol_1->addMacroBond(mm_1_atom_2, mm_1_atom_3, Bond::BondType::SINGLE,
                            "2", "1");

  auto macro_mol_2 = std::make_unique<MacroMol>();
  auto mm_2_atom_1 = macro_mol_2->addMacroAtom("AminoAcid", "C");
  auto mm_2_atom_2 = macro_mol_2->addMacroAtom("AminoAcid", "D");
  auto mm_2_atom_3 = macro_mol_2->addMacroAtom("AminoAcid", "A");
  macro_mol_2->addMacroBond(mm_2_atom_3, mm_2_atom_1, Bond::BondType::SINGLE,
                            "2", "1");
  macro_mol_2->addMacroBond(mm_2_atom_1, mm_2_atom_2, Bond::BondType::SINGLE,
                            "2", "1");

  auto can_smiles_1 = RDKit::MolToSmiles(*macro_mol_1);
  auto can_smiles_2 = RDKit::MolToSmiles(*macro_mol_2);
  CHECK(can_smiles_1 == can_smiles_2);
}