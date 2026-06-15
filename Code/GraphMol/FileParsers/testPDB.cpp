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
#include <GraphMol/FileParsers/MACROMolUtils.h>
#include <GraphMol/FileParsers/SCSRUtils.h>
#include <catch2/catch_all.hpp>
#include <GraphMol/FileParsers/MonomerLibrary.h>

using namespace RDKit;

TEST_CASE("testReadPDB") {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/1DNG.pdb";
  std::unique_ptr<RWMol> mol(PDBFileToMol(fname));
  int expectedNumAtoms = 113;
  int expectedNumBonds = 114;
  int numAtoms = mol->getNumAtoms();
  int numBonds = mol->getNumBonds();

  CHECK(numAtoms == expectedNumAtoms);
  CHECK(numBonds == expectedNumBonds);
}

TEST_CASE("testConvertToMacroMol") {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/1DNG.pdb";
  std::unique_ptr<RWMol> mol(PDBFileToMol(fname));

  MonomerLibrary *globalLib = MonomerLibrary::getGlobalLibrary();
  MolToMACROParams molToMACROParams;

  std::unique_ptr<RDKit::MACROMol> macroMol = RDKit::MolToMACROMol(
      *(mol.get()), globalLib->getMACROMolTemplateLib(), molToMACROParams);

  CHECK(macroMol->getNumAtoms() == 15);
  std::string sequence;
  for (const auto &atom : macroMol->atoms()) {
    std::string templateName =
        atom->getProp<std::string>(common_properties::dummyLabel);
    sequence += templateName;
  }
  CHECK(sequence == "QAPAYEEAAEELAKS");
  CHECK(macroMol->getNumBonds() == 14);
}

TEST_CASE("testMacroMolToScsr") {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/FileParsers/test_data/1DNG.pdb";
  std::unique_ptr<RWMol> mol(PDBFileToMol(fname));

  MonomerLibrary *globalLib = MonomerLibrary::getGlobalLibrary();
  MolToMACROParams molToMACROParams;

  std::unique_ptr<RDKit::MACROMol> macroMol = RDKit::MolToMACROMol(
      *(mol.get()), globalLib->getMACROMolTemplateLib(), molToMACROParams);

  std::string fOutName =
      rdbase + "/build/Testing/Temporary/macro_mol_to_scsr.scsrout.mol";

  MACROMolToSCSRMolFile(*(macroMol.get()), fOutName);

  RDKit::v2::FileParsers::MolFileParserParams pp;
  pp.sanitize = true;
  pp.removeHs = false;
  pp.strictParsing = true;

  RDKit::SCSRBaseHbondOptions scsrBaseHbondOptions;

  auto molReadBackIn =
      RDKit::MACROMolFromSCSRFile(fOutName, pp, scsrBaseHbondOptions);

  CHECK(molReadBackIn->getNumAtoms() == 15);
  std::string sequence;
  for (const auto &atom : molReadBackIn->atoms()) {
    std::string templateName =
        atom->getProp<std::string>(common_properties::dummyLabel);
    sequence += templateName;
  }
  CHECK(sequence == "QAPAYEEAAEELAKS");
  CHECK(molReadBackIn->getNumBonds() == 14);
}