//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD licenseFget
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>

#include <GraphMol/RWMol.h>
#include <GraphMol/MonomerMol/MonomerLibrary.h>
#include <GraphMol/MonomerMol/MonomerMol.h>
#include <GraphMol/MonomerMol/Conversions.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("MonomerLibrary") {
  SECTION("GlobalLibraryBasics") {
    // The global library is available by default
    //CHECK(MonomerLibrary::isUsingGlobalLibrary() == true);

    // Get a library with the globals loaded
    MonomerLibrary *globalLib = MonomerLibrary::getGlobalLibrary();

    // Check that built-in monomers are available
    CHECK(globalLib->hasMonomer("A", "PEPTIDE"));
    CHECK(globalLib->hasMonomer("G", "PEPTIDE"));
    CHECK(globalLib->hasMonomer("C", "PEPTIDE"));

    // Check that non-existent monomers return false
    CHECK(globalLib->hasMonomer("XYZ", "PEPTIDE") == false);
    CHECK(globalLib->hasMonomer("A", "RNA") == false);

    // Get monomer data
    auto alanineData = globalLib->getMonomerData("A", "PEPTIDE");
    REQUIRE(alanineData.has_value());
    CHECK(alanineData->find("C[C@H]") != std::string::npos);  // Alanine SMILES

    // Get monomer info from PDB code
    auto alaInfo = globalLib->getMonomerInfo("ALA", "PEPTIDE");
    REQUIRE(alaInfo.has_value());
    CHECK(std::get<0>(*alaInfo) == "A");  // symbol
    CHECK(std::get<2>(*alaInfo) == "PEPTIDE");  // class

    // Get PDB code from symbol
    auto pdbCode = globalLib->getPdbCode("A", "PEPTIDE");
    REQUIRE(pdbCode.has_value());
    CHECK(*pdbCode == "ALA");
  }

  SECTION("GlobalLibraryWithMonomerMol") {
    // When no custom library is set, MonomerMol uses the global library
    MonomerMol mol(true);  // use global library
    CHECK(mol.hasLocalTemplates() == false);

    // getGlobalLibrary() returns the global library
    MonomerLibrary &lib(*mol.getGlobalLibrary());
    CHECK(lib.hasMonomer("A", "PEPTIDE"));

    // Build a simple peptide using the global library
    mol.addMonomer("A", 1, "PEPTIDE", "PEPTIDE1");
    mol.addMonomer("G");
    mol.addMonomer("C");
    mol.addConnection(0, 1, BACKBONE_LINKAGE);
    mol.addConnection(1, 2, BACKBONE_LINKAGE);

    CHECK(mol.getNumAtoms() == 3);
    CHECK(mol.getNumBonds() == 2);

    // Convert to atomistic - uses global library
    auto atomistic = toAtomistic(mol);
    CHECK(atomistic->getNumAtoms() > 0);
  }

  SECTION("CustomLibraryViaConstructor") {
    // Create a custom library with built-ins plus a non-standard monomer
    auto customLib = std::make_shared<MonomerLibrary>(); 

    // Add a custom monomer (fictitious amino acid "X")
    customLib->addMonomerFromSmiles(
        "CC(C)(C)[C@H](N[H:1])C(=O)[OH:2]",  // tert-butyl glycine
        "X",
        "PEPTIDE",
        "TBG"
    );
    CHECK(customLib->hasMonomer("X", "PEPTIDE"));

    // Create MonomerMol with the custom library
    MonomerMol mol(customLib.get());
    mol.addGlobalLibrary();

    auto &internalLib = mol.getMonomerLibrary();
    internalLib.addMonomerFromSmiles(
        "CC(C)[C@H](N[H:1])C(=O)[OH:2]",  // i-propyl glycine
        "XX",
        "PEPTIDE",
        "IPG"
    );

    CHECK(mol.hasCustomLibrary() == true);
    CHECK(mol.usingGlobalLibrary() == true);
    CHECK(mol.hasLocalTemplates() == true);

    // The custom monomer is accessible
    auto lib = mol.getCustomLibrary();
    CHECK(lib->hasMonomer("X", "PEPTIDE"));
    CHECK(internalLib.hasMonomer("XX", "PEPTIDE"));

    // Build a peptide with the custom monomer
    mol.addMonomer("A", 1, "PEPTIDE", "PEPTIDE1");
    mol.addMonomer("X");  // custom monomer
    mol.addMonomer("XX");  // customer monomer from internal library
    mol.addMonomer("G");
    mol.addConnection(0, 1, BACKBONE_LINKAGE);
    mol.addConnection(1, 2, BACKBONE_LINKAGE);
    mol.addConnection(2, 3, BACKBONE_LINKAGE);

    CHECK(mol.getNumAtoms() == 4);

    // Convert to atomistic - uses custom library
    auto atomistic = toAtomistic(mol);
    CHECK(atomistic->getNumAtoms() > 0);

    // Verify the custom monomer was expanded correctly
    std::string smiles = MolToSmiles(*atomistic);
    CHECK(smiles.find("C(C)(C)") != std::string::npos);  // tert-butyl group
  }

  SECTION("CustomLibraryViaSetMonomerLibrary") {
    // Create a MonomerMol first (uses global by default)
    MonomerMol mol;
    CHECK(mol.hasLocalTemplates() == false);

    // Create and set a custom library
    auto customLib = std::make_shared<MonomerLibrary>();  // load built-ins
    customLib->addMonomerFromSmiles(
        "NCCC[C@H](N[H:1])C(=O)[OH:2]",  // ornithine-like
        "Z",
        "PEPTIDE",
        "ORN"
    );


    mol.setCustomMonomerLibrary(customLib.get());
    mol.addGlobalLibrary();

    CHECK(mol.hasLocalTemplates() == false);
    CHECK(mol.hasCustomLibrary() == true);

    // Now the custom monomer is available
    CHECK(mol.getCustomLibrary()->hasMonomer("Z", "PEPTIDE"));

    // Build peptide with custom monomer
    mol.addMonomer("Z", 1, "PEPTIDE", "PEPTIDE1");
    mol.addMonomer("A");
    mol.addConnection(0, 1, BACKBONE_LINKAGE);

    auto atomistic = toAtomistic(mol);
    CHECK(atomistic->getNumAtoms() > 0);
  }

  SECTION("ClearCustomLibrary") {
    auto customLib = std::make_shared<MonomerLibrary>();
    customLib->addMonomerFromSmiles("CC", "Y", "PEPTIDE");

    MonomerMol mol(customLib.get());
    CHECK(mol.hasCustomLibrary() == true);
    CHECK(mol.getCustomLibrary()->getMACROMolTemplateLib().getNumTemplates() > 0);

    // Clear 
    mol.getCustomLibrary()->getMACROMolTemplateLib().clearTemplateLib();
    CHECK(mol.hasCustomLibrary() == true);
    CHECK(mol.getCustomLibrary()->getMACROMolTemplateLib().getNumTemplates() == 0);

    // Now uses global library again
    mol.addGlobalLibrary();
    CHECK(mol.usingGlobalLibrary());
  }

  SECTION("SharedLibraryBetweenMolecules") {
    // Multiple MonomerMols can share the same custom library
    auto externalLib = std::make_shared<MonomerLibrary>();
    externalLib->addMonomerFromSmiles("CC", "Q1", "PEPTIDE", "QQ1");

    MonomerMol mol1(externalLib.get());
    MonomerMol mol2(externalLib.get());

    CHECK(mol1.hasLocalTemplates() == false);
    CHECK(mol2.hasLocalTemplates() == false);

    // Adding to the shared library affects both
    externalLib->addMonomerFromSmiles("CCC", "Q2", "PEPTIDE", "QQ2");
    CHECK(mol1.getCustomLibrary()->hasMonomer("Q2", "PEPTIDE") ==  true);
    CHECK(mol2.getCustomLibrary()->hasMonomer("Q2", "PEPTIDE") ==  true);


    MonomerMol mol3(externalLib.get());
    CHECK(mol3.hasLocalTemplates() == false);
    CHECK(mol3.getCustomLibrary()->hasMonomer("Q2", "PEPTIDE") ==  true);

  }

  SECTION("CopyPreservesLibrary") {
    auto customLib = std::make_shared<MonomerLibrary>();
    customLib->addMonomerFromSmiles("CC", "W1", "PEPTIDE");

    MonomerMol mol1(customLib.get());
    mol1.addMonomer("A", 1, "PEPTIDE", "PEPTIDE1");
    
    // Copy constructor preserves the library
    MonomerMol mol2(mol1);

    customLib->addMonomerFromSmiles("CCCC", "W2", "PEPTIDE");

    CHECK(mol2.hasLocalTemplates() == false);
    CHECK(mol2.getCustomLibrary()->getMACROMolTemplateLib().getNumTemplates() == mol1.getCustomLibrary()->getMACROMolTemplateLib().getNumTemplates());
    CHECK(mol2.getCustomLibrary()->hasMonomer("W1", "PEPTIDE")== true);
    CHECK(mol2.getCustomLibrary()->hasMonomer("W2", "PEPTIDE")== true);
    
    mol2.getMonomerLibrary().copyMonomerLib(*(customLib.get()));
    CHECK(mol2.getMonomerLibrary().hasMonomer("W2", "PEPTIDE")== true);

    // Copy assignment also preserves
    // MonomerMol mol3;
    // mol3 = mol1;
    // CHECK(mol3.hasLocalTemplates() == true);
  }

  SECTION("MovePreservesLibrary") {
    auto customLib = std::make_shared<MonomerLibrary>();
    customLib->addMonomerFromSmiles("CC", "M1", "PEPTIDE");

    MonomerMol mol1(customLib.get());
    mol1.addMonomer("A", 1, "PEPTIDE", "PEPTIDE1");

    // Move constructor transfers the library
    MonomerMol mol2(std::move(mol1));
    CHECK(mol2.hasLocalTemplates() == false);
    CHECK(mol2.getCustomLibrary()->hasMonomer("M1", "PEPTIDE"));

    MonomerMol mol3 = std::move(mol2);
    CHECK(mol3.hasLocalTemplates() == false);
    CHECK(mol3.getCustomLibrary()->hasMonomer("M1", "PEPTIDE"));

  }

  SECTION("GetMonomer") {
    auto lib = MonomerLibrary::getGlobalLibrary();

    // Get a parsed molecule for alanine
    auto alaMol = lib->getMonomer("A", "PEPTIDE");
    auto test = alaMol != nullptr;
    CHECK(test);
    //CHECK(alaMol->getNumAtoms() > 0);

    // Second call should return the same cached molecule
    // auto alaMol2 = lib->getMonomer("A", "PEPTIDE");
    // CHECK(alaMol == alaMol2);

    // // Non-existent monomer returns nullptr
    // auto notFound = lib->getMonomer("NOTEXIST", "PEPTIDE");
    // CHECK(notFound == nullptr);
  }

  SECTION("AddMonomerWithMol") {
    auto customLib = std::make_shared<MonomerLibrary>();

    // Pre-parse a molecule
    auto mol = std::unique_ptr<RWMol>(SmilesToMol("CC(N)C(=O)O"));
    auto origAtomCount = mol->getNumAtoms();

    // Add with pre-parsed mol (no original data needed)
    customLib->addMonomer(mol, "CC(N)C(=O)O",  "TST" ,"PEPTIDE");

    // getMol returns the same pre-parsed molecule
    auto retrieved = customLib->getMonomer("TST", "PEPTIDE");
    CHECK(retrieved->getNumAtoms() == origAtomCount);
  }

  SECTION("AddMonomerFromSDF") {
    auto customLib = std::make_shared<MonomerLibrary>();

    // Simple alanine-like structure in SDF format
    std::string sdfData = R"(
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  2  0  0  0  0  0  0
    2.2500    1.2990    0.0000 N   0  0  0  0  0  1  0  0  0  0  0  0
    2.2500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  2  0
M  END
)";

    customLib->addMonomerFromSDF(sdfData, "SDF1", "PEPTIDE", "SD1");

    CHECK(customLib->hasMonomer("SDF1", "PEPTIDE"));

    // getMol should parse the SDF and return a molecule
    auto mol = customLib->getMonomer("SDF1", "PEPTIDE");
    REQUIRE(mol != nullptr);
    CHECK(mol->getNumAtoms() == 4);
    CHECK(mol->getNumBonds() == 3);
  }

  // SECTION("GlobalLibraryConfiguration") {
  //   // Save original state
  //   bool originalState = MonomerLibrary::isUsingGlobalLibrary();

  //   // Can toggle global library mode
  //   MonomerLibrary::useGlobalLibrary(false);
  //   CHECK(MonomerLibrary::isUsingGlobalLibrary() == false);

  //   MonomerLibrary::useGlobalLibrary(true);
  //   CHECK(MonomerLibrary::isUsingGlobalLibrary() == true);

  //   // Restore original state
  //   MonomerLibrary::useGlobalLibrary(originalState);
  // }

  SECTION("EmptyLibrary") {
    // Create an empty library (default behavior)
    auto emptyLib = std::make_shared<MonomerLibrary>();

    // Should not have any built-in monomers
    CHECK(emptyLib->hasMonomer("A", "PEPTIDE") == false);
    CHECK(emptyLib->hasMonomer("G", "PEPTIDE") == false);
    CHECK(emptyLib->hasMonomer("C", "PEPTIDE") == false);

    // Can still add custom monomers
    emptyLib->addMonomerFromSmiles(
        "CC[C@H](N[H:1])C(=O)[OH:2]",
        "CUSTOM",
        "PEPTIDE",
        "CUS"
    );
    CHECK(emptyLib->hasMonomer("CUSTOM", "PEPTIDE"));

    // Use with MonomerMol
    MonomerMol mol(emptyLib.get());
    mol.addMonomer("CUSTOM", 1, "PEPTIDE", "PEPTIDE1");
    CHECK(mol.getNumAtoms() == 1);

    // Convert to atomistic works with custom monomer
    auto atomistic = toAtomistic(mol);
    CHECK(atomistic->getNumAtoms() > 0);
  }

  SECTION("LibraryWithBuiltins") {
    // Explicitly create library with built-ins
    auto libWithBuiltins(MonomerLibrary::getGlobalLibrary());

    // Should have built-in monomers
    CHECK(libWithBuiltins->hasMonomer("A", "PEPTIDE"));
    CHECK(libWithBuiltins->hasMonomer("G", "PEPTIDE"));
  }
}
