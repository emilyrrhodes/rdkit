//
//  Copyright (C) 2026 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>

#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MACROMol.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MACROMolUtils.h>
#include <GraphMol/FileParsers/SCSRUtils.h>
#include <GraphMol/FileParsers/MonomerTemplateLib.h>

using namespace RDKit;

namespace {
std::string macromolsDataDir() {
  std::string rdbase = getenv("RDBASE");
  return rdbase + "/Code/GraphMol/FileParsers/test_data/macromols/";
}

void checkPdbToScsr(const std::string &pdbFileName, unsigned nAtoms,
                    unsigned nBonds, unsigned nTemplates,
                    const std::vector<std::string> &expectedTemplateNames) {
  std::string fName = macromolsDataDir() + pdbFileName;

  // read the PDB into an atomistic mol
  std::unique_ptr<RWMol> mol(PDBFileToMol(fName));
  REQUIRE(mol);

  // translate to a MACROMol using the global monomer template library
  auto macroMol = MolToMACROMol(*mol, getGlobalMonomerTemplateLib());
  REQUIRE(macroMol);

  // every residue must resolve to a templated macro atom; a leftover atomistic
  // atom (CLASS-less) here would mean a residue was mis-templated
  CHECK(macroMol->getNumAtoms() == nAtoms);
  CHECK(macroMol->getNumBonds() == nBonds);

  // templates aren't copied locally; the global lib is referenced externally
  CHECK(macroMol->getNumTemplates() == 0);
  CHECK(macroMol->getNumExternalTemplateLibs() == 1);

  // write the MACROMol to SCSR
  std::string scsr = MACROMolToSCSRMolBlock(*macroMol);
  CHECK(!scsr.empty());
  CHECK(scsr.find("BEGIN TEMPLATE") != std::string::npos);
  CHECK(scsr.find("V30") != std::string::npos);

  // the expected monomer template definitions must appear in the TEMPLATE block
  for (const auto &name : expectedTemplateNames) {
    INFO("expected template " << name << " in SCSR output");
    CHECK(scsr.find(name) != std::string::npos);
  }

  // the write used a throwaway copy, so the caller's MACROMol is unchanged
  CHECK(macroMol->getNumTemplates() == 0);
  CHECK(macroMol->getNumExternalTemplateLibs() == 1);

  // round-trip: templates from the SCSR TEMPLATE block land in the local lib
  RDKit::v2::FileParsers::MolFileParserParams pp;
  pp.sanitize = true;
  pp.removeHs = false;
  pp.strictParsing = true;

  auto rt = MACROMolFromSCSRBlock(scsr, pp);
  REQUIRE(rt);
  CHECK(rt->getNumAtoms() == nAtoms);
  CHECK(rt->getNumBonds() == nBonds);
  CHECK(rt->getNumTemplates() == nTemplates);
}
}  // namespace

TEST_CASE("pdbToMacroMolToScsr") {
  SECTION("1dng") {
    // single peptide chain of 15 residues:
    //   GLN ALA PRO ALA TYR GLU GLU ALA ALA GLU GLU LEU ALA LYS SER
    // -> 15 templated macro atoms, 14 backbone bonds, 8 distinct templates.
    checkPdbToScsr(
        "1dng.pdb", 15, 14, 8,
        {"ALA", "PRO", "TYR", "GLU", "LEU", "LYS", "SER", "GLN"});
  }
}
