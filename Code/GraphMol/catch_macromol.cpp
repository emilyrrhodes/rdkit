//
//  Copyright (C) 2026 RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/MACROMol.h>

using namespace RDKit;

TEST_CASE("MACROMol addMacroAtom with numbering preserves original and exposes current", "[macromol]") {
    MACROMol mol;
    auto idx = mol.addMacroAtom("AA", "ALA", "A", 1);
    CHECK(mol.getNumAtoms() == 1);

    // original numbering is set and frozen (no public setter exists)
    auto origChain = mol.getOriginalChainId(idx);
    auto origResNum = mol.getOriginalResidueNumber(idx);
    REQUIRE(origChain.has_value());
    REQUIRE(origResNum.has_value());
    CHECK(*origChain == "A");
    CHECK(*origResNum == 1);

    // current numbering starts equal to original
    CHECK(mol.getCurrentChainId(idx) == "A");
    CHECK(mol.getCurrentResidueNumber(idx) == 1);

    // current numbering can be updated (e.g. by assignChains())
    mol.setCurrentChainId(idx, "B");
    mol.setCurrentResidueNumber(idx, 2);
    CHECK(mol.getCurrentChainId(idx) == "B");
    CHECK(mol.getCurrentResidueNumber(idx) == 2);

    // original is unchanged
    CHECK(*mol.getOriginalChainId(idx) == "A");
    CHECK(*mol.getOriginalResidueNumber(idx) == 1);
}

TEST_CASE("MACROMol addMacroAtom without numbering has no original", "[macromol]") {
    MACROMol mol;
    auto idx = mol.addMacroAtom("AA", "ALA");
    CHECK_FALSE(mol.getOriginalChainId(idx).has_value());
    CHECK_FALSE(mol.getOriginalResidueNumber(idx).has_value());
}
