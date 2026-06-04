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
#include <GraphMol/MonomerInfo.h>

using namespace RDKit;

TEST_CASE("MACROMol addMacroAtom attaches monomer info", "[macromol]") {
    MACROMol mol;
    auto idx = mol.addMacroAtom("AA", "ALA", "A", 1);

    CHECK(mol.getNumAtoms() == 1);

    auto* atom = mol.getAtomWithIdx(idx);
    REQUIRE(atom != nullptr);

    auto* stored = static_cast<const AtomMonomerInfo*>(atom->getMonomerInfo());
    REQUIRE(stored != nullptr);
    CHECK(stored->getResidueName() == "ALA");
    CHECK(stored->getResidueNumber() == 1);
    CHECK(stored->getMonomerClass() == "AA");
    CHECK(stored->getChainId() == "A");
}
