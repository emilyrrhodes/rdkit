/* -------------------------------------------------------------------------
 * Declares RDKit atomistic ROMol <-> MACROMol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#pragma once

#include <RDGeneral/export.h>

#include <memory>

namespace RDKit {
class ROMol;
class RWMol;
class MACROMol;

/**
 * Converts an atomistic ROMol into a MACROMol.
 *
 * @param atomistic_mol Atomistic molecule to convert to MACROMol
 * @return MACROMol
 */
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::MACROMol> toMonomeric(
    const ROMol &atomistic_mol);

/**
 * Build an atomistic molecule from a MACROMol using monomers
 * from the molecule's MonomerLibrary.
 *
 * @param macro_mol Monomeric molecule to convert to atomistic
 * @return Atomistic molecule
 */
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> toAtomistic(
    const MACROMol &macro_mol);

RDKIT_FILEPARSERS_EXPORT bool hasPdbResidueInfo(const ROMol &mol);

}  // namespace RDKit
