/* -------------------------------------------------------------------------
 * Declares RDKit atomistic ROMol <-> MacroMol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#pragma once

#include <RDGeneral/export.h>

#include <memory>

namespace RDKit {
class ROMol;
class RWMol;
class MacroMol;
class MacroMolTemplateLib;
namespace v2 {
namespace FileParsers {
struct MolFileParserParams;
}
}  // namespace v2

enum class MacroMolTemplateNames {
  AsEntered,
  UseFirstName,
  UseSecondName,
  All,
};

struct RDKIT_FILEPARSERS_EXPORT MolToMacroMolParams {
  MacroMolTemplateNames macroTemplateNames = MacroMolTemplateNames::UseFirstName;
};

struct RDKIT_FILEPARSERS_EXPORT MolFromMacroMolParams {
  bool includeLeavingGroups = true;
  bool outputSgroups = true;
  MacroMolTemplateNames macroTemplateNames = MacroMolTemplateNames::AsEntered;
};

/**
 * Converts an atomistic ROMol into a MacroMol.
 *
 * @param atomistic_mol Atomistic molecule to convert to MacroMol
 * @param templates Templates used to identify macro atoms
 * @param molToMacroMolParams Conversion options
 * @return MacroMol
 */
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<MacroMol> MolToMacroMol(
    const ROMol &atomistic_mol, const MacroMolTemplateLib &templates,
    MolToMacroMolParams molToMacroMolParams = MolToMacroMolParams());

/**
 * Build an atomistic molecule from a MacroMol using templates from the
 * molecule's MacroMolTemplateLib (falling back to the global library).
 *
 * @param monomer_mol Monomeric molecule to convert to atomistic
 * @return Atomistic molecule
 */
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMacroMol(
    const MacroMol &monomer_mol);

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolFromMacroMol(
    const MacroMol *macroMol,
    const v2::FileParsers::MolFileParserParams &molFileParserParams,
    const MolFromMacroMolParams &molFromMacroMolParams =
        MolFromMacroMolParams());

}  // namespace RDKit
