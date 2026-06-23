//
//  Copyright (C) 2026 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MONOMERTEMPLATES_H
#define RD_MONOMERTEMPLATES_H

#include <RDGeneral/export.h>

#include <memory>
#include <string>

namespace RDKit {

class MacroMolTemplate;
class MacroMolTemplateLib;

//! Build a single amino-acid MacroMolTemplate from a map-numbered SMILES.
/*!
    Each map-numbered atom in \c smiles is treated as a single-atom leaving
    group; the map number is its connection label and its single non-mapped
    neighbor is the attachment atom.  The returned template carries a main
    "SUP"/CLASS="AA" sgroup over the core atoms (with one attachment
    point per leaving atom) and one "SUP"/CLASS="LGRP" sgroup per leaving atom.
*/
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<MacroMolTemplate> buildAminoAcidTemplate(
    const std::string &symbol, const std::string &smiles);

//! Load the builtin amino-acid templates into \c lib.
RDKIT_FILEPARSERS_EXPORT void loadBuiltinMacroMolTemplates(
    MacroMolTemplateLib &lib);

}  // namespace RDKit

#endif
