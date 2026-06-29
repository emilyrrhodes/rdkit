//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MOLTOMACROMOL_H
#define RD_MOLTOMACROMOL_H

#include <GraphMol/MacroMol.h>
#include <GraphMol/MacroMolTemplate.h>
#include <RDGeneral/export.h>

#include <memory>

namespace RDKit {

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<MacroMol> MolToMacroMol(
    const ROMol &mol, const MacroMolTemplateLibrary &templates);

}  // namespace RDKit

#endif
