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

//! Builds a MacroMol by collapsing template matches in \c mol into macro atoms.
/*!
  Each entry in \c templates is substructure-matched against \c mol; a matched
  region becomes a single macro atom and bonds crossing the region boundary
  become macro bonds via the templates' attachment points. Atoms not covered by
  any template are copied through as regular atoms, and regular bonds between
  regular atoms retain their full metadata (direction, stereo, query state,
  properties).

  \b Note: assignment is greedy — templates are tried largest-main-group-first
  and a match is rejected if it overlaps atoms already claimed. This does not
  guarantee a maximal or globally optimal cover when monomers can overlap. The
  result is deterministic (matches are claimed in a stable order).

  Each atom is assumed to carry at most one attachment point; a template atom
  with more than one attachment point is unsupported and throws. Attachment
  point ids must be numeric; a non-numeric id throws. Template matches that
  would require more than one external bond per attachment atom are rejected
  and fall back to plain atoms.
*/
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<MacroMol> MolToMacroMol(
    const ROMol &mol, const MacroMolTemplateLibrary &templates);

}  // namespace RDKit

#endif
