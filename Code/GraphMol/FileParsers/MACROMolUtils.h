//
//  Copyright (C) 2010-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MACROMOLRUTILS_H
#define RD_MACROMOLRUTILS_H

#include <string>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MACROMol.h>

#include <string_view>

namespace RDKit {
class RWMol;
class Conformer;

enum class MACROTemplateNames {
  AsEntered,      //<! use the name of the temlate as entered in the MACRO Mol
  UseFirstName,   //<!Use the first name in the template
                  // def (For AA, the 3 letter code
  UseSecondName,  //<!use the second name in the tempate def (
                  // For AA, the 1 letter code)
  All             //<! use all names in the template def
};

struct RDKIT_FILEPARSERS_EXPORT MolFromMACROMolParams {
  bool includeLeavingGroups =
      true; /**< when true, leaving groups on atoms that are not exo-bonded
               are retained.  When false, no leaving groups are retained */
  bool outputSgroups = true;
  MACROTemplateNames macroTemplateNames = MACROTemplateNames::All;
};

enum class MACROUseTemplateName {
  UseFirstName,   //<!Use the first name in the template
                  // def (For AA, the 3 letter code
  UseSecondName,  //<!use the second name in the tempate def (
                  // For AA, the 1 letter code)
};

struct RDKIT_FILEPARSERS_EXPORT MolToMACROParams {
  MACROUseTemplateName macroUseTemplateName =
      MACROUseTemplateName::UseFirstName;
};

RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::RWMol> MolFromMACROMol(
    const MACROMol *macroMol,
    const RDKit::v2::FileParsers::MolFileParserParams &molFileParserParams,
    const RDKit::MolFromMACROMolParams &molFromMACROMolParams);

RDKIT_FILEPARSERS_EXPORT void MACROMolToSCSRMolFile(
    RDKit::MACROMol &macroMol, const std::string &fName,
    const RDKit::MolWriterParams &params, int confId);

//! Build a MACROMol from an atomistic mol using the supplied template library.
/*!
    \warning The resulting MACROMol does NOT copy \c templates into its own local
    library; instead it retains a reference to \c templates as an external
    template library (see MACROMol::addTemplateLibrary).  Therefore \c templates
    MUST outlive the returned MACROMol.  The intended caller passes the
    process-lifetime singleton from getGlobalMonomerTemplateLib(); callers that
    pass a stack/temporary library are responsible for keeping it alive for as
    long as the MACROMol is used.
*/
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RDKit::MACROMol> MolToMACROMol(
    const ROMol &mol, RDKit::MACROMolTemplateLib &templates,
    MolToMACROParams molToMACROMolParams = MolToMACROParams());

//! \overload
/*!
    \warning Same lifetime contract as the overload above: \c templates is
    referenced (not copied) by \c res and must outlive it.
*/
RDKIT_FILEPARSERS_EXPORT void MolToMACROMol(MACROMol *res,
    const ROMol &mol, RDKit::MACROMolTemplateLib &templates,
    MolToMACROParams molToMACROMolParams = MolToMACROParams());

RDKIT_FILEPARSERS_EXPORT void ensureLocalTemplatesForScsr(MACROMol &macroMol);

}  // namespace RDKit

#endif
