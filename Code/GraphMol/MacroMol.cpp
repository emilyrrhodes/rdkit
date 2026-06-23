//
//  Copyright (C) 2024 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include "FileParsers/FileParsers.h"
#include "FileParsers/FileParserUtils.h"

#include "MacroMol.h"
#include "Atom.h"
#include "MonomerInfo.h"

#include <optional>

namespace RDKit {

std::string monomerClassToString(MonomerClass monomerClass) {
  switch (monomerClass) {
    case MonomerClass::AA:
      return "AA";
    case MonomerClass::NA:
      return "NA";
    case MonomerClass::CHEM:
      return "CHEM";
    case MonomerClass::OTHER:
      return "OTHER";
  }
  return "OTHER";  // unreachable, but silences -Wreturn-type
}

unsigned int RDKit::MacroMol::addMacroAtom(
    MonomerClass monomerClass, std::string templateName,
    std::optional<unsigned int> residueNumber,
    std::optional<std::string> chainId) {
  auto className = monomerClassToString(monomerClass);
  auto atom = new Atom(0);

  atom->setProp(common_properties::dummyLabel, templateName);
  atom->setProp(common_properties::molAtomClass, className);

  auto *monomer_info = new ::RDKit::AtomMonomerInfo();
  monomer_info->setResidueName(templateName);
  monomer_info->setMonomerClass(className);
  if (residueNumber) {
    monomer_info->setResidueNumber(*residueNumber);
  }
  if (chainId) {
    monomer_info->setChainId(*chainId);
  }
  atom->setMonomerInfo(monomer_info);

  return this->addAtom(atom, false, true);
}

void RDKit::MacroMol::addMacroBond(unsigned int fromAtomIdx,
                                   unsigned int toAtomIdx,
                                   int fromConnectionPoint,
                                   int toConnectionPoint,
                                   std::optional<Bond::BondType> bondType) {
  const auto resolvedBondType = bondType.value_or(Bond::BondType::SINGLE);
  auto bondIdx = this->addBond(fromAtomIdx, toAtomIdx, resolvedBondType) - 1;
  auto bond = this->getBondWithIdx(bondIdx);
  this->setBondBookmark(bond, bondIdx);
  bond->setProp(common_properties::_MolFileBondAttachPt2, toConnectionPoint);
  bond->setProp(common_properties::_MolFileBondAttachPt1, fromConnectionPoint);
}

const MacroMolTemplate *RDKit::MacroMol::getTemplate(
    unsigned int atomIdx) const {
  const auto atom = this->getAtomWithIdx(atomIdx);

  std::string atomClass;
  std::string dummyLabel = "";
  if (!atom->getPropIfPresent(common_properties::dummyLabel, dummyLabel) ||
      dummyLabel == "" ||
      !atom->getPropIfPresent(common_properties::molAtomClass, atomClass) ||
      atomClass == "") {
    return nullptr;  // ordinary atom, not a macro atom
  }

  auto templatePtr = d_templateLibrary->find(atomClass, dummyLabel);
  if (templatePtr == nullptr) {
    // Fall back to the global builtin library.  This is intentional: it is how
    // hand-built MacroMols (which carry no templates of their own) resolve the
    // standard builtin monomers.  Templates are keyed by (class, name), so a
    // local template always takes precedence over a builtin of the same key.
    templatePtr =
        MacroMolTemplateLib::getGlobalLibrary().find(atomClass, dummyLabel);
  }
  if (templatePtr == nullptr) {
    std::ostringstream errout;
    errout << "Template not found for atom " << dummyLabel;
    throw RDKit::FileParseException(errout.str());
  }
  return templatePtr;
}
}  // namespace RDKit
