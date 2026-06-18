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

MacroMol::MacroMol(bool useGlobalLibrary) {
  if (useGlobalLibrary) {
    d_globalLib = &MacroMolTemplateLib::getGlobalLibrary();
  }
}

unsigned int RDKit::MacroMol::addMacroAtom(
    std::string className, std::string templateName,
    std::optional<unsigned int> residueNumber,
    std::optional<std::string> chainId) {
  auto atom = new Atom(0);
  atom->setAtomicNum(0);

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
                                   Bond::BondType bondType,
                                   std::string fromConnectionPoint,
                                   std::string toConnectionPoint) {
  auto bondIdx = this->addBond(fromAtomIdx, toAtomIdx, bondType) - 1;
  auto bond = this->getBondWithIdx(bondIdx);
  this->setBondBookmark(bond, bondIdx);

  if (!toConnectionPoint.empty()) {
    bond->setProp(common_properties::_MolFileBondAttachPt2, toConnectionPoint);
  }

  if (!fromConnectionPoint.empty()) {
    bond->setProp(common_properties::_MolFileBondAttachPt1,
                  fromConnectionPoint);
  }
}

MacroMolTemplate *MacroMol::getMutableTemplate(unsigned int atomIdx) {
  const auto atom = this->getAtomWithIdx(atomIdx);

  std::string atomClass;
  std::string dummyLabel = "";
  if (!atom->getPropIfPresent(common_properties::dummyLabel, dummyLabel) ||
      dummyLabel == "" ||
      !atom->getPropIfPresent(common_properties::molAtomClass, atomClass) ||
      atomClass == "") {
    return nullptr;  // ordinary atom, not a macro atom
  }

  auto templatePtr = d_templateLibrary.findMutable(atomClass, dummyLabel);
  if (templatePtr == nullptr) {
    std::ostringstream errout;
    errout << "Template not found for atom " << dummyLabel;
    throw RDKit::FileParseException(errout.str());
  }
  return templatePtr;
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

  auto templatePtr = d_templateLibrary.find(atomClass, dummyLabel);
  if (templatePtr == nullptr) {
    const MacroMolTemplateLib *globalLib =
        d_globalLib ? d_globalLib : &MacroMolTemplateLib::getGlobalLibrary();
    templatePtr = globalLib->find(atomClass, dummyLabel);
  }
  if (templatePtr == nullptr) {
    std::ostringstream errout;
    errout << "Template not found for atom " << dummyLabel;
    throw RDKit::FileParseException(errout.str());
  }
  return templatePtr;
}
}  // namespace RDKit
