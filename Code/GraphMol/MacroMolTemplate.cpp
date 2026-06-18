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

#include "MacroMolTemplate.h"
#include "Atom.h"

namespace RDKit {

void RDKit::MacroMolTemplate::findMainSgroupForTemplate(
    std::string className, std::string templateName) const {
  d_mainSgroupIdx = UINT_MAX;

  auto sgroups = RDKit::getSubstanceGroups(*this);

  for (unsigned int sgroupIdx = 0; sgroupIdx < sgroups.size(); ++sgroupIdx) {
    auto &sgroupToTest = sgroups[sgroupIdx];
    std::string sup;
    std::string sgroupAtomClass;
    if (sgroupToTest.getPropIfPresent("TYPE", sup) && sup == "SUP" &&
        sgroupToTest.getPropIfPresent("CLASS", sgroupAtomClass) &&
        sgroupAtomClass == className) {
      d_mainSgroupIdx = sgroupIdx;
      break;
    }
  }

  if (d_mainSgroupIdx == UINT_MAX) {
    std::ostringstream errout;
    errout << "No main sgroup found for template: " << className << "/"
           << templateName;
    throw RDKit::FileParseException(errout.str());
  }
}

void RDKit::MacroMolTemplate::init(
    std::string className, std::vector<std::string> templateNames,
    std::vector<std::pair<std::string, std::string>> templateAttrs) {
  PRECONDITION(className.empty() == false, "no className for template");
  PRECONDITION(templateNames.size() > 0, "no template names for template");

  this->setProp(RDKit::common_properties::molAtomClass, className);
  this->setProp(RDKit::common_properties::templateNames, templateNames);

  for (auto templateAttr : templateAttrs) {
    this->setProp(templateAttr.first, templateAttr.second);
  }
  d_mainSgroupIdx = UINT_MAX;
  this->updatePropertyCache(false);
}

MacroMolTemplate::MacroMolTemplate(
    std::unique_ptr<RWMol> &mol, std::string className,
    std::vector<std::string> templateNames,
    std::vector<std::pair<std::string, std::string>> templateAttrs)
    : RWMol(std::move(*mol)) {
  init(className, templateNames, templateAttrs);
}

MacroMolTemplate::MacroMolTemplate(const MacroMolTemplate &other)
    : RWMol(other) {
  d_mainSgroupIdx = UINT_MAX;
}

MacroMolTemplate::MacroMolTemplate(
    std::unique_ptr<RWMol> &mol, std::string className,
    std::string templateName,
    std::vector<std::pair<std::string, std::string>> templateAttrs)
    : RWMol(std::move(*mol)) {
  PRECONDITION(!templateName.empty(), "no name for template");

  std::vector<std::string> templateNames(1, templateName);
  MacroMolTemplate::init(className, templateNames, templateAttrs);
}

void MacroMolTemplate::initMainSgroupIdx() const {
  std::call_once(d_mainSgroupIdxOnceFlag, [this]() {
    std::string className = "";
    std::vector<std::string> templateNames;
    if (!getPropIfPresent(RDKit::common_properties::molAtomClass, className) ||
        !getPropIfPresent(RDKit::common_properties::templateNames,
                          templateNames)) {
      std::ostringstream errout;
      errout << "Template molecule is missing required properties: "
             << RDKit::common_properties::molAtomClass << " and/or "
             << RDKit::common_properties::templateNames;
      throw RDKit::FileParseException(errout.str());
    }

    findMainSgroupForTemplate(className, templateNames[0]);
  });
}

RDKit::SubstanceGroup *MacroMolTemplate::getMainSgroup() {
  initMainSgroupIdx();
  return &RDKit::getSubstanceGroups(*this)[d_mainSgroupIdx];
}

const RDKit::SubstanceGroup *MacroMolTemplate::getMainSgroup() const {
  initMainSgroupIdx();
  return &RDKit::getSubstanceGroups(*this)[d_mainSgroupIdx];
}

void MacroMolTemplateLib::addTemplate(
    std::unique_ptr<MacroMolTemplate> &templateMolToAdd) {
  PRECONDITION(templateMolToAdd, "bad template molecule");

  d_templates.push_back(std::move(templateMolToAdd));

  auto templateMol = d_templates.back().get();
  std::string templateClass;
  std::vector<std::string> templateNames;
  templateMol->getPropIfPresent<std::string>(common_properties::molAtomClass,
                                             templateClass);
  templateMol->getPropIfPresent<std::vector<std::string>>(
      common_properties::templateNames, templateNames);

  for (auto templateName : templateNames) {
    MacroMolTemplateKey key(std::pair(templateClass, templateName));

    d_keyToIndex[key] = d_templates.size() - 1;
  }
}

void MacroMolTemplateLib::copyTemplateLib(
    const MacroMolTemplateLib &libToCopy) {
  // clear any existing templates in this library
  d_templates.clear();
  d_keyToIndex.clear();

  for (auto &templateToCopy : libToCopy) {
    // make a copy of the template
    auto templateCopy = std::unique_ptr<MacroMolTemplate>(
        new MacroMolTemplate(*(templateToCopy.get())));
    this->addTemplate(templateCopy);
  }
}
}  // namespace RDKit
