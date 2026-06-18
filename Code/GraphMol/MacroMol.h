//
//  Copyright (C) 2024 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MACROMOL_H
#define RD_MACROMOL_H

#include "MacroMolTemplate.h"

#include <optional>

namespace RDKit {

class RDKIT_GRAPHMOL_EXPORT MacroMol : public RWMol {
 private:
  friend class SCSRUtils;
  // elements (MacroMolTemplate items) of the library are owned by the library
  MacroMolTemplateLib d_templateLibrary;
  // pointer to the global library if loaded (Owned globally)
  const MacroMolTemplateLib *d_globalLib = nullptr;

 protected:
  // the following are used in SCSR parsing, where the templates are still being
  // contructed, and should not be used by other callers
  MacroMolTemplate *getMutableTemplate(unsigned int atomIdx);

 public:
  MacroMol() = default;
  explicit MacroMol(bool useGlobalLibrary);
  MacroMol(const MacroMol &other) : RWMol((RWMol)other) {
    d_templateLibrary.copyTemplateLib(*other.getTemplateLibrary());
    d_globalLib = other.d_globalLib;  // share the same global; do not clone
  }
  MacroMol(MacroMol &&other) noexcept = delete;
  MacroMol &operator=(MacroMol &&other) noexcept = delete;

  MacroMol &operator=(const MacroMol &) = delete;  // disable assignment

  MacroMol(std::unique_ptr<RWMol> &rwMol) : RWMol(std::move(*(rwMol))) {
    rwMol = nullptr;
  }

  ~MacroMol() {}

  const MacroMolTemplate *getTemplate(unsigned int atomIdx) const;

  const MacroMolTemplateLib *getTemplateLibrary() const {
    return &d_templateLibrary;
  }

  // the following adds a template to the internal libraty for this MacroMol
  void addTemplate(std::unique_ptr<MacroMolTemplate> &templateMol) {
    PRECONDITION(templateMol, "bad template molecule");
    d_templateLibrary.addTemplate(templateMol);
  }

  unsigned int size() const { return d_templateLibrary.size(); }

  unsigned int addMacroAtom(
      std::string className, std::string templateName,
      std::optional<unsigned int> residueNumber = std::nullopt,
      std::optional<std::string> chainId = std::nullopt);

  void addMacroBond(unsigned int fromAtomIdx, unsigned int toAtomIdx,
                    Bond::BondType bondType, std::string fromConnectionPoint,
                    std::string toConnectionPoint);
};
typedef boost::shared_ptr<MacroMol> MacroMol_SPTR;
}  // namespace RDKit

#endif
