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
#include <string>

namespace RDKit {

enum class MonomerClass {
  AA,
  NA,
  CHEM,
  OTHER
};

RDKIT_GRAPHMOL_EXPORT std::string monomerClassToString(
    MonomerClass monomerClass);

class RDKIT_GRAPHMOL_EXPORT MacroMol : public RWMol {
 private:
  // elements (MacroMolTemplate items) of the library are owned by the library
  MacroMolTemplateLib d_templateLibrary;

 public:
  MacroMol() = default;
  MacroMol(const MacroMol &other) : RWMol((RWMol)other) {
    d_templateLibrary.copyTemplateLib(*other.getTemplateLibrary());
  }
  MacroMol(MacroMol &&other) noexcept = delete;
  MacroMol &operator=(MacroMol &&other) noexcept = delete;

  MacroMol &operator=(const MacroMol &) = delete;  // disable assignment

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

  unsigned int addMacroAtom(
      MonomerClass monomerClass, std::string templateName,
      std::optional<unsigned int> residueNumber = std::nullopt,
      std::optional<std::string> chainId = std::nullopt);

  void addMacroBond(unsigned int fromAtomIdx, unsigned int toAtomIdx,
                    Bond::BondType bondType, int fromConnectionPoint,
                    int toConnectionPoint);
};
}  // namespace RDKit

#endif
