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

#include <memory>
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
  // The template library (and the MacroMolTemplate items it owns) is shared
  // between copies; templates are immutable once added, so copies can safely
  // share the same library.  addTemplate() copies-on-write to preserve value
  // semantics when the library is shared.
  std::shared_ptr<MacroMolTemplateLib> d_templateLibrary =
      std::make_shared<MacroMolTemplateLib>();

 public:
  MacroMol() = default;
  MacroMol(const MacroMol &other)
      : RWMol((RWMol)other), d_templateLibrary(other.d_templateLibrary) {}
  MacroMol(MacroMol &&other) noexcept = delete;
  MacroMol &operator=(MacroMol &&other) noexcept = delete;

  MacroMol &operator=(const MacroMol &) = delete;  // disable assignment

  const MacroMolTemplate *getTemplate(unsigned int atomIdx) const;

  const MacroMolTemplateLib *getTemplateLibrary() const {
    return d_templateLibrary.get();
  }

  // the following adds a template to the internal library for this MacroMol
  void addTemplate(std::unique_ptr<MacroMolTemplate> &templateMol) {
    PRECONDITION(templateMol, "bad template molecule");
    // copy-on-write: if the library is shared with another MacroMol, take a
    // private copy before mutating it.
    if (d_templateLibrary.use_count() > 1) {
      auto privateLib = std::make_shared<MacroMolTemplateLib>();
      privateLib->copyTemplateLib(*d_templateLibrary);
      d_templateLibrary = std::move(privateLib);
    }
    d_templateLibrary->addTemplate(templateMol);
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
