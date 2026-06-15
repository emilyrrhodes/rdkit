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

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include "FileParsers/FileParsers.h"
#include "FileParsers/FileParserUtils.h"
#include <mutex>
#include <unordered_map>


namespace RDKit {

// acts as a view to a polymer chain
struct Chain {
  std::vector<unsigned int> atoms;
  std::vector<unsigned int> bonds;
  std::string annotation;
};

const std::string LINKAGE{"attachmentPoints"};
const std::string EXTRA_LINKAGE{"extraAttachmentPoints"};
const std::string ATOM_LABEL{"atomLabel"};

// Some default linkage options
const std::string BRANCH_LINKAGE{"R3-R1"};
const std::string BACKBONE_LINKAGE{"R2-R1"};
const std::string CROSS_LINKAGE{"R3-R3"};
const std::string HYDROGEN_LINKAGE{"pair-pair"};

// Monomer properties stored on atoms
const std::string BRANCH_MONOMER{"isBranchMonomer"};
const std::string SMILES_MONOMER{"isSmilesMonomer"};

// Substance group property to indicate an annotation on a chain
const std::string ANNOTATION{"annotation"};

class RDKIT_GRAPHMOL_EXPORT MACROMolTemplate : public RDKit::RWMol {
 public:
  MACROMolTemplate(std::unique_ptr<RWMol> &mol, std::string className,
                   std::vector<std::string> templateNames,
                   std::vector<std::pair<std::string,std::string>> templateAttrs);

  MACROMolTemplate(std::unique_ptr<RWMol> &mol, std::string className,
                   std::string templateName,
                   std::vector<std::pair<std::string,std::string>> templateAttrs);

  MACROMolTemplate() = delete;
  MACROMolTemplate(const MACROMolTemplate &other);
  MACROMolTemplate(MACROMolTemplate &&other) noexcept = delete;
  MACROMolTemplate &operator=(MACROMolTemplate &&other) noexcept = delete;
  MACROMolTemplate &operator=(const MACROMolTemplate &) = delete; 
  ~MACROMolTemplate() {}

  RDKit::SubstanceGroup *getMainSgroup();
  const RDKit::SubstanceGroup *getMainSgroup() const;

 private:
  void init(std::string className,
            std::vector<std::string> templateNames,
            std::vector<std::pair<std::string,std::string>> templateAttrs);

  void findMainSgroupForTemplate(std::string className,
                                 std::string templateName) const;
  void initMainSgroupIdx() const;
  mutable unsigned int d_mainSgroupIdx;
  mutable std::once_flag d_mainSgroupIdxOnceFlag;
};

class RDKIT_GRAPHMOL_EXPORT MACROMolTemplateLib {
 private:
  // All templates in the library are owned by the library. Consumers get
  // const access; the library mutates them only during construction.
  std::vector<std::unique_ptr<MACROMolTemplate>> d_templates;

  //! Key is (monomer_class, symbol)
  using MACROMolTemplateKey = std::pair<std::string, std::string>;

  //! Hash function for MACROMolTemplateKey
  struct MACROMolTemplateKeyHash {
    std::size_t operator()(const MACROMolTemplateKey &key) const {
      std::size_t h1 = std::hash<std::string>{}(key.first);
      std::size_t h2 = std::hash<std::string>{}(key.second);
      return h1 ^ (h2 << 1);
    }
  };

  std::unordered_map<MACROMolTemplateKey, unsigned int, MACROMolTemplateKeyHash>
      d_keyToIndex;

 public:
  MACROMolTemplateLib() = default;
  MACROMolTemplateLib(const MACROMolTemplateLib &other) = delete;
  MACROMolTemplateLib(MACROMolTemplateLib &&other) noexcept = delete;
  MACROMolTemplateLib &operator=(MACROMolTemplateLib &&other) noexcept = delete;
  MACROMolTemplateLib &operator=(const MACROMolTemplateLib &) = delete;
  ~MACROMolTemplateLib() {}

  std::vector<std::unique_ptr<MACROMolTemplate>>::const_iterator begin() const {
    return d_templates.begin();
  }
  std::vector<std::unique_ptr<MACROMolTemplate>>::const_iterator end() const {
    return d_templates.end();
  }

  void addTemplate(std::unique_ptr<MACROMolTemplate> &templateMol);
  void copyTemplateLib(const MACROMolTemplateLib &libToCopy);

  //! Look up a template by (class, name). Templates are owned by the library
  //! and consumers only get const access, so this is only available on a
  //! const library; the non-const overload throws to catch lookups attempted
  //! while the library is still being mutated during construction.
  const MACROMolTemplate *getTemplate(const std::string &templateClass,
                                      const std::string &templateName) const {
    return find(templateClass, templateName);
  }
  MACROMolTemplate *getTemplate(const std::string & /*templateClass*/,
                                const std::string & /*templateName*/) {
    throw ValueErrorException(
        "getTemplate is not supported on non-const MACROMolTemplateLib");
  }

  unsigned int size() const { return d_templates.size(); }

  const RDKit::MACROMolTemplate *find(const std::string &templateClass,
                                      const std::string &templateName) const {
    auto iter = d_keyToIndex.find({templateClass, templateName});
    if (iter == d_keyToIndex.end()) {
      return nullptr;
    }
    return d_templates.at(iter->second).get();
  }

  bool contains(const std::string &templateClass,
                const std::string &templateName) const {
    return d_keyToIndex.contains({templateClass, templateName});
  }

  bool doesLibHaveCoords() const {
    for (auto const &macroTemplate : d_templates) {
      if (macroTemplate->getNumConformers() == 0) {
        return false;
      }
    }
    return true;
  }
};

class RDKIT_GRAPHMOL_EXPORT MACROMol : public RWMol {
 private:
  // elements (MACROMolTemplate items) of the library are owned by the library
  MACROMolTemplateLib d_templateLibrary;
  // pointer to the global library if loaded (Owned globally)
  const MACROMolTemplateLib *d_globalLib = nullptr;

 public:
  MACROMol() = default;

  //! Construct optionally attaching the global monomer template library.
  /*!
    \param useGlobalLibrary  if true, this molecule consults the global
                             MonomerLibrary's template library

    Defined out-of-line in the FileParsers layer (MonomerMol.cpp) because it
    depends on MonomerLibrary, which is not visible from GraphMol.
  */
  explicit MACROMol(bool useGlobalLibrary);

  MACROMol(const MACROMol &other) : RWMol((RWMol)other) {
    d_templateLibrary.copyTemplateLib(*other.getTemplateLibrary());
    d_globalLib = other.d_globalLib;  // share the same global; do not clone
  }
  MACROMol(MACROMol &&other) noexcept = delete;
  MACROMol &operator=(MACROMol &&other) noexcept = delete;

  MACROMol &operator=(const MACROMol &) = delete;  // disable assignment

  MACROMol(std::unique_ptr<RWMol> &rwMol)
      : RWMol(std::move(*(rwMol))) {
    rwMol = nullptr;
  }

  ~MACROMol() {}

  const MACROMolTemplate *getTemplate(unsigned int atomIdx) const;

  // provides non-const access to the template library for this MACROMol (local
  // library)
  MACROMolTemplateLib *getTemplateLibrary() { return &d_templateLibrary; }

  const MACROMolTemplateLib *getTemplateLibrary() const {
    return &d_templateLibrary;
  }

  // attach the shared, read-only global library this molecule will consult.
  // borrowed, not owned — caller guarantees it outlives this molecule.
  void setGlobalLibrary(const MACROMolTemplateLib *lib) {
    PRECONDITION(lib, "null template library");
    d_globalLib = lib;
  }

  // the following adds a template to the internal libraty for this MACROMol
  void addTemplate(std::unique_ptr<MACROMolTemplate> &templateMol) {
    PRECONDITION(templateMol, "bad template molecule");
    d_templateLibrary.addTemplate(templateMol);
  }

  unsigned int size() const { return d_templateLibrary.size(); }

  unsigned int addMacroAtom(std::string className, std::string templateName);

  void addMacroBond(unsigned int fromAtomIdx, unsigned int toAtomIdx,
                    Bond::BondType bondType, std::string fromConnectionPoint,
                    std::string toConnectionPoint);

  // ---- Query Operations ----

  [[nodiscard]] Chain getPolymer(std::string_view polymer_id) const;

  [[nodiscard]] std::vector<std::string> getPolymerIds() const;

  // ---- Chain Assignment ----

  // Discards existing chains and reassigns monomers to sequential chains where
  // monomers are reordered based on connectivity.
  void assignChains();
};
typedef boost::shared_ptr<MACROMol> MACROMol_SPTR;
typedef boost::shared_ptr<MACROMolTemplate> MACROMolTemplate_SPTR;
}  // namespace RDKit

#endif
