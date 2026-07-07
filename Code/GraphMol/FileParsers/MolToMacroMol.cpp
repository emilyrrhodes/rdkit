//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/FileParsers/MolToMacroMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MacroMol.h>
#include <GraphMol/MacroMolTemplate.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/Invariant.h>

#include <algorithm>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace RDKit {

namespace {

struct AttachmentPointSpec {
  unsigned int sourceAtomIdx;
  int attachmentIdx;
};

struct MacroAtomSpec {
  MonomerClass monomerClass = MonomerClass::OTHER;
  std::string symbol;
};

struct AtomSpec {
  std::vector<unsigned int> sourceAtomIndices;
  std::optional<MacroAtomSpec> macroAtom;
  std::vector<AttachmentPointSpec> attachmentPoints;

  std::optional<int> getAttachmentIdx(unsigned int sourceAtomIdx) const {
    auto it = std::find_if(attachmentPoints.begin(), attachmentPoints.end(),
                           [sourceAtomIdx](const AttachmentPointSpec &spec) {
                             return spec.sourceAtomIdx == sourceAtomIdx;
                           });
    if (it == attachmentPoints.end()) {
      return std::nullopt;
    }
    return it->attachmentIdx;
  }
};

struct BondSpec {
  unsigned int sourceBondIdx;
  unsigned int beginOutputAtomIdx;
  unsigned int endOutputAtomIdx;
  std::optional<int> beginAttachPt;
  std::optional<int> endAttachPt;
};

struct MacroMolSpec {
  explicit MacroMolSpec(unsigned int numSourceAtoms)
      : sourceAtomToOutputAtom(numSourceAtoms) {}

  std::vector<AtomSpec> atoms;
  std::vector<BondSpec> bonds;
  std::vector<std::optional<unsigned int>> sourceAtomToOutputAtom;

  unsigned int addAtom(AtomSpec atomSpec) {
    std::sort(atomSpec.sourceAtomIndices.begin(),
              atomSpec.sourceAtomIndices.end());

    const unsigned int outputAtomIdx = atoms.size();
    atoms.push_back(std::move(atomSpec));
    for (const auto sourceAtomIdx : atoms.back().sourceAtomIndices) {
      sourceAtomToOutputAtom[sourceAtomIdx] = outputAtomIdx;
    }
    return outputAtomIdx;
  }

  bool hasOutputAtomForSourceAtom(unsigned int sourceAtomIdx) const {
    return sourceAtomIdx < sourceAtomToOutputAtom.size() &&
           sourceAtomToOutputAtom[sourceAtomIdx].has_value();
  }

  unsigned int getOutputAtomIdxBySourceAtom(unsigned int sourceAtomIdx) const {
    CHECK_INVARIANT(sourceAtomIdx < sourceAtomToOutputAtom.size(),
                    "source atom index out of range");
    const auto outputAtomIdx = sourceAtomToOutputAtom[sourceAtomIdx];
    CHECK_INVARIANT(outputAtomIdx, "source atom has no output atom");
    return *outputAtomIdx;
  }
};

struct MainGroupQuery {
  std::unique_ptr<RWMol> query;
  std::vector<unsigned int> queryToTemplateAtom;
};

MainGroupQuery makeMainGroupQuery(const MacroMolTemplate &templ,
                                  const SubstanceGroup &mainSgroup) {
  const auto &mainAtoms = mainSgroup.getAtoms();
  std::vector<unsigned int> sortedMainAtoms(mainAtoms.begin(), mainAtoms.end());
  std::sort(sortedMainAtoms.begin(), sortedMainAtoms.end());

  const std::unordered_set<unsigned int> mainAtomSet(sortedMainAtoms.begin(),
                                                     sortedMainAtoms.end());
  auto query = std::make_unique<RWMol>(templ);
  query->beginBatchEdit();
  for (auto atom : query->atoms()) {
    if (mainAtomSet.find(atom->getIdx()) == mainAtomSet.end()) {
      query->removeAtom(atom);
    }
  }
  query->commitBatchEdit();
  query->updatePropertyCache(false);

  return MainGroupQuery{std::move(query), std::move(sortedMainAtoms)};
}

std::vector<MatchVectType> findTemplateMatches(const ROMol &mol,
                                               const RWMol &query) {
  SubstructMatchParameters params;
  params.recursionPossible = false;
  params.useChirality = true;
  params.useQueryQueryMatches = false;
  params.uniquify = false;
  params.maxMatches = 0;  // no limit
  auto matches = SubstructMatch(mol, query, params);

  // Claim atoms in a deterministic order so the result does not depend on the
  // enumeration order of symmetric matches. Keep all automorphic mappings so
  // attachment-point validation can choose a valid orientation.
  auto matchKey = [](const MatchVectType &match) {
    std::vector<int> sourceAtoms;
    std::vector<std::pair<int, int>> queryToSource;
    sourceAtoms.reserve(match.size());
    queryToSource.reserve(match.size());
    for (const auto &pair : match) {
      sourceAtoms.push_back(pair.second);
      queryToSource.push_back(pair);
    }
    std::sort(sourceAtoms.begin(), sourceAtoms.end());
    std::sort(queryToSource.begin(), queryToSource.end());
    return std::make_pair(sourceAtoms, queryToSource);
  };
  std::sort(matches.begin(), matches.end(),
            [&matchKey](const MatchVectType &lhs, const MatchVectType &rhs) {
              return matchKey(lhs) < matchKey(rhs);
            });
  return matches;
}

// A match is only usable if every bond that crosses the monomer boundary is
// accounted for by an attachment point. The design invariant is that each atom
// has at most one attachment point, so an attach atom may have at most one
// external bond and a non-attach atom may have none.
bool matchViolatesAttachmentInvariant(const ROMol &mol,
                                      const AtomSpec &atomSpec) {
  std::vector<unsigned int> sortedSourceAtoms(atomSpec.sourceAtomIndices);
  std::sort(sortedSourceAtoms.begin(), sortedSourceAtoms.end());

  for (const auto sourceAtomIdx : atomSpec.sourceAtomIndices) {
    unsigned int externalBondCount = 0;
    const auto *sourceAtom = mol.getAtomWithIdx(sourceAtomIdx);
    for (const auto *bond : mol.atomBonds(sourceAtom)) {
      const auto otherSourceAtomIdx = bond->getOtherAtomIdx(sourceAtomIdx);
      if (!std::binary_search(sortedSourceAtoms.begin(),
                              sortedSourceAtoms.end(), otherSourceAtomIdx)) {
        ++externalBondCount;
      }
    }

    const unsigned int allowedExternalBonds =
        atomSpec.getAttachmentIdx(sourceAtomIdx) ? 1 : 0;
    if (externalBondCount > allowedExternalBonds) {
      return true;
    }
  }
  return false;
}

// Attachment-point ids are required to be numeric (they name the integer attach
// point used by the MacroMol bonds). A template whose id is non-numeric is
// malformed for this API; callers must convert before use.
int parseAttachmentId(const std::string &id) {
  try {
    size_t consumed = 0;
    const int value = std::stoi(id, &consumed);
    if (consumed == id.size()) {
      return value;
    }
  } catch (const std::exception &) {
    // fall through to the throw below
  }
  throw ValueErrorException(
      "MolToMacroMol requires numeric attachment-point ids; got '" + id + "'");
}

std::optional<AtomSpec> buildAtomSpecForMatch(
    const MatchVectType &match, const ROMol &mol, const MacroMolTemplate &templ,
    const MainGroupQuery &mainQuery, const MacroMolSpec &macroMolSpec) {
  const auto *mainSgroup = templ.getMainSgroup();

  // Reject matches that overlap atoms already claimed by an earlier template.
  // Templates are processed largest-first, so the earlier (larger) claim wins.
  for (const auto &matchPair : match) {
    if (macroMolSpec.hasOutputAtomForSourceAtom(matchPair.second)) {
      return std::nullopt;
    }
  }

  AtomSpec atomSpec;
  std::unordered_map<unsigned int, unsigned int> templateAtomToSourceAtom;
  for (const auto &matchPair : match) {
    const auto queryAtomIdx = static_cast<unsigned int>(matchPair.first);
    const auto sourceAtomIdx = static_cast<unsigned int>(matchPair.second);
    atomSpec.sourceAtomIndices.push_back(sourceAtomIdx);

    const auto templateAtomIdx = mainQuery.queryToTemplateAtom[queryAtomIdx];
    templateAtomToSourceAtom[templateAtomIdx] = sourceAtomIdx;
  }

  std::unordered_set<unsigned int> attachAtomsSeen;
  for (const auto &attachPoint : mainSgroup->getAttachPoints()) {
    const auto it = templateAtomToSourceAtom.find(attachPoint.aIdx);
    if (it == templateAtomToSourceAtom.end()) {
      // The attachment point references an atom outside the main group: a
      // malformed template. Reject the match rather than corrupt the output.
      return std::nullopt;
    }
    const auto sourceAtomIdx = it->second;
    if (!attachAtomsSeen.insert(sourceAtomIdx).second) {
      throw ValueErrorException(
          "MolToMacroMol does not support template atoms carrying more than "
          "one attachment point");
    }

    AttachmentPointSpec attachPtSpec;
    attachPtSpec.sourceAtomIdx = sourceAtomIdx;
    attachPtSpec.attachmentIdx = parseAttachmentId(attachPoint.id);
    atomSpec.attachmentPoints.push_back(attachPtSpec);
  }

  if (matchViolatesAttachmentInvariant(mol, atomSpec)) {
    return std::nullopt;
  }

  return atomSpec;
}

void recordAcceptedAtomSpec(const MacroMolEntry &entry, AtomSpec atomSpec,
                            MacroMolSpec &macroMolSpec) {
  atomSpec.macroAtom =
      MacroAtomSpec{monomerClassFromString(entry.monomerClass), entry.symbol};
  macroMolSpec.addAtom(std::move(atomSpec));
}

void addTemplateMatchesToSpec(const ROMol &mol, const MacroMolEntry &entry,
                              MacroMolSpec &macroMolSpec) {
  if (!entry.molTemplate) {
    return;
  }

  const auto &templ = *entry.molTemplate;
  const auto *mainSgroup = templ.getMainSgroup();
  if (!mainSgroup) {
    return;
  }

  const auto mainQuery = makeMainGroupQuery(templ, *mainSgroup);
  for (const auto &match : findTemplateMatches(mol, *mainQuery.query)) {
    auto atomSpec =
        buildAtomSpecForMatch(match, mol, templ, mainQuery, macroMolSpec);
    if (atomSpec) {
      recordAcceptedAtomSpec(entry, std::move(*atomSpec), macroMolSpec);
    }
  }
}

void addUnmatchedAtomsToSpec(const ROMol &mol, MacroMolSpec &macroMolSpec) {
  for (const auto *atom : mol.atoms()) {
    const auto sourceAtomIdx = atom->getIdx();
    if (macroMolSpec.hasOutputAtomForSourceAtom(sourceAtomIdx)) {
      continue;
    }

    AtomSpec atomSpec;
    atomSpec.sourceAtomIndices.push_back(sourceAtomIdx);
    macroMolSpec.addAtom(std::move(atomSpec));
  }
}

void groupSourceAtoms(const ROMol &mol,
                      const MacroMolTemplateLibrary &templates,
                      MacroMolSpec &macroMolSpec) {
  // templates.entries() is ordered largest-main-group-first, so larger monomers
  // claim their atoms before smaller ones.
  for (const auto &entry : templates.entries()) {
    addTemplateMatchesToSpec(mol, *entry, macroMolSpec);
  }
  // Must run after all template matches so it only picks up genuine leftovers.
  addUnmatchedAtomsToSpec(mol, macroMolSpec);
}

void deriveBonds(const ROMol &mol, MacroMolSpec &macroMolSpec) {
  for (const auto *bond : mol.bonds()) {
    const auto sourceBeginAtomIdx = bond->getBeginAtomIdx();
    const auto sourceEndAtomIdx = bond->getEndAtomIdx();
    const auto beginOutputAtomIdx =
        macroMolSpec.getOutputAtomIdxBySourceAtom(sourceBeginAtomIdx);
    const auto endOutputAtomIdx =
        macroMolSpec.getOutputAtomIdxBySourceAtom(sourceEndAtomIdx);
    if (beginOutputAtomIdx == endOutputAtomIdx) {
      // Internal to a single output (macro) atom; not a bond in the MacroMol.
      continue;
    }

    macroMolSpec.bonds.push_back(
        BondSpec{bond->getIdx(), beginOutputAtomIdx, endOutputAtomIdx,
                 macroMolSpec.atoms[beginOutputAtomIdx].getAttachmentIdx(
                     sourceBeginAtomIdx),
                 macroMolSpec.atoms[endOutputAtomIdx].getAttachmentIdx(
                     sourceEndAtomIdx)});
  }
}

// Output atoms are created in ascending order of their smallest source atom
// index, giving a stable result-atom numbering.
std::vector<unsigned int> outputAtomOrder(const MacroMolSpec &macroMolSpec) {
  std::vector<unsigned int> order(macroMolSpec.atoms.size());
  std::iota(order.begin(), order.end(), 0u);
  std::sort(order.begin(), order.end(),
            [&macroMolSpec](unsigned int lhs, unsigned int rhs) {
              return macroMolSpec.atoms[lhs].sourceAtomIndices.front() <
                     macroMolSpec.atoms[rhs].sourceAtomIndices.front();
            });
  return order;
}

// Copies a regular (non-macro) source bond into the MacroMol, preserving all of
// its metadata (direction, stereo, query state, aromatic/conjugation flags,
// custom properties) and remapping its atoms to the output numbering.
void copyRegularBond(MacroMol &macroMol, const Bond &sourceBond,
                     unsigned int resultBeginAtomIdx,
                     unsigned int resultEndAtomIdx,
                     const MacroMolSpec &macroMolSpec,
                     const std::vector<unsigned int> &resultAtomByOutputIdx) {
  std::unique_ptr<Bond> newBond(sourceBond.copy());
  newBond->setBeginAtomIdx(resultBeginAtomIdx);
  newBond->setEndAtomIdx(resultEndAtomIdx);
  for (auto &stereoAtom : newBond->getStereoAtoms()) {
    if (stereoAtom >= 0) {
      stereoAtom = static_cast<int>(
          resultAtomByOutputIdx[macroMolSpec.getOutputAtomIdxBySourceAtom(
              static_cast<unsigned int>(stereoAtom))]);
    }
  }
  macroMol.addBond(newBond.get());
}

std::unique_ptr<MacroMol> materialize(const ROMol &mol,
                                      const MacroMolSpec &macroMolSpec) {
  auto macroMol = std::make_unique<MacroMol>();

  std::vector<unsigned int> resultAtomByOutputIdx(macroMolSpec.atoms.size());
  for (const auto outputIdx : outputAtomOrder(macroMolSpec)) {
    const auto &atomSpec = macroMolSpec.atoms[outputIdx];
    unsigned int resultAtomIdx;
    if (atomSpec.macroAtom) {
      resultAtomIdx = macroMol->addMacroAtom(atomSpec.macroAtom->monomerClass,
                                             atomSpec.macroAtom->symbol);
    } else {
      const auto *sourceAtom =
          mol.getAtomWithIdx(atomSpec.sourceAtomIndices.front());
      auto newAtom = std::unique_ptr<Atom>(sourceAtom->copy());
      resultAtomIdx = macroMol->addAtom(newAtom.get(), true, true);
      newAtom.release();
    }
    resultAtomByOutputIdx[outputIdx] = resultAtomIdx;
  }

  for (const auto &bondSpec : macroMolSpec.bonds) {
    const auto resultBeginAtomIdx =
        resultAtomByOutputIdx[bondSpec.beginOutputAtomIdx];
    const auto resultEndAtomIdx =
        resultAtomByOutputIdx[bondSpec.endOutputAtomIdx];
    const auto *sourceBond = mol.getBondWithIdx(bondSpec.sourceBondIdx);
    const auto bondType = sourceBond->getBondType();

    if (bondSpec.beginAttachPt && bondSpec.endAttachPt) {
      macroMol->addMacroBond(resultBeginAtomIdx, resultEndAtomIdx,
                             *bondSpec.beginAttachPt, *bondSpec.endAttachPt,
                             bondType);
    } else if (bondSpec.beginAttachPt) {
      macroMol->addMacroAtomToAtomBond(resultBeginAtomIdx, resultEndAtomIdx,
                                       *bondSpec.beginAttachPt, bondType);
    } else if (bondSpec.endAttachPt) {
      macroMol->addAtomToMacroAtomBond(resultBeginAtomIdx, resultEndAtomIdx,
                                       *bondSpec.endAttachPt, bondType);
    } else {
      // Both endpoints are regular atoms: preserve the full source bond.
      copyRegularBond(*macroMol, *sourceBond, resultBeginAtomIdx,
                      resultEndAtomIdx, macroMolSpec, resultAtomByOutputIdx);
    }
  }

  return macroMol;
}

}  // namespace

std::unique_ptr<MacroMol> MolToMacroMol(
    const ROMol &mol, const MacroMolTemplateLibrary &templates) {
  MacroMolSpec macroMolSpec(mol.getNumAtoms());
  groupSourceAtoms(mol, templates, macroMolSpec);
  deriveBonds(mol, macroMolSpec);
  return materialize(mol, macroMolSpec);
}

}  // namespace RDKit
