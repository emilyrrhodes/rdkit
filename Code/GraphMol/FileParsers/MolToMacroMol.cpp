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
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/Invariant.h>

#include <algorithm>
#include <array>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace RDKit {

namespace {

const std::string ORIG_ATOM_ID = "origAtomId";

struct DirectedBond {
  unsigned int from;
  unsigned int to;

  bool operator<(const DirectedBond &other) const {
    if (from != other.from) {
      return from < other.from;
    }
    return to < other.to;
  }
};

struct GroupBond {
  unsigned int beginSourceAtomIdx;
  unsigned int endSourceAtomIdx;
  std::optional<int> beginAttachPt;
  std::optional<int> endAttachPt;
  Bond::BondType bondType;
};

struct MacroAtomSpec {
  MonomerClass monomerClass = MonomerClass::OTHER;
  std::string name;
};

struct AtomGroup {
  std::vector<unsigned int> sourceAtomIndices;
  std::optional<MacroAtomSpec> macroAtom;
  std::map<DirectedBond, int> sourceBondToAttachPt;
};

struct GroupingPlan {
  std::map<unsigned int, AtomGroup> groupsByStartAtom;
  std::vector<std::optional<unsigned int>> sourceAtomToGroupStartAtom;
  std::vector<GroupBond> bonds;
};

struct AcceptedMatch {
  std::vector<unsigned int> sourceAtomIndices;
  std::map<DirectedBond, int> attachPoints;
};

const std::array<std::pair<const char *, MonomerClass>, 4> monomerClasses = {{
    {"AA", MonomerClass::AA},
    {"NA", MonomerClass::NA},
    {"CHEM", MonomerClass::CHEM},
    {"OTHER", MonomerClass::OTHER},
}};

int parseAttachPointId(const std::string &attachPointId) {
  if (attachPointId.empty()) {
    throw FileParseException("empty macro attachment point id");
  }
  try {
    return std::stoi(attachPointId);
  } catch (const std::exception &) {
    throw FileParseException("non-numeric macro attachment point id: " +
                             attachPointId);
  }
}

MonomerClass stringToMonomerClass(const std::string &monomerClass) {
  for (const auto &[name, value] : monomerClasses) {
    if (monomerClass == name) {
      return value;
    }
  }
  return MonomerClass::OTHER;
}

std::unique_ptr<RWMol> makeMainGroupQuery(const MacroMolTemplate &templ,
                                          const SubstanceGroup &mainSgroup) {
  auto query = std::make_unique<RWMol>(templ);
  for (auto atom : query->atoms()) {
    atom->setProp(ORIG_ATOM_ID, atom->getIdx());
  }

  const auto &mainAtoms = mainSgroup.getAtoms();
  const std::unordered_set<unsigned int> mainAtomSet(mainAtoms.begin(),
                                                     mainAtoms.end());
  query->beginBatchEdit();
  for (auto atom : query->atoms()) {
    if (mainAtomSet.find(atom->getIdx()) == mainAtomSet.end()) {
      query->removeAtom(atom);
    }
  }
  query->commitBatchEdit();
  query->updatePropertyCache(false);
  return query;
}

std::vector<MatchVectType> findTemplateMatches(const ROMol &mol,
                                               const RWMol &query) {
  SubstructMatchParameters params;
  params.recursionPossible = false;
  params.useChirality = true;
  params.useQueryQueryMatches = false;
  params.maxMatches = 0;
  return SubstructMatch(mol, query, params);
}

bool leavingGroupMatchesSourceAtom(const ROMol &mol,
                                      const MacroMolTemplate &templ,
                                      unsigned int sourceAtomIdx,
                                      const SubstanceGroup::AttachPoint &sap) {
  if (sap.lvIdx < 0 ||
      static_cast<unsigned int>(sap.lvIdx) >= templ.getNumAtoms()) {
    return false;
  }

  const auto *sourceAtom = mol.getAtomWithIdx(sourceAtomIdx);
  const auto *templateAtom =
      templ.getAtomWithIdx(static_cast<unsigned int>(sap.lvIdx));
  return sourceAtom->getDegree() == 1 && templateAtom->getDegree() == 1 &&
         sourceAtom->getAtomicNum() == templateAtom->getAtomicNum() &&
         sourceAtom->getTotalNumHs() == templateAtom->getTotalNumHs();
}

bool handleExternalBond(
    const ROMol &mol, const MacroMolTemplate &templ,
    unsigned int sourceAtomIdx, unsigned int externalAtomIdx,
    unsigned int templateAtomIdx,
    const std::vector<SubstanceGroup::AttachPoint> &attachPoints,
    std::vector<bool> &attachPointUsed, const GroupingPlan &plan,
    AcceptedMatch &acceptedMatch) {
  for (unsigned int attachPointIdx = 0; attachPointIdx < attachPoints.size();
       ++attachPointIdx) {
    const auto &attachPoint = attachPoints[attachPointIdx];
    if (attachPointUsed[attachPointIdx] ||
        attachPoint.aIdx != templateAtomIdx) {
      continue;
    }
    if (leavingGroupMatchesSourceAtom(mol, templ, externalAtomIdx,
                                         attachPoint)) {
      attachPointUsed[attachPointIdx] = true;
      if (plan.sourceAtomToGroupStartAtom[externalAtomIdx]) {
        return false;
      }
      acceptedMatch.sourceAtomIndices.push_back(externalAtomIdx);
      return true;
    }
    attachPointUsed[attachPointIdx] = true;
    acceptedMatch.attachPoints[{sourceAtomIdx, externalAtomIdx}] =
        parseAttachPointId(attachPoint.id);
    return true;
  }

  return false;
}

std::optional<AcceptedMatch> validateMatch(const MatchVectType &match,
                                           const ROMol &mol,
                                           const MacroMolTemplate &templ,
                                           const RWMol &query,
                                           const GroupingPlan &plan) {
  const auto *mainSgroup = templ.getMainSgroup();
  if (!mainSgroup) {
    return std::nullopt;
  }

  AcceptedMatch acceptedMatch;
  std::map<unsigned int, unsigned int> sourceToQueryAtom;
  for (const auto &matchPair : match) {
    const auto queryAtomIdx = matchPair.first;
    const auto sourceAtomIdx = matchPair.second;
    if (plan.sourceAtomToGroupStartAtom[sourceAtomIdx]) {
      return std::nullopt;
    }
    sourceToQueryAtom[sourceAtomIdx] = queryAtomIdx;
    acceptedMatch.sourceAtomIndices.push_back(sourceAtomIdx);
  }

  const auto attachPoints = mainSgroup->getAttachPoints();
  std::vector<bool> attachPointUsed(attachPoints.size(), false);

  for (const auto &sourceAndQueryAtom : sourceToQueryAtom) {
    const auto sourceAtomIdx = sourceAndQueryAtom.first;
    const auto queryAtomIdx = sourceAndQueryAtom.second;
    const auto *sourceAtom = mol.getAtomWithIdx(sourceAtomIdx);
    const auto *queryAtom = query.getAtomWithIdx(queryAtomIdx);
    const auto templateAtomIdx = queryAtom->getProp<unsigned int>(ORIG_ATOM_ID);

    for (const auto *bond : mol.atomBonds(sourceAtom)) {
      const auto neighborIdx = bond->getOtherAtomIdx(sourceAtomIdx);
      auto neighborQueryIt = sourceToQueryAtom.find(neighborIdx);
      if (neighborQueryIt != sourceToQueryAtom.end()) {
        if (!query.getBondBetweenAtoms(queryAtomIdx, neighborQueryIt->second)) {
          return std::nullopt;
        }
        continue;
      }

      if (!handleExternalBond(mol, templ, sourceAtomIdx, neighborIdx,
                              templateAtomIdx, attachPoints, attachPointUsed,
                              plan, acceptedMatch)) {
        return std::nullopt;
      }
    }
  }

  return acceptedMatch;
}

void recordAcceptedMatchAtoms(const MacroMolEntry &entry,
                              const AcceptedMatch &acceptedMatch,
                              GroupingPlan &plan) {
  auto sourceAtomIndices = acceptedMatch.sourceAtomIndices;
  std::sort(sourceAtomIndices.begin(), sourceAtomIndices.end());
  sourceAtomIndices.erase(
      std::unique(sourceAtomIndices.begin(), sourceAtomIndices.end()),
      sourceAtomIndices.end());
  CHECK_INVARIANT(!sourceAtomIndices.empty(),
                  "macro atom group has no source atoms");

  const auto groupStartAtom = sourceAtomIndices.front();
  AtomGroup group;
  group.sourceAtomIndices = sourceAtomIndices;
  group.macroAtom =
      MacroAtomSpec{stringToMonomerClass(entry.monomerClass), entry.symbol};
  group.sourceBondToAttachPt = acceptedMatch.attachPoints;
  const auto insertResult =
      plan.groupsByStartAtom.emplace(groupStartAtom, std::move(group));
  CHECK_INVARIANT(insertResult.second, "source atom group already exists");

  for (const auto sourceAtomIdx : sourceAtomIndices) {
    auto &groupStart = plan.sourceAtomToGroupStartAtom[sourceAtomIdx];
    CHECK_INVARIANT(!groupStart || *groupStart == groupStartAtom,
                    "source atom assigned to a different atom group");
    groupStart = groupStartAtom;
  }
}

void addTemplateMatchesToPlan(const ROMol &mol, const MacroMolEntry &entry,
                              GroupingPlan &plan) {
  if (!entry.molTemplate) {
    return;
  }

  const auto &templ = *entry.molTemplate;
  const auto *mainSgroup = templ.getMainSgroup();
  if (!mainSgroup) {
    return;
  }

  auto query = makeMainGroupQuery(templ, *mainSgroup);
  for (const auto &match : findTemplateMatches(mol, *query)) {
    auto acceptedMatch = validateMatch(match, mol, templ, *query, plan);
    if (acceptedMatch) {
      recordAcceptedMatchAtoms(entry, *acceptedMatch, plan);
    }
  }
}

void addUnmatchedAtomsToPlan(const ROMol &mol, GroupingPlan &plan) {
  for (const auto *atom : mol.atoms()) {
    const auto sourceAtomIdx = atom->getIdx();
    if (plan.sourceAtomToGroupStartAtom[sourceAtomIdx]) {
      continue;
    }

    AtomGroup group;
    group.sourceAtomIndices.push_back(sourceAtomIdx);
    const auto insertResult =
        plan.groupsByStartAtom.emplace(sourceAtomIdx, std::move(group));
    CHECK_INVARIANT(insertResult.second, "source atom group already exists");
    plan.sourceAtomToGroupStartAtom[sourceAtomIdx] = sourceAtomIdx;
  }
}

void addBondsToPlan(const ROMol &mol, GroupingPlan &plan) {
  for (const auto *bond : mol.bonds()) {
    const auto sourceBeginAtomIdx = bond->getBeginAtomIdx();
    const auto sourceEndAtomIdx = bond->getEndAtomIdx();
    const auto &beginGroupStartAtom =
        plan.sourceAtomToGroupStartAtom[sourceBeginAtomIdx];
    CHECK_INVARIANT(beginGroupStartAtom, "begin source atom has no atom group");
    const auto &endGroupStartAtom =
        plan.sourceAtomToGroupStartAtom[sourceEndAtomIdx];
    CHECK_INVARIANT(endGroupStartAtom, "end source atom has no atom group");
    if (*beginGroupStartAtom == *endGroupStartAtom) {
      continue;
    }

    const auto beginGroupIt = plan.groupsByStartAtom.find(*beginGroupStartAtom);
    CHECK_INVARIANT(beginGroupIt != plan.groupsByStartAtom.end(),
                    "missing begin atom group");
    const auto endGroupIt = plan.groupsByStartAtom.find(*endGroupStartAtom);
    CHECK_INVARIANT(endGroupIt != plan.groupsByStartAtom.end(),
                    "missing end atom group");
    const auto &beginGroup = beginGroupIt->second;
    const auto &endGroup = endGroupIt->second;

    GroupBond groupBond;
    groupBond.beginSourceAtomIdx = *beginGroupStartAtom;
    groupBond.endSourceAtomIdx = *endGroupStartAtom;
    groupBond.bondType = bond->getBondType();

    const auto beginAttachPtIt = beginGroup.sourceBondToAttachPt.find(
        {sourceBeginAtomIdx, sourceEndAtomIdx});
    if (beginAttachPtIt != beginGroup.sourceBondToAttachPt.end()) {
      groupBond.beginAttachPt = beginAttachPtIt->second;
    }
    const auto endAttachPtIt = endGroup.sourceBondToAttachPt.find(
        {sourceEndAtomIdx, sourceBeginAtomIdx});
    if (endAttachPtIt != endGroup.sourceBondToAttachPt.end()) {
      groupBond.endAttachPt = endAttachPtIt->second;
    }

    plan.bonds.push_back(std::move(groupBond));
  }
}

std::unique_ptr<MacroMol> buildMacroMolFromPlan(const ROMol &mol,
                                                const GroupingPlan &plan) {
  auto macroMol = std::make_unique<MacroMol>();
  std::map<unsigned int, unsigned int> resultAtomByGroupStart;

  for (unsigned int sourceAtomIdx = 0; sourceAtomIdx < mol.getNumAtoms();
       ++sourceAtomIdx) {
    const auto &groupStartAtom =
        plan.sourceAtomToGroupStartAtom[sourceAtomIdx];
    CHECK_INVARIANT(groupStartAtom, "source atom has no atom group");
    if (*groupStartAtom != sourceAtomIdx) {
      continue;
    }

    const auto groupIt = plan.groupsByStartAtom.find(*groupStartAtom);
    CHECK_INVARIANT(groupIt != plan.groupsByStartAtom.end(),
                    "missing atom group");
    const auto &group = groupIt->second;
    unsigned int resultAtomIdx;
    if (group.macroAtom) {
      resultAtomIdx = macroMol->addMacroAtom(group.macroAtom->monomerClass,
                                             group.macroAtom->name);
    } else {
      CHECK_INVARIANT(group.sourceAtomIndices.size() == 1,
                      "atomistic groups must contain a single source atom");
      const auto *sourceAtom =
          mol.getAtomWithIdx(group.sourceAtomIndices.front());
      resultAtomIdx = macroMol->addAtom(new Atom(*sourceAtom), true, true);
    }
    resultAtomByGroupStart[*groupStartAtom] = resultAtomIdx;
  }

  for (const auto &groupBond : plan.bonds) {
    const auto resultBeginAtomIdx =
        resultAtomByGroupStart.at(groupBond.beginSourceAtomIdx);
    const auto resultEndAtomIdx =
        resultAtomByGroupStart.at(groupBond.endSourceAtomIdx);
    if (groupBond.beginAttachPt && groupBond.endAttachPt) {
      macroMol->addMacroBond(resultBeginAtomIdx, resultEndAtomIdx,
                             *groupBond.beginAttachPt, *groupBond.endAttachPt,
                             groupBond.bondType);
    } else if (groupBond.beginAttachPt) {
      macroMol->addMacroAtomToAtomBond(resultBeginAtomIdx, resultEndAtomIdx,
                                       *groupBond.beginAttachPt,
                                       groupBond.bondType);
    } else if (groupBond.endAttachPt) {
      macroMol->addAtomToMacroAtomBond(resultBeginAtomIdx, resultEndAtomIdx,
                                       *groupBond.endAttachPt,
                                       groupBond.bondType);
    } else {
      macroMol->addBond(resultBeginAtomIdx, resultEndAtomIdx,
                        groupBond.bondType);
    }
  }

  return macroMol;
}

}  // namespace

std::unique_ptr<MacroMol> MolToMacroMol(
    const ROMol &mol, const MacroMolTemplateLibrary &templates) {
  GroupingPlan plan;
  plan.sourceAtomToGroupStartAtom.resize(mol.getNumAtoms());
  for (const auto &entry : templates.entries()) {
    addTemplateMatchesToPlan(mol, *entry, plan);
  }
  addUnmatchedAtomsToPlan(mol, plan);
  addBondsToPlan(mol, plan);
  return buildMacroMolFromPlan(mol, plan);
}

}  // namespace RDKit
