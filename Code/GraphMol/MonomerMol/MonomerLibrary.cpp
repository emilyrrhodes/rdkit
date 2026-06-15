//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MonomerLibrary.h"

#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {

// Static member definitions
std::unique_ptr<MonomerLibrary> MonomerLibrary::s_globalLibrary;
std::unique_ptr<MACROMolTemplateLib> MonomerLibrary::s_globalMACROMOLTemplateLib
;
std::once_flag MonomerLibrary::s_globalLibraryOnce;

// Built-in monomer definitions (symbol -> SMILES)
namespace {

    /* items are the symbol, monomer_class,  the PDB code, and the smiles */
const std::list<std::tuple<std::string, std::string, std::string, std::string>> builtin_monomer_data = {
    {"A", "PEPTIDE", "ALA" ,"C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CB :1.pdbName. CA "
          ":2.pdbName. N  :3.pdbName. H  :4.pdbName. C  :5.pdbName. O  "
          ":6.pdbName. OXT|"},
    {"C", "PEPTIDE", "CYS", "O=C([C@H](CS[H:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. SG "
          ":5.pdbName. HG :6.pdbName. N  :7.pdbName. H  :8.pdbName. OXT|"},
    {"D", "PEPTIDE", "ASP/ASH", 
     "O=C([C@H](CC(=O)[OH:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
     "OD1:6.pdbName. OD2:7.pdbName. N  :8.pdbName. H  :9.pdbName. OXT|"},
    {"E", "PEPTIDE", "GLU/GLH", "O=C([C@H](CCC(=O)[OH:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
          ":5.pdbName. CD :6.pdbName. OE1:7.pdbName. OE2:8.pdbName. N  "
          ":9.pdbName. H  :10.pdbName. OXT|"},
    {"F", "PEPTIDE", "PHE", 
     "O=C([C@H](Cc1ccccc1)N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C "
     " :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
     "CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. CE2:9.pdbName. "
     "CD2:10.pdbName. N  :11.pdbName. H  :12.pdbName. OXT|"},
    {"G", "PEPTIDE", "GLY", "O=C(CN[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C  "
          ":2.pdbName. CA :3.pdbName. N  :4.pdbName. H  :5.pdbName. OXT|"},
    {"H", "PEPTIDE", "HIS/HIE/HIP/HID/HSE/HSD/HSP", "O=C([C@H](Cc1cnc[nH]1)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
          ":5.pdbName. CD2:6.pdbName. NE2:7.pdbName. CE1:8.pdbName. "
          "ND1:9.pdbName. N  :10.pdbName. H  :11.pdbName. OXT|"},
    {"I", "PEPTIDE", "ILE", 
     "CC[C@H](C)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CD1:1.pdbName. "
     "CG1:2.pdbName. CB :3.pdbName. CG2:4.pdbName. CA :5.pdbName. N  "
     ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|"},
    {"K", "PEPTIDE", "LYS/LYN", "O=C([C@H](CCCCN[H:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
          ":5.pdbName. CD :6.pdbName. CE :7.pdbName. NZ :8.pdbName. "
          "HZ1:9.pdbName. N  :10.pdbName. H  :11.pdbName. OXT|"},
    {"L", "PEPTIDE", "LEU", "CC(C)C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CD1:1.pdbName. "
          "CG :2.pdbName. CD2:3.pdbName. CB :4.pdbName. CA :5.pdbName. N  "
          ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|"},
    {"M", "PEPTIDE", "MET",  "CSCC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CE :1.pdbName. "
          "SD :2.pdbName. CG :3.pdbName. CB :4.pdbName. CA :5.pdbName. N  "
          ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|"},
    {"N", "PEPTIDE", "ASN", 
     "NC(=O)C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. ND2:1.pdbName. CG "
     ":2.pdbName. OD1:3.pdbName. CB :4.pdbName. CA :5.pdbName. N  :6.pdbName. "
     "H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|"},
    {"O", "PEPTIDE", "PYL", "C[C@@H]1CC=N[C@H]1C(=O)NCCCC[C@H](N[H:1])C(=O)[OH:2] "
          "|atomProp:0.pdbName. CB2:1.pdbName. CG2:2.pdbName. CD2:3.pdbName. "
          "CE2:4.pdbName. N2 :5.pdbName. CA2:6.pdbName. C2 :7.pdbName. O2 "
          ":8.pdbName. NZ :9.pdbName. CE :10.pdbName. CD :11.pdbName. CG "
          ":12.pdbName. CB :13.pdbName. CA :14.pdbName. N  :15.pdbName. H  "
          ":16.pdbName. C  :17.pdbName. O  :18.pdbName. OXT|"},
    {"P", "PEPTIDE", "PRO",  "O=C([C@@H]1CCCN1[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C "
          " :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD "
          ":6.pdbName. N  :7.pdbName. H  :8.pdbName. OXT|"},
    {"Q", "PEPTIDE", "GLN",
     "NC(=O)CC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. NE2:1.pdbName. CD "
     ":2.pdbName. OE1:3.pdbName. CG :4.pdbName. CB :5.pdbName. CA :6.pdbName. "
     "N  :7.pdbName. H  :8.pdbName. C  :9.pdbName. O  :10.pdbName. OXT|"},
    {"R", "PEPTIDE", "ARG/ARN", "N=C(N)NCCC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. "
          "NH2:1.pdbName. CZ :2.pdbName. NH1:3.pdbName. NE :4.pdbName. CD "
          ":5.pdbName. CG :6.pdbName. CB :7.pdbName. CA :8.pdbName. N  "
          ":9.pdbName. H  :10.pdbName. C  :11.pdbName. O  :12.pdbName. OXT|"},
    {"S", "PEPTIDE", "SER/SRO",  "O=C([C@H](CO)N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C  "
          ":2.pdbName. CA :3.pdbName. CB :4.pdbName. OG :5.pdbName. N  "
          ":6.pdbName. H  :7.pdbName. OXT|"},
    {"T", "PEPTIDE", "THR/THO", "C[C@@H](O)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. "
          "CG2:1.pdbName. CB :2.pdbName. OG1:3.pdbName. CA :4.pdbName. N  "
          ":5.pdbName. H  :6.pdbName. C  :7.pdbName. O  :8.pdbName. OXT|"},
    {"U", "PEPTIDE", "SEC", "O=C([C@H](C[SeH])N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. "
          "C  :2.pdbName. CA :3.pdbName. CB :4.pdbName.SE  :5.pdbName. N  "
          ":6.pdbName. H  :7.pdbName. OXT|"},
    {"V", "PEPTIDE", "VAL", "CC(C)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CG1:1.pdbName. "
          "CB :2.pdbName. CG2:3.pdbName. CA :4.pdbName. N  :5.pdbName. H  "
          ":6.pdbName. C  :7.pdbName. O  :8.pdbName. OXT|"},
    {"W", "PEPTIDE", "TRP",  "O=C([C@H](Cc1c[nH]c2ccccc12)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
          ":5.pdbName. CD1:6.pdbName. NE1:7.pdbName. CE2:8.pdbName. "
          "CZ2:9.pdbName. CH2:10.pdbName. CZ3:11.pdbName. CE3:12.pdbName. "
          "CD2:13.pdbName. N  :14.pdbName. H  :15.pdbName. OXT|"},
    {"Y", "PEPTIDE", "TYR/TYO:",
     "O=C([C@H](Cc1ccc(O)cc1)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
     "CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. OH :9.pdbName. "
     "CE2:10.pdbName. CD2:11.pdbName. N  :12.pdbName. H  :13.pdbName. OXT|"}
};

// PDB three-letter code to single-letter symbol mapping
// const std::unordered_map<std::string, std::string> pdb_to_symbol({
//     {"ALA", "A"}, // Alanine
//     {"ARG", "R"}, // Arginine
//     {"ASH", "D"}, // Protonated Aspartic
//     {"ASN", "N"}, // Asparagine
//     {"ASP", "D"}, // Aspartic acid
//     {"CYS", "C"}, // Cysteine
//     {"GLN", "Q"}, // Glutamine
//     {"GLU", "E"}, // Glutamic
//     {"GLY", "G"}, // Glycine
//     {"HIS", "H"}, // Histidine
//     {"ILE", "I"}, // Isoleucine
//     {"LEU", "L"}, // Leucine
//     {"LYS", "K"}, // Lysine
//     {"MET", "M"}, // Methionine
//     {"PHE", "F"}, // Phenylalanine
//     {"PRO", "P"}, // Proline
//     {"SER", "S"}, // Serine
//     {"THR", "T"}, // Threonine
//     {"TRP", "W"}, // Tryptophan
//     {"TYR", "Y"}, // Tyrosine
//     {"VAL", "V"}, // Valine
//     {"ARN", "R"}, // Neutral-Arginine
//     {"GLH", "E"}, // Protonated Glutamic
//     {"HID", "H"}, // Histidine (protonated at delta N)
//     {"HIE", "H"}, // Histidine (protonated at epsilon N)
//     {"HIP", "H"}, // Histidine (protonated at both N)
//     {"HSD", "H"}, // Histidine (protonated at delta N, CHARMM name)
//     {"HSE", "H"}, // Histidine (protonated at epsilon N, CHARMM name)
//     {"HSP", "H"}, // Histidine (protonated at both N, CHARMM name)
//     {"LYN", "K"}, // Protonated Lysine
//     {"SRO", "S"}, // Ionized Serine
//     {"THO", "T"}, // Ionized Threonine
//     {"TYO", "Y"}, // Ionized Tyrosine
//     {"SEC", "U"}, // Selenocysteine
//     {"PYL", "O"}, // Pyrrolysine
//     {"XXX", "X"} // Unknown
// });

// Single-letter symbol to canonical PDB code mapping
// const std::unordered_map<std::string, std::string> symbol_to_pdb({
//     {"A", "ALA"}, // Alanine
//     {"R", "ARG"}, // Arginine
//     {"N", "ASN"}, // Asparagine
//     {"D", "ASP"}, // Aspartic acid
//     {"C", "CYS"}, // Cysteine
//     {"Q", "GLN"}, // Glutamine
//     {"E", "GLU"}, // Glutamic acid
//     {"G", "GLY"}, // Glycine
//     {"H", "HIS"}, // Histidine
//     {"I", "ILE"}, // Isoleucine
//     {"L", "LEU"}, // Leucine
//     {"K", "LYS"}, // Lysine
//     {"M", "MET"}, // Methionine
//     {"F", "PHE"}, // Phenylalanine
//     {"P", "PRO"}, // Proline
//     {"S", "SER"}, // Serine
//     {"T", "THR"}, // Threonine
//     {"W", "TRP"}, // Tryptophan
//     {"Y", "TYR"}, // Tyrosine
//     {"V", "VAL"}, // Valine
//     {"U", "SEC"}, // Selenocysteine
//     {"O", "PYL"}  // Pyrrolysine
// });

} // anonymous namespace

// MonomerLibrary implementation
// MonomerLibrary::MonomerLibrary (bool loadBuiltins) {
//     if (loadBuiltins) {
//         //loadBuiltinDefinitions();

//         auto &globalLib = MonomerLibrary::getGlobalLibrary();


//         for (auto globalTemplate : globalLib) {
//             this->addTemplate(globalTemplate);   // NOT owned by this lib, but still owned by the global lib
//         }
//     }
// }

// MonomerLibrary::MonomerLibrary ([[maybe_unused]] std::string_view database_path) {
//     // TODO: Load from database/JSON file
//     //loadBuiltinDefinitions();

//     auto &globalLib = MonomerLibrary::getGlobalLibrary();


//     for (auto globalTemplate : globalLib) {
//         this->addTemplate(globalTemplate);   // NOT owned by this lib, but still owned by the global lib
//     }
// }

MonomerLibrary::~MonomerLibrary() = default;

void MonomerLibrary::loadBuiltinDefinitions() {
    // Load all built-in peptide monomers. default PEPTIDE to begin
    for (const auto& [symbol, monomerClass, pdbCodes, data] : builtin_monomer_data) {


        addMonomerFromSmiles(data, symbol, monomerClass, pdbCodes);
        //addMonomerFromSmiles(data, symbol, "PEPTIDE", symbol_to_pdb.count(symbol) ? symbol_to_pdb.at(symbol) : "");
                                          
            
        // MonomerEntry entry{
        //     symbol,
        //     data,
        //     std::shared_ptr<ROMol>(SmilesToMol(data, 0, false)),
        //     "PEPTIDE",
        //     symbol_to_pdb.count(symbol) ? symbol_to_pdb.at(symbol) : ""
        // };
        // d_monomers[key] = std::move(entry);

        // Add to PDB code index

        // MonomerKey key{symbol, "PEPTIDE"};

        // if (symbol_to_pdb.count(symbol)) {
        //     d_pdbCodeIndex[symbol_to_pdb.at(symbol)] = key;
        // }
    }

    // // Add reverse mappings from PDB codes (for variants like HID, HIE, etc.)
    // for (const auto& [pdb_code, symbol] : pdb_to_symbol) {
    //     if (d_pdbCodeIndex.find(pdb_code) == d_pdbCodeIndex.end()) {
    //         // Map variant PDB codes to the canonical symbol
    //         d_pdbCodeIndex[pdb_code] = {symbol, "PEPTIDE"};
    //     }
    // }
 
}

std::optional<std::string> MonomerLibrary::getMonomerData(
    const std::string& monomer_id,
    const std::string& monomer_class) const {
    // MonomerKey key{monomer_id, monomer_class};
    // auto it = d_monomers.find(key);
    // if (it != d_monomers.end()) {
    //     return it->second.original_data;
    // }
    // return std::nullopt;

     auto templateMol = this->getMACROMolTemplateLib().getTemplate(monomer_class, monomer_id);
     std::string originalData;
      
    if (templateMol->getPropIfPresent("OriginalData", originalData)) {
        return originalData;
    }
    return std::nullopt;

}

MonomerLibrary::monomer_info_t
MonomerLibrary::getMonomerInfo(const std::string& symbol, const std::string monomer_class) const {
    // auto it = d_pdbCodeIndex.find(pdb_code);
    // if (it != d_pdbCodeIndex.end()) {
    //     // const auto& [symbol, monomer_class] = it->second;
        // auto monomer_it = d_monomers.find(it->second);
        // if (monomer_it != d_monomers.end()) {
        //     return std::make_tuple(
        //         symbol,
        //         monomer_it->second.original_data,
        //         monomer_class
        //     );
        // }

        // auto symbol = it->second.first;
        // auto monomer_class = it->second.second;
        
        std::string originalData= "";
        std::vector<std::string> templateNames;
        auto templateMol = this->getMACROMolTemplateLib().getTemplate(monomer_class, symbol);
        if(templateMol ) {

            templateMol->getPropIfPresent("OriginalData", originalData);

            if(templateMol->getPropIfPresent(common_properties::templateNames, templateNames)) {
                return std::make_tuple(
                    templateNames[0],
                    originalData,
                    monomer_class
                );
            }
        }
    return std::nullopt;
    
}

std::optional<std::string>
MonomerLibrary::getPdbCode(const std::string& symbol,
                           const std::string& monomer_class) const {
    // MonomerKey key{symbol, monomer_class};
    // auto it = d_monomers.find(key);
    // if (it != d_monomers.end() && !it->second.pdb_code.empty()) {
    //     return it->second.pdb_code;
    // }

    auto pdbCodes = getPdbCodes(symbol, monomer_class);
    
    if (pdbCodes && pdbCodes->size() > 0) {
        return pdbCodes->at(0);  // just the first one
    }
    
    return std::nullopt;
}


std::optional<std::vector<std::string>> 
MonomerLibrary::getPdbCodes(const std::string& symbol,
                           const std::string& monomer_class) const {
  
    auto templateMol = this->getMACROMolTemplateLib().getTemplate(monomer_class, symbol);
    std::vector<std::string >pdbCodes;
    std::string pdbCodesStr;
    if (templateMol->getPropIfPresent("PDBCodes", pdbCodesStr)){
        boost::algorithm::split(pdbCodes, pdbCodesStr, boost::algorithm::is_any_of("/"));
        return pdbCodes;  // just the first one
    }
    return std::nullopt;
}

//std::shared_ptr<ROMol> MonomerLibrary::getMonomer(
RDKit::MACROMolTemplate *MonomerLibrary::getMonomer(
    const std::string& monomer_id,
    const std::string& monomer_class) const {

    // MonomerKey key{monomer_id, monomer_class};
    // auto it = d_monomers.find(key);
    // if (it != d_monomers.end()) {
    //     return it->second.mol;
    // }
    // return nullptr;

    if (!this->hasMonomer(monomer_id, monomer_class )) {
        return nullptr;
    }

    return this->d_macroMolTemplateLib.getTemplate(monomer_class, monomer_id);
    
}


void MonomerLibrary::addMonomer(std::unique_ptr<RWMol> &templateMol,
                                const std::string& data,
                                const std::string& symbol,
                                const std::string& monomer_class,
                                const std::string& pdb_codes) {
                            
    std::vector<std::pair<std::string, std::string>> otherAttrs;
    std::pair<std::string, std::string>  pairToAdd = {"OriginalData", data};
    otherAttrs.emplace_back(pairToAdd);
    if (!pdb_codes.empty()) {
        pairToAdd = {"PDBCodes", pdb_codes};
        otherAttrs.emplace_back(pairToAdd);
    }

    //  mols defining monomers may be in rgroup form like
    // *N[C@H](C(=O)O)S* |$_R1;;;;;;;_R3$| or use atom map numbers like
    // [*:1]N[C@H](C(=O)O)S[*:3]. Translate the RGroup to atom map
    // numbers
    for (auto atom : templateMol->atoms()) {
        std::string rgroup_label;
        if (atom->getPropIfPresent(RDKit::common_properties::atomLabel,
                                    rgroup_label) &&
            rgroup_label.find("_R") == 0) {
            auto rgroup_num = std::stoi(rgroup_label.substr(2));
            atom->setAtomMapNum(rgroup_num);
            atom->clearProp(RDKit::common_properties::atomLabel);
        }
    }

    //the definition mols can also have Rgroups as the leaving atom, or can have the attachment point
    // of can have the H or OH as the leaving group with map numbers
    // or have the map numbers on the atom to attach to (no R-groups or *)
    // in smiles parlance, this is either 
    // [*:1]N[C@H](C(=O)O)S[*:3]
    // or
    // [NH:1[C@H](C(=O)O)[S:3]


    // add the main Sgroup for the template - this is all atoms in the data

    std::vector<std::unique_ptr<SubstanceGroup>> newSgroups;
    const std::string typ = "SUP";
    newSgroups.emplace_back(new SubstanceGroup((ROMol *)templateMol.get(), typ));
    auto mainSgroup = newSgroups.back().get();
    mainSgroup->setProp("LABEL", symbol);
    mainSgroup->setProp("CLASS", monomer_class);
    mainSgroup->clearAttachPoints(); 


    // go through the mol and find the attachments points.  They are the atoms with atom map nums.
    // add the connection point and leaving group for each one

    templateMol->updatePropertyCache(false); 
    auto atomCount = templateMol->getNumAtoms();
    for (unsigned int atomIdx = 0 ; atomIdx < atomCount ; ++atomIdx) {
        auto atom = templateMol->getAtomWithIdx(atomIdx);
        unsigned int map_num;
        if (atom->getPropIfPresent(RDKit::common_properties::molAtomMapNumber, map_num)) {
            // if the atom is an R-group (* atom), Change it to be the leaving group atom (H or OH).
            // the one and only atom it is attached to is the  attachment point atom

            std::string leavingGroupId = "R" + std::to_string(map_num);
            atom->clearProp(RDKit::common_properties::molAtomMapNumber);
            
            newSgroups.emplace_back(new SubstanceGroup((ROMol *)templateMol.get(), typ));
            auto lvgSgroup = newSgroups.back().get();
            lvgSgroup->setProp("LABEL", leavingGroupId);
            lvgSgroup->setProp("CLASS", "LGRP");

            // see if this atom IS the leaving group.  It is if it in a R-atom (atomic num 0) or it is an H or OH
            auto atomicNumber  = atom->getAtomicNum();
            if (atomicNumber == 0 ||  atomicNumber  == 1 || (atomicNumber == 8 &&  atom->getDegree() == 1)) {
                // find the attached atom

                unsigned int otherAtomIdx = UINT_MAX;
                for (auto bnd : templateMol->atomBonds(atom)) {
                    // Should be exactly one iteration -- an Rgroup point can
                    // only have one neighbor
                    if (otherAtomIdx != UINT_MAX) {
                        throw std::runtime_error(
                            "Invalid attachment point at index " +
                            std::to_string(atom->getIdx()));
                    }
                    auto otherAtomIdx = bnd->getOtherAtomIdx(atom->getIdx());

                    std::string leavingGroupName;
                    auto otherAtomicNumber = templateMol->getAtomWithIdx(otherAtomIdx)->getAtomicNum();
                    if (atomicNumber == 1 // already an H
                        || otherAtomicNumber == 8 /*oxygen*/ || otherAtomicNumber == 7 /*nitrogen*/ || otherAtomicNumber == 16 /*sulfur*/ ) {
                        atom->setAtomicNum(1);  // change to H atom
                        leavingGroupName = "H";
                    } else {
                        atom->setAtomicNum(8);  // change to oxygen atom
                        leavingGroupName = "OH";
                    }
                
                    mainSgroup->addAttachPoint(otherAtomIdx, atom->getIdx(),  leavingGroupId);
                    
                    mainSgroup->addParentAtomWithIdx(otherAtomIdx);  // add to the main Sgroup

                    // add a SUP substance group for the leaving group

                    lvgSgroup->addAtomWithIdx(atom->getIdx());
                    lvgSgroup->addBondWithIdx(bnd->getIdx());
                    mainSgroup->addBondWithIdx(bnd->getIdx());  // add bond to the main Sgroup
                }
            } else {

                // one new atom  and one new SAP will be added for each attachment point

                std::string leavingGroupName;
                auto oldAtomicNumber = atom->getAtomicNum();
                
                unsigned int newAtomIdx;

                if (oldAtomicNumber == 8 /*oxygen*/ || oldAtomicNumber == 7 /*nitrogen*/ || oldAtomicNumber == 16 /*sulfur*/ ) {
                    newAtomIdx = templateMol->addAtom(1);  // add a hydrogen
                    std::string leavingGroupName = "H";
                } else {
                    newAtomIdx = templateMol->addAtom(8);  // add Oxygen atom
                    std::string leavingGroupName = "OH";
                }

                unsigned int newBondIdx = templateMol->addBond(atomIdx, newAtomIdx, Bond::BondType::SINGLE);

                mainSgroup->addAttachPoint(atom->getIdx(), newAtomIdx,  leavingGroupId);
                mainSgroup->addAtomWithIdx(atom->getIdx());  // add to the main Sgroup

                // add a SUP substance group for each attachment point

                lvgSgroup->addAtomWithIdx(newAtomIdx);
                lvgSgroup->addBondWithIdx(newBondIdx);

                mainSgroup->addBondWithIdx(newBondIdx);  // add bond to the main Sgroup
            }
        } else {
            mainSgroup->addAtomWithIdx(atomIdx);  // add to the main Sgroup
        }
    }

    if (newSgroups.size() > 0) {
        for (auto &sg : newSgroups) {
          addSubstanceGroup(*templateMol, *sg.get());
        }
      }

    std::vector<std::string> templateNames;
    templateNames.push_back(symbol);

    if (!pdb_codes.empty()) {
       
        std::vector<std::string> pdbList;
        boost::algorithm::split(pdbList, pdb_codes, boost::algorithm::is_any_of("/"));
        for (auto pdbCode : pdbList) {
            templateNames.push_back(pdbCode);
        }
    
    }
    auto newTemplate = std::unique_ptr<MACROMolTemplate>(new MACROMolTemplate(templateMol, monomer_class, templateNames, otherAttrs));
    d_macroMolTemplateLib.addTemplate(newTemplate);


    // MonomerEntry entry{
    //     symbol,
    //     data,
    //     std::shared_ptr<ROMol>(SmilesToMol(data, 0, false)),
    //     "PEPTIDE",
    //     symbol_to_pdb.count(symbol) ? symbol_to_pdb.at(symbol) : ""
    // };
    // d_monomers[key] = std::move(entry);

    // Add to PDB code index

    // MonomerKey key{symbol, monomer_class};

    // if (symbol_to_pdb.count(symbol)) {
    //     d_pdbCodeIndex[symbol_to_pdb.at(symbol)] = key;
    // }


    //d_monomers[key] = std::move(entry);

    // if (!pdb_codes.empty()) {
    //     d_pdbCodeIndex[pdb_code] = key;
    // }
}


void MonomerLibrary::addMonomerFromSmiles(const std::string& smiles,                                        
                                          const std::string& symbol,
                                          const std::string& monomer_class,
                                          const std::string& pdb_codes) {
    //MonomerKey key{symbol, monomer_class};
    // MonomerEntry entry{
    //     symbol,
    //     smiles,
    //     std::shared_ptr<ROMol>(SmilesToMol(smiles, 0, false)),
    //     monomer_class,
    //     pdb_code
    // };
    auto templateMol(std::unique_ptr<RDKit::RWMol>(SmilesToMol(smiles, 0, false)));
   
    addMonomer(templateMol, smiles, symbol, monomer_class, pdb_codes);
}

void MonomerLibrary::addMonomerFromSDF(const std::string& sdf_data,
                                       const std::string& symbol,
                                       const std::string& monomer_class,
                                       const std::string& pdb_codes) {
    // MonomerKey key{symbol, monomer_class};
    // MonomerEntry entry{
    //     symbol,
    //     sdf_data,
    //     std::shared_ptr<ROMol>(MolBlockToMol(sdf_data, false, false)),
    //     monomer_class,
    //     pdb_code
    // };
    // d_monomers[key] = std::move(entry);

    // if (!pdb_code.empty()) {
    //     d_pdbCodeIndex[pdb_code] = key;
    // }

    auto templateMol(std::unique_ptr<RDKit::RWMol>(MolBlockToMol(sdf_data, false, false)));
    addMonomer(templateMol, sdf_data, symbol, monomer_class, pdb_codes);
}

// void MonomerLibrary::addMonomer(std::shared_ptr<ROMol> mol,
//                                 const std::string& symbol,
//                                 const std::string& monomer_class,
//                                 const std::string& pdb_code) {
//     MonomerKey key{symbol, monomer_class};
//     MonomerEntry entry{
//         symbol,
//         "",  // no original data
//         std::move(mol),
//         monomer_class,
//         pdb_code
//     };
//     d_monomers[key] = std::move(entry);

//     if (!pdb_code.empty()) {
//         d_pdbCodeIndex[pdb_code] = key;
//     }
// }

bool MonomerLibrary::hasMonomer(const std::string& symbol,
                                const std::string& monomer_class) const {
    return d_macroMolTemplateLib.libContains(monomer_class, symbol);
}

// Static methods for global library management

MonomerLibrary *MonomerLibrary::getGlobalLibrary(){
    std::call_once(s_globalLibraryOnce, []() {
        s_globalMACROMOLTemplateLib = std::make_unique<MACROMolTemplateLib>();

        s_globalLibrary = std::make_unique<MonomerLibrary>(*(s_globalMACROMOLTemplateLib.get()));
        s_globalLibrary->loadBuiltinDefinitions();

    });

    return s_globalLibrary.get();

}

} // namespace RDKit
