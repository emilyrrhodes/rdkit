//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#pragma once

#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <GraphMol/MACROMol.h>

#include <RDGeneral/export.h>

namespace RDKit {

class ROMol;

//! Represents a monomer definition in the library
struct RDKIT_MONOMERMOL_EXPORT MonomerEntry {
    std::string symbol;           // e.g., "A" for Alanine
    std::string original_data;    // Original definition (SMILES, SDF, etc.)
    std::shared_ptr<ROMol> mol;   // Parsed molecule
    std::string monomer_class;    // "PEPTIDE", "RNA", etc.
    std::string pdb_code;         // e.g., "ALA" (optional)
};

//! Library of monomer definitions for use with MonomerMol
/*!
    MonomerLibrary provides access to monomer definitions. The monomer items
    can be owned by an external lib (global libraary) or owned by this library

    Access is typically via MonomerMol::getMonomerLibrary().
*/
class RDKIT_MONOMERMOL_EXPORT MonomerLibrary{
 public:
    //! Return type for getMonomerInfo: {symbol, original_data, monomer_class}
    using monomer_info_t = std::optional<std::tuple<std::string, std::string, std::string>>;


    //! Constructor
    /*!
    */
    explicit MonomerLibrary() 
            : d_localMacroMolTemplateLib(std::unique_ptr<RDKit::MACROMolTemplateLib>(new RDKit::MACROMolTemplateLib())),
            d_macroMolTemplateLib(*(d_localMacroMolTemplateLib.get())) {
    }

    //! Constructor with path to database/json file (future use)
    //explicit MonomerLibrary(std::string_view database_path);
    explicit MonomerLibrary(RDKit::MACROMolTemplateLib &libInit) : d_macroMolTemplateLib(libInit) {};
   
    ~MonomerLibrary();

    // Non-copyable but movable
    MonomerLibrary(const MonomerLibrary&) = delete;
    MonomerLibrary& operator=(const MonomerLibrary&) = delete;
    MonomerLibrary(MonomerLibrary&&) = default;
    MonomerLibrary& operator=(MonomerLibrary&&) = default;

    void copyMonomerLib(const MonomerLibrary &libToCopy){
        d_macroMolTemplateLib.copyTemplateLib(libToCopy.getMACROMolTemplateLib());
    }

    // --- Query operations ---

    //! Get original data (SMILES/SDF/etc) for a monomer by its symbol and class
    [[nodiscard]] std::optional<std::string> getMonomerData(
        const std::string& monomer_id,
        const std::string& monomer_class) const;

    //! Get full monomer info from a three-letter PDB code
    /*!
        \param pdb_code The three-letter PDB residue code (e.g., "ALA")
        \return tuple of {symbol, original_data, monomer_class} or nullopt if not found
    */
    [[nodiscard]] monomer_info_t
    getMonomerInfo(const std::string& pdb_code, const std::string monomer_class) const;

    //! Get PDB code for a monomer symbol
    [[nodiscard]] std::optional<std::string>
    getPdbCode(const std::string& symbol, const std::string& monomer_class) const;

    [[nodiscard]] std::optional<std::vector<std::string>>
    getPdbCodes(const std::string& symbol, const std::string& monomer_class) const;

    MACROMolTemplateLib &getMACROMolTemplateLib() {
        return d_macroMolTemplateLib;
    }
    const MACROMolTemplateLib &getMACROMolTemplateLib() const {
        return d_macroMolTemplateLib;
    }

    //! Get parsed molecule for a monomer
    //[[nodiscard]] std::shared_ptr<ROMol> getMonomer(
    [[nodiscard]] MACROMolTemplate *getMonomer(
        const std::string& monomer_id,
        const std::string& monomer_class) const;

    // --- Mutation operations (for instance libraries) ---

    //! Add a monomer from a SMILES string
    void addMonomerFromSmiles(const std::string& smiles,
                              const std::string& symbol,
                              const std::string& monomer_class,
                              const std::string& pdb_code = "");

    //! Add a monomer from an SDF/molblock string
    void addMonomerFromSDF(const std::string& sdf_data,
                           const std::string& symbol,
                           const std::string& monomer_class,
                           const std::string& pdb_code = "");

    //! Add a monomer with a pre-parsed molecule
    void addMonomer(std::unique_ptr<RWMol> &mol,
                    const std::string& data,
                    const std::string& symbol,
                    const std::string& monomer_class,
                    const std::string& pdb_code = "");

    //! Check if a monomer exists
    [[nodiscard]] bool hasMonomer(const std::string& symbol,
                                  const std::string& monomer_class) const;

    // --- Global library configuration ---

    //! Get the global singleton instance

    static MonomerLibrary *getGlobalLibrary();

 private:
    //! Load built-in monomer definitions
    //! TODO: this may load in a preset json or sqlite DB
    void loadBuiltinDefinitions();

    //! Static state for global mode
    //static bool s_useGlobalLibrary
    static std::unique_ptr<MonomerLibrary> s_globalLibrary;
    static std::unique_ptr<MACROMolTemplateLib> s_globalMACROMOLTemplateLib;
    static std::once_flag s_globalLibraryOnce;


    std::unique_ptr<RDKit::MACROMolTemplateLib> d_localMacroMolTemplateLib; // used for standalone MACROMolTemplateLib
    RDKit::MACROMolTemplateLib &d_macroMolTemplateLib;
    

};

} // namespace RDKit
