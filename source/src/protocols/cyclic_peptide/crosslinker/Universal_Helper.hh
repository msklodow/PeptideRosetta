// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/TMA_Helper.hh
/// @brief A derived class of the CrosslinkerMoverHelper base class, used to set up
/// [...] cross-linker.
/// @author Vikram K. Mulligan (vmullig@uw.edu), MaxHoffmann ()

#ifndef INCLUDED_protocols_cyclic_peptide_crosslinker_Universal_Helper_hh
#define INCLUDED_protocols_cyclic_peptide_crosslinker_Universal_Helper_hh

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/Universal_Helper.fwd.hh>
#include <protocols/cyclic_peptide/crosslinker/CrosslinkerMoverHelper.hh>

// Protocol headers

// Core headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

// Basic/Utility headers
#include <core/types.hh>

#include <string>

namespace protocols {
    namespace cyclic_peptide {
        namespace crosslinker {

/// @brief A derived class of the CrosslinkerMoverHelper base class, used to set up
/// the trimesic acid (TMA) cross-linker.
            class Universal_Helper : public CrosslinkerMoverHelper {

            public: //Constructors

                /// @brief Default constructor
                Universal_Helper();

                /// @brief Copy constructor
                Universal_Helper(Universal_Helper const &src);

                /// @brief Destructor (important for properly forward-declaring smart-pointer members)
                ~Universal_Helper() override;


            public: // public methods

                /// @brief Provide an opportunity to provide a citation for this crosslinker type.
                /// @details The base class implementation does nothing.  This override indicates that this helper is
                /// published in Dang, Wu, Mulligan et al. 2017.
                void
                provide_citation_info(
                        basic::citation_manager::CitationCollectionList & citations
                ) const override;

                /// @brief Given a pose and a selection of exactly three residues, add the TBMB linker,
                /// align it crudely to the selected residues, and set up covalent bonds.
                void add_linker_asymmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const override;

                /// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
                /// @details Called by add_linker_asymmetric().  Version for asymmetric poses.
                void add_linker_bonds_asymmetric(core::pose::Pose &pose, utility::vector1< core::Size > const &res_indices, core::Size const linker_index ) const override;

                /// @brief Given a pose and a selection of exactly three residues, add the TBMB linker,
                /// align it crudely to the selected residues, and set up covalent bonds.
                /// @details Version for symmetric poses.
                void add_linker_symmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const override;

                /// @brief Given a pose and a linker, add bonds between the TBMB linker and the residues that coordinate the linker.
                /// @details Called by add_linker_symmetric().  Version for symmetric poses.
                void add_linker_bonds_symmetric(core::pose::Pose &pose, core::Size const res1, core::Size const linker_index1, core::Size const linker_index2 ) const override;

                /// @brief Given a selection of exactly three residues that have already been connected to a 1,3,5-tris(bromomethyl)benzene crosslinker,
                /// add constraints for the crosslinker.
                void add_linker_constraints_asymmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const override;

                /// @brief Given a selection of exactly three residues that have already been connected to a 1,3,5-tris(bromomethyl)benzene crosslinker,
                /// add constraints for the crosslinker.  This version is for symmetric poses.
                void add_linker_constraints_symmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, bool const linker_was_added) const override;

                /// @brief Given indices of three cysteine residues that are already linked to a TBMB, get the index
                /// of the TBMB residue.
                /// @details Throws an error if the three cysteines are not all linked to the same TBMB residue.
                core::Size get_linker_index_asymmetric( core::pose::Pose const &pose, utility::vector1< core::Size > const & res_indices ) const override;

                /// @brief Given indices of three cysteine residues that are already linked to pieces of a linker, get
                /// of the indices of the symmetric pieces of the linker.
                /// @details Throws an error if a residue is not linked to something.
                void get_linker_indices_symmetric( core::pose::Pose const &pose, utility::vector1< core::Size > const & res_indices, utility::vector1< core::Size > & linker_indices ) const override;

                /// @brief Given a pose with residues selected to be linked by a 1,3,5-tris(bromomethyl)benzene crosslinker,
                /// determine whether the residues are too far apart.
                /// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
                /// @note Higher values of the filter multiplier make it more permissive.
                bool filter_by_sidechain_distance_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const &filter_multiplier ) const override;

                /// @brief Given a pose with residues selected to be linked by a 1,3,5-tris(bromomethyl)benzene crosslinker,
                /// determine whether the residues are too far apart.  This version is for symmetric poses.
                /// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
                /// @note Higher values of the filter multiplier make it more permissive.
                bool filter_by_sidechain_distance_symmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const &filter_multiplier ) const override;

                /// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
                /// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
                /// @note Higher values of the filter multiplier make it more permissive.
                bool filter_by_constraints_energy_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const &filter_multiplier) const override;

                /// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  This version is for symmetric poses.
                /// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
                /// @note Higher values of the filter multiplier make it more permissive.
                bool filter_by_constraints_energy_symmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, bool const linker_was_added, core::Real const &filter_multiplier) const override;

                /// @brief Does this CrosslinkerMoverHelper add a residue for the linker?
                /// @details Yes, it does.
                bool helper_adds_linker_residue() const override { return true; }

            private: // private methods
                /// @brief Read in the params specified in the provided universal_linker_params.txt file.
                void read_params();

                /// @brief Assert if selected residues match allowed types.
                void assert_sidechain_type(const core::pose::Pose &pose, utility::vector1< core::Size > const &res_indices) const;

            private: // data
                std::string linker_name;
                int num_res{};
                std::vector<std::string> AA_res;
                std::vector<std::vector<std::string>> residue_connection_chains;
            };
        }
    }
}
#endif
