// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/Universal_Helper.cc
/// @brief A derived class of the CrosslinkerMoverHelper base class, used to set up
/// [...] cross-linker.
/// @author Vikram K. Mulligan (vmullig@uw.edu), Maximilian Hoffmann (m.hoffmann200398@gmail.com)

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/Universal_Helper.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AA.hh>

// Protocols headers
#include <protocols/simple_moves/DeclareBond.hh>
#include <protocols/cyclic_peptide/CreateDistanceConstraint.hh>
#include <protocols/cyclic_peptide/CreateTorsionConstraint.hh>
#include <protocols/simple_moves/ModifyVariantTypeMover.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <utility/exit.hh>

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <map>
#include <vector>


static basic::Tracer TR( "protocols.cyclic_peptide.crosslinker.Universal_Helper" );

namespace protocols {
    namespace cyclic_peptide {
        namespace crosslinker {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief constructor
            Universal_Helper::Universal_Helper(){
                read_params();
            }

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
////////////////////////////////////////////////////////////////////////////////
            Universal_Helper::Universal_Helper(Universal_Helper const &/*src*/ ) = default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer
/// members)
////////////////////////////////////////////////////////////////////////////////
            Universal_Helper::~Universal_Helper() = default;

//////////////////////
/// Public Methods ///
//////////////////////

/// @brief Provide an opportunity to provide a citation for this crosslinker type.
/// @details The base class implementation does nothing.  This override indicates that this helper is
/// unpublished.
            void
            Universal_Helper::provide_citation_info(
                    basic::citation_manager::CitationCollectionList &citations
            ) const {
                citations.add(
                        utility::pointer::make_shared<basic::citation_manager::UnpublishedModuleInfo>(
                                "Universal_Helper", basic::citation_manager::CitedModuleType::CrosslinkerMoverHelper,
                                "Maximilian Hoffmann / Vikram K. Mulligan",
                                "Institute for Drug Discovery, University of Leipzig / Systems Biology Group, Center for Computational Biology, Flatiron Institute",
                                "m.hoffmann200398@gmail.com / vmulligan@flatironinstitute.org",
                                "Added more universal support for different linker types. So far just covalent bonds are implemented."
                        )
                );
            }

/// @brief Given a pose and a selection of residues, add the linker,
/// align it crudely to the selected residues, and set up covalent bonds.
            void
            Universal_Helper::add_linker_asymmetric(
                    core::pose::Pose &pose,
                    core::select::residue_selector::ResidueSubset const & selection
            ) const {
                utility::vector1< core::Size > res_indices;
                CrosslinkerMoverHelper::get_sidechain_indices( selection, res_indices );

//              assert number residues
                runtime_assert_string_msg(int(res_indices.size()) == num_res,
                                          "Error in protocols::cyclic_peptide::crosslinker::Universal_Helper::add_linker_asymmetric(): The numbers of residues in the provided XML and universal_linker_params.txt files are not matching.");
//              assert sidechain types
                assert_sidechain_type(pose, res_indices);

//              assert if for every residue a connection chain was provided
                runtime_assert_string_msg(int(residue_connection_chains.size()) == num_res,
                                          "Error in protocols::cyclic_peptide::crosslinker::Universal_Helper::add_linker_asymmetric(): The numbers of residues in the provided XML does not match the number of provided 'residue_connection_chains' in universal_linker_params.txt.");

//              add sidechain conjugation type
                protocols::simple_moves::ModifyVariantTypeMover mut1;
                core::select::residue_selector::ResidueIndexSelectorOP index_selector(
                        utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >()
                );
                for ( auto const index : res_indices ) {
                    index_selector->append_index( index );
                }
                mut1.set_residue_selector( index_selector );
                mut1.set_additional_type_to_add( "SIDECHAIN_CONJUGATION" ); // TODO other conjugation type, specified by provided .txt params file?
                mut1.apply(pose);

                core::chemical::ResidueTypeSetCOP standard_residues( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
                core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( standard_residues->name_map(linker_name) ) );
                core::pose::Pose linker_pose;
                linker_pose.append_residue_by_jump(*new_rsd, 1);

//              create alignment map
                std::map <core::id::AtomID, core::id::AtomID> alignment_atoms;
                int iter_count = 1;
//              TODO add individual atom names for alignment map, specified py provided .txt param file
                for ( auto const index : res_indices ) {
                    alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index("V" + std::to_string(iter_count)), 1) ] = core::id::AtomID( pose.residue_type(index).atom_index(residue_connection_chains[iter_count-1][0]), index);
                    alignment_atoms[ core::id::AtomID( new_rsd->type().atom_index(residue_connection_chains[iter_count-1][1]), 1) ] = core::id::AtomID( pose.residue_type(index).atom_index("V1"), index);
                    iter_count++;
                }

                //Align the linker pose to the original pose:
                core::scoring::superimpose_pose( linker_pose, pose, alignment_atoms );

                //Merge the poses:
                pose.append_residue_by_jump( linker_pose.residue(1), res_indices[1] );
                core::Size const linker_res( pose.total_residue() );

                //Declare covalent bonds:
                add_linker_bonds_asymmetric( pose, res_indices, linker_res );

                core::chemical::DeleteAtom delete_virtuals_res("V1");
                core::chemical::MutableResidueType mutable_residue_linker(pose.residue_type(linker_res));
                iter_count = 1;
                for ( auto const index : res_indices ) {
                    core::chemical::MutableResidueType mutable_residue(pose.residue_type(index));
                    delete_virtuals_res.apply(mutable_residue);

                    core::chemical::DeleteAtom delete_virtuals_linker("V" + std::to_string(iter_count));
                    delete_virtuals_linker.apply(mutable_residue_linker);
                    iter_count++;

                }
            }

/// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
/// @details Called by add_linker_asymmetric().  Version for asymmetric poses.
            void
            Universal_Helper::add_linker_bonds_asymmetric(
                    core::pose::Pose &pose,
                    utility::vector1< core::Size > const &res_indices,
                    core::Size const linker_index
            ) const {
//              TODO add individual atom names to set bonds, specified py provided .txt param file
                int iter_count = 1;
                for ( auto const index : res_indices ) {
                    protocols::simple_moves::DeclareBond bond;
                    bond.set( linker_index, residue_connection_chains[iter_count-1][1], index, residue_connection_chains[iter_count-1][0], false );
                    bond.apply(pose);
                    iter_count++;
                }
            }


/// @brief no support for symmetric poses
            void
            Universal_Helper::add_linker_symmetric(
                    core::pose::Pose &pose,
                    core::select::residue_selector::ResidueSubset const & selection
            ) const {
                static_cast<void>(pose);
                static_cast<void>(selection);

                throw std::runtime_error("No support for symmetric poses!");
            }

/// @brief no support for symmetric poses
            void
            Universal_Helper::add_linker_bonds_symmetric(
                    core::pose::Pose &pose,
                    core::Size const res1,
                    core::Size const linker_index1,
                    core::Size const linker_index2
            ) const {

                static_cast<void>(pose);
                static_cast<void>(res1);
                static_cast<void>(linker_index1);
                static_cast<void>(linker_index2);

                throw std::runtime_error("No support for symmetric poses!");
            }


/// @brief Given a selection of residues that have already been connected the linker,
/// add constraints for the linker.
            void
            Universal_Helper::add_linker_constraints_asymmetric(
                    core::pose::Pose &pose,
                    core::select::residue_selector::ResidueSubset const & selection
            ) const {

                //Get indices of residues
                utility::vector1< core::Size > res_indices;
                CrosslinkerMoverHelper::get_sidechain_indices( selection, res_indices );
                core::Size const linker_index( get_linker_index_asymmetric( pose, res_indices) );

//              assert number residues
                runtime_assert_string_msg(int(res_indices.size()) == num_res,
                                          "Error in protocols::cyclic_peptide::crosslinker::Universal_Helper::add_linker_asymmetric(): The numbers of residues in the provided XML and universal_linker_params.txt files are not matching.");
                //Set up distance constraints:
                { //Begin scope
                    std::string const dist_cst_string("HARMONIC 0.0 0.01");
                    protocols::cyclic_peptide::CreateDistanceConstraint dist_csts;
                    utility::vector1< core::Size > res1(num_res * 2), res2(num_res * 2);
////                    len(res1) = 2n
////                    res1[1...n] contains residue indices 1...n
////                    res1[n+1...2n] contains residue indices 1...n
////
////                    len(atom1) = 2n
////                    atom1[1...n] contains names of the proteins 1...n residue atoms connected to the linker
////                    atom1[n+1...2n] contains name of virtual atom of protein residues, which is for all "V1"
////
////                    len(res2) = 2n
////                    res2[1...2n] contains linker index
////
////                    len(atom2) = 2n
////                    atom2[1...n] contains names of linker atoms 1...n connected to the protein
////                    atom1[n+1...2n] contains name of the linkers virtual atoms connected to the protein ("V1"..."Vn")

//                  TODO add individual atom names specified py provided .txt params file
                    utility::vector1< std::string > atom1(num_res * 2), atom2(num_res * 2);
                    utility::vector1< std::string > cst_fxn(num_res * 2);
                    for (int index = 1; index <= num_res; ++index) {
                        res1[index] = res_indices[index];
                        res1[index + num_res] = res_indices[index];
                        atom1[index] = residue_connection_chains[index-1][0];
                        atom1[index + num_res] = "V1";

                        res2[index] = linker_index;
                        res2[index + num_res] = linker_index;
                        atom2[index] = "V" + std::to_string(index);
                        atom2[index + num_res] = residue_connection_chains[index-1][1];

                        cst_fxn[index] = dist_cst_string;
                        cst_fxn[index + num_res] = dist_cst_string;

                    }

                    if ( TR.Debug.visible() ) {
                        TR.Debug << "R1\tA1\tR2\tA2\tFUNC\n";
                        for (int i = 1; i <= num_res * 2; ++i) {
                            TR.Debug << res1[i] << "\t";
                            TR.Debug << atom1[i] << "\t";
                            TR.Debug << res2[i] << "\t";
                            TR.Debug << atom2[i] << "\t";
                            TR.Debug << cst_fxn[i] << "\n";
                        }
                        TR.Debug << std::endl;
                    }
                    dist_csts.set( res1, atom1, res2, atom2, cst_fxn );
                    dist_csts.apply(pose);
                } //End of scope

//              TODO add individual atom names, torsion constraints to be calculated, specified by provided .txt params file
                //Set up torsion constraints:
                { //Begin scope
                    protocols::cyclic_peptide::CreateTorsionConstraint tors_csts;
                    utility::vector1 < core::Size > res1(num_res * 3), res2(num_res * 3), res3(num_res * 3), res4(num_res * 3);
                    utility::vector1 < std::string > atom1(num_res * 3), atom2(num_res * 3), atom3(num_res * 3), atom4(num_res * 3);
                    utility::vector1 < std::string > cst_fxns(num_res * 3);

                    int i; int ii;
                    for (int index = 1; index <= num_res; ++index) {
                        //CM#index - prot_res_conn_atom - CB - CA   should be three-well potential:
                        res1[index] = linker_index; res2[index] = res_indices[index]; res3[index] = res_indices[index]; res4[index] = res_indices[index];
                        atom1[index] = residue_connection_chains[index-1][1]; atom2[index] = residue_connection_chains[index-1][0]; atom3[index] = "CB"; atom4[index] = "CA";
                        cst_fxns[index] = "AMBERPERIODIC 0 3 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)

                        //C#index*2 - C#index*2-1 - CM# - prot_res_conn_atom    should be two-well potential (above/below plane).
                        i = index + num_res;
                        res1[i] = linker_index; res2[i] = linker_index; res3[i] = linker_index; res4[i] = res_indices[index];
                        atom1[i] = residue_connection_chains[index-1][3]; atom2[i] = residue_connection_chains[index-1][2]; atom3[i] = residue_connection_chains[index-1][1]; atom4[i] = residue_connection_chains[index-1][0];
                        cst_fxns[i] = "AMBERPERIODIC 0 2 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)

                        //C#index*2-1 - CM#index - prot_res_conn_atom - CB    should be three-well potential:
                        ii = i + num_res;
                        res1[ii] = linker_index; res2[ii] = linker_index; res3[ii] = res_indices[index]; res4[ii] = res_indices[index];
                        atom1[ii] = residue_connection_chains[index-1][2]; atom2[ii] = residue_connection_chains[index-1][1]; atom3[ii] = residue_connection_chains[index-1][0]; atom4[ii] = "CB";
                        cst_fxns[ii] = "AMBERPERIODIC 0 3 2"; //AMBERPERIODIC has maximum at x=x0; params are x0 (offset), N (periodicity), K (amplitude)
                    }

                    if ( TR.Debug.visible() ) {
                        TR.Debug << "R1\tA1\tR2\tA2\tR3\tA3\tR4\tA4\tFUNC\n";
                        for (int c = 1; c <= num_res * 3; ++c) {
                            TR.Debug << res1[c] << "\t";
                            TR.Debug << atom1[c] << "\t";
                            TR.Debug << res2[c] << "\t";
                            TR.Debug << atom2[c] << "\t";
                            TR.Debug << res3[c] << "\t";
                            TR.Debug << atom3[c] << "\t";
                            TR.Debug << res4[c] << "\t";
                            TR.Debug << atom4[c] << "\t";
                            TR.Debug << cst_fxns[c] << "\n";
                        }
                        TR.Debug << std::endl;
                    }

                    tors_csts.set( res1, atom1, res2, atom2, res3, atom3, res4, atom4, cst_fxns );
                    tors_csts.apply( pose );
                } //End of scope

            }

/// @brief no support for symmetric poses
            void
            Universal_Helper::add_linker_constraints_symmetric(
                    core::pose::Pose &pose,
                    core::select::residue_selector::ResidueSubset const & selection,
                    bool const linker_was_added
            ) const {

                static_cast<void>(pose);
                static_cast<void>(selection);
                static_cast<void>(linker_was_added);

                throw std::runtime_error("No support for symmetric poses!");

            }

/// @brief Given indices of the residues that are already connected to a linker, get the index
/// of the linker residue.
/// @details Throws an error if the residues not all connected to the same linker residue.
            core::Size
            Universal_Helper::get_linker_index_asymmetric(
                    core::pose::Pose const &pose,
                    utility::vector1< core::Size > const & res_indices
            ) const {
                std::string const errmsg( "Error in protocols::cyclic_peptide::crosslinker::Universal_Helper::get_linker_index_asymmetric(): " );
                runtime_assert_string_msg( int(res_indices.size()) == num_res, errmsg + "The wrong number of residues was passed to this function.  A vector of exactly three residues is expected." );

                utility::vector1<core::Size> nconns(num_res);

                core::Size linker_index = 0;

                int iter_counter = 1;
                for (auto const index : res_indices) {
                    nconns[iter_counter] = pose.residue(index).n_possible_residue_connections();

                    if (iter_counter == 1) {
                        linker_index = pose.residue(index).residue_connection_partner(nconns[1]);
                    }

                    runtime_assert_string_msg( pose.residue(index).residue_connection_partner(nconns[iter_counter]) == linker_index,
                                               errmsg + "The latest residue (which is assumed to be the linker) connected to the  " + std::to_string(iter_counter) + ". selected protein residue is NOT the same as for the 1. selected protein residue!" );

                    iter_counter++;

                }
                return linker_index;
            }

/// @brief no support for symmetric poses
            void
            Universal_Helper::get_linker_indices_symmetric(
                    core::pose::Pose const &pose,
                    utility::vector1< core::Size > const & res_indices,
                    utility::vector1< core::Size > & linker_indices
            ) const {

                static_cast<void>(pose);
                static_cast<void>(res_indices);
                static_cast<void>(linker_indices);

                throw std::runtime_error("No support for symmetric poses!");

            }

/// @brief Given a pose with residues selected to be connected by the linker,
/// determine whether the residues are too far apart.
/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
            bool
            Universal_Helper::filter_by_sidechain_distance_asymmetric(
                    core::pose::Pose const &pose,
                    core::select::residue_selector::ResidueSubset const & selection,
                    core::Real const &filter_multiplier
            ) const {

                core::Real const hardcoded_cutoff( 12.0 * filter_multiplier );
                core::Real const hardcoded_cutoff_sq( hardcoded_cutoff*hardcoded_cutoff );

                //Get indices of residues
                utility::vector1< core::Size > res_indices;
                CrosslinkerMoverHelper::get_sidechain_indices( selection, res_indices );

                assert_sidechain_type(pose, res_indices);
                runtime_assert_string_msg(int(res_indices.size()) == num_res,
                                          "Error in protocols::cyclic_peptide::crosslinker::Universal_Helper::add_linker_asymmetric(): The numbers of residues in the provided XML and universal_linker_params.txt files are not matching.");

                core::Real sqdist;
                std::string D;
//              loop over residue sqdist (1,2), ..., (n-1,n)
                for (int i = 1; i < num_res; ++i) {
                    sqdist = pose.residue(res_indices[i]).xyz("CB").distance_squared(pose.residue(res_indices[i + 1]).xyz("CB") );
                    if ( TR.Debug.visible() ) {
                        D = "D: ";
                        D.insert(1, std::to_string(i));
                        TR.Debug << D << sqrt(sqdist) << " MAX: " << hardcoded_cutoff << std::endl;
                    }
                    if (sqdist > hardcoded_cutoff_sq ) return true;
                }
//              calculate residue sqdist (1,n)
                sqdist = pose.residue(res_indices[1]).xyz("CB").distance_squared(pose.residue(res_indices[num_res]).xyz("CB") );
                if ( TR.Debug.visible() ) {
                    D = "D: ";
                    D.insert(1, std::to_string(num_res));
                    TR.Debug << D << sqrt(sqdist) << " MAX: " << hardcoded_cutoff << std::endl;
                }
                if ( sqdist > hardcoded_cutoff_sq ) return true;

                return false;
            }

/// @brief no support of symmetric poses
            bool
            Universal_Helper::filter_by_sidechain_distance_symmetric(
                    core::pose::Pose const &pose,
                    core::select::residue_selector::ResidueSubset const & selection,
                    core::Real const &filter_multiplier
            ) const {

                static_cast<void>(pose);
                static_cast<void>(selection);
                static_cast<void>(filter_multiplier);

                throw std::runtime_error("No support for symmetric poses!");

                return false;

            }

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
            bool
            Universal_Helper::filter_by_constraints_energy_asymmetric(
                    core::pose::Pose const &pose,
                    core::select::residue_selector::ResidueSubset const & selection,
                    core::Real const &filter_multiplier
            ) const {
                return filter_by_constraints_energy( pose, selection, false, false, filter_multiplier );
            }

/// @brief No support for symmetric poses
            bool
            Universal_Helper::filter_by_constraints_energy_symmetric(
                    core::pose::Pose const &pose,
                    core::select::residue_selector::ResidueSubset const & selection,
                    bool const linker_was_added,
                    core::Real const &filter_multiplier
            ) const {

                throw std::runtime_error("No support for symmetric poses!");

                return filter_by_constraints_energy( pose, selection, true, linker_was_added, filter_multiplier );
            }

/********************************************
PRIVATE FUNCTIONS
*********************************************/
/// @brief Read in the params specified in the provided universal_linker_params.txt file.
            void
            Universal_Helper::read_params() {
                // specify if argument can have multiple values
                std::map<std::string , int> multi_ins_map;
                multi_ins_map["linker_name"] = false;
                multi_ins_map["num_res"] = false;
                multi_ins_map["AA_res"] = true;
                multi_ins_map["residue_connection_chains"] = true;

                std::vector<std::string> reads;

                char delimiter = '=';
                char val_separator = ',';

                // Specify the file name
                std::string file_name = "universal_linker_params.txt";

                // Create an input file stream
                std::ifstream infile(file_name);

                // Check if the file was successfully opened
                if (!infile.is_open()) {
                    throw std::runtime_error("Error: Could not open universal_linker_params.txt");
                }
                std::string input;
                while (std::getline(infile, input)) {
                    // Process the line, e.g., print it to the console
                    if (input.find(delimiter) == std::string::npos) {
                        continue;
                    }
                    input.erase(std::remove_if(input.begin(), input.end(), isspace), input.end());
                    if (input.at(0) == '#') {
                        continue;
                    }
                    std::string key = input.substr(0, input.find(delimiter));
                    input.erase(0, key.length() + 1);
                    if (input.empty()) {
                        throw std::runtime_error("No input argument given for key '" + key + "'");
                    }
                    if (input.back() == val_separator) {
                        input.pop_back();
                    }
                    if (multi_ins_map.find(key) != multi_ins_map.end()) {
                        int n = int(std::count(input.begin(), input.end(), val_separator)) + 1;
                        if ((n > 1) && (multi_ins_map[key] == false)) {
                            throw std::runtime_error("Too many input arguments for key '" + key + "'. Can take only one argument.");
                        } else if (n > 1) {
                            // read multiple values
                            int pos;
                            while ((pos = int(input.find(val_separator))) != int(std::string::npos)) {
                                reads.push_back(input.substr(0, pos));
                                input.erase(0, pos + 1);
                            }
                            reads.push_back(input);
                        } else if (n == 1) {
                            // read single value
                            reads.push_back(input);
                        }
                        // assign value
                        if (key == "linker_name") {
                            linker_name = reads.back();
                        } else if (key == "num_res") {
                            num_res = std::stoi(reads.back());
                        } else if (key == "AA_res") {
                            for (const std::string& value: reads) {
                                AA_res.push_back(value);
                            }
                        } else if (key == "residue_connection_chains") {
                            residue_connection_chains.push_back(reads);
                        }

                        reads.clear();
                    } else {
                        throw std::runtime_error("Unknown key '" + key + "'");
                    }
                }
                // Close the file when done
                infile.close();
            }


            /// @brief Assert if selected residues match allowed types.
            void
            Universal_Helper::assert_sidechain_type(
                    const core::pose::Pose &pose,
                    const utility::vector1<core::Size> &res_indices
            ) const {
                if (AA_res[0] != "ALLOW_ALL") {
                    for ( auto const index : res_indices ) {
                        std::string const &basetype_name(pose.residue_type(index).base_name());
                        runtime_assert_string_msg(std::find(AA_res.begin(), AA_res.end(), basetype_name) != AA_res.end(),
                                                  "Error in protocols::cyclic_peptide::crosslinker::Universal_Helper::add_linker_asymmetric(): Selected residue Type does not match one of the allowed types.");
                    }
                }
            }

        } //crosslinker
    } //protocols
} //cyclic_peptide
