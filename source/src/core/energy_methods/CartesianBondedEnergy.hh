// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/CartesianBondedEnergy.hh
/// @brief
/// @author Frank DiMaio
/// @author Andrew Leaver-Fay optimized the code a bit

#ifndef INCLUDED_core_energy_methods_CartesianBondedEnergy_hh
#define INCLUDED_core_energy_methods_CartesianBondedEnergy_hh

// Unit headers
#include <core/energy_methods/CartesianBondedEnergy.fwd.hh>
#include <core/scoring/methods/CartBondedParameters.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/RestypeDestructionEvent.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// boost
#include <boost/unordered_map.hpp>
#include <boost/tuple/tuple.hpp>

// C++ headers
#include <iosfwd>
#include <map>

#include <utility/VirtualBase.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>

#if defined MULTI_THREADED

#include <mutex>
#include <utility/thread/ReadWriteMutex.hh>

#endif
//#include <map>


typedef boost::tuples::tuple< std::string, std::string, std::string, std::string, std::string > atm_name_quad;
typedef boost::tuples::tuple< std::string, std::string, std::string, std::string > atm_name_triple;
typedef boost::tuples::tuple< std::string, std::string, std::string > atm_name_pair;
typedef boost::tuples::tuple< std::string, std::string > atm_name_single;

namespace boost {
namespace tuples {

std::size_t hash_value(atm_name_quad const& e);
std::size_t hash_value(atm_name_triple const& e);
std::size_t hash_value(atm_name_pair const& e);
std::size_t hash_value(atm_name_single const& e);
bool operator==(atm_name_quad const& a,atm_name_quad const& b);
bool operator==(atm_name_triple const& a,atm_name_triple const& b);
bool operator==(atm_name_pair const& a,atm_name_pair const& b);
bool operator==(atm_name_single const& a,atm_name_single const& b);

}
}

namespace core {
namespace energy_methods {


typedef utility::vector1< std::pair< atm_name_quad, core::scoring::methods::CartBondedParametersCOP > > torsionparam_vector;

class ResidueCartBondedParameters : public utility::VirtualBase {
public:
	typedef utility::fixedsizearray1< Size, 2 > Size2;
	typedef utility::fixedsizearray1< Size, 3 > Size3;
	typedef utility::fixedsizearray1< Size, 4 > Size4;
	typedef std::pair< Size2, core::scoring::methods::CartBondedParametersCOP > length_parameter;
	typedef std::pair< Size3, core::scoring::methods::CartBondedParametersCOP > angle_parameter;
	typedef std::pair< Size4, core::scoring::methods::CartBondedParametersCOP > torsion_parameter;

public:
	ResidueCartBondedParameters();
	~ResidueCartBondedParameters() override;

	void add_length_parameter(  Size2 atom_inds, core::scoring::methods::CartBondedParametersCOP );
	void add_angle_parameter(   Size3 atom_inds, core::scoring::methods::CartBondedParametersCOP );
	void add_torsion_parameter( Size4 atom_inds, core::scoring::methods::CartBondedParametersCOP );
	void add_improper_parameter( Size4 atom_inds, core::scoring::methods::CartBondedParametersCOP );
	void add_bbdep_length_parameter(  Size2 atom_inds, core::scoring::methods::CartBondedParametersCOP );
	void add_bbdep_angle_parameter(   Size3 atom_inds, core::scoring::methods::CartBondedParametersCOP );
	void add_lower_connect_angle_params( Size3 atom_inds, core::scoring::methods::CartBondedParametersCOP );
	void add_upper_connect_angle_params( Size3 atom_inds, core::scoring::methods::CartBondedParametersCOP );

	void bb_N_index( Size index );
	void bb_CA_index( Size index );
	void bb_C_index( Size index );
	void bb_O_index( Size index );
	void bb_H_index( Size index );
	void pro_CD_index( Size index );

	void ca_cprev_n_h_interres_improper_params( core::scoring::methods::CartBondedParametersCOP );
	void oprev_cprev_n_h_interres_improper_params( core::scoring::methods::CartBondedParametersCOP );
	void ca_nnext_c_o_interres_improper_params( core::scoring::methods::CartBondedParametersCOP );
	void pro_cd_cprev_n_ca_interres_improper_params( core::scoring::methods::CartBondedParametersCOP );
	void cprev_n_bond_length_params( core::scoring::methods::CartBondedParametersCOP );

	utility::vector1< length_parameter > const &
	length_parameters() const {
		return length_params_;
	}

	utility::vector1< angle_parameter > const &
	angle_parameters() const {
		return angle_params_;
	}

	utility::vector1< torsion_parameter > const &
	torsion_parameters() const {
		return torsion_params_;
	}

	/// @brief Exactly the same as proper torsion parameters, but parceled out
	/// into their own section so that debugging information can be given for
	/// these torsions in particular.
	utility::vector1< torsion_parameter > const &
	improper_parameters() const {
		return improper_params_;
	}

	/// @brief just the list of length parameters that are dependent on phi and psi; used for calculating dE/dphi and dE/dpsi
	utility::vector1< length_parameter > const &
	bbdep_length_parameters() const {
		return bbdep_length_params_;
	}

	/// @brief just the list of angle parameters that are dependent on phi and psi; used for calculating dE/dphi and dE/dpsi
	utility::vector1< angle_parameter > const &
	bbdep_angle_parameters() const {
		return bbdep_angle_params_;
	}

	utility::vector1< angle_parameter > const &
	lower_connect_angle_params() const {
		return lower_connect_angle_params_;
	}

	utility::vector1< angle_parameter > const &
	upper_connect_angle_params() const {
		return upper_connect_angle_params_;
	}

	Size bb_N_index()  const { return bb_N_index_;  }
	Size bb_CA_index() const { return bb_CA_index_; }
	Size bb_C_index()  const { return bb_C_index_;  }
	Size bb_O_index()  const { return bb_O_index_;  }
	Size bb_H_index()  const { return bb_H_index_;  }
	Size pro_CD_index( )  const { return pro_CD_index_;  }

	core::scoring::methods::CartBondedParametersCOP
	ca_cprev_n_h_interres_improper_params() const {
		return ca_cprev_n_h_interres_improper_params_;
	}

	core::scoring::methods::CartBondedParametersCOP
	oprev_cprev_n_h_interres_improper_params() const {
		return oprev_cprev_n_h_interres_improper_params_;
	}

	core::scoring::methods::CartBondedParametersCOP
	ca_nnext_c_o_interres_improper_params() const {
		return ca_nnext_c_o_interres_improper_params_;
	}

	core::scoring::methods::CartBondedParametersCOP
	pro_cd_cprev_n_ca_interres_improper_params() const {
		return pro_cd_cprev_n_ca_interres_improper_params_;
	}

	core::scoring::methods::CartBondedParametersCOP
	cprev_n_bond_length_params() const {
		return cprev_n_bond_length_params_;
	}


private:
	utility::vector1< length_parameter  > length_params_;
	utility::vector1< angle_parameter   > angle_params_;
	utility::vector1< torsion_parameter > torsion_params_;
	utility::vector1< torsion_parameter > improper_params_;
	utility::vector1< length_parameter  > bbdep_length_params_;
	utility::vector1< angle_parameter   > bbdep_angle_params_;

	/// For amino acids only: if they have a lower connection,
	/// then what are the angle parameters for Cprev-at1-at2 for all
	/// atoms at2 bonded to lower-connect-atom at1?
	utility::vector1< angle_parameter   > lower_connect_angle_params_;

	/// For amino acids only: if they have an upper connection,
	/// then what are the angle parameters for Nnext-at1-at2 for all
	/// atoms at2 bonded to upper-connect-atom at1?
	utility::vector1< angle_parameter   > upper_connect_angle_params_;


	Size bb_N_index_;
	Size bb_CA_index_;
	Size bb_C_index_;
	Size bb_O_index_;
	Size bb_H_index_;
	Size pro_CD_index_;

	core::scoring::methods::CartBondedParametersCOP ca_cprev_n_h_interres_improper_params_;
	core::scoring::methods::CartBondedParametersCOP oprev_cprev_n_h_interres_improper_params_;
	core::scoring::methods::CartBondedParametersCOP ca_nnext_c_o_interres_improper_params_;
	core::scoring::methods::CartBondedParametersCOP pro_cd_cprev_n_ca_interres_improper_params_;
	core::scoring::methods::CartBondedParametersCOP cprev_n_bond_length_params_;

};

////////////////////
//fpd
//  Database stores all ideal parameters
class IdealParametersDatabase  : public utility::VirtualBase {
public:
	// The per-residue maping for the torsions
	typedef boost::unordered_map< atm_name_quad, core::scoring::methods::CartBondedParametersOP > TorsionsIndepSubmap;

	IdealParametersDatabase(Real k_len, Real k_ang, Real k_tors, Real k_tors_prot, Real k_tors_improper);

	// If you enable copy and assignment operators, you need to account for
	// Detaching old/attaching new destruction observer to the ResidueType
	IdealParametersDatabase( IdealParametersDatabase const & ) = delete;
	IdealParametersDatabase & operator=( IdealParametersDatabase const & ) = delete;

	~IdealParametersDatabase() override;

	core::scoring::methods::CartBondedParametersCOP
	lookup_improper(
		core::chemical::ResidueType const & rsd_type,
		std::string const & atm1_name,
		std::string const & atm2_name,
		std::string const & atm3_name,
		std::string const & atm4_name
	);

	/// @brief Get the improper torsion constraints for the particular residue.
	TorsionsIndepSubmap
	generate_impropers_map_res(
		core::chemical::ResidueType const & restype
	);

	// needs both names (for keying off databases) and indices (for building those not found from ideal)
	core::scoring::methods::CartBondedParametersCOP
	lookup_angle(
		core::chemical::ResidueType const & rsd_type,
		bool pre_proline,
		std::string const & atm1_name,
		std::string const & atm2_name,
		std::string const & atm3_name,
		int atm1idx,
		int atm2idx,
		int atm3idx
	);

	// needs both names (for keying off databases) and indices (for building those not found from ideal)
	core::scoring::methods::CartBondedParametersCOP
	lookup_length(
		core::chemical::ResidueType const & rsd_type,
		bool pre_proline,
		std::string const & atm1_name,
		std::string const & atm2_name,
		int atm1idx,
		int atm2idx
	);

	// old-style interface to database
	void
	lookup_torsion_legacy( core::chemical::ResidueType const & restype,
		int atm1, int atm2, int atm3, int atm4, Real &Kphi, Real &phi0, Real &phi_step ) const;

	// old-style interface to database
	void
	lookup_angle_legacy( core::pose::Pose const & pose, core::conformation::Residue const & res,
		int atm1, int atm2, int atm3, Real &Ktheta, Real &d0) const;

	// old-style interface to database
	void
	lookup_length_legacy( core::pose::Pose const & pose, core::conformation::Residue const & res, int atm1, int atm2, Real &Kd, Real &d0 ) const;

	bool bbdep_bond_params() const { return bbdep_bond_params_; }
	bool bbdep_bond_devs() const { return bbdep_bond_devs_; }

	/// @brief Return a list of all the bond lengths, bond angles, and bond torsions
	/// for a single residue type.  This list is constructed lazily as required.
	ResidueCartBondedParameters const &
	parameters_for_restype(
		core::chemical::ResidueType const & restype,
		bool prepro
	);

	// getters
	Real k_length() const { return k_length_; }
	Real k_angle() const { return k_angle_; }
	Real k_torsion() const { return k_torsion_; }
	Real k_torsion_proton() const { return k_torsion_proton_; }
	Real k_torsion_improper() const { return k_torsion_improper_; }

	void restype_destruction_observer( core::chemical::RestypeDestructionEvent const & event );

private:
	void init(Real k_len, Real k_ang, Real k_tors, Real k_tors_prot, Real k_tors_improper);

	// helper functions: find the ideal values by constructing from Rosetta's params file
	void
	lookup_bondangle_buildideal(
		core::chemical::ResidueType const & restype,
		int atm1,
		int atm2,
		int atm3,
		Real &Ktheta,
		Real &theta0
	) const;

	void
	lookup_bondlength_buildideal(
		core::chemical::ResidueType const & restype,
		int atm1,
		int atm2,
		Real &Kd,
		Real &d0
	) const;

	// read bb indep tables
	void read_length_database( std::string const & infile, bool const symmetrize_gly );
	void read_angle_database(std::string const & infile, bool const symmetrize_gly);
	void read_torsion_database(std::string const & infile, bool const symmetrize_gly);
	void read_improper_database(std::string const & infile, bool const symmetrize_gly);
	void add_impropers_from_stream(std::istream & instream);

	// another helper function: read backbone dependent db files
	void
	read_bbdep_table(
		std::string const &filename,
		boost::unordered_map< atm_name_single, core::scoring::methods::CartBondedParametersOP > &bondlengths,
		boost::unordered_map< atm_name_pair, core::scoring::methods::CartBondedParametersOP > &bondangles,
		std::string const &res,
		bool const symmetrize_table
	);

	/// @brief Generate cached parameters for the given ResidueType
	ResidueCartBondedParametersCOP
	create_parameters_for_restype(
		core::chemical::ResidueType const & restype,
		bool prepro
	);

private:
	/// @brief Symmetrize the glycine backbone-dependent table.
	/// @details Only called if the score::symmetric_gly_tables option is used.  Intended for design
	/// with glyceine in a mixed D/L context (in which there should be no preference for a left-handed
	/// conformation over a right).
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	//fd : make this private since it should only be called during initialization
	void symmetrize_tables( ObjexxFCL::FArray2D<core::Real> &table ) const;


	// TODO: Make the following statement actually true.
	// All of these values (except *_restype_data_) are set in the constructor and
	// are read-only, meaning they don't need to be mutex protected.

	// defaults (they should be rarely used as everything should be in the DB now)
	Real k_length_, k_angle_, k_torsion_, k_torsion_proton_, k_torsion_improper_;

	// options
	bool bbdep_bond_params_, bbdep_bond_devs_;

#if defined MULTI_THREADED
	// The following two mutexes are separate because in the generation of the
	// ResidueCartBondedParametersCOP, we access the data maps, and recursive mutex locking is hairy at best

	// mutex protecting shared access to bondlengths_indep_, bondangles_indep_
	// (These are the only data maps that are altered after initial loading.)
	utility::thread::ReadWriteMutex data_map_mutex_;
	// mutex protecting shared access to the residue-pointer indexed data
	utility::thread::ReadWriteMutex restype_db_mutex_;
#endif

	// backbone-independent parameters (keyed on atom names)
	boost::unordered_map< atm_name_pair, core::scoring::methods::CartBondedParametersOP > bondlengths_indep_;
	boost::unordered_map< atm_name_triple, core::scoring::methods::CartBondedParametersOP > bondangles_indep_;
	boost::unordered_multimap< atm_name_quad, core::scoring::methods::CartBondedParametersOP > torsions_indep_;
	boost::unordered_map< atm_name_quad,  core::scoring::methods::CartBondedParametersOP > impropers_indep_;

	// backbone-dependent parameter sets
	boost::unordered_map< atm_name_pair, core::scoring::methods::CartBondedParametersOP >
		bondangles_bbdep_def_, bondangles_bbdep_pro_, bondangles_bbdep_valile_, bondangles_bbdep_prepro_, bondangles_bbdep_gly_;
	boost::unordered_map< atm_name_single, core::scoring::methods::CartBondedParametersOP >
		bondlengths_bbdep_def_, bondlengths_bbdep_pro_, bondlengths_bbdep_valile_, bondlengths_bbdep_prepro_, bondlengths_bbdep_gly_;

	// per residue-type data
	// The raw pointers here are registered in the destruction observer for their respective ResidueType
	std::map< chemical::ResidueType const *, ResidueCartBondedParametersCOP > prepro_restype_data_;
	std::map< chemical::ResidueType const *, ResidueCartBondedParametersCOP > nonprepro_restype_data_;

};


///////////////////
///
/// the energy method
class CartesianBondedEnergy : public core::scoring::methods::ContextIndependentLRTwoBodyEnergy {
public:
	typedef core::scoring::methods::ContextIndependentLRTwoBodyEnergy  parent;

public:
	CartesianBondedEnergy() = delete; // Need options to initialize this structure

	CartesianBondedEnergy( core::scoring::methods::EnergyMethodOptions const & options );

	CartesianBondedEnergy( CartesianBondedEnergy const & src );

	~CartesianBondedEnergy() override;

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void
	setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & sfxn ) const override;

	/// @brief Idealize the virtual NV atom of every proline in the pose. This
	///prevents innacurate pro-close scores when switching between cartesian
	///and non-cartesian score functions.
	void
	idealize_proline_nvs(
		pose::Pose & pose
	) const;

	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const override;


	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &
	) const override;

	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		core::scoring::ResSingleMinimizationData const & res_data_cache,
		pose::Pose const & pose,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & atom_derivs
	) const override;

	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData const & min_data,
		pose::Pose const &,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const override;

	void
	eval_residue_pair_derivatives_sorted(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData const & min_data,
		pose::Pose const &,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const;

	// dof (bbdep) derivatives
	Real
	eval_intraresidue_dof_derivative(
		conformation::Residue const & rsd,
		core::scoring::ResSingleMinimizationData const & min_data,
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights
	) const override;

	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override { return true; }

	bool
	defines_intrares_dof_derivatives( pose::Pose const & ) const override { return true; }

	//fpd  use the new minimizer interface
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override { return false; }

	virtual
	Distance
	atomic_interaction_cutoff() const;

	void indicate_required_context_graphs( utility::vector1< bool > & ) const override;

	core::scoring::methods::LongRangeEnergyType
	long_range_type() const override;

private:

	//////////////////////////////
	/// Score evaluation methods
	//////////////////////////////

	/// @brief here, rsd1.seqpos() < rsd2.seqpos()
	void
	residue_pair_energy_sorted(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sf,
		core::scoring::EnergyMap & emap
	) const;

	/// Methods for intra-residue energies

	//@brief because of separate pre-proline distributions
	//    this function must be called from residue_pair_energy
	//    rather than defined in intrares energy
	void
	eval_singleres_energy(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & rsdparams,
		Real phi,
		Real psi,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;

	/// @brief helper function to handle intrares torsions
	void
	eval_singleres_torsion_energies(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;


	/// @brief helper function to handle intrares ring torsions & angles
	void
	eval_singleres_ring_energies(
		conformation::Residue const & rsd,
		//ResidueCartBondedParameters const & resparams,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;

	/// @brief helper function to handle intrares bond improper torsions
	void
	eval_singleres_improper_energies(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real const phi,
		Real const psi,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;

	/// @brief helper function to handle intrares bond angles
	void
	eval_singleres_angle_energies(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real const phi,
		Real const psi,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;

	/// @brief helper function to handle intrares bond lengths
	void
	eval_singleres_length_energies(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real const phi,
		Real const psi,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;

	/// Methods for inter-residue energies

	/// @brief Evaluate all the inter-residue components
	void
	eval_residue_pair_energies(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi1,
		Real psi1,
		Real phi2,
		Real psi2,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;

	/// @brief Evaluate all inter
	void
	eval_interresidue_angle_energies_two_from_rsd1(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi1,
		Real psi1,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;

	void
	eval_interresidue_angle_energies_two_from_rsd2(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi2,
		Real psi2,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;

	void
	eval_interresidue_bond_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi2,
		Real psi2,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;


	void
	eval_interresidue_improper_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;

	void
	eval_interresidue_ring_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;

	/////////////////////////////////
	/// Derivative evaluation methods
	/////////////////////////////////

	// single-residue derivatives

	/// @brief evaluate all intra-residue derivatives
	void
	eval_singleres_derivatives(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real phi,
		Real psi,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r_atom_derivs
	) const;

	/// @brief evaluate intra-residue ring derivatives
	void
	eval_singleres_ring_derivatives(
		conformation::Residue const & rsd,
		//ResidueCartBondedParameters const & resparams,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r_atom_derivs
	) const;


	/// @brief evaluate intra-residue improper torsion derivatives
	void
	eval_singleres_improper_derivatives(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real const phi,
		Real const psi,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r_atom_derivs
	) const;

	/// @brief evaluate intra-residue angle derivatives
	void
	eval_singleres_angle_derivatives(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real const phi,
		Real const psi,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r_atom_derivs
	) const;

	/// @brief evaluate intra-residue bond-length derivatives
	void
	eval_singleres_length_derivatives(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real const phi,
		Real const psi,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r_atom_derivs
	) const;


	/// @brief evaluate intra-residue torsion derivatives
	void
	eval_singleres_torsion_derivatives(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r_atom_derivs
	) const;

	// residue-pair derivatives

	/// @brief evaluate inter-residue angle derivatives where
	/// two of the atoms defining the angle are from rsd1
	void
	eval_interresidue_angle_derivs_two_from_rsd1(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi1,
		Real psi1,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const;

	/// @brief evaluate inter-residue angle derivatives where
	/// two of the atoms defining the angle are from rsd2
	void
	eval_interresidue_angle_derivs_two_from_rsd2(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi2,
		Real psi2,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const;

	/// @brief evaluate inter-residue bond-length derivatives
	void
	eval_interresidue_bond_length_derivs(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi2,
		Real psi2,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const;

	/// @brief evaluate inter-residue improper torsion derivatives
	void
	eval_interresidue_improper_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & res1params,
		ResidueCartBondedParameters const & res2params,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const;

	void
	eval_interresidue_ring_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const;

	////////////////////////////////////////////////////////////////////
	/// Common to evaluating the score for torsions, angles, and lengths
	////////////////////////////////////////////////////////////////////


	/// @brief Evaluate either the harmonic or linearized-harmonic energy
	/// given by either:
	/// score = 0.5 * K * (val-val0)^2
	/// or
	/// score = 0.5 * K * (val-val0)^2 if std::abs(val-val0) < 1
	///       = 0.5 * K * std::abs(val-val0) otherwise
	Real eval_score( Real val, Real K, Real val0 ) const;

	/// @brief Evaluate the derivative for a

private:

	/// @brief Accessor function which gives access to the appropriate IdealParametersDatabase
	static
	IdealParametersDatabaseOP
	get_db(Real k_len, Real k_ang, Real k_tors, Real k_tors_prot, Real k_tors_improper);

private:
	// the ideal parameter database used by this method.
	IdealParametersDatabaseOP db_;

	// option
	bool linear_bonded_potential_;

	core::Size version() const override;

	std::string pro_nv_;
	bool skip_cutpoints_; // skip evaluation if rsd1/rsd2 are defined as cutpoint

};

} // namespace energy_methods
} // namespace core


#endif // INCLUDED_core_energy_methods_CartesianBondedEnergy_HH
