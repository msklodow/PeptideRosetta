// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rotamers/SingleBasicRotamerLibrary.hh
/// @brief  "Basic" RotamerLibrary class
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_pack_rotamers_SingleBasicRotamerLibrary_hh
#define INCLUDED_core_pack_rotamers_SingleBasicRotamerLibrary_hh

// Unit Headers
#include <core/pack/rotamers/SingleBasicRotamerLibrary.fwd.hh>

// Package Headers
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

// Project Headers
#include <core/conformation/Residue.fwd.hh>

namespace core {
namespace pack {
namespace rotamers {

/// @brief A simple Rotamer library, which serves as a default for ResidueTypes which don't have some other
/// more specific rotamer library.
/// @details In practice just diversifies the proton chi records.
class SingleBasicRotamerLibrary : public rotamers::SingleResidueRotamerLibrary
{
public:
	SingleBasicRotamerLibrary();

	~SingleBasicRotamerLibrary() override;

	/// @brief Adheres to the contract from SingleBasicRotamerLibrary
	void
	rotamer_energy_deriv(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::id::TorsionID const & tor_id,
		TorsionEnergy & tderiv
	) const override;

	/// @brief Adheres to the contract from SingleBasicRotamerLibrary
	Real
	rotamer_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		TorsionEnergy & energy
	) const override;

	std::set< id::PartialAtomID >
	atoms_w_dof_derivatives(
		conformation::Residue const & rsd,
		pose::Pose const & pose
	) const override;

	Real
	best_rotamer_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		bool curr_rotamer_only
	) const override;

	void
	assign_random_rotamer_with_bias(
		conformation::Residue const &,// rsd,
		pose::Pose const & /*pose*/,
		numeric::random::RandomGenerator &,// RG,
		dunbrack::ChiVector &,// new_chi_angles,
		bool //perturb_from_rotamer_center
	) const override;

	/// @brief Adheres to the contract from SingleBasicRotamerLibrary
	void
	fill_rotamer_vector(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		pack::task::PackerTask const & task,
		utility::graph::GraphCOP,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const& existing_residue,
		utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
		bool buried,
		rotamers::RotamerVector & rotamers
	) const override;

	/// @brief Adheres to the contract from SingleBasicRotamerLibrary
	void
	write_to_file( utility::io::ozstream &out ) const override;

private:

}; // SingleBasicRotamerLibrary


} // namespace rotamers
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_pack_rotamers_SingleBasicRotamerLibrary_hh
