// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/WaterAdductIntraEnergy.hh
/// @brief  Energetic offset/cost for placing a water adduct on an amino or nucleic acid
/// @author Jim Havranek

// Unit headers
#include <core/energy_methods/WaterAdductIntraEnergy.hh>
#include <core/energy_methods/WaterAdductIntraEnergyCreator.hh>

// Package headers
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



/// @details This must return a fresh instance of the WaterAdductIntraEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
WaterAdductIntraEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< WaterAdductIntraEnergy >();
}

core::scoring::ScoreTypes
WaterAdductIntraEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( h2o_intra );
	return sts;
}


/// ctor
WaterAdductIntraEnergy::WaterAdductIntraEnergy() :
	parent( utility::pointer::make_shared< WaterAdductIntraEnergyCreator >() )
{}

/// clone
core::scoring::methods::EnergyMethodOP
WaterAdductIntraEnergy::clone() const
{
	return utility::pointer::make_shared< WaterAdductIntraEnergy >();
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////


void
WaterAdductIntraEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & ,
	core::scoring::EnergyMap & emap
) const
{
	// Sum over all waters
	for ( int atm = 1, atme = rsd.natoms() ; atm <= atme ; ++atm ) {
		if ( rsd.atom_type( atm ).is_h2o() ) emap[ core::scoring::h2o_intra ] += 1.0;
	}
}


Real
WaterAdductIntraEnergy::eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const &, //  tor_id
	pose::Pose const &, // pose
	core::scoring::ScoreFunction const &, //sfxn
	core::scoring::EnergyMap const & // weights
) const
{
	return 0.0;
}

/// @brief WaterAdductIntraEnergy is context independent; indicates that no
/// context graphs are required
void
WaterAdductIntraEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
WaterAdductIntraEnergy::version() const
{
	return 1; // Initial versioning
}


} // scoring
} // core

