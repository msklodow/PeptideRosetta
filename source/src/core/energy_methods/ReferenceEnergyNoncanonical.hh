// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/ReferenceEnergyNoncanonical.hh
/// @brief  Reference energy method declaration
/// @author Fang-Chieh Chou (fcchou@stanford.edu)



#ifndef INCLUDED_core_energy_methods_ReferenceEnergyNoncanonical_hh
#define INCLUDED_core_energy_methods_ReferenceEnergyNoncanonical_hh

// Unit headers
#include <core/energy_methods/ReferenceEnergyNoncanonical.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace energy_methods {


class ReferenceEnergyNoncanonical : public core::scoring::methods::ContextIndependentOneBodyEnergy
{
public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy parent;

public:


	ReferenceEnergyNoncanonical();


	ReferenceEnergyNoncanonical( utility::vector1< Real > const & weight_list );

	///
	~ReferenceEnergyNoncanonical() override;

	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const override;


	virtual
	Real
	eval_dof_derivative(
		id::DOF_ID const & dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights
	) const;

	/// @brief DunbrackEnergy is context independent; indicates that no
	/// context graphs are required
	void indicate_required_context_graphs( utility::vector1< bool > & ) const override;

	void init_res_list();

private:
	utility::vector1< Real > weights_;
	utility::vector1< std::string > res_list_;
	core::Size version() const override;

};

} // scoring
} // core


#endif
