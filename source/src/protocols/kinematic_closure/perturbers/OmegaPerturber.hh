// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_kinematic_closure_perturbers_OmegaPerturber_HH
#define INCLUDED_protocols_kinematic_closure_perturbers_OmegaPerturber_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/perturbers/OmegaPerturber.fwd.hh>

namespace protocols {
namespace kinematic_closure {
namespace perturbers {

class OmegaPerturber : public Perturber {

public:

	/// @copydoc Perturber::get_name()
	std::string get_name() const override { return "OmegaPerturber"; }

	/// @copydoc Perturber::perturb_subset()
	void perturb_subset(
		Pose const & pose,
		IndexList const & residues,
		ClosureProblemOP problem) override;

};

}
}
}

#endif
