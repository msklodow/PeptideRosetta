// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SetReturningPackRotamersMover.hh
/// @brief A PackRotamers mover which returns a set of packed poses for when ndruns/nloop is in use.
/// @author Ron Jacak

#ifndef INCLUDED_protocols_minimization_packing_SetReturningPackRotamersMover_hh
#define INCLUDED_protocols_minimization_packing_SetReturningPackRotamersMover_hh

// Unit headers
#include <protocols/minimization_packing/SetReturningPackRotamersMover.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh> // needed b/c we are extending from PackRotamersMover

// Project headers
#include <core/types.hh>
#include <core/pack/task/PackerTask.fwd.hh> // for PT COP
#include <core/scoring/ScoreFunction.fwd.hh> // for SF COP

#include <core/pose/Pose.hh>

#include <utility/vector1.hh>




namespace protocols {
namespace minimization_packing {


class SetReturningPackRotamersMover : public protocols::minimization_packing::PackRotamersMover {

public:
	SetReturningPackRotamersMover( core::Size ndruns );
	// custom constructor
	SetReturningPackRotamersMover( core::scoring::ScoreFunctionCOP scorefxn, core::pack::task::PackerTaskCOP task, core::Size ndruns );

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;
	void get_repacked_poses( utility::vector1< core::pose::Pose > & v );
	void output_repacked_poses( std::string filename_prefix );

private:
	utility::vector1< core::pose::Pose > repacked_poses_;
	core::Size ndruns_;

};

} // moves
} // protocols


#endif
