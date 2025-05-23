// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/AddSapMathConstraintMover.fwd.hh
/// @brief Mover that adds the SapMathConstraint to the pose
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_AddSapMathConstraintMover_fwd_hh
#define INCLUDED_protocols_simple_moves_AddSapMathConstraintMover_fwd_hh
// Utility Headers

#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace simple_moves {

class AddSapMathConstraintMover;
typedef utility::pointer::shared_ptr< AddSapMathConstraintMover > AddSapMathConstraintMoverOP;
typedef utility::pointer::shared_ptr< AddSapMathConstraintMover const > AddSapMathConstraintMoverCOP;

} // simple_moves
} // protocols


#endif
