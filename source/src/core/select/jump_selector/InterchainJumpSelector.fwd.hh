// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/InterchainJumpSelector.fwd.hh
/// @brief  Forward declaration of a class that combines JumpSelector logic
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_jump_selector_InterchainJumpSelector_FWD_HH
#define INCLUDED_core_select_jump_selector_InterchainJumpSelector_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace select {
namespace jump_selector {

class InterchainJumpSelector;

typedef utility::pointer::shared_ptr< InterchainJumpSelector > InterchainJumpSelectorOP;
typedef utility::pointer::shared_ptr< InterchainJumpSelector const > InterchainJumpSelectorCOP;

} //namespace jump_selector
} //namespace select
} //namespace core


#endif
