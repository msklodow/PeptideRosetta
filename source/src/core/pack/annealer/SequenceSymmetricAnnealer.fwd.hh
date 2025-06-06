// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/annealer/SequenceSymmetricAnnealer.fwd.hh
/// @author Jack Maguire

#ifndef INCLUDED_core_pack_annealer_SequenceSymmetricAnnealer_fwd_hh
#define INCLUDED_core_pack_annealer_SequenceSymmetricAnnealer_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace annealer {

class SequenceSymmetricAnnealer;

typedef utility::pointer::shared_ptr< SequenceSymmetricAnnealer > SequenceSymmetricAnnealerOP;

}//end namespace annealer
}//end namespace pack
}//end namespace core


#endif

