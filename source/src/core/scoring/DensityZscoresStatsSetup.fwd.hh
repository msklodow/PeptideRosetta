// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/DensityZscoresStatsSetup.fwd.hh
/// @brief Forward declarations for a helper class that stores the setup information for Density Z-scores
/// @details
/// @author Gabriella Reggiano

#ifndef INCLUDED_core_scoring_DensityZscoresStatsSetup_fwd_hh
#define INCLUDED_core_scoring_DensityZscoresStatsSetup_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

class DensityZscoresStatsSetup;

typedef utility::pointer::shared_ptr< DensityZscoresStatsSetup > DensityZscoresStatsSetupOP;
typedef utility::pointer::shared_ptr< DensityZscoresStatsSetup const > DensityZscoresStatsSetupCOP;

} // scoring
} // core


#endif
