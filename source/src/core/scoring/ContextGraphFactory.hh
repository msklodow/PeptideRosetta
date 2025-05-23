// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ContextGraphFactory.hh
/// @brief  Context graph factory class declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_ContextGraphFactory_hh
#define INCLUDED_core_scoring_ContextGraphFactory_hh

// Unit Headers
#include <core/scoring/ContextGraphFactory.fwd.hh>

// PackageHeaders
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>

namespace core {
namespace scoring {

class ContextGraphFactory
{
public:
	static
	ContextGraphOP
	create_context_graph( ContextGraphType type );
};

} // scoring
} // core

#endif

