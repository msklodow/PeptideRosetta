// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/P_AA_pp_Energy.hh
/// @brief  Probability of observing an amino acid, given its phi/psi energy method forward declaration
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_energy_methods_P_AA_pp_Energy_fwd_hh
#define INCLUDED_core_energy_methods_P_AA_pp_Energy_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace energy_methods {


class P_AA_pp_Energy;

typedef utility::pointer::shared_ptr< P_AA_pp_Energy > P_AA_pp_EnergyOP;

} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
