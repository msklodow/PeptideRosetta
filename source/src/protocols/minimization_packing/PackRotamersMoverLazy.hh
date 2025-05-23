// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Hetu Kamisetty

#ifndef INCLUDED_protocols_minimization_packing_PackRotamersMoverLazy_hh
#define INCLUDED_protocols_minimization_packing_PackRotamersMoverLazy_hh

// Unit headers
#include <protocols/minimization_packing/PackRotamersMoverLazy.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>

//#ifdef __clang__
//#endif
#include <utility/vector0.hh>


namespace protocols {
namespace minimization_packing {

/// @brief a mover that packs the side-chains using a rotamer library in a lazy fashion.
// specifically. if given an ig, it will use it instead of trying to make its own. this allows
// user to call multiple times with the same ig, but different arguments (eg. rot_to_pack)
class PackRotamersMoverLazy: public protocols::minimization_packing::PackRotamersMover{

public:
	/// @brief default constructor

	PackRotamersMoverLazy(
		ScoreFunctionCOP scorefxn,
		PackerTaskCOP task = nullptr,
		core::Size nloop = 1
	);
	PackRotamersMoverLazy();

	/// @brief constructor with typename
	// Undefinede, commenting out to fix PyRosetta build  PackRotamersMoverLazy( std::string const & );

	// destructor (important for properly forward-declaring smart-pointer members)
	~PackRotamersMoverLazy() override;

	// methods
	virtual void call_setup(Pose & pose);
	virtual void apply_to_rotpack( Pose & pose,  utility::vector0< int > rot_to_pack);
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &
	) override;

private:

};

// note: it is better to create new files, instead of adding additional classes here

} // moves
} // protocols

#endif
