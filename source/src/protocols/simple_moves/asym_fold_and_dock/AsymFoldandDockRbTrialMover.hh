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
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_simple_moves_asym_fold_and_dock_AsymFoldandDockRbTrialMover_hh
#define INCLUDED_protocols_simple_moves_asym_fold_and_dock_AsymFoldandDockRbTrialMover_hh

// Unit headers
#include <protocols/simple_moves/asym_fold_and_dock/AsymFoldandDockRbTrialMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>


// Utility Headers

namespace protocols {
namespace simple_moves {
namespace asym_fold_and_dock {
///////////////////////////////////////////////////////////////////////////////
class AsymFoldandDockRbTrialMover : public moves::Mover
{
public:

	// default constructor
	AsymFoldandDockRbTrialMover();

	AsymFoldandDockRbTrialMover( core::scoring::ScoreFunctionCOP scorefxn );

	AsymFoldandDockRbTrialMover( core::scoring::ScoreFunctionCOP scorefxn, bool smooth_move );

	AsymFoldandDockRbTrialMover(
		core::scoring::ScoreFunctionCOP scorefxn,
		bool smooth_move,
		core::Real rot_mag,
		core::Real trans_mag
	);

	~AsymFoldandDockRbTrialMover() override{}

	void apply( core::pose::Pose & pose ) override;
	//void find_new_jump_residue( core::pose::Pose & pose );
	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	std::string get_name() const override;

	static
	std::string mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::scoring::ScoreFunctionCOP scorefxn_;
	bool smooth_move_;
	core::Real rot_mag_;
	core::Real trans_mag_;
	core::Size rigid_body_cycles_;
	bool mc_filter_;
};

} // asym_fold_and_dock
}
} // rosetta

#endif
