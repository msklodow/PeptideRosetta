// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/bin_transitions/PerturbByBins.cc
/// @brief  Headers for the PerturbByBins mover.  This mover takes a stretch of backbone and perturbs its mainchain torsions based
/// on the probabilities of transitions from one torsion bin to another.
/// @details Bin transitions are read from database files.  The algorithm is: set the first residue based on the probability of a residue
/// being in a bin.  Set subsequent residues based on the probability of a residue being in a bin given that the previous residue is in
/// a particular bin.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_bin_transitions_PerturbByBins_HH_
#define INCLUDED_protocols_simple_moves_bin_transitions_PerturbByBins_HH_

#include <protocols/simple_moves/bin_transitions/PerturbByBins.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Bin transition calculator headers:
#include <core/scoring/bin_transitions/BinTransitionCalculator.fwd.hh>

//parsing
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>



// Utility headers

// C++ headers

// Unit headers

namespace protocols {
namespace simple_moves {
namespace bin_transitions {

/// @brief A mover to set mainchain torsions by bin transition probabilities
///
class PerturbByBins : public protocols::moves::Mover
{
private:
	typedef protocols::moves::Mover moverclass;
public:
	/// @brief Default constructor.
	///
	PerturbByBins();

	/// @brief Copy constructor.
	///
	PerturbByBins( PerturbByBins const &src );

	/// @brief Destructor.
	///
	~PerturbByBins() override;


	/// @brief Clone -- i.e. create a new object copying this one and return an owning pointer to the copy.
	///
	protocols::moves::MoverOP clone() const override {
		return (utility::pointer::make_shared< protocols::simple_moves::bin_transitions::PerturbByBins >( *this ) );
	}

	/// @brief Get a new instance of this mover (NOT copying).
	///
	protocols::moves::MoverOP fresh_instance() const override {
		return utility::pointer::make_shared< PerturbByBins >();
	}

	/// @brief Apply the mover to a pose.
	///
	void apply( core::pose::Pose & pose ) override;

	/// @brief Parse XML for RosettaScripts.
	///
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	/******************************
	** SETTERS                  **
	******************************/

	/// @brief Set the bin transition probability file.
	/// @details Also, loads the object.
	void set_binfile_and_load( std::string const &name );

	/// @brief Set the residue ranges.  If set to (0,0), the
	/// start and end of the pose are used as the range bounds.
	void set_residue_range( core::Size const start, core::Size const end );

	/// @brief Set the number of repeats.
	/// @details A value of 1 means that a single residue in the range will randomly be selected and
	/// flipped to another bin (with probabilities based on its neighbours and the bin transition
	/// probabilities.  Higher values mean that this operation will be repeated.  If set to 0, defaults
	/// to 1.
	void set_repeats( core::Size const repeats_in );

	/// @brief Sets whether the residue that is being perturbed can stay within its own bin (in which case new mainchain
	/// torsions are drawn from within the bin), or whether it must jump to a different bin.  True means it must jump.
	void set_must_switch_bins( bool const val );

	/******************************
	** GETTERS                  **
	******************************/

	/// @brief Return the number of iterations at apply time.
	/// @details A value of 1 means that a single residue in the range will randomly be selected and
	/// flipped to another bin (with probabilities based on its neighbours and the bin transition
	/// probabilities.  Higher values mean that this operation will be repeated.  If set to 0, defaults
	/// to 1.
	core::Size repeats() const { return repeats_; }

	/// @brief Can the residue that is being perturbed stay within its own bin (in which case new mainchain
	/// torsions are drawn from within the bin), or must it jump to a different bin?
	bool must_switch_bins() const { return must_switch_bins_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Provide citations to the passed CitationCollectionList
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const override;

private:
	///////////////////////
	// PRIVATE VARIABLES //
	///////////////////////

	/// @brief Start of residue range.
	///
	core::Size start_res_ = 0;

	/// @brief End of residue range.
	///
	core::Size end_res_ = 0;

	/// @brief Bin transition probability data file.
	/// @details Defaults to "ABBA".
	std::string binfile_ = "ABBA";

	/// @brief Owning pointer to the BinTransitionCalculator object used by this mover.
	core::scoring::bin_transitions::BinTransitionCalculatorCOP bin_transition_calculator_;

	/// @brief How many times will this be applied?
	/// @details A value of 1 means that a single residue in the range will randomly be selected and
	/// flipped to another bin (with probabilities based on its neighbours and the bin transition
	/// probabilities.  Higher values mean that this operation will be repeated.  Defaults to 1.
	core::Size repeats_ = 1;

	/// @brief Can the residue that is being perturbed stay within its own bin (in which case new mainchain
	/// torsions are drawn from within the bin), or must it jump to a different bin?  Default false (i.e.
	/// need not switch bins).
	bool must_switch_bins_ = false;

	///////////////////////
	// PRIVATE FUNCTIONS //
	///////////////////////


};

} // bin_transitions
} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_bin_transitions_PerturbByBins_HH_
