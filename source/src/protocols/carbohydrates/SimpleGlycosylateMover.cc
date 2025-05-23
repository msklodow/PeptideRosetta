// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/SimpleGlycosylateMover.cc
/// @brief A mover for glycosylation of common biological glycosylations.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/carbohydrates/SimpleGlycosylateMover.hh>
#include <protocols/carbohydrates/SimpleGlycosylateMoverCreator.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/PDBInfo.hh>

#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/carbohydrates/database_io.hh>

#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/select/residue_selector/ResidueSelector.hh> // AUTO IWYU For ResidueSelector

static basic::Tracer TR( "protocols.carbohydrates.SimpleGlycosylateMover" );


namespace protocols {
namespace carbohydrates {
using namespace core::chemical::carbohydrates;
using namespace core::pose::carbohydrates;
using namespace core::kinematics;
using namespace protocols::rosetta_scripts;

SimpleGlycosylateMover::SimpleGlycosylateMover():
	protocols::moves::Mover( "SimpleGlycosylateMover" ),
	strip_existing_glycans_( true ),
	ref_pose_name_( "" ),
	idealize_glycosylation_( false )

{

}

SimpleGlycosylateMover::~SimpleGlycosylateMover()= default;

SimpleGlycosylateMover::SimpleGlycosylateMover( SimpleGlycosylateMover const & src ):
	protocols::moves::Mover( src ),
	glycosylations_( src.glycosylations_ ),
	glycosylation_weights_( src.glycosylation_weights_ ),
	parsed_positions_( src.parsed_positions_ ),
	positions_( src.positions_ ),
	parsed_atom_names_( src.parsed_atom_names_ ),
	atom_names_( src.atom_names_ ),
	strip_existing_glycans_( src.strip_existing_glycans_ ),
	ref_pose_name_( src.ref_pose_name_ ),
	idealize_glycosylation_( src.idealize_glycosylation_)
{
	if ( src.selector_ ) {
		selector_ = src.selector_->clone();
	}
}

void
SimpleGlycosylateMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data
) {

	using namespace core::kinematics;

	glycosylations_.clear();
	glycosylation_weights_.clear();
	parsed_positions_.clear();
	positions_.clear();
	parsed_atom_names_.clear();
	atom_names_.clear();

	if ( tag->hasOption("glycosylation") ) {
		glycosylations_.push_back( tag->getOption< std::string >("glycosylation"));

	} else if ( tag->hasOption("glycosylations") ) {
		glycosylations_ = utility::string_split_multi_delim( tag->getOption< std::string >("glycosylations"), ",'`~+*&|;. ");
	} else {
		utility_exit_with_message("Must pass either glycosylation or glycosylations!");
	}

	if ( tag->hasOption("position") ) {
		parsed_positions_.push_back( tag->getOption< std::string >("position") );
	} else if ( tag->hasOption("positions") ) {
		parsed_positions_ = utility::string_split_multi_delim( tag->getOption< std::string >("positions"), ",'`~+*&|;. ");
	} else if  ( tag->hasOption( "residue_selector") ) {
		selector_ = parse_residue_selector( tag, data );
	} else {
		utility_exit_with_message(" Must pass either position or positions");
	}

	if ( parsed_positions_.size() > 0 && selector_ ) {
		utility_exit_with_message(" Cannot set position(s) and residue_selector! ");
	}

	if ( tag->hasOption("atom_name") ) {
		parsed_atom_names_.push_back( tag->getOption< std::string >("atom_name") );
	} else if ( tag->hasOption("atom_names") ) {
		parsed_atom_names_ =  utility::string_split_multi_delim( tag->getOption< std::string >("atom_names"), ",'`~+*&|;. ");
	}

	if ( parsed_atom_names_.size() > 0 && selector_ ) {
		utility_exit_with_message(" Cannot set atom_name(s) and residue_selector! ");
	}


	if ( tag->hasOption("weights") ) {
		utility::vector1< std::string > weights = utility::string_split_multi_delim( tag->getOption< std::string >("weights"), ",'`~+*&|;. ");
		for ( core::Size i = 1; i <= weights.size(); ++i ) {
			glycosylation_weights_.push_back( utility::string2Real( weights[ i ]));
		}
	}

	//Convert positions to proper resnums.

	strip_existing_glycans_ = tag->getOption<bool>("strip_existing", strip_existing_glycans_);

	ref_pose_name_ = tag->getOption< std::string >("ref_pose_name", ref_pose_name_);
	idealize_glycosylation_ = tag->getOption< bool >("idealize_glycosylation", idealize_glycosylation_);
}

protocols::moves::MoverOP
SimpleGlycosylateMover::clone() const{
	return utility::pointer::make_shared< SimpleGlycosylateMover >( *this );
}

/*
SimpleGlycosylateMover & SimpleGlycosylateMoveroperator=( SimpleGlycosylateMover const & src){
return SimpleGlycosylateMover( src );
}
*/


moves::MoverOP
SimpleGlycosylateMover::fresh_instance() const
{
	return utility::pointer::make_shared< SimpleGlycosylateMover >();
}


void
SimpleGlycosylateMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

std::ostream &operator<< (std::ostream &os, SimpleGlycosylateMover const &mover)
{
	mover.show(os);
	return os;
}

void
SimpleGlycosylateMover::set_position(core::Size position){
	positions_.clear();
	positions_.push_back( position );

}

void
SimpleGlycosylateMover::set_positions(const utility::vector1<bool> &positions){
	positions_.clear();
	if ( positions.size() == 0 ) {
		utility_exit_with_message(" Positions empty! ");
	}

	for ( core::Size i = 1; i <= positions.size(); ++i ) {
		if ( positions[ i ] ) {
			positions_.push_back( i );
		}
	}
}

void
SimpleGlycosylateMover::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector ){
	positions_.clear();
	selector_ = selector->clone();
}


void
SimpleGlycosylateMover::set_glycosylation(const std::string & iupac_or_common_string){
	glycosylations_.clear();
	glycosylations_.push_back(iupac_or_common_string);
}

void
SimpleGlycosylateMover::set_glycosylations(const utility::vector1<std::string> &iupac_strings){
	glycosylations_ = iupac_strings;
}

void
SimpleGlycosylateMover::set_glycosylation_weights(const utility::vector1<core::Real> &weights){
	glycosylation_weights_ = weights;
}

void
SimpleGlycosylateMover::set_strip_existing_glycans(bool strip_existing){
	strip_existing_glycans_ = strip_existing;
}

void
SimpleGlycosylateMover::remove_index( utility::vector1< core::Size > & current_vector, core::Size resnum) const{

	//Not the fastest way to do this, but it should work for now.
	utility::vector1< core::Size > vector_copy = current_vector;
	current_vector.clear();
	for ( core::Size i = 1; i <= vector_copy.size(); ++i ) {
		core::Size item = vector_copy[ i ];
		if ( item != resnum ) {
			current_vector.push_back( vector_copy[ i ]);
		} else {
			continue;
		}
	}
}

utility::vector1< std::string >
SimpleGlycosylateMover::setup_and_load_iupac_sequences() const {

	using namespace core::chemical::carbohydrates;

	std::map< std::string, std::string > const & short_names = core::chemical::carbohydrates::CarbohydrateInfoManager::get_short_name_to_iupac_strings_map();

	std::string common_names_db = "chemical/carbohydrates/common_glycans/";
	utility::vector1< std::string > full_iupac_glycans;

	for ( core::Size i = 1; i <= glycosylations_.size(); ++i ) {
		std::string set_glycan = glycosylations_[ i ];

		std::map< std::string, std::string >::const_iterator it;
		it  = short_names.find( set_glycan );



		// Check for iupac extension and attempt to load the sequence from the iupac file.
		if ( it != short_names.end() ) {
			full_iupac_glycans.push_back( it->second );
			continue;
		} else if ( set_glycan.find(".iupac") != std::string::npos ) {
			// Check for short name
			std::string iupac_file_path =basic::database::find_database_path( common_names_db, set_glycan );
			std::string full_glycan = read_glycan_sequence_file( iupac_file_path );
			full_iupac_glycans.push_back( full_glycan );

			continue;

		} else {
			// Otherwise, it should be a full iupac glycan name.  Push it back.
			full_iupac_glycans.push_back( set_glycan );
			continue;
		}
	}

	return full_iupac_glycans;

}


void
SimpleGlycosylateMover::apply( core::pose::Pose& pose ){

	using namespace numeric::random;
	using namespace core::pose::carbohydrates;
	using namespace core::chemical::carbohydrates;

	//Since we may be deleting residues, we need to add a reference pose.

	//Convert parsed positions.
	if ( parsed_positions_.size() > 0 && parsed_atom_names_.size() == 0 ) {
		positions_.clear();
		for ( core::Size i = 1; i <= parsed_positions_.size(); ++i ) {
			core::Size resnum = core::pose::parse_resnum( parsed_positions_[ i ], pose);
			if ( ref_pose_name_ != "" ) {
				resnum = pose.corresponding_residue_in_current( resnum, ref_pose_name_);
			}
			positions_.push_back(resnum);
		}
	} else if ( parsed_positions_.size() > 0 && parsed_atom_names_.size() > 0 ) {
		positions_.clear();
		atom_names_.clear();
		for ( core::Size i = 1; i <= parsed_positions_.size(); ++i ) {
			core::Size resnum = core::pose::parse_resnum( parsed_positions_[ i ], pose);
			if ( ref_pose_name_ != "" ) {
				resnum = pose.corresponding_residue_in_current( resnum, ref_pose_name_);
			}
			positions_.push_back(resnum);
			std::string linker_atom = parsed_atom_names_ [ i ];
			// it would be a good idea to check that the resnum has that atom
			atom_names_[ resnum ] = linker_atom;
		}
	} else if ( selector_ ) {
		if ( ref_pose_name_ != "" ) {
			utility_exit_with_message("Using a ResidueSelector with a set ref pose name is currently unsupported!");
		}
		utility::vector1< bool > subset = selector_->apply( pose );
		set_positions( subset);
		if ( positions_.size() == 0 ) {
			TR << "Empty subset passed! This is probably part of a larger protocol.  Continuing." << std::endl;
			return;
		}
	}

	if ( positions_.size() == 0 ) {
		utility_exit_with_message(" Position(s) need to be set! ");
	}
	if ( glycosylations_.size() == 0 ) {
		utility_exit_with_message(" Glycosylation(s) need to be set! ");
	}
	if ( glycosylation_weights_.size() > 0 && glycosylation_weights_.size() != glycosylations_.size() ) {
		utility_exit_with_message("Number of weights must equal number of glycosylations!");
	}

	if ( atom_names_.size() > 0 && positions_.size() != atom_names_.size() ) {
		utility_exit_with_message("Number of atom_names must equal number of positions, one atom_name for each position in the list.");
	}

	std::string ref_pose_name = "simple_glycosylate_mover";
	pose.reference_pose_from_current( ref_pose_name , true /* override */); //Make a refpose as we may be deleting/adding/etc.

	utility::vector1< std::string > glycosylations =  setup_and_load_iupac_sequences();
	WeightedSampler sampler;


	//Go through positions randomly until all are tried.
	// Each round we remove the residue from local positions.
	utility::vector1< core::Size > local_positions = positions_;
	for ( core::Size round = 1; round <= positions_.size(); ++round ) {


		core::Size local_p = numeric::random::random_range(1, local_positions.size());
		core::Size old_resnum = local_positions[ local_p ];
		core::Size resnum = pose.corresponding_residue_in_current( old_resnum , ref_pose_name);
		if ( local_positions.size() > 1 ) {
			remove_index(local_positions, resnum); //Why can't we have a method that takes an index and deletes it?
		}

		//Does the position already have a glycan attached to it?

		if ( strip_existing_glycans_ == false ) {
			utility_exit_with_message(" Glycan extension not currently implemented!");
		} else {
			delete_carbohydrate_branch( pose, resnum ); //Delete any carbohydrates currently attached.
		}

		std::string glycosylation;

		if ( glycosylations.size() == 1 ) {
			glycosylation = glycosylations[ 1 ];
		} else {
			if ( glycosylation_weights_.size() > 0 ) {
				sampler.weights( glycosylation_weights_ );
				glycosylation = glycosylations[ sampler.random_sample(numeric::random::rg()) ];
			} else {
				glycosylation = glycosylations[ random_range( 1, glycosylations.size() ) ];
			}
		}

		if ( atom_names_.size() == 0 ) { // the base case, where canonical amino acids are glycosylated and the link atom is defined by consensus
			//Glycosylate the pose.
			TR << "Glycosylating at " << old_resnum << " : " << resnum << " " << glycosylation << std::endl;

			glycosylate_pose( pose, resnum, glycosylation, idealize_glycosylation_ /* idealize linkages - Seems to be a bug here!*/);
		}

		// If defining the atom of the residue to be chemically-linked
		// do we check that the residue contains the atom here? or is that done later
		if ( atom_names_.size() >= 1 ) {
			glycosylate_pose( pose, resnum, atom_names_[ resnum ], glycosylation, idealize_glycosylation_);
		}




	}

	//Fix PDBINFO.  We are adding Glycans, we should retain any PDB Info set.
	pose.pdb_info()->obsolete(false);


	/* This doesn't work for O-linked glycosylations
	//Check positions are ASN.  Optionally mutate to ASN and motif?
	simple_moves::MutateResidue mutate = simple_moves::MutateResidue();
	mutate.set_preserve_atom_coords( true );
	for (core::Size i = 1; i <= positions_.size(); ++i ){
	core::Size resnum = positions_[ i ];
	if (! pose.residue( resnum ).aa() == core::chemical::aa_asn ){
	TR << "Mutating "<< resnum << " to ASN for glycosylation. " << std::endl;
	mutate.set_target( resnum );
	mutate.set_res_name( core::chemical::aa_asn );
	mutate.apply( pose );
	}
	}
	*/

}


// Citation Management
//
/// @brief Provide the citation.
void
SimpleGlycosylateMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	using namespace basic::citation_manager;

	citations.add(
		utility::pointer::make_shared< UnpublishedModuleInfo >(
		mover_name(),
		CitedModuleType::Mover,
		"Jared Adolf-Bryfogle",
		"The Scripps Research Institute, La Jolla, CA",
		"jadolfbr@gmail.com"
		)
	);
}


/////////////// Creator ///////////////

std::string SimpleGlycosylateMover::get_name() const {
	return mover_name();
}

std::string SimpleGlycosylateMover::mover_name() {
	return "SimpleGlycosylateMover";
}

void SimpleGlycosylateMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;

	XMLSchemaRestriction delimited_real_list;
	//std::string real_list_regex = real_regex_pattern() + "([,'`~\\+\\*&\\|;\\. ]" + real_regex_pattern() + ")*";
	// AMW: temp -- can't have &
	std::string real_list_regex = real_regex_pattern() + "([,'`~\\+\\*\\|;\\. ]" + real_regex_pattern() + ")*";
	delimited_real_list.name( "delimited_real_list" );
	delimited_real_list.base_type( xs_string );
	delimited_real_list.add_restriction( xsr_pattern, real_list_regex );
	xsd.add_top_level_element( delimited_real_list );


	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "glycosylation", xs_string, "String specifying the glycosylation to add to this pose (IUPAC Glycan String) (or a file name containing the glycosylation with a .iupac name)" )
		+ XMLSchemaAttribute( "glycosylations", xs_string, "String or file name specifying multiple possible glycosylations to add to this pose (IUPAC Glycan String) (or a file name containing the glycosylation with a .iupac name) (will be sampled randomly" )
		+ XMLSchemaAttribute( "position", xsct_refpose_enabled_residue_number, "Position to add glycosylation" )
		+ XMLSchemaAttribute( "positions", xsct_refpose_enabled_residue_number_cslist, "Positions to add glycosylations" )
		+ XMLSchemaAttribute( "atom_name", xs_string, "Atom to which glycan is chemically linked" )
		+ XMLSchemaAttribute( "weights", "delimited_real_list", "Sampling weights corresponding to the provided set of glycans" )
		+ XMLSchemaAttribute( "strip_existing", xsct_rosetta_bool, "Strip existing glycosylations from the pose" )
		+ XMLSchemaAttribute( "ref_pose_name", xs_string, "Name of saved reference pose" )
		+ XMLSchemaAttribute( "idealize_glycosylation", xsct_rosetta_bool, "Idealize glycosylations on the pose?" );

	rosetta_scripts::attributes_for_parse_residue_selector( attlist );

	XMLSchemaSimpleSubelementList subelements;

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(),
		"Authors: Jared Adolf-Bryfogle (jadolfbr@gmail.com) and Jason W. Labonte (JWLabonte@jhu.edu)\n"
		"Mover to add specified glycosylations to a pose in the specified positions", attlist, subelements );
}

std::string SimpleGlycosylateMoverCreator::keyname() const {
	return SimpleGlycosylateMover::mover_name();
}

protocols::moves::MoverOP
SimpleGlycosylateMoverCreator::create_mover() const {
	return utility::pointer::make_shared< SimpleGlycosylateMover >();
}

void SimpleGlycosylateMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SimpleGlycosylateMover::provide_xml_schema( xsd );
}


} //protocols
} //carbohydrates
