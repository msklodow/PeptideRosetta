// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/SaveResfileToDiskFilterCreator.hh
/// @brief  FilterCreator for the SaveResfileToDiskFilter
/// @author Neil King (neilking@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_SaveResfileToDiskFilterCreator_hh
#define INCLUDED_protocols_simple_filters_SaveResfileToDiskFilterCreator_hh


// Package Headers
#include <protocols/filters/FilterCreator.hh>

// Utility Headers

// c++ headers
#include <string>

namespace protocols {
namespace simple_filters {

class SaveResfileToDiskFilterCreator : public protocols::filters::FilterCreator
{
public:
	protocols::filters::FilterOP create_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};


} //namespace simple_filters
} //namespace protocols

#endif
