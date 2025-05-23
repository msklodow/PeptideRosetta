// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/SidechainRmsdFilter.hh
/// @brief A filter based on automorphic sidechain RMSD
/// @author Noah Ollikainen

#ifndef INCLUDED_protocols_simple_filters_SidechainRmsdFilter_hh
#define INCLUDED_protocols_simple_filters_SidechainRmsdFilter_hh

#include <protocols/simple_filters/SidechainRmsdFilter.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.hh>

namespace protocols {
namespace simple_filters {

class SidechainRmsdFilter : public filters::Filter
{
public:
	SidechainRmsdFilter();
	SidechainRmsdFilter( std::string const & res1, std::string const & res2, core::Real const rmsd_threshold );
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const override;
	filters::FilterOP fresh_instance() const override;

	~SidechainRmsdFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string res1_, res2_;
	core::Real rmsd_threshold_;
	core::pose::PoseCOP reference_pose_;
	bool include_backbone_;
};

}
}

#endif
