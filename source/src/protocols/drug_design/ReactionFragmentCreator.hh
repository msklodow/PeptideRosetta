// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/drug_design/ReactionFragmentCreator.hh
/// @brief  Class for instantiating a ReactionFragment
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_ReactionFragmentCreator_HH
#define INCLUDED_protocols_drug_design_ReactionFragmentCreator_HH

// Package headers
#include <protocols/chemistries/ChemistryCreator.hh>
#include <protocols/chemistries/Chemistry.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace drug_design {

class ReactionFragmentCreator : public protocols::chemistries::ChemistryCreator {
public:
	protocols::chemistries::ChemistryOP create_chemistry() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};


} //namespace drug_design
} //namespace protocols


#endif
