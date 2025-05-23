// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ResLvlTaskOperationFactory.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResLvlTaskOperationFactory_hh
#define INCLUDED_core_pack_task_operation_ResLvlTaskOperationFactory_hh

// Unit Headers
#include <core/pack/task/operation/ResLvlTaskOperationFactory.fwd.hh>

// Package Headers
#include <core/pack/task/operation/ResLvlTaskOperation.fwd.hh>
#include <core/pack/task/operation/ResLvlTaskOperationCreator.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <map>
#include <string>

#ifdef PYROSETTA
#include <utility/tag/Tag.hh>
#endif

namespace core {
namespace pack {
namespace task {
namespace operation {

// singleton class
class ResLvlTaskOperationFactory : public utility::SingletonBase< ResLvlTaskOperationFactory >
{
public:
	friend class utility::SingletonBase< ResLvlTaskOperationFactory >;
public:
	typedef std::map< std::string, ResLvlTaskOperationCreatorOP > RLTOC_Map;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;

public:
	void factory_register( ResLvlTaskOperationCreatorOP );

	/// @brief add a prototype, using its default type name as the map key
	void add_creator( ResLvlTaskOperationCreatorOP );
	bool has_type( std::string const & ) const;

	/// @brief Get the XML schema for a given ResLvlTaskOperation.
	/// @details Throws an error if the ResLvlTaskOperation is unknown to Rosetta.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void provide_xml_schema( std::string const &task_operation_name, utility::tag::XMLSchemaDefinition & xsd ) const;

	/// @brief Get a human-readable listing of the citations for a given residue level task operation, by taskop name.
	/// @details Returns an empty string if there are no citations.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	std::string
	get_citation_humanreadable(
		std::string const & taskop_name
	) const;

	/// @brief return new ResLvlTaskOperation by key lookup in rlto_map_
	ResLvlTaskOperationOP newRLTO( std::string const & ) const;

	/// @brief return new ResLvlTaskOperation by key lookup in rlto_map_ (with the
	/// new ResLvlTaskOperation parsing the input Tag)
	ResLvlTaskOperationOP newRLTO(
		std::string const & class_keyname,
		TagCOP tag
	) const;

	/// @brief The %ResLvlTaskOperationFactory is the point of entry for the definition of the XML Schemas
	/// for every ResLvlTaskOperation that may be instantiated from a file.  It is  responsible for defining
	/// an xs:group named "res_lvl_task_op" listing each of the residue-level-task-operation-complex
	/// types that may be initialized using the %ResLvlTaskOperationFactory and to iterate across each of the
	/// ResLvlTaskOperationCreator it contains asking them for the XML schema of the ResLvlTaskOperation they
	/// are responsible for creating.
	void define_res_lvl_task_op_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	/// @brief The name given to the XML schema group of all ResLvlTaskOperation s.
	static std::string res_lvl_task_op_xml_schema_group_name();

	/// @brief Access the creator map (for unit testing).
	RLTOC_Map const & creator_map() const;

private:

	ResLvlTaskOperationFactory();
	virtual ~ResLvlTaskOperationFactory();

	RLTOC_Map rltoc_map_;

};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
