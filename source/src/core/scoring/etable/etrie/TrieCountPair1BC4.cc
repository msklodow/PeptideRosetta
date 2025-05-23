// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/etable/etrie/TrieCountPair1BC4.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/etable/etrie/TrieCountPair1BC4.hh>

// Package Headers
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/methods/MMLJEnergyInter.hh>

#include <core/scoring/etable/etrie/CountPairData_1_1.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.hh>

#include <core/scoring/trie/trie_vs_trie.hh>
#include <core/scoring/trie/trie_vs_path.hh>

#include <core/scoring/lkball/lkbtrie/LKBAtom.fwd.hh>
#include <core/scoring/lkball/lkbtrie/LKBTrieEvaluator.hh>

// STL Headers

#include <core/scoring/etable/etrie/EtableAtom.hh>
#include <core/scoring/trie/RotamerTrie.hh>
#include <utility/vector1.hh>

#include <core/scoring/elec/electrie/ElecTrieEvaluator.hh> // AUTO IWYU For ElecTrieEvaluator
#include <core/scoring/lkball/lkbtrie/LKBAtom.hh> // AUTO IWYU For LKBAtom


namespace core {
namespace scoring {
namespace etable {
namespace etrie {

using namespace trie;

///---------- TYPE RESOLUTION FUNCTIONS ----------///


void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::TableLookupEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

///////////////////////// EtableEnergy -- analytic evaluation //////////////////////////////////

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::AnalyticEtableEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

//XRW_E_T1

// HBONDS
void
TrieCountPair1BC4::resolve_trie_vs_trie(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	hbonds::HBondEnergy const & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > &,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const )
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with HBondEnergy" );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	hbonds::HBondEnergy const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > &,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const )
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with HBondEnergy" );
}


/// Hack Elec E
void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_1 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_2 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairData_1_3 > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< elec::electrie::ElecAtom, CountPairDataGeneric > const & trie2,
	elec::electrie::ElecTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

////////////////////////////////// lkball //////////////////////////////////
void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_1 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_2 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairData_1_3 > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void
TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< lkball::lkbtrie::LKBAtom, CountPairDataGeneric > const & trie2,
	lkball::lkbtrie::LKBTrieEvaluator const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

/////////////////////////////// MMLJEnergyInter ////////////////////////////

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const cached_data //Can be nullptr
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector, cached_data );
}


/////////////////////////////// VDW Energy ////////////////////////////

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_trie(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*pair_energy_table*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & /*temp_table*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_trie reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_1 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_2 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairData_1_3 > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

void TrieCountPair1BC4::resolve_trie_vs_path(
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie1*/,
	trie::RotamerTrie< vdwaals::VDWAtom, CountPairDataGeneric > const & /*trie2*/,
	vdwaals::VDWTrieEvaluator const & /*sfxn*/,
	utility::vector1< core::PackerEnergy > & /*pair_energy_vector*/,
	utility::vector1< core::PackerEnergy > & /*temp_vector*/,
	core::scoring::trie::TrieVsTrieCachedDataContainerBase const * const /*cached_data*/
)
{
	utility_exit_with_message( "etable::etrie::TrieCountPair1BC4::resolve_trie_vs_path reached with VDW_Energy" );
}

} // namespace trie
} // namespace etable
} // namespace scoring
} // namespace core

