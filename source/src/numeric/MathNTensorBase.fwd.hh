// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/MathNTensorBase.fwd.hh
/// @brief  Owning pointer declarations for MathNTensorBase template class.
/// @details Note that if you're using an owning pointer to a MathNTensor, you need to either know at pointer declaration time
/// the dimensionality of the MathNTensor, or you need to be using an owning pointer to the base class (MathNTensorBase).
/// @author Vikram K. Mulligan (vmullig@uw.edu).


#ifndef INCLUDED_numeric_MathNTensorBase_fwd_hh
#define INCLUDED_numeric_MathNTensorBase_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <numeric/types.hh>

namespace numeric {

// Forward declaration
template< class T > class MathNTensorBase;

// Owning pointer
template< class T >
using MathNTensorBaseOP = utility::pointer::shared_ptr< MathNTensorBase< T > >; //Vikram is reluctantly using a C++11 feature.  Mrph.

// Const-access owning pointer
template< class T >
using MathNTensorBaseCOP = utility::pointer::shared_ptr< MathNTensorBase< T > const >;

template< class T >
MathNTensorBaseOP< T > deep_copy( MathNTensorBase< T > const & source );

} // namespace numeric

#endif // INCLUDED_numeric_MathNTensorBase_fwd_hh

