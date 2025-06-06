// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/RosettaTensorflowSessionContainer.tmpl.hh
/// @brief A container for Rosetta Tensorflow sessions, allowing sessions to be loaded once and stored in the global Tensorflow session manager.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_tensorflow_manager_RosettaTensorflowSessionContainer_tmpl_hh
#define INCLUDED_basic_tensorflow_manager_RosettaTensorflowSessionContainer_tmpl_hh

#ifdef USE_TENSORFLOW

#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.hh>

#include <basic/Tracer.hh>

#include <cstring> //memcpy

#ifdef USE_TENSORFLOW
#include <core/types.hh>
#endif

namespace basic {
namespace tensorflow_manager {

static basic::Tracer TRimpl( "basic.tensorflow_manager.RosettaTensorflowSessionContainer" );

/// @brief Given an input tensor and a place to stick the outputs, run the Tensorflow session.
/// @details This method is provided so that no one needs to handle the TF_Session object directly (and its
/// creation and destruction can be handled safely by the RosettaTensorflowSessionContainer).
/// @param[in]	input_name		The name of the input tensor in the Tensorflow session graph.
/// @param[in]	output_name		The name of the output tensor in the Tensorflow session graph.
/// @param[in]	input_tensor	The tensor of inputs.
/// @param[out]	ouput_tensor	The place to stick the outputs, overwritten by this operation.
/// @param[out]	runtime			The time for the actual Tensorflow session evaluation, in microseconds.  The contents of runtime will be overwritten by tihs operation.
/// @note There is no tracer output produced by this operation.  If you wish to write out runtime information, do something
/// with the runtime output variable.
template< typename T1, typename T2 >
void
RosettaTensorflowSessionContainer::run_session(
	std::string const & input_name,
	std::string const & output_name,
	RosettaTensorflowTensorContainer<T1> const & input_tensor,
	RosettaTensorflowTensorContainer<T2> & output_tensor,
	std::chrono::duration< double, std::micro > & runtime
) const {
	//Debug-mode sanity checks:
	debug_assert( status_ != nullptr );
	debug_assert( graph_ != nullptr );
	debug_assert( session_ != nullptr );

	//Get the locations of the input to the graph and the output from the graph:
	TF_Operation* graph_op_in( TF_GraphOperationByName( graph_, input_name.c_str() ) );
	TF_Operation* graph_op_out( TF_GraphOperationByName( graph_, output_name.c_str() ) );
	debug_assert( graph_op_in != nullptr );
	debug_assert( graph_op_out != nullptr );
	TF_Output graph_inputs( TF_Output{ graph_op_in, 0 } );
	TF_Output graph_outputs( TF_Output{ graph_op_out, 0 } );

	debug_assert( input_tensor.tensor_ != nullptr );

	//Actually run the network given the input data.  (Note that the RosettaTensorflowSessionContainer is the sole class allowed to access the private
	//member data of the RosettaTensorflowTensorContainer class.):
	std::chrono::time_point<std::chrono::system_clock> const starttime( ROSETTA_TENSORFLOW_CLOCK::now() );
	TF_SessionRun( session_, nullptr, &graph_inputs, &(input_tensor.tensor_), 1, &graph_outputs, &(output_tensor.tensor_), 1, nullptr, 0, nullptr, status_);
	std::chrono::time_point<std::chrono::system_clock> const endtime( ROSETTA_TENSORFLOW_CLOCK::now() );
	runtime = endtime - starttime;

	//Check that the run was okay:
	if ( TF_GetCode(status_) != TF_OK ) {
		utility_exit_with_message( "Unable to evaluate TensorFlow session: " + std::string( TF_Message(status_) ) );
	}

	//Ensure that the output tensor is in a good state for reading out the output:
	output_tensor.update_tensor_data_pointer();
}

/// @brief Given a vector of input tensors and vector of outputs, run the Tensorflow session.
/// @details This method is provided so that no one needs to handle the TF_Session object directly (and its
/// creation and destruction can be handled safely by the RosettaTensorflowSessionContainer).
/// @param[in]	input_name				The name of the input tensor in the Tensorflow session graph.
/// @param[in]	output_name				The name of the output tensor in the Tensorflow session graph.
/// @param[in]	input_tensor_vector		The vector of tensor of inputs.
/// @param[out]	output_tensor_vector	The vector of places to stick the outputs, overwritten by this operation.
/// @param[out]	runtime					The time for the actual Tensorflow session evaluation, in microseconds.  The contents of runtime will be overwritten by tihs operation.
/// @note There is no tracer output produced by this operation.  If you wish to write out runtime information, do something
/// with the runtime output variable.  Also note that the tensor types must be the same in the input and output vectors.
template< typename T1, typename T2 >
void
RosettaTensorflowSessionContainer::multirun_session(
	std::string const & input_name,
	std::string const & output_name,
	utility::vector1< RosettaTensorflowTensorContainer<T1> > const & input_tensor_vector,
	utility::vector1< RosettaTensorflowTensorContainer<T2> > & output_tensor_vector,
	std::chrono::duration< double, std::micro > & runtime
) const {
	RosettaTensorflowTensorContainer< T1 > input_tensors;
	RosettaTensorflowTensorContainer< T1 >::combine_tensors( input_tensor_vector, input_tensors );
	RosettaTensorflowTensorContainer< T2 > output_tensors;
	RosettaTensorflowTensorContainer< T2 >::combine_tensors( output_tensor_vector, output_tensors );

	run_session( input_name, output_name, input_tensors, output_tensors, runtime );

	//Copy data back into original output:
	RosettaTensorflowTensorContainer< T2 >::split_combined_tensors( output_tensors, output_tensor_vector );
}


/// @brief Equivalent of run_session() but for networks with multiple heads (multiple input tensors per run).
/// @details This method is provided so that no one needs to handle the TF_Session object directly (and its
/// creation and destruction can be handled safely by the RosettaTensorflowSessionContainer).
/// @param[in]	input_names		The name of the input tensors in the Tensorflow session graph.
/// @param[in]	output_name		The name of the output tensor in the Tensorflow session graph.
/// @param[in]	input_tensor	The tensors of inputs, must be in the same order as input_names.
/// @param[out]	ouput_tensor	The place to stick the outputs, overwritten by this operation.
/// @param[out]	runtime			The time for the actual Tensorflow session evaluation, in microseconds.
/// The contents of runtime will be overwritten by tihs operation.
/// @note There is no tracer output produced by this operation.  If you wish to write out runtime information,
/// do something with the runtime output variable.
template< typename T1, typename T2 >
void
RosettaTensorflowSessionContainer::run_multiinput_session(
 	utility::vector1< std::string > const & input_names,
 	std::string const & output_name,
 	utility::vector1< RosettaTensorflowTensorContainer<T1> > const & input_tensors,
 	RosettaTensorflowTensorContainer<T2> & output_tensor,
 	std::chrono::duration< double, std::micro > & runtime
) const {
	//Debug-mode sanity checks:
	debug_assert( status_ != nullptr );
 	debug_assert( graph_ != nullptr );
 	debug_assert( session_ != nullptr );
	debug_assert( ! input_names.empty() );
 	debug_assert( input_names.size() == input_tensors.size() );

	platform::Size const ninputs = input_names.size();
	platform::Size const noutputs = 1;

	std::unique_ptr<TF_Output []> inputs(new TF_Output [ninputs]);
	std::unique_ptr<TF_Tensor * []> input_values(new TF_Tensor * [ninputs]);
	std::unique_ptr<TF_Output[]> outputs(new TF_Output[noutputs]);

	for( platform::Size ii = 1; ii <= input_tensors.size(); ++ii ){
		input_values.get()[ ii - 1 ] = input_tensors[ ii ].tensor_;
		debug_assert( input_tensors[ ii ].tensor_ != nullptr );

		TF_Operation * const op = TF_GraphOperationByName( graph_, input_names[ ii ].c_str() );
		inputs.get()[ ii - 1 ] = TF_Output{ op, 0 };
		debug_assert( op != nullptr );
	}

	{
		TF_Operation * const out_op = TF_GraphOperationByName( graph_, output_name.c_str() );
		outputs.get()[0] = TF_Output{ out_op, 0 };
		debug_assert( out_op != nullptr);
	}

	//Actually run the network given the input data.  (Note that the RosettaTensorflowSessionContainer is the sole class allowed to access the private
	//member data of the RosettaTensorflowTensorContainer class.):
 	std::chrono::time_point<std::chrono::system_clock> const starttime( ROSETTA_TENSORFLOW_CLOCK::now() );
	TF_SessionRun( session_, nullptr, inputs.get(), input_values.get(), ninputs,
		outputs.get(), &(output_tensor.tensor_), noutputs, nullptr, 0,
		nullptr, status_ );
 	std::chrono::time_point<std::chrono::system_clock> const endtime( ROSETTA_TENSORFLOW_CLOCK::now() );
 	runtime = endtime - starttime;

	//Check that the run was okay:
 	if ( TF_GetCode(status_) != TF_OK ) {
 		utility_exit_with_message( "Unable to evaluate TensorFlow session: " + std::string( TF_Message(status_) ) );
 	}


	//Ensure that the output tensor is in a good state for reading out the output:
	output_tensor.update_tensor_data_pointer();
}

/// @brief Equivalent of run_session() but for networks with multiple heads (multiple input tensors per run) and
/// multiple tails (multiple output tensors per run).
/// @details This method is provided so that no one needs to handle the TF_Session object directly (and its
/// creation and destruction can be handled safely by the RosettaTensorflowSessionContainer).
/// @param[in]	input_names		The names of the input tensors in the Tensorflow session graph.
/// @param[in]	output_names	The names of the output tensor in the Tensorflow session graph.
/// @param[in]	input_tensor	The tensors of inputs.  Must be in the same order as input_names.
/// @param[out]	ouput_tensor	The places to stick the outputs, overwritten by this operation.   Must be in the same order as output_names.
/// @param[out]	runtime			The time for the actual Tensorflow session evaluation, in microseconds.  The contents of runtime will be
/// overwritten by tihs operation.
/// @note There is no tracer output produced by this operation.  If you wish to write out runtime information, do something
/// with the runtime output variable.
/// @author Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitute.org).
template< typename T1, typename T2 >
void
RosettaTensorflowSessionContainer::run_multiinput_multioutput_session(
 	utility::vector1< std::string > const & input_names,
 	utility::vector1< std::string > const & output_names,
 	utility::vector1< RosettaTensorflowTensorContainer<T1> > const & input_tensors,
 	utility::vector1< RosettaTensorflowTensorContainer<T2> > & output_tensors,
 	std::chrono::duration< double, std::micro > & runtime
) const {
	//Debug-mode sanity checks:
	debug_assert( status_ != nullptr );
 	debug_assert( graph_ != nullptr );
 	debug_assert( session_ != nullptr );
	debug_assert( ! input_names.empty() );
 	debug_assert( input_names.size() == input_tensors.size() );
	debug_assert( !output_names.empty() );
	debug_assert( output_names.size() == output_tensors.size() );

	platform::Size const ninputs = input_names.size();
	platform::Size const noutputs = output_names.size();

	std::unique_ptr<TF_Output []> inputs(new TF_Output [ninputs]);
	std::unique_ptr<TF_Tensor * []> input_values(new TF_Tensor * [ninputs]);
	std::unique_ptr<TF_Output []> outputs(new TF_Output[noutputs]);
	std::unique_ptr<TF_Tensor * []> output_values(new TF_Tensor * [noutputs]);

	for( platform::Size ii = 1; ii <= ninputs; ++ii ){
		debug_assert( input_tensors[ ii ].tensor_ != nullptr );
		input_values.get()[ ii - 1 ] = input_tensors[ ii ].tensor_;

		TF_Operation * const op = TF_GraphOperationByName( graph_, input_names[ ii ].c_str() );
		inputs.get()[ ii - 1 ] = TF_Output{ op, 0 };
		debug_assert( op != nullptr );
	}

	for( platform::Size jj = 1; jj <= noutputs; ++jj ) {
		// debug_assert( output_tensors[ jj ].tensor_ != nullptr );
		// output_values.get()[ jj - 1 ] = output_tensors[ jj ].tensor_;

		TF_Operation * const out_op = TF_GraphOperationByName( graph_, output_names[ jj ].c_str() );
		outputs.get()[jj - 1] = TF_Output{ out_op, 0 };
		debug_assert( out_op != nullptr);
	}

	//Actually run the network given the input data.  (Note that the RosettaTensorflowSessionContainer is the sole class allowed to access the private
	//member data of the RosettaTensorflowTensorContainer class.):
 	std::chrono::time_point<std::chrono::system_clock> const starttime( ROSETTA_TENSORFLOW_CLOCK::now() );
	TF_SessionRun( session_, nullptr, inputs.get(), input_values.get(), ninputs,
		outputs.get(), output_values.get(), noutputs, nullptr, 0,
		nullptr, status_ );
 	std::chrono::time_point<std::chrono::system_clock> const endtime( ROSETTA_TENSORFLOW_CLOCK::now() );
 	runtime = endtime - starttime;

	//Check that the run was okay:
 	if ( TF_GetCode(status_) != TF_OK ) {
 		utility_exit_with_message( "Unable to evaluate TensorFlow session: " + std::string( TF_Message(status_) ) );
 	}

	//Ensure that the output tensors are in a good state for reading out the output:
	for( core::Size i(1); i<=noutputs; ++i ) {
		if(output_tensors[i].tensor_ != nullptr) {
			TF_DeleteTensor( output_tensors[i].tensor_ );
		}
		output_tensors[i].tensor_ = output_values[i-1];
		output_tensors[i].update_tensor_data_pointer();
	}
}

template< typename T1, typename T2 >
void
RosettaTensorflowSessionContainer::multirun_multiinput_session(
	utility::vector1< std::string > const & input_names,
	std::string const & output_name,
	utility::vector1< utility::vector1< RosettaTensorflowTensorContainer<T1> > > const & input_tensors_vector,
	utility::vector1< RosettaTensorflowTensorContainer<T2> > & output_tensor_vector,
	std::chrono::duration< double, std::micro > & runtime
) const {
	if( input_tensors_vector.empty() ) return;
	runtime_assert( input_tensors_vector[1].size() == output_tensor_vector.size() );

#ifndef NDEBUG
	for( platform::Size ii = 2; ii <= input_tensors_vector.size(); ++ii ){
		debug_assert( input_tensors_vector[ ii ].size() == input_tensors_vector[ 1 ].size() );
	}
#endif

	utility::vector1< RosettaTensorflowTensorContainer< T1 > > combined_input_tensors( input_tensors_vector.size() );
	for( platform::Size ii = 1; ii <= input_tensors_vector.size(); ++ii ){
		RosettaTensorflowTensorContainer< T1 >::combine_tensors( input_tensors_vector[ ii ], combined_input_tensors[ ii ] );
	}

	RosettaTensorflowTensorContainer< T2 > output_tensors;
	RosettaTensorflowTensorContainer< T2 >::combine_tensors( output_tensor_vector, output_tensors );

	run_multiinput_session( input_names, output_name, combined_input_tensors, output_tensors, runtime );

	//Copy data back into original output:
	RosettaTensorflowTensorContainer< T2 >::split_combined_tensors( output_tensors, output_tensor_vector );
}


} //namespace tensorflow_manager
} //namespace basic

#endif //USE_TENSORFLOW
#endif //INCLUDED_basic_tensorflow_manager_RosettaTensorflowSessionContainer_tmpl_hh
