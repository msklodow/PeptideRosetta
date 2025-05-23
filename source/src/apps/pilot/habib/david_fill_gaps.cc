// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk + dj

#include <iostream>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>




// Utility Headers
#include <string>


//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>

#include <fstream> // AUTO IWYU For ifstream


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR( "apps.pilot.david_recompute_score_and_rmsd.main" );


OPT_KEY( String, fa_file )

//set to store pdb info keys
std::vector<std::string> surface;
std::set <std::string> interface;


/// General testing code
int
main( int argc, char * argv [] )
{
	try {

		std::vector<std::string> surface;
		std::set <std::string> interface;
		NEW_OPT ( fa_file, "File name for fasta data","");

		using namespace core;
		using namespace core::scoring;

		devel::init(argc, argv);
		pose::Pose pose;
		std::string const input_pdb_name ( basic::options::start_file() );
		core::import_pose::pose_from_file( pose, input_pdb_name , core::import_pose::PDB_file);

		std::string const ffilename = option[ fa_file ];
		std::string seq;
		if ( ffilename != "" ) {
			std::ifstream ifs(ffilename.c_str(), std::ifstream::in);
			if ( !ifs.is_open() ) {
				std::cout<< "Error opening fasta file "<<ffilename<<std::endl;
				return -100;
			}
			while ( ifs.good() ) {
				std::string intres;
				ifs >> intres;
				if ( intres.length() >0 ) {
					if ( intres[0]!= '>' ) {
						seq.append(intres);
					}
				}
			}
		}

		if ( seq.length()==0 ) {
			std::cout<< "Fasta file contains no sequence: "<<ffilename<<std::endl;
			return -101;

		}

		int last = -100;
		core::chemical::ResidueTypeSet const & rsd_set( *pose.residue_type_set_for_pose() );

		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( last == -100 ) {
				last = pose.pdb_info()->number(i);
				std::cout<<i<<" "<< pose.pdb_info()->number(i)<<" "<<pose.residue(i).name1()<<"\n" ;
				continue;
			}
			if ( last > pose.pdb_info()->number(i)+1 ) {
				for ( int j = last-1; j > pose.pdb_info()->number(i); j-- ) {
					std::cout<<seq[j-1]<<"\n";
					// The representative type should have no/minimal variants
					core::chemical::ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name1( seq[j-1] ) );
					core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );
					pose.append_polymer_residue_after_seqpos(*new_rsd, i, false );
					pose.pdb_info()->number(i+2, pose.pdb_info()->number(i)-1);
					pose.pdb_info()->chain(i+2, pose.pdb_info()->chain(i));
				}

			} else if ( last < pose.pdb_info()->number(i)-1 ) {
				for ( int j = last+1; j < pose.pdb_info()->number(i); j++ ) {
					std::cout<<seq[j-1]<<"\n";
					// The representative type should have no/minimal variants
					core::chemical::ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name1( seq[j-1] ) );
					core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );
					pose.append_polymer_residue_after_seqpos(*new_rsd, i-1, false );
					pose.pdb_info()->number(i, pose.pdb_info()->number(i-1)+1);
					pose.pdb_info()->chain(i, pose.pdb_info()->chain(i-1));
					i++;
				}


			}
			std::cout<<i<<" "<< pose.pdb_info()->number(i)<<" "<<pose.residue(i).name1()<<"\n" ;
			last = pose.pdb_info()->number(i);
		}

		for ( Size i = 1; i <= pose.size(); ++i ) {
			std::cout<<pose.pdb_info()->number(i)<<"\n";
		}
		for ( int i = 1; i <pose.pdb_info()->number(1); ++i ) {
			std::cout<<" ";
		}
		for ( Size i = 1; i <= pose.size(); ++i ) {
			std::cout<<pose.residue(i).name1();
		}
		std::cout<<"\n"<<seq<<"\n";

		pose.pdb_info()->obsolete(false);
		pose.dump_pdb("test.pdb");

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;

}


