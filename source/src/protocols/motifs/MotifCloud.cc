// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief aginsparg, jcheng, ccummings, sthyme; MotifCloud object, which is a saveable spatial projection of motifs about residues in a PDB system


//#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/io/pdb/pdb_writer.hh> // pose_from_pdb
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/dna/setup.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/motifs/LigandMotifSearch.hh>
#include <protocols/motifs/LigandMotifSearch.fwd.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ChemicalManager.hh>
#include <protocols/motifs/Motif.hh>
#include <core/chemical/residue_io.hh>
#include <protocols/motifs/MotifHit.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <protocols/motifs/BuildPosition.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <core/chemical/GlobalResidueTypeSet.hh>

//#include <protocols/motifs/MotifLigandPacker.hh>
#include <protocols/dna/PDBOutput.hh>
#include <protocols/dna/util.hh>
//#include <protocols/dna/DnaInterfacePacker.hh>
//#include <protocols/dna/DnaInterfaceFinder.hh>
#include <protocols/motifs/motif_utils.hh>
using namespace protocols::dna;

#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh> // Need since refactor

#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>

// Utility Headers
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh> // file_exists
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
using utility::vector1;
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
using utility::string_split;

// c++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <queue>
#include <functional>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>
#include <basic/options/keys/protein_grid.OptionKeys.gen.hh>

//this may solve atom clash issue
#include <core/pose/xyzStripeHashPose.hh>
#include <numeric/geometry/hashing/xyzStripeHash.hh>
#include <numeric/geometry/hashing/xyzStripeHash.fwd.hh>

#include <numeric/xyzVector.hh>
#include <protocols/ligand_docking/InterfaceScoreCalculator.hh>
#include <protocols/ligand_docking/ligand_scores.hh>

#include <protocols/qsar/scoring_grid/AtrGrid.hh>
#include <protocols/qsar/scoring_grid/RepGrid.hh>
#include <protocols/qsar/scoring_grid/VdwGrid.hh>
#include <protocols/qsar/scoring_grid/HbaGrid.hh>
#include <protocols/qsar/scoring_grid/HbdGrid.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>

#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/LigandArea.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>
#include <protocols/ligand_docking/HighResDocker.hh>
#include <protocols/ligand_docking/ligand_dock_impl.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/chemical/PoseResidueTypeSet.hh>

#include <ObjexxFCL/FArray1D.hh>

#include <protocols/motifs/LigandDiscoverySearch.hh>

#include <protocols/motifs/IdentifyLigandMotifs.hh>

#include <utility/Binary_Util.hh>

#include <protocols/protein_grid/ProteinGrid.hh>

// Time profiling header
#include <time.h>

//for data type debugging
#include <typeinfo>

using namespace core;
using namespace basic;
using namespace chemical;
using namespace pack;
using namespace task;
using namespace scoring;
using namespace options;
using namespace OptionKeys;

//declare global tracer for class
static basic::Tracer ms_tr( "protocols.motifs.MotifCloud" );

// destructor
MotifCloud::~MotifCloud() = default;

// @brief create a new motifcloud from a given pose and motifs library that are hashed based on amino acid residue and atoms on the ligand side
MotifCloud::MotifCloud(core::pose::PoseOP in_pose, protocols::motifs::MotifCOPs motif_library)
{
	//initialize options/key variables
	initialize();

	//seed the working_pose_ from the in_pose and motif_library_ from motif_library
	working_pose_ = in_pose;
	motif_library_ = motif_library;

	//wrap the matrix around the pose
	wrap_matrix_around_pose();

	//project motifs from protein residues
	project_motifs_from_residues();
}

// @brief fill out a motifcloud from data read in from a file, where we had previously saved a motifcloud before
//DO THIS LATER AFTER ACTUALLY MAKING A MOTIFCLOUD
//MotifGrid::MotifMatrix(whatever data is needed to make the motifmatrix from how we save it externally)
//{
//
//}

// @brief blanket function to initialize relevant variables and options upon creation of a 
void MotifCloud::initialize()
{
	// atom type set operations
	generate_fa_standard_atom_type_set();
	set_num_atom_types();
}

// @brief seed the atom_type_set_ with the fa_standard library
void MotifCloud::generate_fa_standard_atom_type_set()
{
	//make and set a fa_standard atom type set
	atom_type_set_ = core::chemical::ChemicalManager::get_instance()->atom_type_set( "fa_standard" );
}

// @brief seed the num_atom_types_ value based on the number of atom types in the atom_type_set_
void MotifCloud::set_num_atom_types()
{
	num_atom_types_ = atom_type_set_->n_atomtypes();
}

// @brief critical function that builds/rebuilds (overwrites existing) the motif_matrix_ around the pose in working_pose_
//called on creation of the object
//this will scale the size of the first 3 dimensions of the matrix about the working pose and provide integer shifts to effectively translate coordinates in the pose to voxelized indices in the matrix
//all values in the 4th dimension (counts of each fa_standard atom type) will be seeded to 0 in each voxel in the 3d space
void MotifCloud::wrap_matrix_around_pose()
{
	//run through all atoms to derive a range of dimensions to contain the protein in a 3D  space
	//since we can't have negative indices, we need to normalize the coordinate values so that everything is positive
	//derive constant values based  on the most negative values in each dimension, and then add that constant to all coordinates
	//to be safe, values need to be seeded with an initial value, or errors could be thrown when deriving shift values

	int smallest_x = 1;
	int smallest_y = 1;
	int smallest_z = 1;

	int largest_x = 1;
	int largest_y = 1;
	int largest_z = 1;

	//create a list of coordinates of each atom to hold and work with to fill the protein_representation_matrix
	//can't seem to make a vector of xyzVector objects, so will need to just make a custome 2D vector  to  hold the data
	utility::vector1<numeric::xyzVector<int>> atom_coordinates;

	//determine largest and smallest x,y,z  values to determine dimensions of matrix
	for ( core::Size res_num = 1; res_num <= working_pose_->size(); ++res_num ) {
		for ( core::Size atom_num = 1; atom_num <= working_pose_->residue(res_num).natoms(); ++atom_num ) {
			//get the x,y,z data of the atom, rounded to the closest value
			numeric::xyzVector<int> atom_xyz;
			//floor the coordinates down for a constant negative directional shift
			atom_xyz.x() = std::floor(working_pose_->residue(res_num).xyz(atom_num).x());
			atom_xyz.y() = std::floor(working_pose_->residue(res_num).xyz(atom_num).y());
			atom_xyz.z() = std::floor(working_pose_->residue(res_num).xyz(atom_num).z());

			//safe handling for the first atom encountered to be set as the smallest and largest values
			if ( res_num == 1 && atom_num == 1 ) {
				smallest_x = atom_xyz.x();
				smallest_y = atom_xyz.y();
				smallest_z = atom_xyz.z();
				largest_x = atom_xyz.x();
				largest_y = atom_xyz.y();
				largest_z = atom_xyz.z();
				continue;
			}

			//determine if any of the values  are the smallest
			if ( smallest_x > atom_xyz.x() ) {
				smallest_x = atom_xyz.x();
			}
			if ( smallest_y > atom_xyz.y() ) {
				smallest_y = atom_xyz.y();
			}
			if ( smallest_z > atom_xyz.z() ) {
				smallest_z = atom_xyz.z();
			}

			//determine if any  are the largest
			if ( largest_x < atom_xyz.x() ) {
				largest_x = atom_xyz.x();
			}
			if ( largest_y < atom_xyz.y() ) {
				largest_y = atom_xyz.y();
			}
			if ( largest_z < atom_xyz.z() ) {
				largest_z = atom_xyz.z();
			}

			atom_coordinates.push_back(atom_xyz);

		}
	}

	//set (or reset) the shift and bound vectors to be 3 values of 0 before setting them
	//this is necessary at least for the initial setting when the object is constructed
	reset_xyz_vectors();

	//take negative values of the smallest values and then add 1 to derive the constants
	//the logic here should apply, whether the smallest value is positive or negative
	//for the smallest value in the system to be indexed to 1, you add the negative of itself + 1; this shift would be applied to all other atom coordinates
	xyz_shift_[1] = std::floor(((smallest_x * -1) + 1));
	xyz_shift_[2] = std::floor(((smallest_y * -1) + 1));
	xyz_shift_[3] = std::floor(((smallest_z * -1) + 1));

	//apply shift values to largest to get boundaries
	xyz_bound_[1] = std::floor((xyz_shift_[1] + largest_x));
	xyz_bound_[2] = std::floor((xyz_shift_[2] + largest_y));
	xyz_bound_[3] = std::floor((xyz_shift_[3] + largest_z));

	//create 3D matrix to roughly represent 3D coordinate space of protein
	ms_tr.Debug << "Creating protein clash coordinate matrix. Dimensions of matrix (in angstroms) are " << xyz_bound_[1] << "," << xyz_bound_[2] << "," << xyz_bound_[3] << std::endl;
	ms_tr.Debug << "Shift from from original coordinates (in angstroms) are: " << xyz_shift_[1] << "," << xyz_shift_[2] << "," << xyz_shift_[3] << std::endl;

	//wipe the current contents of the protein_matrix_ and reset fullness values
	//create an empty dummy vector, and then assign protein_matrix_ with it
	MotifMatrix dummy_matrix;
	motif_matrix_ = dummy_matrix;

	//matrix_volume_ = get_grid_volume();

	for ( core::Size x = 1; x <= xyz_bound_[1]; ++x ) {
		//make a 3D  matrix
		utility::vector1<utility::vector1<utility::vector1<core::Size>>> sub_matrix;

		for ( core::Size y = 1; y <= xyz_bound_[2]; ++y ) {
			//make a 2D  matrix
			utility::vector1<utility::vector1<core::Size>> sub_sub_matrix;

			//CONTINUE HERE! ADD 4th dimension of atom type vector
			for ( core::Size z = 1; z <= xyz_bound_[3]; ++z ) {
				//make a 1D matrix, length = the number of atom types in the atom type set + 1 (the last index is intended to be a sum so we don't have to bother calculating a sum in runtime for determining how dense a voxel is in occupation) seed with 0
				utility::vector1<core::Size> sub_sub_sub_matrix(num_atom_types_ + 1,  0);
				//push 1D  matrix into 2D
				sub_sub_matrix.push_back(sub_sub_sub_matrix);
			}
			//push the 2D matrix into 3D
			sub_matrix.push_back(sub_sub_matrix);
		}
		//push a 3D  matrix into the 2D matrix
		motif_matrix_.push_back(sub_matrix);
	}
}

// @brief reset (or set) the xyz shift and bound matrices to be 3 values of zeroes
void MotifCloud::reset_xyz_vectors()
{
	xyz_shift_ = utility::vector1<int>(3, 0);
	xyz_bound_ = utility::vector1<int>(3, 0);
}

// @brief function to take a motifs list and project it about all atoms in the working_pose_
void MotifCloud::project_motifs_from_residues()
{
	//iterate over each residue in the pose
	for ( core::Size resi_pos = 1; resi_pos <= working_pose_->size(); ++resi_pos ) {

		//very likely want to integrate code that allows a user to only make a cloud on specific residue indices
		//integrate a residue index gate if the user does specify wanting to only investigate select residues for a motif cloud

		//get the residue name
		std::string res_name = working_pose_->residue(resi_pos).name3();

		// if the map does not already contain a motif sublibrary for this residue type,
		// create it and store it
		auto iter = motif_library_map_by_residue_.find( res_name );

		//if the motif_library_map_by_residue_ map doesn't already have motifs for the workign residue type, attempt to derive them and populate the map
		//if not in the library map
		if( iter == motif_library_map_by_residue_.end() )
		{
			//call get_motif_sublibrary_by_aa and add the library to the map
			motif_library_map_by_residue_[res_name] = get_motif_sublibrary_by_aa( motif_library_, res_name );


			// update the iterator now that the map has been updated
			iter = motif_library_map_by_residue_.find( res_name );
		}

		//safety/sanity check if we somehow fail to get motifs (I guess if like a bad list is included or something)
		if ( iter == motif_library_map_by_residue_.end() ) {
			ms_tr.Warning << "Could not find or create motif sub-library for residue " << res_name << " at position " << resi_pos << std::endl;
			continue;
		}		

		//iterate over each motif for the residue and project atoms relative to the residue's position
		for ( auto motifcop : motif_library_map_by_residue_[ res_name ] ) {
			//collect residue name from motifcop
			std::string motif_residue_name(motifcop->restype_name1());

			//place the motif atoms in a 3D space and derive their coordinates so that we may add them to the motifmatrix
		}

	}
}