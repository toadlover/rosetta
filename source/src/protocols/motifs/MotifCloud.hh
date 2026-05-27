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


#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueTypeFinder.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.fwd.hh>
#include <protocols/motifs/LigandMotifSearch.fwd.hh>
#include <protocols/motifs/MotifLibrary.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <protocols/motifs/Motif.fwd.hh>
#include <protocols/motifs/MotifHit.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.fwd.hh>
#include <protocols/motifs/BuildPosition.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/minimization_packing/MinMover.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/pose/PDBInfo.fwd.hh>

#include <core/chemical/GlobalResidueTypeSet.fwd.hh>

#include <protocols/dna/PDBOutput.fwd.hh>

#include <basic/prof.fwd.hh>
#include <basic/Tracer.fwd.hh>



// Utility Headers
#include <utility/io/ozstream.fwd.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/excn/Exceptions.fwd.hh>

// c++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <queue>
#include <functional>
#include <map>

// option key includes

#include <numeric/xyzVector.fwd.hh>
#include <protocols/ligand_docking/InterfaceScoreCalculator.fwd.hh>

#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <protocols/ligand_docking/LigandArea.fwd.hh>
#include <protocols/ligand_docking/InterfaceBuilder.fwd.hh>
#include <protocols/ligand_docking/MoveMapBuilder.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/ligand_docking/HighResDocker.fwd.hh>
#include <protocols/ligand_docking/LigandDockProtocol.fwd.hh>
#include <core/scoring/func/HarmonicFunc.fwd.hh>
#include <core/chemical/PoseResidueTypeSet.fwd.hh>

#include <ObjexxFCL/FArray1D.fwd.hh>

#include <protocols/protein_grid/ProteinGrid.fwd.hh>

// Time profiling header
#include <time.h>

//for data type debugging
#include <typeinfo>

class MotifCloud
{
public:

	/// @brief condensing the data type name for the main 4D vector (or matrix) that represents the motif cloud
	// currently laid out in that the matrix internal data can only be of Size.
	// 4D vector has the first 3 dimensinos represent 1 angstrom vectors projected about a pose; 4th dimension is a vector1 with length equal to the number of atom types in fa_standard with unsigned ints to represent counts of each atom type found in each voxel (default to 0)
	typedef utility::vector1<utility::vector1<utility::vector1<utility::vector1<core::Size>>>> MotifMatrix;

	// destructor
	~MotifCloud();

	//constructors:

	// @brief create a new motifcloud from a given pose and motifs library that are hashed based on amino acid residue and atoms on the ligand side
	MotifCloud(core::pose::PoseOP in_pose, protocols::motifs::MotifCOPs motif_library);

	// @brief fill out a motifcloud from data read in from a file, where we had previously saved a motifcloud before
	//DO THIS LATER AFTER ACTUALLY MAKING A MOTIFCLOUD
	//MotifMatrix(whatever data is needed to make the motifmatrix from how we save it externally);



private:

	//Functions==========================================================================================

	// @brief default constructor
	//will need to use class functions to seed values for input pose and other potential input data
	//should only use parameterized or copy constructor
	MotifCloud();

	// @brief blanket function to initialize relevant variables and options upon creation of a 
	void initialize();

	// @brief seed the atom_type_set_ with the fa_standard library
	void generate_fa_standard_atom_type_set();

	// @brief seed the num_atom_types_ value based on the number of atom types in the atom_type_set_
	void set_num_atom_types();

	// @brief critical function that builds/rebuilds (overwrites existing) the motif_matrix_ around the pose in working_pose_
	//called on creation of the object
	//this will scale the size of the first 3 dimensions of the matrix about the working pose and provide integer shifts to effectively translate coordinates in the pose to voxelized indices in the matrix
	//all values in the 4th dimension (counts of each fa_standard atom type) will be seeded to 0 in each voxel in the 3d space
	void wrap_matrix_around_pose();

	// @brief reset (or set) the xyz shift and bound matrices to be 3 values of zeroes
	void reset_xyz_vectors();

	//Member Variables===================================================================================

	/// @brief condensing the data type name for the main 4D vector (or matrix) that represents the motif cloud
	// currently laid out in that the matrix internal data can only be of Size.
	// 4D vector has the first 3 dimensinos represent 1 angstrom vectors projected about a pose; 4th dimension is a vector1 with length equal to the number of atom types in fa_standard with unsigned ints to represent counts of each atom type found in each voxel (default to 0)
	MotifMatrix motif_matrix_;

	// @brief the corresponding pose that the motif_matrix_ is wrapped around
	//it is possible that this will go unused if the user loads in a motif_matrix, as the original pose isn't explicitly necessary to represent the motif matrix if we have already projected motifs
	core::pose::PoseOP working_pose_;

	// @brief the full input motifs list provided by the user
	protocols::motifs::MotifCOPs motif_library_;

	// @brief a map of smaller motif libraries separated by protein residue for faster processing at the cost of memory overhead
	//will be filled out when iterating over residues to project motifs from, and populated using the motif_utils get_motif_sublibrary_by_aa() function. Keys will be strings for residue 3 letter codes
	std::map< std::string, protocols::motifs::MotifCOPs > motif_library_map_by_residue_;

	// @brief vector to hold the values that atom coordinates shift by in the x,y,z directions
	utility::vector1<int> xyz_shift_;

	// @brief vector to hold the xyz boundaries of the matrix (not the pose coordinates, but the boundaries that are shifted and potentially scaled within the matrix)
	utility::vector1<core::Size> xyz_bound_;

	// @brief atom type set to derive atom type ints from to pudate the cloud
	core::chemical::AtomTypeSetCOP atom_type_set_;

	// @brief the number of atom types in the atom_type_set_
	core::Size num_atom_types_;

};