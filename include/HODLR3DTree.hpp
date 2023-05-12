#ifndef _HODLR3DTree_hpp__
#define _HODLR3DTree_hpp__
#include <cstdlib>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "definitions.hpp"
#include "kernel.hpp"
#include "ACA.hpp"

class HODLR3DBox {
public:
	int boxNumber;
	int parentNumber;
	int childrenNumbers[8];
	int neighborNumbers[26];//other than self
	std::vector<int> interactionList;
	pts3D center;
	Eigen::VectorXd charges, potential;
	Eigen::MatrixXd L[216], R[216]; //upper bound of neighbors for FMM3D is 27; upper bound of IL for FMM3D is 216=27*8;
	std::vector<int> row_basis[216], col_basis[216];
  std::vector<int> chargeLocations;
	HODLR3DBox ();
};

class HODLR3DTree {
public:
	userkernel* K;
	int nLevels;			//	Number of levels in the tree.
	int N;					//	Number of particles.
	int cubeRootN;					//	Number of particles.
	double L;				//	Semi-length of the simulation box.
	double smallestBoxSize;	//	This is L/2.0^(nLevels).

	std::vector<int> nBoxesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<double> boxRadius;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<std::vector<HODLR3DBox> > tree;	//	The tree storing all the information.

	double ACA_epsilon;
	int nParticlesInLeafAlong1D;
  int nParticlesInLeaf;
  std::vector<double> Nodes1D;
	std::vector<pts3D> Nodes;
  std::vector<pts3D> gridPoints; //all particles in domain
	int TOL_POW;
	double assTime, matVecTime;

// public:
HODLR3DTree(userkernel* K, int cubeRootN, int nLevels, int nParticlesInLeafAlong1D, double L, int TOL_POW);

  void set_Standard_Cheb_Nodes();

	void set_Uniform_Nodes();

  void shift_Nodes(double radius, pts3D center, std::vector<pts3D> &particle_loc);

	void createTree();

	//	Assigns the interactions for child0 of a box
	void assign_Child0_Interaction(int j, int k);

	//	Assigns the interactions for child1 of a box
	void assign_Child1_Interaction(int j, int k);

	//	Assigns the interactions for child2 of a box
	void assign_Child2_Interaction(int j, int k);

	//	Assigns the interactions for child3 of a box
	void assign_Child3_Interaction(int j, int k);

	//	Assigns the interactions for child4 of a box
	void assign_Child4_Interaction(int j, int k);

	//	Assigns the interactions for child5 of a box
	void assign_Child5_Interaction(int j, int k);

	//	Assigns the interactions for child6 of a box
	void assign_Child6_Interaction(int j, int k);

	//	Assigns the interactions for child7 of a box
	void assign_Child7_Interaction(int j, int k);

	//	Assigns the interactions for the children of a box
	void assign_Box_Interactions(int j, int k);

	//	Assigns the interactions for the children all boxes at a given level
	void assign_Level_Interactions(int j);

	//	Assigns the interactions for the children all boxes in the tree
	void assign_Tree_Interactions();

	void assign_Center_Location();

	void assignChargeLocations();

	void assignNonLeafChargeLocations();

	void getUVtTree();

	void getUVtInstance(int j, int k, int ki, Eigen::MatrixXd& L,Eigen::MatrixXd& R, std::vector<int>& row_indices, std::vector<int>& col_indices);

	void assignCharges(Eigen::VectorXd &charges);

	void evaluateFarField();

	void evaluate_NearField();

	void findMemory(double &sum);

	void collectPotential(Eigen::VectorXd &potential);

	void reorder(Eigen::VectorXd &potential);

	int getMaxRank();

};

#endif
