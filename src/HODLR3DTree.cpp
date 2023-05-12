#include "HODLR3DTree.hpp"

	HODLR3DBox::HODLR3DBox () {
		boxNumber		=	-1;
		parentNumber	=	-1;
		for (int l=0; l<8; ++l) {
			childrenNumbers[l]	=	-1;
		}
		for (int l=0; l<26; ++l) {
			neighborNumbers[l]	=	-1;
		}
	}

HODLR3DTree::HODLR3DTree(userkernel* K, int N, int nLevels, int nParticlesInLeafAlong1D, double L, int TOL_POW) {
		this->K					=	K;
		this->N	=	N;
		this->nLevels		=	nLevels;
		this->L					=	L;
    this->nParticlesInLeafAlong1D = nParticlesInLeafAlong1D;
    this->nParticlesInLeaf = nParticlesInLeafAlong1D*nParticlesInLeafAlong1D*nParticlesInLeafAlong1D;
		this->TOL_POW = TOL_POW;
    nBoxesPerLevel.push_back(1);
		boxRadius.push_back(L);
		for (int k=1; k<=nLevels; ++k) {
			nBoxesPerLevel.push_back(8*nBoxesPerLevel[k-1]);
			boxRadius.push_back(0.5*boxRadius[k-1]);
		}
		this->smallestBoxSize	=	boxRadius[nLevels];
		K->a					=	smallestBoxSize;
		this->assTime = 0.0;
		this->matVecTime = 0.0;
	}

  // void HODLR3DTree::set_Standard_Cheb_Nodes() {
	// 	for (int k=0; k<nParticlesInLeafAlong1D; ++k) {
	// 		Nodes1D.push_back(-cos((k+0.5)/nParticlesInLeafAlong1D*PI));
	// 	}
	// 	pts3D temp1;
	// 	for (int j=0; j<nParticlesInLeafAlong1D; ++j) {
	// 		for (int k=0; k<nParticlesInLeafAlong1D; ++k) {
	// 			for (int i=0; i<nParticlesInLeafAlong1D; ++i) {
	// 				temp1.x	=	Nodes1D[k];
	// 				temp1.y	=	Nodes1D[j];
	// 				temp1.z	=	Nodes1D[i];
	// 				Nodes.push_back(temp1);
	// 			}
	// 		}
	// 	}
	// }

	// void HODLR3DTree::set_Uniform_Nodes() {
	// 	for (int k=0; k<cubeRootN; ++k) {
	// 		Nodes1D.push_back(-L+2.0*L*(k+1.0)/(cubeRootN+1.0));
	// 	}
	// 	pts3D temp1;
	// 	for (int j=0; j<cubeRootN; ++j) {
	// 		for (int k=0; k<cubeRootN; ++k) {
	// 			for (int i=0; i<cubeRootN; ++i) {
	// 				temp1.x	=	Nodes1D[k];
	// 				temp1.y	=	Nodes1D[j];
	// 				temp1.z	=	Nodes1D[i];
	// 				K->particles.push_back(temp1);
	// 			}
	// 		}
	// 	}
	// }

  void HODLR3DTree::shift_Nodes(double radius, pts3D center, std::vector<pts3D> &particle_loc) {
    for (size_t i = 0; i < Nodes.size(); i++) {
      pts3D temp;
      temp.x = Nodes[i].x*radius + center.x;
      temp.y = Nodes[i].y*radius + center.y;
			temp.z = Nodes[i].z*radius + center.z;
      particle_loc.push_back(temp);
    }
  }

	void HODLR3DTree::createTree() {
		//	First create root and add to tree
		HODLR3DBox root;
		root.boxNumber		=	0;
		root.parentNumber	=	-1;
		#pragma omp parallel for
		for (int l=0; l<8; ++l) {
			root.childrenNumbers[l]	=	l;
		}
		#pragma omp parallel for
		for (int l=0; l<26; ++l) {
			root.neighborNumbers[l]	=	-1;
		}
		std::vector<HODLR3DBox> rootLevel;
		rootLevel.push_back(root);
		tree.push_back(rootLevel);

		for (int j=1; j<=nLevels; ++j) {
			std::vector<HODLR3DBox> level;
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				HODLR3DBox box;
				box.boxNumber		=	k;
				box.parentNumber	=	k/8;
				for (int l=0; l<8; ++l) {
					box.childrenNumbers[l]	=	8*k+l;
				}
				level.push_back(box);
			}
			tree.push_back(level);
		}
	}

	//	Assigns the interactions for child0 of a box
	void HODLR3DTree::assign_Child0_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	8*k;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC].neighborNumbers[13]	=	8*k+1;
			tree[nL][nC].neighborNumbers[16]	=	8*k+2;
			tree[nL][nC].neighborNumbers[15]	=	8*k+3;
			tree[nL][nC].neighborNumbers[21]	=	8*k+4;
			tree[nL][nC].neighborNumbers[22]	=	8*k+5;
			tree[nL][nC].interactionList.push_back(8*k+6);
			tree[nL][nC].neighborNumbers[24]	=	8*k+7;
		}

		//	Assign children of parent's first neighbor
		{
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
			}
		}

		//	Assign children of parent's third neighbor
		{
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN!=-1) {
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[5];
				// tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[5];
				// tree[nL][nC].neighborNumbers[8]	=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's ninth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[9];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[9]	=	tree[j][nN].childrenNumbers[2];
				// tree[nL][nC].neighborNumbers[17]=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's tenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[10];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[11]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[10] =	tree[j][nN].childrenNumbers[3];
				// tree[nL][nC].neighborNumbers[19] =	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].neighborNumbers[18] =	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
			}
		}

		//	Assign children of parent's eleventh neighbor
		{
			nN	=	tree[j][k].neighborNumbers[11];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's twelveth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[12];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[12]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[14] =	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[20] =	tree[j][nN].childrenNumbers[5];
				// tree[nL][nC].neighborNumbers[23] =	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		for (size_t n = 13; n <= 25; n++) {
			//	Assign children of parent's nth neighbor
			{
				if (n==17 || n==19 || n==23 || n==25) {}
				else {
					nN	=	tree[j][k].neighborNumbers[n];
					if (nN != -1) {
						for (size_t i = 0; i < 8; i++) {
							tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
						}
					}
				}
			}
		}

	}

	//	Assigns the interactions for child1 of a box
	void HODLR3DTree::assign_Child1_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	8*k+1;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC].neighborNumbers[12]	=	8*k+0;
			tree[nL][nC].neighborNumbers[15]	=	8*k+2;
			tree[nL][nC].neighborNumbers[14]	=	8*k+3;
			tree[nL][nC].neighborNumbers[20]	=	8*k+4;
			tree[nL][nC].neighborNumbers[21]	=	8*k+5;
			tree[nL][nC].neighborNumbers[24]	=	8*k+6;
			// tree[nL][nC].neighborNumbers[23]	=	8*k+7;
			tree[nL][nC].interactionList.push_back(8*k+7);
		}

		//	Assign children of parent's first neighbor
		{
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[6];
				// tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
			}
		}

		//	Assign children of parent's third neighbor
		{
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN!=-1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[6];
				// tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[4];
				// tree[nL][nC].neighborNumbers[8]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's ninth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[9];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's tenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[10];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[10]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[9]	  =	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].neighborNumbers[18]	=	tree[j][nN].childrenNumbers[6];
				// tree[nL][nC].neighborNumbers[17]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
			}
		}

		//	Assign children of parent's eleventh neighbor
		{
			nN	=	tree[j][k].neighborNumbers[11];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[11]	=	tree[j][nN].childrenNumbers[3];
				// tree[nL][nC].neighborNumbers[19]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
			}
		}

		//	Assign children of parent's twelveth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[12];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's thirteenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[13];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[13]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[16]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].neighborNumbers[22]	=	tree[j][nN].childrenNumbers[4];
				// tree[nL][nC].neighborNumbers[25]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
			}
		}

		for (size_t n = 14; n <= 25; n++) {
			//	Assign children of parent's nth neighbor
			{
				if (n==17 || n==19 || n==23 || n==25) {}
				else {
					nN	=	tree[j][k].neighborNumbers[n];
					if (nN != -1) {
						for (size_t i = 0; i < 8; i++) {
							tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
						}
					}
				}
			}
		}

	}

	//	Assigns the interactions for child2 of a box
	void HODLR3DTree::assign_Child2_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	8*k+2;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC].neighborNumbers[9]	=	8*k+0;
			tree[nL][nC].neighborNumbers[10]	=	8*k+1;
			tree[nL][nC].neighborNumbers[12]	=	8*k+3;
			// tree[nL][nC].neighborNumbers[17]	=	8*k+4;
			tree[nL][nC].interactionList.push_back(8*k+4);
			tree[nL][nC].neighborNumbers[18]	=	8*k+5;
			tree[nL][nC].neighborNumbers[21]	=	8*k+6;
			tree[nL][nC].neighborNumbers[20]	=	8*k+7;
		}

		//	Assign children of parent's first neighbor
		{
			nN	=	tree[j][k].neighborNumbers[1];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's third neighbor
		{
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN!=-1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's ninth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[9];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's tenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[10];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's eleventh neighbor
		{
			nN	=	tree[j][k].neighborNumbers[11];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's twelveth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[12];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's thirteenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[13];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[11]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[13]	=	tree[j][nN].childrenNumbers[3];
				// tree[nL][nC].neighborNumbers[19]	=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].neighborNumbers[22]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
			}
		}

		//	Assign children of parent's fourteenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[14];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's fifteenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[15];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[14]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[15]	=	tree[j][nN].childrenNumbers[1];
				// tree[nL][nC].neighborNumbers[23]	=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].neighborNumbers[24]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's sixteenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[16];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[16]	=	tree[j][nN].childrenNumbers[0];
				// tree[nL][nC].neighborNumbers[25]	=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		for (size_t n = 17; n <= 25; n++) {
			//	Assign children of parent's nth neighbor
			{
				if (n==17 || n==19 || n==23 || n==25) {}
				else {
					nN	=	tree[j][k].neighborNumbers[n];
					if (nN != -1) {
						for (size_t i = 0; i < 8; i++) {
							tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
						}
					}
				}
			}
		}

	}

	//	Assigns the interactions for child3 of a box
	void HODLR3DTree::assign_Child3_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	8*k+3;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC].neighborNumbers[10]	=	8*k+0;
			tree[nL][nC].neighborNumbers[11]	=	8*k+1;
			tree[nL][nC].neighborNumbers[13]	=	8*k+2;
			tree[nL][nC].neighborNumbers[18]	=	8*k+4;
			// tree[nL][nC].neighborNumbers[19]	=	8*k+5;
			tree[nL][nC].interactionList.push_back(8*k+5);
			tree[nL][nC].neighborNumbers[22]	=	8*k+6;
			tree[nL][nC].neighborNumbers[21]	=	8*k+7;
		}

		for (size_t n = 1; n <= 1; n++) {
			//	Assign children of parent's nth neighbor
			{
				nN	=	tree[j][k].neighborNumbers[n];
				if (nN != -1) {
					for (size_t i = 0; i < 8; i++) {
						tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
					}
				}
			}
		}

		//	Assign children of parent's third neighbor
		{
			nN	=	tree[j][k].neighborNumbers[3];
			if (nN!=-1) {
				// tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[4];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[4];
				// tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[5];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			nN	=	tree[j][k].neighborNumbers[7];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[4];
				// tree[nL][nC].neighborNumbers[8]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		for (size_t n = 9; n <= 11; n++) {
			//	Assign children of parent's nth neighbor
			{
				nN	=	tree[j][k].neighborNumbers[n];
				if (nN != -1) {
					for (size_t i = 0; i < 8; i++) {
						tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
					}
				}
			}
		}

		//	Assign children of parent's twelveth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[12];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[9]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[12]	=	tree[j][nN].childrenNumbers[2];
				// tree[nL][nC].neighborNumbers[17]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].neighborNumbers[20]	=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's thirteenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[13];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's fourteenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[14];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[14]	=	tree[j][nN].childrenNumbers[1];
				// tree[nL][nC].neighborNumbers[23]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's fifteenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[15];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[15]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[16]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[24]	=	tree[j][nN].childrenNumbers[4];
				// tree[nL][nC].neighborNumbers[25]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		for (size_t n = 16; n <= 25; n++) {
			//	Assign children of parent's nth neighbor
			{
				if (n==17 || n==19 || n==23 || n==25) {}
				else {
					nN	=	tree[j][k].neighborNumbers[n];
					if (nN != -1) {
						for (size_t i = 0; i < 8; i++) {
							tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
						}
					}
				}
			}
		}

	}

	//	Assigns the interactions for child4 of a box
	void HODLR3DTree::assign_Child4_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	8*k+4;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC].neighborNumbers[4]	=	8*k+0;
			tree[nL][nC].neighborNumbers[5]	=	8*k+1;
			// tree[nL][nC].neighborNumbers[8]	=	8*k+2;
			tree[nL][nC].interactionList.push_back(8*k+2);
			tree[nL][nC].neighborNumbers[7]	=	8*k+3;
			tree[nL][nC].neighborNumbers[13]	=	8*k+5;
			tree[nL][nC].neighborNumbers[16]	=	8*k+6;
			tree[nL][nC].neighborNumbers[15]	=	8*k+7;
		}

		for (size_t n = 1; n <= 8; n++) {
			//	Assign children of parent's nth neighbor
			{
				if (n==2 || n==6 || n==8) {}
				else {
					nN	=	tree[j][k].neighborNumbers[n];
					if (nN != -1) {
						for (size_t i = 0; i < 8; i++) {
							tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
						}
					}
				}
			}
		}

		//	Assign children of parent's nineth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[9];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].neighborNumbers[9]	=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's tenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[10];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].neighborNumbers[11]	=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].neighborNumbers[10]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
			}
		}

		//	Assign children of parent's eleventh neighbor
		{
			nN	=	tree[j][k].neighborNumbers[11];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's twelveth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[12];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[1];
				// tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].neighborNumbers[12]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].neighborNumbers[14]	=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		for (size_t n = 13; n <= 16; n++) {
			//	Assign children of parent's nth neighbor
			{
				nN	=	tree[j][k].neighborNumbers[n];
				if (nN != -1) {
					for (size_t i = 0; i < 8; i++) {
						tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
					}
				}
			}
		}

		//	Assign children of parent's eighteenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[18];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[19]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].neighborNumbers[18]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's twentieth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[20];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[20]	=	tree[j][nN].childrenNumbers[1];
				// tree[nL][nC].neighborNumbers[23]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's twentyfirst neighbor
		{
			nN	=	tree[j][k].neighborNumbers[21];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[21]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[22]	=	tree[j][nN].childrenNumbers[1];
				// tree[nL][nC].neighborNumbers[25]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].neighborNumbers[24]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}


		for (size_t n = 22; n <= 25; n++) {
			//	Assign children of parent's nth neighbor
			{
				if (n==23 || n==25) {}
				else {
					nN	=	tree[j][k].neighborNumbers[n];
					if (nN != -1) {
						for (size_t i = 0; i < 8; i++) {
							tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
						}
					}
				}
			}
		}

	}

	//	Assigns the interactions for child5 of a box
	void HODLR3DTree::assign_Child5_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	8*k+5;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC].neighborNumbers[3]	=	8*k+0;
			tree[nL][nC].neighborNumbers[4]	=	8*k+1;
			tree[nL][nC].neighborNumbers[7]	=	8*k+2;
			// tree[nL][nC].neighborNumbers[6]	=	8*k+3;
			tree[nL][nC].interactionList.push_back(8*k+3);
			tree[nL][nC].neighborNumbers[12]	=	8*k+4;
			tree[nL][nC].neighborNumbers[15]	=	8*k+6;
			tree[nL][nC].neighborNumbers[14]	=	8*k+7;
		}

		for (size_t n = 1; n <= 8; n++) {
			//	Assign children of parent's nth neighbor
			{
				if (n==2 ||n==6 || n==8) {}
				else {
					nN	=	tree[j][k].neighborNumbers[n];
					if (nN != -1) {
						for (size_t i = 0; i < 8; i++) {
							tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
						}
					}
				}
			}
		}

		//	Assign children of parent's nineth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[9];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's tenth neighbor
		{
			nN	=	tree[j][k].neighborNumbers[10];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[2];
				// tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].neighborNumbers[10]	=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].neighborNumbers[9]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
			}
		}

		//	Assign children of parent's eleventh neighbor
		{
			nN	=	tree[j][k].neighborNumbers[11];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].neighborNumbers[11]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
			}
		}

		//	Assign children of parent's 12th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[12];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's 13 neighbor
		{
			nN	=	tree[j][k].neighborNumbers[13];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				// tree[nL][nC].neighborNumbers[8]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].neighborNumbers[13]	=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].neighborNumbers[16]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
			}
		}

		for (size_t n = 14; n <= 16; n++) {
			//	Assign children of parent's nth neighbor
			{
				nN	=	tree[j][k].neighborNumbers[n];
				if (nN != -1) {
					for (size_t i = 0; i < 8; i++) {
						tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
					}
				}
			}
		}

		//	Assign children of parent's 18th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[18];
			if (nN!=-1) {
				tree[nL][nC].neighborNumbers[18]	=	tree[j][nN].childrenNumbers[2];
				// tree[nL][nC].neighborNumbers[17]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 20th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[20];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's 21st neighbor
		{
			nN	=	tree[j][k].neighborNumbers[21];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[20]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[21]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[24]	=	tree[j][nN].childrenNumbers[2];
				// tree[nL][nC].neighborNumbers[23]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 22nd neighbor
		{
			nN	=	tree[j][k].neighborNumbers[22];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[22]	=	tree[j][nN].childrenNumbers[0];
				// tree[nL][nC].neighborNumbers[25]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		for (size_t n = 24; n <= 24; n++) {
			//	Assign children of parent's nth neighbor
			{
				nN	=	tree[j][k].neighborNumbers[n];
				if (nN != -1) {
					for (size_t i = 0; i < 8; i++) {
						tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
					}
				}
			}
		}

	}

	//	Assigns the interactions for child6 of a box
	void HODLR3DTree::assign_Child6_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	8*k+6;
		int nN;

		//	Assign siblings
		{
			// tree[nL][nC].neighborNumbers[0]	=	8*k+0;
			tree[nL][nC].interactionList.push_back(8*k+0);
			tree[nL][nC].neighborNumbers[1]	=	8*k+1;
			tree[nL][nC].neighborNumbers[4]	=	8*k+2;
			tree[nL][nC].neighborNumbers[3]	=	8*k+3;
			tree[nL][nC].neighborNumbers[9]	=	8*k+4;
			tree[nL][nC].neighborNumbers[10]	=	8*k+5;
			tree[nL][nC].neighborNumbers[12]	=	8*k+7;
		}

		for (size_t n = 1; n <= 12; n++) {
			//	Assign children of parent's nth neighbor
			{
				if (n==2 || n==6 || n==8) {}
				else {
					nN	=	tree[j][k].neighborNumbers[n];
					if (nN != -1) {
						for (size_t i = 0; i < 8; i++) {
							tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
						}
					}
				}
			}
		}

		//	Assign children of parent's 13th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[13];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].neighborNumbers[11]	=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].neighborNumbers[13]	=	tree[j][nN].childrenNumbers[7];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
			}
		}

		//	Assign children of parent's 14th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[14];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's 15th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[15];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[14]=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].neighborNumbers[15]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 16th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[16];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[8]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].neighborNumbers[16]	=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		for (size_t n = 18; n <= 20; n++) {
			//	Assign children of parent's nth neighbor
			{
				if (n == 19) {}
				else {
					nN	=	tree[j][k].neighborNumbers[n];
					if (nN != -1) {
						for (size_t i = 0; i < 8; i++) {
							tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
						}
					}
				}
			}
		}

		//	Assign children of parent's 21st neighbor
		{
			nN	=	tree[j][k].neighborNumbers[21];
			if (nN!=-1) {
				// tree[nL][nC].neighborNumbers[17]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].neighborNumbers[18]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[21]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[20]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 22nd neighbor
		{
			nN	=	tree[j][k].neighborNumbers[22];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[19]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].neighborNumbers[22]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 24th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[24];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[23]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].neighborNumbers[24]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}
	}

	//	Assigns the interactions for child7 of a box
	void HODLR3DTree::assign_Child7_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	8*k+7;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC].neighborNumbers[1]	=	8*k+0;
			// tree[nL][nC].neighborNumbers[2]	=	8*k+1;
			tree[nL][nC].interactionList.push_back(8*k+1);
			tree[nL][nC].neighborNumbers[5]	=	8*k+2;
			tree[nL][nC].neighborNumbers[4]	=	8*k+3;
			tree[nL][nC].neighborNumbers[10]=	8*k+4;
			tree[nL][nC].neighborNumbers[11]=	8*k+5;
			tree[nL][nC].neighborNumbers[13]=	8*k+6;
		}

		for (size_t n = 1; n <= 11; n++) {
			//	Assign children of parent's nth neighbor
			{
				if (n==2 || n==6 || n==8) {}
				else {
					nN	=	tree[j][k].neighborNumbers[n];
					if (nN != -1) {
						for (size_t i = 0; i < 8; i++) {
							tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
						}
					}
				}
			}
		}

		//	Assign children of parent's 12th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[12];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[9]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].neighborNumbers[12]	=	tree[j][nN].childrenNumbers[6];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 13th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[13];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's 14th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[14];
			if (nN != -1) {
				// tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].neighborNumbers[14]	=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 15th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[15];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[0];
				// tree[nL][nC].neighborNumbers[8] =	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].neighborNumbers[15]=	tree[j][nN].childrenNumbers[4];
				tree[nL][nC].neighborNumbers[16]=	tree[j][nN].childrenNumbers[5];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		for (size_t n = 16; n <= 18; n++) {
			//	Assign children of parent's nth neighbor
			{
				if (n==17) {}
				else {
					nN	=	tree[j][k].neighborNumbers[n];
					if (nN != -1) {
						for (size_t i = 0; i < 8; i++) {
							tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
						}
					}
				}
			}
		}

		//	Assign children of parent's 20th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[20];
			if (nN!=-1) {
				// tree[nL][nC].neighborNumbers[17]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].neighborNumbers[20]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[0]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 21st neighbor
		{
			nN	=	tree[j][k].neighborNumbers[21];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[18]	=	tree[j][nN].childrenNumbers[0];
				// tree[nL][nC].neighborNumbers[19]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].neighborNumbers[22]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[21]	=	tree[j][nN].childrenNumbers[3];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 22nd neighbor
		{
			nN	=	tree[j][k].neighborNumbers[22];
			if (nN != -1) {
				for (size_t i = 0; i < 8; i++) {
					tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's 24th neighbor
		{
			nN	=	tree[j][k].neighborNumbers[24];
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[24]	=	tree[j][nN].childrenNumbers[0];
				// tree[nL][nC].neighborNumbers[25]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[1]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[2]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[3]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[4]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[5]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[6]);
				tree[nL][nC].interactionList.push_back(tree[j][nN].childrenNumbers[7]);
			}
		}
	}

	//	Assigns the interactions for the children of a box
	void HODLR3DTree::assign_Box_Interactions(int j, int k) {
		assign_Child0_Interaction(j,k);
		assign_Child1_Interaction(j,k);
		assign_Child2_Interaction(j,k);
		assign_Child3_Interaction(j,k);
		assign_Child4_Interaction(j,k);
		assign_Child5_Interaction(j,k);
		assign_Child6_Interaction(j,k);
		assign_Child7_Interaction(j,k);
	}

	//	Assigns the interactions for the children all boxes at a given level
	void HODLR3DTree::assign_Level_Interactions(int j) {
		#pragma omp parallel for
		for (int k=0; k<nBoxesPerLevel[j]; ++k) {
			assign_Box_Interactions(j,k);
		}
	}

	//	Assigns the interactions for the children all boxes in the tree
	void HODLR3DTree::assign_Tree_Interactions() {
		for (int j=0; j<nLevels; ++j) {
			assign_Level_Interactions(j);
		}
	}

	void HODLR3DTree::assign_Center_Location() {
		int J;
		tree[0][0].center.x	=	0.0;
		tree[0][0].center.y	=	0.0;
		tree[0][0].center.z	=	0.0;
		for (int j=0; j<nLevels; ++j) {
			J	=	j+1;
			double shift	=	0.5*boxRadius[j];
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				tree[J][8*k].center.x		=	tree[j][k].center.x-shift;
				tree[J][8*k+1].center.x	=	tree[j][k].center.x+shift;
				tree[J][8*k+2].center.x	=	tree[j][k].center.x+shift;
				tree[J][8*k+3].center.x	=	tree[j][k].center.x-shift;
				tree[J][8*k+4].center.x	=	tree[j][k].center.x-shift;
				tree[J][8*k+5].center.x	=	tree[j][k].center.x+shift;
				tree[J][8*k+6].center.x	=	tree[j][k].center.x+shift;
				tree[J][8*k+7].center.x	=	tree[j][k].center.x-shift;

				tree[J][8*k].center.y		=	tree[j][k].center.y-shift;
				tree[J][8*k+1].center.y	=	tree[j][k].center.y-shift;
				tree[J][8*k+2].center.y	=	tree[j][k].center.y+shift;
				tree[J][8*k+3].center.y	=	tree[j][k].center.y+shift;
				tree[J][8*k+4].center.y	=	tree[j][k].center.y-shift;
				tree[J][8*k+5].center.y	=	tree[j][k].center.y-shift;
				tree[J][8*k+6].center.y	=	tree[j][k].center.y+shift;
				tree[J][8*k+7].center.y	=	tree[j][k].center.y+shift;

				tree[J][8*k].center.z		=	tree[j][k].center.z-shift;
				tree[J][8*k+1].center.z	=	tree[j][k].center.z-shift;
				tree[J][8*k+2].center.z	=	tree[j][k].center.z-shift;
				tree[J][8*k+3].center.z	=	tree[j][k].center.z-shift;
				tree[J][8*k+4].center.z	=	tree[j][k].center.z+shift;
				tree[J][8*k+5].center.z	=	tree[j][k].center.z+shift;
				tree[J][8*k+6].center.z	=	tree[j][k].center.z+shift;
				tree[J][8*k+7].center.z	=	tree[j][k].center.z+shift;
			}
		}
	}

	void HODLR3DTree::assignChargeLocations() {
		// K->particles = Nodes;//object of base class FMM_Matrix
		for (size_t i = 0; i < N; i++) {
			tree[0][0].chargeLocations.push_back(i);
		}
		for (size_t j = 0; j < nLevels; j++) { //assign particles to its children
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				int J = j+1;
				int Kp = 8*k;
				for (size_t i = 0; i < tree[j][k].chargeLocations.size(); i++) {
					int index = tree[j][k].chargeLocations[i];
					if (K->particles[index].z <= tree[j][k].center.z) { //children 0,1,2,3
						if (K->particles[index].x <= tree[j][k].center.x) { //children 0,3
							if (K->particles[index].y <= tree[j][k].center.y) { //child 0
								tree[J][Kp].chargeLocations.push_back(index);
							}
							else { //child 3
								tree[J][Kp+3].chargeLocations.push_back(index);
							}
						}
						else { //children 1,2
							if (K->particles[index].y <= tree[j][k].center.y) { //child 1
								tree[J][Kp+1].chargeLocations.push_back(index);
							}
							else { //child 2
								tree[J][Kp+2].chargeLocations.push_back(index);
							}
						}
					}
					else {//children 4,5,6,7
						if (K->particles[index].x <= tree[j][k].center.x) { //children 4,7
							if (K->particles[index].y <= tree[j][k].center.y) { //child 4
								tree[J][Kp+4].chargeLocations.push_back(index);
							}
							else { //child 7
								tree[J][Kp+7].chargeLocations.push_back(index);
							}
						}
						else { //children 5,6
							if (K->particles[index].y <= tree[j][k].center.y) { //child 5
								tree[J][Kp+5].chargeLocations.push_back(index);
							}
							else { //child 6
								tree[J][Kp+6].chargeLocations.push_back(index);
							}
						}
					}
				}
			}
		}
	}

	void HODLR3DTree::assignNonLeafChargeLocations() {
		for (int j = nLevels-1; j >= 1; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				tree[j][k].chargeLocations.clear();
				for (size_t c = 0; c < 8; c++) {
					tree[j][k].chargeLocations.insert(tree[j][k].chargeLocations.end(), tree[j+1][8*k+c].chargeLocations.begin(), tree[j+1][8*k+c].chargeLocations.end());
				}
			}
		}
	}

	void HODLR3DTree::getUVtTree() {
		for (size_t j = 1; j <= nLevels; j++) {
			#pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				for (size_t i = 0; i < tree[j][k].interactionList.size(); i++) {
					int ki = tree[j][k].interactionList[i];
					getUVtInstance(j,k,ki,tree[j][k].L[i],tree[j][k].R[i], tree[j][k].row_basis[i],tree[j][k].col_basis[i]);
				}
			}
		}
	}

	void HODLR3DTree::getUVtInstance(int j, int k, int ki, Eigen::MatrixXd& L,Eigen::MatrixXd& R, std::vector<int>& row_indices, std::vector<int>& col_indices) {
		int computed_rank;
		std::vector<int> row_indices_local,col_indices_local;
		LowRank* LR		=	new LowRank(K, TOL_POW, tree[j][k].chargeLocations, tree[j][ki].chargeLocations);
		LR->ACA_only_nodesCUR(row_indices_local, col_indices_local, computed_rank, L, R);
		for (int r = 0; r < computed_rank; r++) {
			row_indices.push_back(tree[j][k].chargeLocations[row_indices_local[r]]);
		}
		for (int c = 0; c < computed_rank; c++) {
			col_indices.push_back(tree[j][ki].chargeLocations[col_indices_local[c]]);
		}
		delete LR;
	}

	void HODLR3DTree::assignCharges(Eigen::VectorXd &charges) {
		int start = 0;
		for (size_t j = 1; j <= nLevels; j++) {
			#pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				tree[j][k].charges = Eigen::VectorXd::Zero(tree[j][k].chargeLocations.size());
				for (size_t i = 0; i < tree[j][k].chargeLocations.size(); i++) {
					tree[j][k].charges[i] = charges[tree[j][k].chargeLocations[i]];
				}
			}
		}
	}

	void HODLR3DTree::evaluateFarField() {
		for (size_t j = 1; j <= nLevels; j++) {
			double assTimeLevel = 0.0;
			double matVecTimeLevel = 0.0;

			#pragma omp parallel
			{
				double start, end;
				double assTimeThread = 0.0;
				double matVecTimeThread = 0.0;

				#pragma omp for nowait
				for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
					tree[j][k].potential = Eigen::VectorXd::Zero(tree[j][k].charges.size());
					for (size_t i = 0; i < tree[j][k].interactionList.size(); i++) {
						if (tree[j][k].col_basis[i].size() > 0) {
							int ki = tree[j][k].interactionList[i];
							start = omp_get_wtime();
							Eigen::MatrixXd Ac = K->getMatrix(tree[j][k].chargeLocations,tree[j][k].col_basis[i]);
							Eigen::MatrixXd Ar = K->getMatrix(tree[j][k].row_basis[i],tree[j][ki].chargeLocations);
							end = omp_get_wtime();
							assTimeThread += end-start;

							start = omp_get_wtime();
							Eigen::VectorXd t0 = Ar*tree[j][ki].charges;
							Eigen::VectorXd t1 = tree[j][k].L[i].triangularView<Eigen::Lower>().solve(t0);
							Eigen::VectorXd t2 = tree[j][k].R[i].triangularView<Eigen::Upper>().solve(t1);
							tree[j][k].potential += Ac*t2;
							end = omp_get_wtime();
							matVecTimeThread += end-start;
						}
					}
				}
				#pragma omp critical
				{
					if (matVecTimeLevel < matVecTimeThread) {
						matVecTimeLevel = matVecTimeThread;
					}
					if (assTimeLevel < assTimeThread) {
						assTimeLevel = assTimeThread;
					}
				}
			}
			assTime += assTimeLevel;
			matVecTime += matVecTimeLevel;
		}
	}

	void HODLR3DTree::evaluate_NearField() { // evaluating at chargeLocations
		// #pragma omp parallel for reduction(+:assTime,matVecTime)
		double assTimeLevel = 0.0;
		double matVecTimeLevel = 0.0;
		#pragma omp parallel
		{
			double start, end;
			double assTimeThread = 0.0;
			double matVecTimeThread = 0.0;
			#pragma omp for nowait
			for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
				for (size_t n = 0; n < 26; n++) {
					int nn = tree[nLevels][k].neighborNumbers[n];
					if(nn != -1) {
						start = omp_get_wtime();
						Eigen::MatrixXd R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][nn].chargeLocations);
						end = omp_get_wtime();
						assTimeThread += end-start;

						start = omp_get_wtime();
						tree[nLevels][k].potential += R*tree[nLevels][nn].charges;
						end = omp_get_wtime();
						matVecTimeThread += end-start;
					}
				}
				start = omp_get_wtime();
				Eigen::MatrixXd R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][k].chargeLocations);
				end = omp_get_wtime();
				assTimeThread += end-start;

				start = omp_get_wtime();
				tree[nLevels][k].potential += R*tree[nLevels][k].charges; //self Interaction
				end = omp_get_wtime();
				matVecTimeThread += end-start;
			}
			#pragma omp critical
			{
				if (matVecTimeLevel < matVecTimeThread) {
					matVecTimeLevel = matVecTimeThread;
				}
				if (assTimeLevel < assTimeThread) {
					assTimeLevel = assTimeThread;
				}
			}
		}
		assTime += assTimeLevel;
		matVecTime += matVecTimeLevel;
	}

	void HODLR3DTree::findMemory(double &sum) {
		sum = 0;
		for (size_t j = 1; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				for (size_t i = 0; i < tree[j][k].interactionList.size(); i++) {
					int ki = tree[j][k].interactionList[i];
					sum += tree[j][k].L[i].cols()*tree[j][k].chargeLocations.size();
					sum += tree[j][k].L[i].cols()*tree[j][k].L[i].cols();
					sum += tree[j][ki].chargeLocations.size()*tree[j][k].L[i].cols();
				}
			}
		}
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			sum += tree[nLevels][k].chargeLocations.size()*tree[nLevels][k].chargeLocations.size(); //self
			for (size_t n = 0; n < 26; n++) {
				int nn = tree[nLevels][k].neighborNumbers[n];
				if(nn != -1) {
					sum += tree[nLevels][k].chargeLocations.size()*tree[nLevels][nn].chargeLocations.size();
				}
			}
		}
	}


	void HODLR3DTree::collectPotential(Eigen::VectorXd &potential) {
		potential = Eigen::VectorXd::Zero(N);
		for (size_t j = 1; j <= nLevels; j++) {
			int start = 0;
			Eigen::VectorXd potentialTemp = Eigen::VectorXd::Zero(N);
			// #pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) { //using the fact that all the leaves have same number of particles
				potentialTemp.segment(start, tree[j][k].potential.size()) = tree[j][k].potential;
				start += tree[j][k].potential.size();
			}
			potential = potential + potentialTemp;
		}
	}

	void HODLR3DTree::reorder(Eigen::VectorXd &potential) {
		Eigen::VectorXd potentialTemp = potential;
		int start = 0;
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			for (size_t i = 0; i < tree[nLevels][k].chargeLocations.size(); i++) {
				int index = tree[nLevels][k].chargeLocations[i];
				potential(index) = potentialTemp(start);
				start++;
			}
		}
	}

	int HODLR3DTree::getMaxRank() {
		int maxRank = 0;
		for (size_t j = 1; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				for (size_t i = 0; i < tree[j][k].interactionList.size(); i++) {
					if (maxRank < tree[j][k].L[i].cols())
						maxRank = tree[j][k].L[i].cols();
				}
			}
		}
		return maxRank;
	}
