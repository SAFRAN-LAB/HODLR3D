#ifndef _FMM3DTreeRAMeff2_HPP__
#define _FMM3DTreeRAMeff2_HPP__
class HODLRBox {
public:
	int boxNumber;
	int parentNumber;
	int childrenNumbers[8];

	HODLRBox () {
		boxNumber		=	-1;
		parentNumber	=	-1;
		for (int l=0; l<8; ++l) {
			childrenNumbers[l]	=	-1;
		}
	}

	pts3D center;
	Eigen::VectorXd charges, potential;
	// Eigen::MatrixXd L[216], R[216]; //upper bound of neighbors for FMM3D is 27; upper bound of IL for FMM3D is 216=27*8;
	Eigen::MatrixXd *L, *R; //upper bound of neighbors for FMM3D is 27; upper bound of IL for FMM3D is 216=27*8;
	std::vector<int> row_basis[216], col_basis[216];
	// std::vector<Eigen::MatrixXd> U;
	// std::vector<Eigen::MatrixXd> Vt;
  std::vector<int> chargeLocations;
};

template <typename kerneltype>
class HODLRTree {
public:
	kerneltype* K;
	int nLevels;			//	Number of levels in the tree.
	int N;					//	Number of particles.
	int cubeRootN;					//	Number of particles.
	double L;				//	Semi-length of the simulation box.
	double smallestBoxSize;	//	This is L/2.0^(nLevels).

	std::vector<int> nBoxesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<double> boxRadius;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<std::vector<HODLRBox> > tree;	//	The tree storing all the information.

	double ACA_epsilon;
	int nParticlesInLeafAlong1D;
  int nParticlesInLeaf;
  std::vector<double> Nodes1D;
	// std::vector<pts3D> Nodes;
  std::vector<pts3D> gridPoints; //all particles in domain
	int TOL_POW;
	double assTime, matVecTime;

// public:
	HODLRTree(kerneltype* K, int cubeRootN, int nLevels, int nParticlesInLeafAlong1D, double L, int TOL_POW) {
		this->K					=	K;
		this->cubeRootN	=	cubeRootN;
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
		this->N					=	cubeRootN*cubeRootN*cubeRootN;
		this->assTime = 0.0;
		this->matVecTime = 0.0;
	}

  // void set_Standard_Cheb_Nodes() {
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
	//
	void set_Uniform_Nodes() {
		for (int k=0; k<cubeRootN; ++k) {
			Nodes1D.push_back(-L+2.0*L*(k+1.0)/(cubeRootN+1.0));
		}
		pts3D temp1;
		for (int j=0; j<cubeRootN; ++j) {
			for (int k=0; k<cubeRootN; ++k) {
				for (int i=0; i<cubeRootN; ++i) {
					temp1.x	=	Nodes1D[k];
					temp1.y	=	Nodes1D[j];
					temp1.z	=	Nodes1D[i];
					K->particles.push_back(temp1);
				}
			}
		}
	}

  // void shift_Nodes(double radius, pts3D center, std::vector<pts3D> &particle_loc) {
  //   for (size_t i = 0; i < Nodes.size(); i++) {
  //     pts3D temp;
  //     temp.x = Nodes[i].x*radius + center.x;
  //     temp.y = Nodes[i].y*radius + center.y;
	// 		temp.z = Nodes[i].z*radius + center.z;
  //     particle_loc.push_back(temp);
  //   }
  // }

	void createTree() {
		//	First create root and add to tree
		HODLRBox root;
		root.boxNumber		=	0;
		root.parentNumber	=	-1;
		root.L = new Eigen::MatrixXd[216];
		root.R = new Eigen::MatrixXd[216];
		#pragma omp parallel for
		for (int l=0; l<8; ++l) {
			root.childrenNumbers[l]	=	l;
		}
		std::vector<HODLRBox> rootLevel;
		rootLevel.push_back(root);
		tree.push_back(rootLevel);

		for (int j=1; j<=nLevels; ++j) {
			std::vector<HODLRBox> level;
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				HODLRBox box;
				box.boxNumber		=	k;
				box.parentNumber	=	k/8;
				box.L = new Eigen::MatrixXd[216];
				box.R = new Eigen::MatrixXd[216];
				for (int l=0; l<8; ++l) {
					box.childrenNumbers[l]	=	8*k+l;
				}
				level.push_back(box);
			}
			tree.push_back(level);
		}
	}

	void assign_Center_Location() {
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

	void assignChargeLocations() {
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

	void check() {
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			for (size_t i = 0; i < tree[nLevels][k].chargeLocations.size(); i++) {
				std::cout << tree[nLevels][k].chargeLocations[i] << "," << std::endl;
			}
			std::cout << "-------" << std::endl;
		}
	}

  // void assignLeafChargeLocations() {
	// 	for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
  //     int startIndex = gridPoints.size();
  //     for (size_t i = 0; i < nParticlesInLeaf; i++) {
  //       tree[nLevels][k].chargeLocations.push_back(startIndex+i);
  //     }
  //     shift_Nodes(boxRadius[nLevels], tree[nLevels][k].center, gridPoints);
  //   }
  //   K->particles_X = gridPoints;//object of base class FMM_Matrix
  //   K->particles_Y = gridPoints;
  // }

	void assignNonLeafChargeLocations() {
		for (int j = nLevels-1; j >= 1; j--) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				tree[j][k].chargeLocations.clear();
				for (size_t c = 0; c < 8; c++) {
					tree[j][k].chargeLocations.insert(tree[j][k].chargeLocations.end(), tree[j+1][8*k+c].chargeLocations.begin(), tree[j+1][8*k+c].chargeLocations.end());
				}
				// std::cout << "tree[j][k].charges: " << tree[j][k].charges.size() << std::endl;
			}
		}
	}

	void getUVtTree() {
		// #pragma omp parallel for
		for (size_t j = 1; j <= nLevels; j++) {
			#pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				int Kp = k/8;
				// #pragma omp parallel for
				for (size_t c = 0; c < 8; c++) {
					int ki = 8*Kp+c;
					if (ki != k) {
						getUVtInstance(j,k,ki,tree[j][k].L[c],tree[j][k].R[c], tree[j][k].row_basis[c],tree[j][k].col_basis[c]);
					}
				}
			}
		}
	}

	void getUVtInstance(int j, int k, int ki, Eigen::MatrixXd& L,Eigen::MatrixXd& R, std::vector<int>& row_indices, std::vector<int>& col_indices) {
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

	void assignCharges(Eigen::VectorXd &charges) {
		int start = 0;
		// std::cout << "charges" << std::endl;
		for (size_t j = 1; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				tree[j][k].charges = Eigen::VectorXd::Zero(tree[j][k].chargeLocations.size());
				for (size_t i = 0; i < tree[j][k].chargeLocations.size(); i++) {
					tree[j][k].charges[i] = charges[tree[j][k].chargeLocations[i]];
					// std::cout << tree[j][k].chargeLocations[i] << std::endl;
				}
				// std::cout << "---------" << std::endl;
			}
		}
	}

	// void assignLeafCharges(Eigen::VectorXd &charges) {
	// 	int start = 0;
	// 	for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
	// 		tree[nLevels][k].charges	=	charges.segment(start, nParticlesInLeaf);
	// 		start += nParticlesInLeaf;
	// 	}
	// }

	// void assignNonLeafCharges() {
	// 	for (int j = nLevels-1; j >= 1; j--) {
	// 		for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
	// 			int sum = 0;
	// 			for (size_t c = 0; c < 8; c++) {
	// 				sum += tree[j+1][8*k+c].charges.size();
	// 			}
	// 			tree[j][k].charges	=	Eigen::VectorXd(sum);
	// 			int start = 0;
	// 			for (size_t c = 0; c < 8; c++) {
	// 				tree[j][k].charges.segment(start, tree[j+1][8*k+c].charges.size()) = tree[j+1][8*k+c].charges;
	// 				start += tree[j+1][8*k+c].charges.size();
	// 			}
	// 			// std::cout << "tree[j][k].charges: " << tree[j][k].charges.size() << std::endl;
	// 		}
	// 	}
	// }

	// void fullMatrix() {
	// 	std::cout << K->getMatrix(0,0,N,N) << std::endl;
	// }


	void evaluateFarField() {
		// #pragma omp parallel for
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
					int kp = k/8;
					for (size_t c = 0; c < 8; c++) {
						int ki = 8*kp+c;
						if (ki != k && tree[j][k].col_basis[c].size() > 0) {
							start = omp_get_wtime();
							// Eigen::MatrixXd R = K->getMatrix(tree[j][k].chargeLocations,tree[j][ki].chargeLocations);
							// tree[j][k].potential += R*tree[j][ki].charges;

							Eigen::MatrixXd Ac = K->getMatrix(tree[j][k].chargeLocations,tree[j][k].col_basis[c]);
							Eigen::MatrixXd Ar = K->getMatrix(tree[j][k].row_basis[c],tree[j][ki].chargeLocations);
							end = omp_get_wtime();
							assTimeThread += end-start;

							start = omp_get_wtime();
							Eigen::VectorXd t0 = Ar*tree[j][ki].charges;
							Eigen::VectorXd t1 = tree[j][k].L[c].triangularView<Eigen::Lower>().solve(t0);
							Eigen::VectorXd t2 = tree[j][k].R[c].triangularView<Eigen::Upper>().solve(t1);
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

	void evaluate_NearField() { // evaluating at chargeLocations
		double assTimeLevel = 0.0;
		double matVecTimeLevel = 0.0;
		#pragma omp parallel
		{
			double start, end;
			double assTimeThread = 0.0;
			double matVecTimeThread = 0.0;
			#pragma omp for nowait
			for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
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

	void collectPotential(Eigen::VectorXd &potential) {
		potential = Eigen::VectorXd::Zero(N);
		// #pragma omp parallel for
		for (size_t j = 1; j <= nLevels; j++) {
			int start = 0;
			Eigen::VectorXd potentialTemp = Eigen::VectorXd::Zero(N);
			// #pragma omp parallel for
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) { //using the fact that all the leaves have same number of particles
				// potentialTemp.segment(pow(8,nLevels-j)*nParticlesInLeaf*k, tree[j][k].potential.size()) = tree[j][k].potential;
				potentialTemp.segment(start, tree[j][k].potential.size()) = tree[j][k].potential;
				start += tree[j][k].potential.size();
			}
			// if(j==nLevels)
			// reorder(potentialTemp);
			potential = potential + potentialTemp;
		}
	}

	void reorder(Eigen::VectorXd &potential) {
		Eigen::VectorXd potentialTemp = potential;
		int start = 0;
		// std::cout << "index: " << std::endl;
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			for (size_t i = 0; i < tree[nLevels][k].chargeLocations.size(); i++) {
				int index = tree[nLevels][k].chargeLocations[i];
				// std::cout << index << std::endl;
				potential(index) = potentialTemp(start);
				start++;
			}
			// std::cout << "---------" << std::endl;
		}
	}

	int getMaxRank() {
		int maxRank = 0;
		for (size_t j = 1; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				int kp = k/8;
				for (size_t c = 0; c < 8; c++) {
					int ki = 8*kp+c;
					if (ki != k) {
						if (maxRank < tree[j][k].L[c].cols())
							maxRank = tree[j][k].L[c].cols();
					}
				}
			}
		}
		return maxRank;
	}

	void findMemory(double &sum) {
		sum = 0;
		for (size_t j = 1; j <= nLevels; j++) {
			for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
				int kp = k/8;
				for (size_t c = 0; c < 8; c++) {
					int ki = 8*kp+c;
					if (ki != k) {
						//sum += tree[j][k].L[c].cols()*tree[j][k].chargeLocations.size();
						//sum += tree[j][k].L[c].cols()*tree[j][k].L[c].cols();
						//sum += tree[j][ki].chargeLocations.size()*tree[j][k].L[c].cols();
						sum += tree[j][k].R[c].rows()*tree[j][k].R[c].cols();
						sum += tree[j][k].L[c].rows()*tree[j][k].L[c].cols();
						sum += tree[j][k].row_basis[c].size();
						sum += tree[j][k].col_basis[c].size();

					}
				}
			}
		}
		for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) { //self
			sum += tree[nLevels][k].chargeLocations.size()*tree[nLevels][k].chargeLocations.size(); //self
			// sum += tree[nLevels][k].chargeLocations.size()*tree[nLevels][k].chargeLocations.size();
		}
	}

};

#endif
