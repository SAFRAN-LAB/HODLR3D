#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp>
#include "definitions.hpp"

#include <map>

#include "kernel.hpp"
#include "ACA.hpp"
#include "HODLR3DTree.hpp"
#include "HODLR3D.hpp"

void set_Uniform_Nodes(int cubeRootN, double L, std::vector<pts3D>& particles) {
	std::vector<double> Nodes1D;
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
				particles.push_back(temp1);
			}
		}
	}
}

double userkernel::getMatrixEntry(const unsigned i, const unsigned j) {
	if (i==j) {
		return 0.0;
	}
	else {
		if (Qchoice == 0)
			return RBF_Logarithm(i, j);
		else if (Qchoice == 1)
			return RBF_Exponential(i, j);
		else if (Qchoice == 2)
			return RBF_Inverse_Quadric(i, j);
		else if (Qchoice == 3)
			return RBF_Quadric(i, j);
		else if (Qchoice == 4)
			return RBF_Inverse_Multi_Quadric(i, j);
		else if (Qchoice == 5)
			return RBF_Gaussian(i, j);
		else if (Qchoice == 6)
			return RBF_Multi_quadric(i, j);
		else if (Qchoice == 7)
			return Laplacian_3D(i, j);
		else if (Qchoice == 8)
			return oneOverR4(i, j);
		else if (Qchoice == 9)
			return Laplacian_2D(i, j);
		else if (Qchoice == 10)
			return RBF_spline(i, j);
		else if (Qchoice == 11)
			return kernel_besselJ(i, j);
		else if (Qchoice == 12)
			return kernel_besselY(i, j);
		else if (Qchoice == 13)
			return Helmholtz_cos(i, j);
		else if (Qchoice == 14)
			return Feynman(i, j);
		else if (Qchoice == 15)
			return Yukawa(i, j);
		}
}

int main(int argc, char* argv[]) {
	int cubeRootN;
	int nParticlesInLeafAlong1D; // assuming the particles are located at tensor product chebyNodes
	double L;
	int TOL_POW;
	int Qchoice;
	if (argc < 5) {
		cubeRootN		=	40;
		nParticlesInLeafAlong1D	=	8; // assuming the particles are located at tensor product chebyNodes
		L			=	1.0;
		TOL_POW = 6;
		Qchoice = 7;
	}
	else {
		cubeRootN		=	atoi(argv[1]);
		nParticlesInLeafAlong1D	=	atoi(argv[2]); // assuming the particles are located at tensor product chebyNodes
		L			=	atof(argv[3]);
		TOL_POW = atoi(argv[4]);
		Qchoice = atoi(argv[5]);
	}
	double start, end;

	int nLevels		=	ceil(3*log(double(cubeRootN)/nParticlesInLeafAlong1D)/log(8));
	std::cout << "nLevels: " << nLevels << std::endl;

	// Defining particles at the tensor product grid of uniformly distributed pow(cubeRootN,3) number of particles in [-L,L]^3.
	std::vector<pts3D> particles;
	set_Uniform_Nodes(cubeRootN, L, particles);
	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	HODLR3D *H = new HODLR3D(nParticlesInLeafAlong1D, L, TOL_POW, Qchoice, particles);
	end		=	omp_get_wtime();
	double timeCreateTree	=	(end-start);

	int N = H->N;
	std::cout << std::endl << "Number of particles is: " << N << std::endl;
	std::cout << std::endl << "Time taken to create the tree is: " << timeCreateTree << std::endl;
	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	H->assemble();
	end		=	omp_get_wtime();

	double timeHODLR3DRepres =	(end-start);

	/////////////////////////////////////////////////////////////////////////
	// For testing the code, b is considered to be a vector of ones and zeros in the following code snippet
	// accordingly true_Ab is calculated later

		Eigen::VectorXd b=Eigen::VectorXd::Zero(N);
	  int n = N/500; //randomly choosing n different indices where b is set to 1, b at the rest of the indices is set to 0
	  srand(time(NULL));
	  std::set<int> s;
	  while(s.size() < n) {
	    int index	=	rand()%N;
	    s.insert(index);
	  }
	  std::set<int>::iterator it;
	  for (it = s.begin(); it != s.end(); it++) {
	    b(*it) = 1.0;
	  }

	Eigen::VectorXd Hb;
	H->MatVecProduct(b, Hb);
	timeHODLR3DRepres += H->A->assTime;
	double timeMatVecProduct = H->A->matVecTime;

	double err;
	Eigen::VectorXd true_Ab = Eigen::VectorXd::Zero(N);
  for (it = s.begin(); it != s.end(); it++) {
		true_Ab = true_Ab + H->mykernel->getCol(*it);
	}
	std::cout << std::endl << "Time taken to construct HODLR3D representation is: " << timeHODLR3DRepres << std::endl;
	std::cout << std::endl << "Time taken to do Mat-Vec product is: " << timeMatVecProduct << std::endl << std::endl;

	err = (true_Ab - Hb).norm()/true_Ab.norm();

	double sum;
	H->A->findMemory(sum);
	std::cout << "Memory in GB: " << sum/8*pow(10,-9) << std::endl << std::endl;
	std::cout << "CR: " << double(sum)/N/N << std::endl << std::endl; //compression rate
	std::cout << "max rank: " << H->A->getMaxRank() << std::endl << std::endl;
	std::cout << "relative forward error: " << err << std::endl << std::endl;

	delete H;
}
