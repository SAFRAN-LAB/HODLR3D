#include "kernel.hpp"
#include "ACA.hpp"
#include "HODLR3DTree.hpp"

int main(int argc, char* argv[]) {
	int cubeRootN		=	atoi(argv[1]);
	int nParticlesInLeafAlong1D	=	atoi(argv[2]); // assuming the particles are located at tensor product chebyNodes
	double L			=	atof(argv[3]);
	int TOL_POW = atoi(argv[4]);
	int Qchoice = atoi(argv[5]);
	double start, end;
	int nLevels		=	ceil(3*log(double(cubeRootN)/nParticlesInLeafAlong1D)/log(8));
	std::cout << "nLevels: " << nLevels << std::endl;

	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	std::vector<pts3D> particles;
	userkernel* mykernel		=	new userkernel(particles, Qchoice);
	HODLR3DTree<userkernel>* A	=	new HODLR3DTree<userkernel>(mykernel, cubeRootN, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);

	A->set_Uniform_Nodes();
	A->createTree();
	A->assign_Tree_Interactions();
	A->assign_Center_Location();
	A->assignChargeLocations();

	A->assignNonLeafChargeLocations();
	end		=	omp_get_wtime();
	double timeCreateTree	=	(end-start);

	// std::cout << std::endl << "Number of particles is: " << A->N << std::endl;
	// std::cout << std::endl << "Time taken to create the tree is: " << timeCreateTree << std::endl;
	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	A->getUVtTree();
	end		=	omp_get_wtime();

	double timeHODLR3DRepres =	(end-start);

	/////////////////////////////////////////////////////////////////////////
	// For testing the code, b is considered to be a vector of ones and zeros in the following code snippet
	// accordingly true_Ab is calculated later
	{
		int N = A->N;
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
	}
	A->assignCharges(b);

	A->evaluateFarField();
	A->evaluate_NearField();
	Eigen::VectorXd AFMM_Ab;
	A->collectPotential(AFMM_Ab);
	A->reorder(AFMM_Ab);
	timeHODLR3DRepres += A->assTime;
	double timeMatVecProduct = A->matVecTime;

	double err;
	Eigen::VectorXd true_Ab = Eigen::VectorXd::Zero(N);
  for (it = s.begin(); it != s.end(); it++) {
		true_Ab = true_Ab + mykernel->getCol(*it);
	}
	std::cout << std::endl << "Time taken to construct HODLR3D representation is: " << timeHODLR3DRepres << std::endl;
	std::cout << std::endl << "Time taken to do Mat-Vec product is: " << timeMatVecProduct << std::endl << std::endl;

	err = (true_Ab - AFMM_Ab).norm()/true_Ab.norm();

	double sum;
	A->findMemory(sum);
	std::cout << "Memory in GB: " << sum/8*pow(10,-9) << std::endl << std::endl;
	std::cout << "CR: " << double(sum)/N/N << std::endl << std::endl; //compression rate
	std::cout << "max rank: " << A->getMaxRank() << std::endl << std::endl;
	std::cout << "relative forward error: " << err << std::endl << std::endl;

	delete A;
	delete mykernel;
}
