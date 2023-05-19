#include "kernel.hpp"
#include "ACA.hpp"
#include "HTree.hpp"

int main(int argc, char* argv[]) {
	int cubeRootN;
	int nParticlesInLeafAlong1D;
	double L;
	int TOL_POW;
	int Qchoice;
	if (argc < 5) {
		cubeRootN		=	20;
		nParticlesInLeafAlong1D	=	6; // assuming the particles are located at tensor product chebyNodes
		L			=	1.0;
		TOL_POW = 7;
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
	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	std::vector<pts3D> particles;
	userkernel* mykernel		=	new userkernel(particles, Qchoice);
	H3DTree<userkernel>* A	=	new H3DTree<userkernel>(mykernel, cubeRootN, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);

	A->set_Uniform_Nodes();
	// A->set_Standard_Cheb_Nodes();
	A->createTree();
	A->assign_Tree_Interactions();
	A->assign_Center_Location();
	A->assignChargeLocations();
	A->assignNonLeafChargeLocations();
	end		=	omp_get_wtime();
	double timeCreateTree	=	(end-start);

	std::cout << std::endl << "Number of particles is: " << A->N << std::endl;
	std::cout << std::endl << "Time taken to create the tree is: " << timeCreateTree << std::endl;
	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	A->getUVtTree();
	end		=	omp_get_wtime();

	double timeHODLR3DRepres =	(end-start);

	/////////////////////////////////////////////////////////////////////////
	int N = A->N;

	Eigen::VectorXd b=Eigen::VectorXd::Zero(N);
  int n = N/500;
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

	A->assignCharges(b);
	A->evaluateFarField();
	A->evaluate_NearField();
	Eigen::VectorXd AFMM_Ab;
	A->collectPotential(AFMM_Ab);
	A->reorder(AFMM_Ab);

	timeHODLR3DRepres += A->assTime;

	double timeMatVecProduct = A->matVecTime;
	std::cout << std::endl << "Time taken to construct HODLR3D representation is: " << timeHODLR3DRepres << std::endl;
	std::cout << std::endl << "Time taken to do Mat-Vec product is: " << timeMatVecProduct << std::endl << std::endl;
	double err;
  std::string fname;
  Eigen::VectorXd true_Ab = Eigen::VectorXd::Zero(N);
  for (it = s.begin(); it != s.end(); it++) {
    true_Ab = true_Ab + mykernel->getCol(*it);
  }

	err = (true_Ab - AFMM_Ab).norm()/true_Ab.norm();

	double sum;
	A->findMemory(sum);
	std::cout << "Memory in GB: " << sum/8*pow(10,-9) << std::endl << std::endl;
	std::cout << "CR: " << double(sum)/N/N << std::endl << std::endl;
	std::cout << "max rank: " << A->getMaxRank() << std::endl << std::endl;
	std::cout << "relative forward error: " << err << std::endl << std::endl;

	delete A;
	delete mykernel;
}
