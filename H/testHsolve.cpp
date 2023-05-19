#include "kernel.hpp"
#include "ACA.hpp"
#include "HTree.hpp"
#include "gmres.hpp"

class HODLR3D {
public:
	H3DTree<userkernel>* A;
	HODLR3D(int cubeRootN, int nParticlesInLeafAlong1D, double L, int TOL_POW, int Qchoice) {
		std::vector<pts3D> particles;
		userkernel* mykernel		=	new userkernel(particles, Qchoice);
		int nLevels		=	ceil(3*log(double(cubeRootN)/nParticlesInLeafAlong1D)/log(8));
		A	=	new H3DTree<userkernel>(mykernel, cubeRootN, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);

		A->set_Uniform_Nodes();
		particles = A->K->particles;
		A->K		=	new userkernel(particles, Qchoice);
		A->createTree();
		A->assign_Tree_Interactions();
		A->assign_Center_Location();
		A->assignChargeLocations();
		A->assignNonLeafChargeLocations();
	}
	void factorize() {
		A->getUVtTree();
	}
	void MatVecProduct(Eigen::VectorXd &b, Eigen::VectorXd &AFMM_Ab) {
		A->assignCharges(b);
		A->evaluateFarField();
		A->evaluate_NearField();
		A->collectPotential(AFMM_Ab);
		A->reorder(AFMM_Ab);
	}
	~HODLR3D() {};
};

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
		Qchoice = 16;
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

	start	=	omp_get_wtime();
	HODLR3D *H = new HODLR3D(cubeRootN, nParticlesInLeafAlong1D, L, TOL_POW, Qchoice);
	end		=	omp_get_wtime();
	double timeCreateTree	=	(end-start);

	std::cout << std::endl << "Number of particles is: " << H->A->N << std::endl;
	std::cout << std::endl << "Time taken to create the tree is: " << timeCreateTree << std::endl << std::endl;
	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	H->factorize();
	end		=	omp_get_wtime();

	double timeHODLR3DRepres =	(end-start);
	/////////////////////////////////////////////////////////////////////////
	int N = H->A->N;
	Eigen::VectorXd x=Eigen::VectorXd::Zero(N);
  int n = N/500; //randomly choosing n different indices where x is set to 1, x at the rest of the indices is set to 0
  srand(time(NULL));
  std::set<int> s;
  while(s.size() < n) {
    int index	=	rand()%N;
    s.insert(index);
  }
  std::set<int>::iterator it;
  for (it = s.begin(); it != s.end(); it++) {
    x(*it) = 1.0;
  }

	Eigen::VectorXd true_Ax = Eigen::VectorXd::Zero(N);
  for (it = s.begin(); it != s.end(); it++) {
		true_Ax = true_Ax + H->A->K->getCol(*it);
	}

	double sum;
	H->A->findMemory(sum);
	///////////////////////////////////// solving system ////////////////////////////////////
	int maxIterations = 400;
	double GMRES_threshold = 1e-10;
	double GMRES_residual;
	Vec H3D_x;
	int noOfIterations;
	std::vector<double> e;
	classGMRES *G = new classGMRES();

	start	=	omp_get_wtime();
	G->gmres(H, true_Ax, maxIterations, GMRES_threshold, H3D_x, GMRES_residual, noOfIterations, e);
	end		=	omp_get_wtime();
	double timeGMRES = end-start;
	std::cout << "time GMRES: " << timeGMRES << std::endl << std::endl;
	std::cout << "GMRES residual err: " << GMRES_residual << std::endl << std::endl;
	std::cout << "GMRES no. of iterations: " << noOfIterations << std::endl << std::endl;
	Vec err = H3D_x-x;
	std::cout << "relative forward error in solution: " << err.norm()/x.norm() << std::endl << std::endl;
	delete H;
}
