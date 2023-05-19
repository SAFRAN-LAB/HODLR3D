#include "kernel.hpp"
#include "ACA.hpp"
#include "FMM3DTreeRAMeff2.hpp"
#include "parallelH3D.hpp" 
#include <omp.h>
#include <mpi.h>

int main(int argc, char* argv[]) {
	MPI_Init(NULL, NULL);
	int myid,numprocs,sz;
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    std::cout << "MPI Code with " << numprocs << " processors.." << std::endl;
	int cubeRootN		=	atoi(argv[1]);
	// assuming the particles are located at tensor product of 1D nodes
	int nParticlesInLeafAlong1D	=	atoi(argv[2]); 
	double L			=	atof(argv[3]);
	int TOL_POW = atoi(argv[4]);
	int Qchoice = atoi(argv[5]);
	double start, end;
	int nLevels		=	ceil(3*log(double(cubeRootN)/nParticlesInLeafAlong1D)/log(8));

	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	std::vector<pts3D> particles;
	std::vector<pts3D> particles_;
	kernel* mykernel		=	new kernel(particles);
	kernel* mykernel_		=	new kernel(particles_);
	FMM3DTree* A	=	new FMM3DTree(mykernel, cubeRootN, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);
	std::cout << " Tree formed.. with " << nLevels <<" levels" << std::endl;
	std::cout << "System setting - " << particles.size() << "," << TOL_POW << std::endl;
	parallelH3Dtree Bb(numprocs, mykernel_, cubeRootN, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);

	Bb.Initialize_parallelH3D();
	std::cout << " Initialised" << std::endl;

	int N = Bb.system_size();
	Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
	int n = 50;
	//std::cout << "Number of rows is: " << n << std::endl;
	// randomly choosing n different indices where b is set to 1, b at the rest of the indices is set to 0
	srand(1999);
	std::set<int> s;
	while (s.size() < n){
		int index = rand() % N;
		s.insert(index);
	}
	std::set<int>::iterator it;
	for (it = s.begin(); it != s.end(); it++){
		b(*it) = 1.0;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	Eigen::VectorXd AFMM_Ab = Bb.mat_vec(b);
	MPI_Barrier(MPI_COMM_WORLD);
	// std::cout << "Collected potential.." << std::endl;
	if (myid == 0){
		double err;
		Eigen::VectorXd true_Ab = Eigen::VectorXd::Zero(N);
		for (it = s.begin(); it != s.end(); it++){
			true_Ab = true_Ab + mykernel_->getCol(*it);
		}
		err = (true_Ab - AFMM_Ab).norm() / true_Ab.norm();
		//std::cout << "Error in Ab.." << true_Ab.norm() << std::endl;
		//std::cout << "Error in AFFMM Ab.." << AFMM_Ab.norm() << std::endl;
		std::cout << "Error in sol.." << err << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	/*A->set_Uniform_Nodes();
	A->createTree();
	A->assign_Tree_Interactions();
	A->assign_Center_Location();
	A->assignChargeLocations();

	// A->assignLeafChargeLocations();
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
	// std::cout << std::endl << "Time taken to construct HODLR3D representation is: " << timeHODLR3DRepres << std::endl;

	/////////////////////////////////////////////////////////////////////////
	int N = A->N;
	Eigen::VectorXd b=Eigen::VectorXd::Zero(N);
	int n = N/500; 
	//randomly choosing n different indices where b is set to 1, b at the rest of the indices is set to 0
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
	double err;
	Eigen::VectorXd true_Ab = Eigen::VectorXd::Zero(N);
  	for (it = s.begin(); it != s.end(); it++) {
		true_Ab = true_Ab + mykernel->getCol(*it);
	}
	err = (true_Ab - AFMM_Ab).norm()/true_Ab.norm();
	double sum;
	A->findMemory(sum);
	std::cout << N << " " << TOL_POW << " " << Qchoice << " " << sum/8*pow(10,-9) << " " << double(sum)/N/N << " " << timeHODLR3DRepres << " " << timeMatVecProduct << " " << A->getMaxRank() << " " << err << std::endl;
	delete A;
	delete A;
	delete mykernel;
	delete mykernel_;*/
	MPI_Finalize();
	return 0;
}

/*
icpx -I"/usr/local/include/include/eigen3" -I"/home/kva/Github/HODLR3D/HODLR3D/" -I"/home/kva/intel/oneapi/mpi/2021.7.0/include" -L"/home/kva/intel/oneapi/mpi/2021.7.0/lib/release" -L"/home/kva/intel/oneapi/mpi/2021.7.0/lib" -Xlinker --enable-new-dtags -Xlinker -rpath -Xlinker "$/home/kva/intel/oneapi/mpi/2021.7.0/lib/release" -Xlinker -rpath -Xlinker "/home/kva/intel/oneapi/mpi/2021.7.0/lib" -lmpicxx -lmpifort -std=c++17 -lmpi -ldl -lrt -lpthread -qopenmp /home/kva/Github/HODLR3D/HODLR3D/testFMM3DRAMeff2.cpp -o hodlr3d
*/