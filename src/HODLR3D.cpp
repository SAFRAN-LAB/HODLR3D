#include "HODLR3D.hpp"

HODLR3D::HODLR3D(int nParticlesInLeafAlong1D, double L, int TOL_POW, int Qchoice, std::vector<pts3D>& particles) {
		mykernel		=	new userkernel(particles, Qchoice);
		N = particles.size();
		double cubeRootN = pow(N,1.0/3.0);
		int nLevels		=	ceil(3*log(double(cubeRootN)/nParticlesInLeafAlong1D)/log(8));
		A	=	new HODLR3DTree(mykernel, N, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);
		A->K		=	new userkernel(particles, Qchoice);
		A->createTree();
		A->assign_Tree_Interactions();
		A->assign_Center_Location();
		A->assignChargeLocations();
		A->assignNonLeafChargeLocations();
	}

	void HODLR3D::assemble() {
		A->getUVtTree();
	}

	void HODLR3D::MatVecProduct(Eigen::VectorXd &b, Eigen::VectorXd &AFMM_Ab) {
		A->assignCharges(b);
		A->evaluateFarField();
		A->evaluate_NearField();
		A->collectPotential(AFMM_Ab);
		A->reorder(AFMM_Ab);
	}

	HODLR3D::~HODLR3D() {
		delete A;
		delete mykernel;
	};
