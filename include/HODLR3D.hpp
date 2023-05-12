#ifndef _HODLR3D_hpp__
#define _HODLR3D_hpp__

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "definitions.hpp"
#include "kernel.hpp"
#include "ACA.hpp"
#include "HODLR3DTree.hpp"

class HODLR3D {
public:
	HODLR3DTree* A;
	userkernel* mykernel;
	int N;
	HODLR3D(int cubeRootN, int nParticlesInLeafAlong1D, double L, int TOL_POW, int Qchoice);
	void factorize();
	void MatVecProduct(Eigen::VectorXd &b, Eigen::VectorXd &AFMM_Ab);

	~HODLR3D();
};
#endif
