#ifndef _gmres_hpp__
#define _gmres_hpp__

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "HODLR3D.hpp"

typedef double Dtype;
typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXd Mat;

class classGMRES {
public:
classGMRES() {};
/// Calculate the Given rotation matrix----%%
void givensRotation(Dtype v1, Dtype v2, Dtype& cs, Dtype& sn);

void applyGivensRotation(Vec& h, Vec cs, Vec sn, int k, Dtype& cs_k, Dtype& sn_k);

void arnoldi(HODLR3D* M, Mat Q, int k, Vec &h, Vec &q);

void gmres(HODLR3D* M, Vec b, int maxIterations, double threshold, Vec& x, double& error, int& noOfIterations, std::vector<double>& e);

~classGMRES();

};
#endif
