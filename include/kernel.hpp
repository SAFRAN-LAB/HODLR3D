#ifndef _kernel_HPP__
#define _kernel_HPP__

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include "definitions.hpp"
#include <boost/math/special_functions/bessel.hpp>

class kernel {
public:
  double a;
  std::vector<pts3D> particles;

	kernel(std::vector<pts3D>& particles);

	virtual double getMatrixEntry(const unsigned i, const unsigned j) {
		std::cout << "virtual getInteraction" << std::endl;
		return 0.0;
	}

	Eigen::VectorXd getRow(const int j, std::vector<int> col_indices);

  Eigen::VectorXd getCol(const int k, std::vector<int> row_indices);

	Eigen::VectorXd getCol(const int k);

  Eigen::MatrixXd getMatrix(std::vector<int> row_indices, std::vector<int> col_indices);

	Eigen::MatrixXd getMatrix(int row_start_index, int col_start_index, int row_end_index, int col_end_index);

  ~kernel() {};
};

class userkernel: public kernel {
public:
	int Qchoice;

	userkernel(std::vector<pts3D> particles, int Qchoice): kernel(particles) {
		this->Qchoice = Qchoice;
	};

	// RBF Logarithm
	double RBF_Logarithm(const unsigned i, const unsigned j);

	// RBF Exponential
	double RBF_Exponential(const unsigned i, const unsigned j);

	// RBF Inverse Quadric
	double RBF_Inverse_Quadric(const unsigned i, const unsigned j);

	// RBF Quadric
	double  RBF_Quadric(const unsigned i, const unsigned j);

	// RBF Inverse Multi-Quadric
	double RBF_Inverse_Multi_Quadric(const unsigned i, const unsigned j);

	// RBF Gaussian`
	double RBF_Gaussian(const unsigned i, const unsigned j);

	// RBF Multi-quadric
	double RBF_Multi_quadric(const unsigned i, const unsigned j);

	// 3D Laplacian
	double Laplacian_3D(const unsigned i, const unsigned j);

	// 1/r^4
	double oneOverR4(const unsigned i, const unsigned j);

	// 2D Laplacian
	double Laplacian_2D(const unsigned i, const unsigned j);

	// RBF r^{2}log(r)
	double RBF_spline(const unsigned i, const unsigned j);

	// RBF J_{0}(r)
	double kernel_besselJ(const unsigned i, const unsigned j);

	// RBF Y_{0}(r)
	double kernel_besselY(const unsigned i, const unsigned j);

	// RBF cos(kr)/r
	double Helmholtz_cos(const unsigned i, const unsigned j);
	// RBF Feynman(r)
	double Feynman(const unsigned i, const unsigned j);

	// RBF Yukawa(r)
	double Yukawa(const unsigned i, const unsigned j);

	double getMatrixEntry(const unsigned i, const unsigned j);

	~userkernel() {};
};

#endif
