//#include <bits/stdc++.h>
#include <iostream>
#include <map>
#include<cfloat>
#include <cmath>
#include <set>
#include <vector>
#include <Eigen/Dense>
// #include <Eigen/src/Core/ArithmeticSequence.h>
//#include <boost/math/special_functions/bessel.hpp>
// #define EIGEN_DONT_PARALLELIZE

const double PI	=	3.1415926535897932384;
struct pts3D {
	double x,y,z;
};

class kernel {
public:
  double a;
  std::vector<pts3D> particles;

	kernel(std::vector<pts3D>& particles) {
			this->particles = particles;
	}

	Eigen::VectorXd getRow(const int j, std::vector<int> col_indices) {
		int n_cols = col_indices.size();
		Eigen::VectorXd row(n_cols);
    #pragma omp parallel for
    for(int k = 0; k < n_cols; k++) {
        row(k) = this->getMatrixEntry(j, col_indices[k]);
    }
    return row;
  }

  Eigen::VectorXd getCol(const int k, std::vector<int> row_indices) {
		int n_rows = row_indices.size();
    Eigen::VectorXd col(n_rows);
    #pragma omp parallel for
    for (int j=0; j<n_rows; ++j) {
			col(j) = this->getMatrixEntry(row_indices[j], k);
    }
    return col;
  }

	Eigen::VectorXd getCol(const int k) {
		int n_rows = particles.size();
		Eigen::VectorXd col(n_rows);
		#pragma omp parallel for
		for (int j=0; j<n_rows; ++j) {
			col(j) = this->getMatrixEntry(j, k);
		}
		return col;
	}

  Eigen::MatrixXd getMatrix(std::vector<int> row_indices, std::vector<int> col_indices) {
		int n_rows = row_indices.size();
		int n_cols = col_indices.size();
    Eigen::MatrixXd mat(n_rows, n_cols);
    // #pragma omp parallel for collapse(2)
		// #pragma omp parallel for
    for (int j=0; j < n_rows; ++j) {
        // #pragma omp parallel for
        for (int k=0; k < n_cols; ++k) {
            mat(j,k) = this->getMatrixEntry(row_indices[j], col_indices[k]);
        }
    }
    return mat;
  }

	Eigen::MatrixXd getMatrix(int row_start_index, int col_start_index, int row_end_index, int col_end_index) {
		Eigen::MatrixXd mat(row_end_index-row_start_index, col_end_index-col_start_index);
		// #pragma omp parallel for
		for (int j=row_start_index; j < row_end_index; ++j) {
				// #pragma omp parallel for
				for (int k=col_start_index; k < col_end_index; ++k) {
						mat(j,k) = this->getMatrixEntry(j, k);
				}
		}
		return mat;
	}
	// 3D Laplacian
	double Laplacian_3D(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R	=	sqrt(R2);
		return 1/R;
	}
	double getMatrixEntry(const unsigned i, const unsigned j) {
		if (i==j) {
			return 0.0;
		}
		else{
			return Laplacian_3D(i, j);
		}
	}

  ~kernel() {};
};

