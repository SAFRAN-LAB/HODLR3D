#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
// #include <Eigen/src/Core/ArithmeticSequence.h>
#include <boost/math/special_functions/bessel.hpp>
// #define EIGEN_DONT_PARALLELIZE
#include "Integral3D.hpp"

// const double PI	=	3.1415926535897932384;
#include <map>
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

	virtual double getMatrixEntry(const unsigned i, const unsigned j) {
		std::cout << "virtual getInteraction" << std::endl;
		return 0.0;
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

  ~kernel() {};
};

class userkernel: public kernel {
public:
	int Qchoice;
	double Kii;
	double h3;
	double chargesFunction(const pts3D r) {
		double q = r.x; //user defined
		return q;
	};
	// #ifdef ONEOVERR
	// userkernel(std::vector<pts3D>& particles_X, std::vector<pts3D>& particles_Y): kernel(particles_X, particles_Y) {
	// };

	// #elif LOGR
	userkernel(std::vector<pts3D> particles, int Qchoice): kernel(particles) {
		this->Qchoice = Qchoice;
		if(Qchoice == 16){
			double h = 1.0/cbrt(double(particles.size()));
			// std::cout << "particles.size(): " << particles.size() << std::endl;
			double *a,*b;
			a = new double[3];
			b = new double[3];
			a[0] = 0;
			a[1] = 0;
			a[2] = 0;

			b[0] = h*0.5;
			b[1] = h*0.5;
			b[2] = h*0.5;
			Kii = triple_integral(a,b);
			Kii += 1; //For Second kind integral equation.
			h3 = 1.0/double(particles.size());
			// std::cout << "Kii: " << Kii << std::endl;
			// std::cout << "h3: " << h3 << std::endl;
		}
	};
	// double getMatrixEntry(const unsigned i, const unsigned j) {
	// 	pts3D r1 = particles_X[i];
	// 	pts3D r2 = particles_X[j];
	// 	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
	// 	double b = 0.01;
	// 	if (R2 < 1e-10) {
	// 		return 0.0;
	// 	}
	// 	else if (R2 < b*b) {
	// 		return 0.5*R2*log(R2)/b/b;
	// 	}
	// 	else {
	// 		return 0.5*log(R2);
	// 	}
	// }
	double IE_CUBE_3D(const unsigned i,const unsigned j){
		// std::cout << "Kii: " << Kii << std::endl;
		// std::cout << "h3: " << h3 << std::endl;
		double res;
		if(i==j)
			res = Kii;//  triple_integral();
		else{
			res = h3 * Laplacian_3D(i,j);
		}
		return res;
	}

	// RBF Logarithm
	double RBF_Logarithm(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R	=	sqrt(R2);
		double b = 1.0;
		return log(1.0+R/b);
	}

	// RBF Exponential
	double RBF_Exponential(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R	=	sqrt(R2);
		double b = 1.0;
		return exp(-R/b);
	}

	// RBF Inverse Quadric
	double RBF_Inverse_Quadric(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R	=	sqrt(R2);
		double b = 1.0;
		return 1.0/(1.0+(R/b)*(R/b));
	}

	// RBF Quadric
	double  RBF_Quadric(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R	=	sqrt(R2);
		double b = 1.0;
		return 1.0+(R/b)*(R/b);
	}

	// RBF Inverse Multi-Quadric
	double RBF_Inverse_Multi_Quadric(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R	=	sqrt(R2);
		double b = 1.0;
		return 1.0/sqrt((R/b)*(R/b)+1);
	}

	// RBF Gaussian`
	double RBF_Gaussian(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R	=	sqrt(R2);
		double b = 1.0;
		return exp(-(R/b)*(R/b));
	}

	// RBF Multi-quadric
	double RBF_Multi_quadric(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R	=	sqrt(R2);
		double b = 1.0;
		return sqrt((R/b)*(R/b)+1);
	}

	// 3D Laplacian
	double Laplacian_3D(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R	=	sqrt(R2);
		return 1/R;
	}

	// 1/r^4
	double oneOverR4(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		// double R	=	sqrt(R2);
		return 1.0/(R2*R2);
	}

	// 2D Laplacian
	double Laplacian_2D(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R = sqrt(R2);
		return log(R);
	}

	// RBF r^{2}log(r)
	double RBF_spline(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R = sqrt(R2);
		return R*R*log(R);
	}

	// RBF J_{0}(r)
	double kernel_besselJ(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R = sqrt(R2);
		// double a = 0.001;
		double kappa = 1.0;
		return boost::math::cyl_bessel_j(0, kappa*R);
	}

	// RBF Y_{0}(r)
	double kernel_besselY(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R = sqrt(R2);
		// double a = 0.001;
		double kappa = 1.0;
		return boost::math::cyl_neumann(0, kappa*R);
	}

	// RBF cos(kr)/r
	double Helmholtz_cos(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R = sqrt(R2);
		// double a = 0.001;
		double kappa = 1.0;
		return cos(kappa*R)/R;
	}

	// RBF Feynman(r)
	double Feynman(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R = sqrt(R2);
		// double a = 0.001;
		// double kappa = 1.0;
		return boost::math::cyl_bessel_k(0.5, R)/sqrt(R);
	}

	// RBF Yukawa(r)
	double Yukawa(const unsigned i, const unsigned j) {
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
		double R = sqrt(R2);
		double m = 50;
		// double kappa = 1.0;
		return exp(-m*R)/R;
	}

	// check(r)
	double check(const unsigned i, const unsigned j) {
		return double(i)/particles.size();
	}

	double getMatrixEntry(const unsigned i, const unsigned j) {
		if (Qchoice != 16) {
			if (i==j) {
				return 0.0;
			}
			else {
				if (Qchoice == 0)
					return RBF_Logarithm(i, j);
				else if (Qchoice == 1)
					return RBF_Exponential(i, j);
				// else if (Qchoice == 2)
				// 	return RBF_Inverse_Quadric(i, j);
				// else if (Qchoice == 3)
				// 	return RBF_Quadric(i, j);
				// else if (Qchoice == 4)
				// 	return RBF_Inverse_Multi_Quadric(i, j);
				// else if (Qchoice == 5)
				// 	return RBF_Gaussian(i, j);
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
			else {
				return IE_CUBE_3D(i,j);
			}
		}

	// #endif
	~userkernel() {};
};
