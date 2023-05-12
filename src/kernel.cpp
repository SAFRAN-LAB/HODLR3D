#include "kernel.hpp"

	kernel::kernel(std::vector<pts3D>& particles) {
			this->particles = particles;
	}

	Eigen::VectorXd kernel::getRow(const int j, std::vector<int> col_indices) {
		int n_cols = col_indices.size();
		Eigen::VectorXd row(n_cols);
    #pragma omp parallel for
    for(int k = 0; k < n_cols; k++) {
        row(k) = this->getMatrixEntry(j, col_indices[k]);
    }
    return row;
  }

  Eigen::VectorXd kernel::getCol(const int k, std::vector<int> row_indices) {
		int n_rows = row_indices.size();
    Eigen::VectorXd col(n_rows);
    #pragma omp parallel for
    for (int j=0; j<n_rows; ++j) {
			col(j) = this->getMatrixEntry(row_indices[j], k);
    }
    return col;
  }

	Eigen::VectorXd kernel::getCol(const int k) {
		int n_rows = particles.size();
		Eigen::VectorXd col(n_rows);
		#pragma omp parallel for
		for (int j=0; j<n_rows; ++j) {
			col(j) = this->getMatrixEntry(j, k);
		}
		return col;
	}

  Eigen::MatrixXd kernel::getMatrix(std::vector<int> row_indices, std::vector<int> col_indices) {
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

	Eigen::MatrixXd kernel::getMatrix(int row_start_index, int col_start_index, int row_end_index, int col_end_index) {
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


		// RBF Logarithm
		double userkernel::RBF_Logarithm(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R	=	sqrt(R2);
			double b = 1.0;
			return log(1.0+R/b);
		}

		// RBF Exponential
		double userkernel::RBF_Exponential(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R	=	sqrt(R2);
			double b = 1.0;
			return exp(-R/b);
		}

		// RBF Inverse Quadric
		double userkernel::RBF_Inverse_Quadric(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R	=	sqrt(R2);
			double b = 1.0;
			return 1.0/(1.0+(R/b)*(R/b));
		}

		// RBF Quadric
		double  userkernel::RBF_Quadric(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R	=	sqrt(R2);
			double b = 1.0;
			return 1.0+(R/b)*(R/b);
		}

		// RBF Inverse Multi-Quadric
		double userkernel::RBF_Inverse_Multi_Quadric(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R	=	sqrt(R2);
			double b = 1.0;
			return 1.0/sqrt((R/b)*(R/b)+1);
		}

		// RBF Gaussian`
		double userkernel::RBF_Gaussian(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R	=	sqrt(R2);
			double b = 1.0;
			return exp(-(R/b)*(R/b));
		}

		// RBF Multi-quadric
		double userkernel::RBF_Multi_quadric(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R	=	sqrt(R2);
			double b = 1.0;
			return sqrt((R/b)*(R/b)+1);
		}

		// 3D Laplacian
		double userkernel::Laplacian_3D(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R	=	sqrt(R2);
			return 1/R;
		}

		// 1/r^4
		double userkernel::oneOverR4(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			return 1.0/(R2*R2);
		}

		// 2D Laplacian
		double userkernel::Laplacian_2D(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R = sqrt(R2);
			return log(R);
		}

		// RBF r^{2}log(r)
		double userkernel::RBF_spline(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R = sqrt(R2);
			return R*R*log(R);
		}

		// RBF J_{0}(r)
		double userkernel::kernel_besselJ(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R = sqrt(R2);
			double kappa = 1.0;
			return boost::math::cyl_bessel_j(0, kappa*R);
		}

		// RBF Y_{0}(r)
		double userkernel::kernel_besselY(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R = sqrt(R2);
			double kappa = 1.0;
			return boost::math::cyl_neumann(0, kappa*R);
		}

		// RBF cos(kr)/r
		double userkernel::Helmholtz_cos(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R = sqrt(R2);
			double kappa = 1.0;
			return cos(kappa*R)/R;
		}

		// RBF Feynman(r)
		double userkernel::Feynman(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R = sqrt(R2);
			return boost::math::cyl_bessel_k(0.5, R)/sqrt(R);
		}

		// RBF Yukawa(r)
		double userkernel::Yukawa(const unsigned i, const unsigned j) {
			pts3D r1 = particles[i];
			pts3D r2 = particles[j];
			double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y) + (r1.z-r2.z)*(r1.z-r2.z);
			double R = sqrt(R2);
			double m = 50;
			return exp(-m*R)/R;
		}
