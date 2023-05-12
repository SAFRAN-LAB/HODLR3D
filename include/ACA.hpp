#ifndef _ACA_HPP__
#define _ACA_HPP__

#include <iostream>
#include <Eigen/Dense>
#include "kernel.hpp"
#include <bits/stdc++.h>

class LowRank {
public:
	userkernel* K;
	double tol_ACA;
	std::vector<int> row_indices;
	std::vector<int> col_indices;

	LowRank(userkernel* K, int tol_pow, std::vector<int>& row_indices, std::vector<int>& col_indices);

	void maxAbsVector(const Eigen::VectorXd& v, const std::set<int>& allowed_indices,
																double max, int& index
															);

	void ACA_only_nodesUV(std::vector<int>& row_bases, std::vector<int>& col_bases, int &computed_rank, Eigen::MatrixXd &U, Eigen::MatrixXd &Vt);


	void ACA_only_nodesCUR(std::vector<int>& row_bases, std::vector<int>& col_bases, int &computed_rank, Eigen::MatrixXd &L, Eigen::MatrixXd &R);


};
#endif
