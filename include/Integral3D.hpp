/*
This header file evaluates numerically the triple integral of Integrand.
triple integral
*/
#ifndef Integral3D_hpp
#define Integral3D_hpp

#include <iostream>
#include <cmath>
#include "definitions.hpp"
#include "Gauss_Legendre_Nodes_and_Weights.hpp"

double Integrand(double x,double y,double z);

double triple_integral(double*& a,double*& b);

#endif
