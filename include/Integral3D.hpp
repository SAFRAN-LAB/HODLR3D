/*
This header file evaluates numerically the triple integral of Integrand.
triple integral
*/
#ifndef Integral3D_hpp
#define Integral3D_hpp

#include <iostream>
#include <cmath>
#include "Gauss_Legendre_Nodes_and_Weights.hpp"


int n_gauss_points = 16;
const double PI = 3.1415926535897932384;


double Integrand(double x,double y,double z)
{
    return 1.0/((4*PI)*sqrt(x*x + y*y + z*z));
}


double triple_integral(double*& a,double*& b)
{
        double result;
        double tx,ty,tz,cx,cy,cz,Lx,Ly,Lz;
        double *nodes_x,*weights_x;
        Gauss_Legendre_Nodes_and_Weights(n_gauss_points, nodes_x, weights_x);
        double *nodes_y,*weights_y;
        Gauss_Legendre_Nodes_and_Weights(n_gauss_points, nodes_y, weights_y);
        double *nodes_z,*weights_z;
        Gauss_Legendre_Nodes_and_Weights(n_gauss_points, nodes_z, weights_z);
        cx = 0.5 * (b[0] - a[0]);
        cy = 0.5 * (b[1] - a[1]);
        cz = 0.5 * (b[2] - a[2]);
        Lx = 0.5 * (b[0] + a[0]);
        Ly = 0.5 * (b[1] + a[1]);
        Lz = 0.5 * (b[2] + a[2]);

        // Gauss-Legendre Quadrature
        for(int i=0; i<n_gauss_points; i++)
        {
            tx = cx*nodes_x[i] + Lx;
            for(int j=0; j<n_gauss_points; j++)
            {
                    ty = cy*nodes_y[j] + Ly;
                    for(int k=0; k<n_gauss_points; k++)
                    {
                            tz = cz*nodes_z[k] + Lz;
                            result += (weights_x[i]*weights_y[j]*weights_z[k] * Integrand(tx,ty,tz));
                    }
            }
        }
        result *= (cx*cy*cz);
        result *= 8; // This is reducing the integral into half domain
        return result;
}
#endif
