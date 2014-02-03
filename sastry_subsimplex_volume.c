// Copyright (C) 2014 Alexey Anikeenko
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <math.h>
#include <stdlib.h>
#include "sastry_subsimplex_volume.h"

#define PI 3.14159265358979323846

/*
 * double subsimplex_void_volume(double x0, y0, z0, rC2, *surf_area)
 *
 * Computes void volume and area in tetrahedron with three orthogonal edges having sphere in one of its vertices.
 * Based on formulas from paper:
 * S. Sastry, D.S. Corti, P.G. Debenedetti, F.H. Stillinger, "Statistical geometry of particle packings. I. Algorithm for exact determination of connectivity,
 * volume, and surface areas of void space in monodisperse and polydisperse sphere packings", Phys. Rev. E, V.56, N5, p. 5524, 1997. doi:10.1103/PhysRevE.56.5524
 * Formula (A8) was improved in work:
 * V.P. Voloshin, A.V. Anikeenko, N.N. Medvedev, A. Geiger, "An algorithm for the calculation of volume and surface of unions of spheres. Application for solvation shells",
 * proceedings of ISVD conference, 2011.

 * Input:
 * ========================
 * coordinates of first vertex (with sphere)  0,  0,  0
 *               second vertex               x0,  0,  0
 *                third vertex               x0, y0,  0
 *               fourth vertex               x0, y0, z0
 * rC2 - squared sphere radius

 * Output:
 * ========================
 * return void volume and area of spherical surface
 *
 */
double subsimplex_void_volume(double x0, double y0, double z0, double rC2, double *surf_area)
{
    double V = 0, S = 0;             // volume occupied by sphere in subsimplex and it's area
    const double Vc = x0*y0*z0/6;

    if (x0 != 0 && y0 != 0 && z0 != 0) {
        double x, a11;
        const double theta = atan(z0/y0);
        const double rE2 = x0*x0 + y0*y0;
        const double rV2 = rE2 + z0*z0;
        const double rE = sqrt(rE2);
        const double rV = sqrt(rV2);
        const double rC = (rC2 > 0) ? sqrt(rC2) : 0;


        a11 = (y0*y0 * rV2 - z0*z0 * x0*x0)/(rE2 * (y0*y0 + z0*z0));   // connection to papers: a11 = (pi/2 + a1)/2
        if (a11 > 1.0) a11 = 0;
        else if (a11 < -1.0) a11 = PI/2;
        else a11 = acos(a11)*0.5;

        if (rC <= x0) 
        {
            x = rC*rC*(theta-a11);
            S = x;
            V = rC*x/3;
        }
        else if (rC <= rE)
        {
            x = rC2*a11;
            S = theta*x0*rC - x;
            V = 0.5*theta*x0*(rC2 - x0*x0/3) - rC*x/3;
        }
        else if (rC < rV)
        {
            double a22, a33, y, x2, y2;
            x2 = rC*(x0/rE);
            y2 = rC*(y0/rE);

            a22 = y0/sqrt(rC2-x0*x0);
            if (a22 > 1.0) a22 = 0;
            else if (a22 < -1.0) a22 = PI;
            else a22 = acos(a22);

            a33 = (x0*x0 - x2*x2 + y2*y2)/(rC2 - x0*x0);
            if (a33 > 1.0) a33 = 0;
            else if (a33 < -1.0) a33 = PI/2;
            else a33 = acos(a33)*0.5;

            y = rC2*(a33 - a11);
            x = x0*(theta-a22);

            S = rC*x + y;
            V = x*(0.5*rC2 - x0*x0/6) + x0*y0*sqrt(rC2 - rE2)/6 + rC*y/3;
        }
        else
        {
            S = 0;
            V = Vc;
        }
    }

    if (surf_area != NULL)
        *surf_area = S;

    return (Vc-V);
}
