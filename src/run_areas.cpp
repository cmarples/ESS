// Run ellipsoid patch area calculation.
// Author: Callum Marples

#include "run_tests.hpp"
#include "surface_area.hpp"
#include "ellipsoid_shape.hpp"
#include "ellipsoid_patches.hpp"
#include "vectors.hpp"
#include <ctime>
#include <cmath>
#include <memory>
#include <vector>
#include <array>

using namespace std;

// For a given set of shapes, determine surface area of all patches in
// a 181 by 360 grid and output to files. This data is then used to 
// get contact distributions per unit area on the same shapes.
void RunAreas() {

    array<vector<double>, 4> Areas;

    double r1 = 1.0;
    double a2 = 3.0;
    double c2 = 1.5;
    double a3 = 3.0;
    double b3 = 2.0;
    double c3 = 1.0;

    EllipsoidShape E1(r1);            // Sphere    (1, 1, 1)
    EllipsoidShape E2(a2, c2);        // Oblate Spheroid  (3, 3, 1.5)
    EllipsoidShape E3(c2, a2);        // Prolate Spheroid  (1.5, 1.5, 3)
    EllipsoidShape E4(a3, b3, c3);    // Ellipsoid (3, 2, 1)

    PolarPatches P(181, 360);
    int noPatches = P.GetNoPatches();

    Areas[0] = EllipsoidPatchAreas(E1.aAxis, E1.bAxis, E1.cAxis, P);
    Areas[1] = EllipsoidPatchAreas(E2.aAxis, E2.bAxis, E2.cAxis, P);
    Areas[2] = EllipsoidPatchAreas(E3.aAxis, E3.bAxis, E3.cAxis, P);
    Areas[3] = EllipsoidPatchAreas(E4.aAxis, E4.bAxis, E4.cAxis, P);

    // Write information to file
    ofstream File("data/patch_areas.csv");
    if (File.is_open()) {
        File << fixed << setprecision(15);
        File << "Surface areas of ellipsoidal patches\n";
        File << "Sphere,Oblate,Prolate,Triaxial\n";
        for (int j=0; j!=noPatches; ++j) {
            File << Areas[0][j] << "," << Areas[1][j] << "," << Areas[2][j] << "," << Areas[3][j] << "\n";
        }
        File.close();
    }

    // Write information to file
    ofstream File2("data/ellipsoid_shapes.csv");
    if (File2.is_open()) {
        File2 << fixed;
        File2 << "Sphere\n";
        File2 << E1.aAxis << "," << E1.bAxis << "," << E1.cAxis << "\n";
        File2 << "Oblate Spheroid\n";
        File2 << E2.aAxis << "," << E2.bAxis << "," << E2.cAxis << "\n";
        File2 << "Prolate Spheroid\n";
        File2 << E3.aAxis << "," << E3.bAxis << "," << E3.cAxis << "\n";
        File2 << "Ellipsoid\n";
        File2 << E4.aAxis << "," << E4.bAxis << "," << E4.cAxis << "\n";
        File2.close();
    }
}
