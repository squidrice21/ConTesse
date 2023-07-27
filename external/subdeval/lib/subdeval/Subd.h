/* Subd.h - Brent Burley, Feb 2005
   Catmull-Clark subdivision implementation.
   Modified by Dylan Lacewell:
   Added limit surface evaluation at any (faceid, u, v) location.
*/

#pragma once

#include <string>
#include <memory>

class SubdInternal;

//static double ROUNDING_THRESHOLD = 1e-10;
static double ROUNDING_THRESHOLD = 1e-5;

class Subd {
 public:
    Subd() = delete;
    Subd(const Subd&);
    Subd(int nverts, const double* verts, int nfaces, const int* nvertsPerFace,
     const int* faceverts);
    ~Subd();
    void subdivide(int levels=1);
    bool eval(int faceid, double u, double v, double* p, double* dPdU, double* dPdV);
    int nverts();
    int nfaces();
    int nfaceverts();
    const double* verts();
    const int* nvertsPerFace();
    const int* faceverts();
    const double* normals();
    const double* limitverts();

 private:
    std::unique_ptr<SubdInternal> impl;
};
