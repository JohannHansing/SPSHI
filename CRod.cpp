#include "headers/CRod.h"

using namespace std;


CRod::CRod(){
    _sign = 0;
    axis = -1;
    coord(0) = 0;
    coord(1) =  0;
    coord(2) =  0;
}


CRod::CRod(const int& ax, const Eigen::Vector3d& initvec, const int& Nspheres, const double& spheredist, const double& boxsize){
    _sign = 1;
    axis = ax;
    coord = initvec;
    assert(coord(axis) == 0. && "Error in CRod init");
    Eigen::Vector3d spherepos = coord;
    spherepos(axis) = -boxsize;
    for (int i=0;i<Nspheres;i++){
        spherepos(2) += spheredist;
        spheres.push_back( CPolySphere(spherepos));
    }
}

