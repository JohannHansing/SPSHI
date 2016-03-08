/*
 * CPolymers.h
 *
 *  Created on: Aug 9, 2013
 *      Author: jh
 */

#ifndef CROD_H_
#define CROD_H_

#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "CPolySphere.h"

using namespace std;


class CRod {//TODO include CPolymers into this here, by setting random sign!
private:
    // ...

public:
    int _sign;
    int axis; // rod parallel to axis 0 = x, axis 1 = y, etc.
    Eigen::Vector3d coord; // Coordinate of rod in 2D plane orthogonal to axis. the coord parallel to axis is always 0. (see initiaion)
    std::vector<CPolySphere> spheres;

    CRod();
    CRod(const int& ax, const Eigen::Vector3d & initvec, const int& Nspheres, const double& spheredist, const double& boxsize);

    void shiftspheres(int crossaxis, double shift){
        for (int i=0; i<spheres.size(); i++){
            spheres[i].pos(crossaxis) += shift;
        }
    }
};



#endif /* CROD_H_ */
