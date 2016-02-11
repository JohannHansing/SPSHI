/*
 * CPolySphere.h
 *
 *  Created on: Mar 11, 2015
 *      Author: johann hansing
 */

#ifndef CPOLYSPHERE_H_
#define CPOLYSPHERE_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <time.h>
#include <Eigen/Dense>

using namespace std;

class CPolySphere {
private:
    double _tracerdistSq;  // distance to tracer - updated during testOverlap()
    Eigen::Vector3d _tracerVec; // Vector pointing from polymerSphere to tracer particle - updated during testOverlap() NONONONO not anymore!!

public:
    CPolySphere();
    CPolySphere( Eigen::Vector3d PolySpherePos );
    Eigen::Vector3d pos;  //fixed length array for position, this allocates less memory than a (dynamic) vector     

};



#endif /* CPOLYSPHERE_H_ */
