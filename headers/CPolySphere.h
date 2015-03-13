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
    Eigen::Vector3d _position;  //fixed length array for position, this allocates less memory than a (dynamic) vector 
    double _tracerdistSq;  // distance to tracer - updated during testOverlap()
    Eigen::Vector3d _tracerVec; // Vector pointing from polymerSphere to tracer particle - updated during testOverlap() NONONONO not anymore!!

public:
    CPolySphere();
    CPolySphere( Eigen::Vector3d PolySpherePos, Eigen::Vector3d tracerPosition );
    Eigen::Vector3d getPosition(){ return _position; }
    double getTracerdistsq(){ return _tracerdistSq; }
    Eigen::Vector3d getTracerVec(){ return _tracerVec; }
    
    void updatePosition( Eigen::Vector3d new_pos );   //Not needed yet for rigid polymers
    void updateTracerVec( Eigen::Vector3d tracerPosition );

};



#endif /* CPOLYSPHERE_H_ */
