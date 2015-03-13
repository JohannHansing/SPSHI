#include "headers/CPolySphere.h"

using namespace std;

CPolySphere::CPolySphere(){
}

CPolySphere::CPolySphere( Eigen::Vector3d PolySpherePos, Eigen::Vector3d tracerPosition ){
    _position = PolySpherePos;
    _tracerVec = tracerPosition - PolySpherePos;
    _tracerdistSq = _tracerVec.transpose() * _tracerVec;  // compute square of distance between tracer and polysphere
}


void CPolySphere::updatePosition( Eigen::Vector3d new_pos ){
    _position = new_pos;
} 

void CPolySphere::updateTracerVec( Eigen::Vector3d tracerPosition ){
    _tracerVec = tracerPosition - _position;
    _tracerdistSq = _tracerVec.transpose() * _tracerVec;
}