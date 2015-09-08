#ifndef CCONFIGURATION_H_
#define CCONFIGURATION_H_

#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <string>
#include <vector>
#include <array>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include <Eigen/Dense>
#include <Eigen/Cholesky>


#include "CPolymers.h"
#include "CPolySphere.h"


class CConfiguration {
    /*Class where all the configuration variables such as potRange etc. and also most functions for the
     * simulation are stored
     */
private:

    //SCALING
    double _timestep;         //This is the RESCALED timestep! timestep = dt * kT / (frictionCoeffcient * particlesize)
    double _mu_sto;
    double _pradius;     //particle size is most relevant for scaling! (default = 1)
    double _boxsize;          // ALWAYS define boxsize through particlesize due to scaling!
    double _epsilon;
	double _hpi_u;
	double _hpi_k;


    //EXPONENTIAL Potential
    double _potRange;         // Avoid values like 10/3 = 3.33333... ALWAYS define the Range of the exp potential through boxsize due to scaling!
    double _potStrength;      // rescaled U_0 for exponential Potential
    double _rodDistance;
    CPolymers _poly;


    //bool Parameters
    bool _ewaldCorr;             // true if the modified exponential potential version is used that is not 3*U_0 at the intersections.
    bool _LJPot;              // if true, then LJ interaction is calculated for steric hindrance
    bool _ranU;
    bool _hpi;
    bool _noLub;

    //COUNTERS AND INIT VALUES
    int _boxnumberXYZ[3];           //counter to calculate the actual position of the particle
    unsigned int _wallcrossings[3]; //counts how many times the particle has crossed the [entry, opposite, side] wall, while
                                            //travelling through the lattice
    int _entryside[3];            //records through which side and in which direction the particle last entered a box. For the face of the cube in the x-z
                                            //plane at y=0, it is entryside[0] = 1, in the x-y plane at z=L it is entryside[1] = -1!
    double _resetpos;
    Eigen::Vector3d _startpos;          //Stores where the particle starting position was. This is needed to calculate the mean square displacement
    Eigen::Vector3d _prevpos;           //Stores previous particle position before particle is moved.
    std::vector<std::vector<std::vector<int> > > _posHistoM;
    int _min, _max;        // parameters for determining up to which order neighboring rods are considered for the potential
    
    string _testcue;

    //Particle parameters
    Eigen::Vector3d _ppos;    //initialize particle position (DO IT LIKE resetpos FOR MOVEPARTICLEFOLLOW/RESET)
    double _upot;
    Eigen::Vector3d _f_mob;   //store mobility and stochastic force
    Eigen::Vector3d _f_sto;
    Eigen::Vector3d _Vdriftdt;
    Eigen::Vector3d _V0dt;
	
	//HI Paramters
	std::vector<CPolySphere> _polySpheres;
	Eigen::MatrixXd _mobilityMatrix;
	Eigen::Matrix3d _tracerMM;     
	Eigen::Matrix3d _resMNoLub;
    Eigen::Matrix3d _RMLub;
	bool _HI;
	double _polyrad;
	int _edgeParticles;
    double _lastcheck[3];
    double _cutoffMMsq;
    double _stericrSq;
    int _n_cellsAlongb;
    // LOAD OF CRAP I THINK! This term is a correction for the tracer displacement at every time step due to the repeated tracer particle. It is needed to ensure that the average velocity of all particles in the system is 0. It equals (1-1/N) where N is the number of particles int the simulation box. For reference, see Durlofsky1987a page 3333 - Fundamental solution for flow in porous media and comparison with the Brinkmann equation.
    //double _pbc_corr; 
	
	// Ewald sum parameters
	int _nmax;
	int _nkmax;
	double _alpha;
	double _r_cutoffsq;
	double _k_cutoff;
    bool _noEwald;
	
	//Lubrication parameters
	int _mmax;
	double _g[3];
    double _cutofflubSq;
    double _V;


    boost::mt19937 *m_igen;                      //generate instance of random number generator "twister".



private:
    void setRanNumberGen(double seed);
    void countWallCrossing(int crossaxis, int exitmarker);
    void calculateExpHPI(const double r, double& U, double& Fr);
    void calculateExpPotential(const double r, double& U, double& Fr);
    void modifyPot(double& U, double& Fr, double dist);
    void calcLJPot(const double r, double &U, double &dU);
    void initPosHisto();
	void initConstMobilityMatrix(bool Ewaldtest);
	Eigen::Matrix3d CholInvertPart (const Eigen::MatrixXd A);
    Eigen::Matrix3d Cholesky3x3(Eigen::Matrix3d mat);
    Eigen::Matrix3d invert3x3 (const Eigen::Matrix3d A);
	Eigen::Matrix3d realSpcSm( const Eigen::Vector3d & rij, const bool self, const double asq );
	Eigen::Matrix3d reciprocalSpcSm( const Eigen::Vector3d & rij, const double asq );
	Eigen::Matrix3d realSpcM(const double & rsq, const Eigen::Vector3d & rij, const double asq);
	Eigen::Matrix3d reciprocalSpcM(const double ksq, const Eigen::Vector3d & kij,  const double asq);
    Eigen::Matrix3d RotnePrager( const Eigen::Vector3d & rij, const double asq);
	
	Eigen::Matrix3d lub2p( Eigen::Vector3d rij, double rsq, unsigned int mmax );
	Eigen::Matrix3d lubricate( const Eigen::Vector3d & rij );
    Eigen::Vector3d midpointScheme(Eigen::Vector3d V0dt, Eigen::Vector3d F);
    void calcTracerMobilityMatrix(bool full);
    
    void report(std::string);

    template<typename T>
    string toString(const T& value){
        ostringstream oss;
        oss << value;
        return oss.str();
    }


public:
    CConfiguration();
    CConfiguration(
            double timestep,  double potRange,  double potStrength,  double boxsize, double rodDistance, const bool ewaldCorr, double psize,
            const bool posHisto, const bool steric, const bool ranU,  bool hpi, double hpi_u, double hpi_k, double polymersize);
    void resetParameters(double timestep, double potRange, double potStrength, double boxsize);
    void updateStartpos();
    void resetposition();
    int makeStep();
    int checkBoxCrossing();
    void calcStochasticForces();
    void calcMobilityForces();
    void saveXYZTraj(string name, const int& move, string a_w);
    void positionHistogram(double block, double possq, double pposprev, int l, int posHisto[]);

    double getPosVariance();
	std::vector<double> getppos();
    double get1DPosVariance(int dim);
    double getUpot(){ return _upot; }
  //  double getDisplacement();
    unsigned int getwallcrossings(int i){ return _wallcrossings[i]; }
    bool checkFirstPassage(double mfpPos, int dim);
    bool testOverlap();
    void moveBack();
    void addHistoValue();
    void printHistoMatrix(string folder);
    void checkDisplacementforMM();
    string getTestCue(){ return _testcue; };





};



#endif /* CCONFIGURATION_H_ */
