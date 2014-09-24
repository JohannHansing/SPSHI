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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp> 


#include "CPolymers.h"



namespace ublas = boost::numeric::ublas;

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
    bool _potMod;             // true if the modified exponential potential version is used that is not 3*U_0 at the intersections.
    bool _LJPot;              // if true, then LJ interaction is calculated for steric hindrance
    bool _ranU;
    bool _hpi;

    //COUNTERS AND INIT VALUES
    int _boxnumberXYZ[3];           //counter to calculate the actual position of the particle
    unsigned int _wallcrossings[3]; //counts how many times the particle has crossed the [entry, opposite, side] wall, while
                                            //travelling through the lattice
    int _entryside[3];            //records through which side and in which direction the particle last entered a box. For the face of the cube in the x-z
                                            //plane at y=0, it is entryside[0] = 1, in the x-y plane at z=L it is entryside[1] = -1!
    double _resetpos;
    double _startpos[3];          //Stores where the particle starting position was. This is needed to calculate the mean square displacement
    double _prevpos[3];           //Stores previous particle position before particle is moved.
    std::vector<std::vector<std::vector<int> > > _posHistoM;

    int _min, _max;        // parameters for determining up to which order neighboring rods are considered for the potential

    //Particle parameters
    double _ppos[3];    //initialize particle position (DO IT LIKE resetpos FOR MOVEPARTICLEFOLLOW/RESET)
    double _upot;
    double _f_mob[3];   //store mobility and stochastic force
    ublas::vector<double> _f_sto;
	
	//HI Paramters
	std::vector<std::array<double, 3> > _epos;
	ublas::symmetric_matrix<double> _mobilityMatrix;
	ublas::symmetric_matrix<double> _tracerMM;     
	ublas::symmetric_matrix<double> _resMNoLub;
	bool _HI;
	double _polyrad;
	int _edgeParticles;
	
	
	// Ewald sum parameters
	int _nmax;
	int _nkmax;
	double _alpha;
	double _r_cutoffsq;
	double _k_cutoff;
	
	//Lubrication parameters
	int _mmax;
	double _g[3];


    boost::mt19937 *m_igen;                      //generate instance of random number generator "twister".



private:
    void setRanNumberGen(double seed);
    void countWallCrossing(int crossaxis, int exitmarker);
    void calculateExpHPI(const double r, double& U, double& Fr);
    void calculateExpPotential(const double r, double& U, double& Fr);
    void calculateExpPotentialMOD(const double r, double& U, double& Fr, int index);
    void modifyPot(double& U, double& Fr, double dist);
    void calcLJPot(const double r, double &U, double &dU);
    void initPosHisto();
	template<class T> 
	bool InvertMatrix(const ublas::matrix<T>& input, ublas::matrix<T>& inverse);
	ublas::matrix<double> RotnePragerDiffRad(const double & r, const double & rsq, const ublas::vector<double> & rij);
	ublas::matrix<double> RotnePrager(const double & r, const double & rsq, const ublas::vector<double> & rij);
	void initConstMobilityMatrix();
	template<class MATRIX>
	bool CholInvertPart (const MATRIX& A, MATRIX& partInv) ;
	ublas::matrix<double> realSpcSm( const ublas::vector<double> & rij, const bool self, const double asq );
	ublas::matrix<double> reciprocalSpcSm( const ublas::vector<double> & rij, const double asq );
	ublas::matrix<double> realSpcM(const double & rsq, const ublas::vector<double> & rij, const bool self, const double asq);
	ublas::matrix<double> reciprocalSpcM(const double ksq, const ublas::vector<double> & kij,  const double asq);
	
	ublas::symmetric_matrix<double> lub2p( ublas::vector<double> rij, double rsq, unsigned int mmax );
	ublas::symmetric_matrix<double> lubricate( const ublas::vector<double> & rij );




public:
    CConfiguration();
    CConfiguration(
            double timestep,  double potRange,  double potStrength,  double boxsize, double rodDistance, const bool potMod, double psize,
            const bool posHisto, const bool steric, const bool ranU,  bool hpi, double hpi_u, double hpi_k, double polymersize);
    void resetParameters(double timestep, double potRange, double potStrength, double boxsize);
    void updateStartpos();
    void resetposition();
    void makeStep();
    void checkBoxCrossing();
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
	void calcTracerMobilityMatrix(bool full);





};



#endif /* CCONFIGURATION_H_ */
