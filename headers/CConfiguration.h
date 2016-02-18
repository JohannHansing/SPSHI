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
#include <boost/timer.hpp>

#include <Eigen/Dense>
#include <Eigen/Cholesky>


#include "CPolymers.h"
#include "CPolySphere.h"
#include "parameter_structs.h"


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
    bool _ranSpheres;
    bool _trueRan;
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
    //std::vector<std::vector<std::vector<int> > > _posHistoM;
    int _posHistoM[100][100][100];
    int _min, _max;        // parameters for determining up to which order neighboring rods are considered for the potential
    int _EwaldTest;

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
    bool _fitRPinv;
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
    std::vector<double> _fXm;
    std::vector<double> _fYm;
    double _cutofflubSq;
    double _V;
    Eigen::Matrix<double, 6, 6> _RP2p;
    std::vector<double> _fitpIs;
    std::vector<double> _fitprrs;


    boost::mt19937 *m_igen;                      //generate instance of random number generator "twister".



private:
    void setRanNumberGen(double seed);
    void countWallCrossing(int crossaxis, int exitmarker);
    void calculateExpHPI(const double r, double& U, double& Fr);
    void calculateExpPotential(const double r, double& U, double& Fr);
    void modifyPot(double& U, double& Fr, double dist);
    void calcLJPot(const double r, double &U, double &dU);
	void initConstMobilityMatrix();
	Eigen::Matrix3d CholInvertPart (const Eigen::MatrixXd A);
    Eigen::Matrix3d Cholesky3x3(Eigen::Matrix3d mat);
    Eigen::Matrix3d invert3x3 (const Eigen::Matrix3d A);
	Eigen::Matrix3d realSpcSm( const Eigen::Vector3d & rij, const bool self, const double asq );
	Eigen::Matrix3d reciprocalSpcSm( const Eigen::Vector3d & rij, const double asq );
	Eigen::Matrix3d realSpcM(const double & rsq, const Eigen::Vector3d & rij, const double asq);
	Eigen::Matrix3d reciprocalSpcM(const double ksq, const Eigen::Vector3d & kij,  const double asq);
    Eigen::Matrix3d RotnePrager( const Eigen::Vector3d & rij, const double asq);

	Eigen::Matrix3d lub2p( Eigen::Vector3d rij, double rsq );
	Eigen::Matrix3d lubricate( const Eigen::Vector3d & rij );
    Eigen::Vector3d midpointScheme(Eigen::Vector3d V0dt, Eigen::Vector3d F);
    void calcTracerMobilityMatrix(bool full);
    void initPolySpheres();

    void report(std::string);

    template<typename T>
    string toString(const T& value){
        ostringstream oss;
        oss << value;
        return oss.str();
    }
    
    



public:
    CConfiguration();
    CConfiguration( double timestep, model_param_desc modelpar, sim_triggers triggers, file_desc files);
    void resetParameters(double timestep, double potRange, double potStrength, double boxsize);
    void updateStartpos();
    int resetposition();
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
    void initPosHisto();
    void addHistoValue();
    void printHistoMatrix(string folder);
    void checkDisplacementforMM();
    string getTestCue(){ return _testcue; };


    Eigen::Vector3d minImage(Eigen::Vector3d rij){
        // returns disctance vector with minimal image convention.
        // For info - Check wikipedia
        double bhalf = _boxsize/2.;
        for (int i=0;i<3;i++){
            if (rij(i) > bhalf)          rij(i) -= _boxsize;
            else if (rij(i) <= -bhalf)   rij(i) += _boxsize;
        }
        return rij;
    }

private:
    void initLubStuff(double particlesize, double polymersize){
        const double lam = _polyrad/_pradius;
        const double c1 = pow(1+lam, -3);

        _g[0] = 2 * pow(lam, 2) * c1;
        _g[1] = lam/5 * ( 1 + 7*lam + lam*lam ) * c1;
        _g[2] = 1/42 * ( 1 + lam*(18 - lam*(29 + lam*(18 + lam)))) * c1;

        double lampow[11];
        for (int i=0;i<11;i++){
            lampow[i] = pow(lam,i);
        }

        double f0 = 1;
        double f1 = 3*lam;
        double f2 = 9*lam;
        double f3 = -4*lam+27*lampow[2]-4*lampow[3];
        double f4 = -24*lam + 81*lampow[2] + 36*lampow[3];
        double f5 = 72*lampow[2] + 243*lampow[3] + 72*lampow[4];
        double f6 = 16*lam + 108*lampow[2]+ 281*lampow[3] + 648*lampow[4]+ 144*lampow[5];
        double f7 = 1;
        double f8 = 576*lampow[2] + 4848*lampow[3]+ 5409*lampow[4] + 4524*lampow[5] + 3888*lampow[6]+ 576*lampow[7];
        double f9 = 1;
        double f10 = 2304*lampow[2] + 20736*lampow[3]+ 42804*lampow[4]+ 115849*lampow[5]+ 76176*lampow[6] + 39264*lampow[7]+ 20736*lampow[8] + 2304*lampow[9];
        // THe next bit looks ugly but I did it like this for a reason: http://stackoverflow.com/questions/259297/how-do-you-copy-the-contents-of-an-array-to-a-stdvector-in-c-without-looping
        const double tmpArrX[] =  {f0,f2,f4,f6,f8,f10};
        unsigned ArrSize = 6;
        _fXm.insert(_fXm.end(), &tmpArrX[0], &tmpArrX[ArrSize]);

        double fy0 = 1;
        double fy2 = 9/4*lam;
        double fy4 = 6*lam+81/16*lampow[2]+ 18*lampow[3];
        double fy6 = 4*lam + 54*lampow[2]+ 1241/64 *lampow[3 ]+ 81*lampow[4] + 72*lampow[5];
        double fy8 = 279*lampow[2] + 4261/8*lampow[3] + 126369/256*lampow[4] - 117/8*lampow[5] + 648*lampow[6] + 288*lampow[7];
        double fy10 = 1152*lampow[2] +7857/4*lampow[3] +98487/16*lampow[4] + 10548393/1024*lampow[5] +67617/8*lampow[6] - 351/2*lampow[7 ]+ 3888*lampow[8] + 1152*lampow[9];
        const double tmpArrY[] =  {fy0,fy2,fy4,fy6,fy8,fy10};
        _fYm.insert(_fYm.end(), &tmpArrY[0], &tmpArrY[ArrSize]);
        for (int i=0; i<ArrSize; i++){
            _fYm[i] /= pow(2*(1+lam),2*i);
            _fXm[i] /= pow(2*(1+lam),2*i);
            // cout << "fYm[i] =" << _fYm[i] << endl;
            // cout << "fXm[i] =" << _fXm[i] << endl;
        }
        
        
        // Here, I create the matrix to store the two-particle Rotne-Prager calculation which needs to be inverted and subtracted from the two-particle resistance matrix (Brady1988)
        _RP2p = Eigen::MatrixXd::Identity(6,6);
        _RP2p.block<3,3>(3,3) = _pradius/_polyrad * Eigen::Matrix3d::Identity();
        
        
        // reading two particle Rotne Prager fit stuff from file
        std::vector<double> fitpIs;
        std::string filename = "fits/fitp" + toString(particlesize) + "a" + toString(polymersize) + ".txt";
        std::ifstream infile(filename.c_str());
        std::string line;
        while (std::getline(infile, line)){
            std::istringstream iss(line);
            double fitpI, fitprr;
            if (!(iss >> fitpI >> fitprr)) { break; } // error

            _fitpIs.push_back(fitpI);
            _fitprrs.push_back(fitprr);
        }
    }

};



#endif /* CCONFIGURATION_H_ */
