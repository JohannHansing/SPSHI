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
#include <boost/filesystem.hpp>
#include <boost/timer.hpp>

#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/IterativeLinearSolvers>

#include "CPolymers.h"
#include "CPolySphere.h"
#include "CRod.h"
#include "parameter_structs.h"


#define ifdebug(x)
#define iftestEwald(x)   
#define iftestLub2p(x)
#define ifdebugPreComp(x)

class CConfiguration {
    /*Class where all the configuration variables such as potRange etc. and also most functions for the
     * simulation are stored
     */
public:
    //STRUCT
    struct sim_triggers _triggers;
    struct file_desc _files;
    struct model_param_desc _modelpar;

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
    bool _ranRod;
    bool _2DLattice;
    bool _HI2;
    bool _preEwald;
    
    //nan-Stuff
    bool _nanfound=false;
    int _nancount = 0;
    Eigen::Matrix3d _resMPrev = Eigen::Matrix3d::Identity();
    
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
    Eigen::Vector3d _vdisp;
    
    //Steric interaction parameters
    double _stericrSq;
    double _epsilonLJ = 1.;
    int _Nrods;// (n_cells+1)^2 would be for first order summation, but here we always do second order summation (n_cells+3)^2
    std::vector<double> _rSq_arr, _ri_arr, _rk_arr;

	//HI Paramters
	std::vector<CPolySphere> _polySpheres;
	Eigen::MatrixXd _mobilityMatrix;
	Eigen::Matrix3d _tracerMM;
	Eigen::Matrix3d _resMNoLub;
    Eigen::Matrix3d _RMLub;
    Eigen::MatrixXd _prevCG;
	bool _HI;
	double _polyrad;
    double _polydiamSq;
    double _Invpolyrad;
	int _edgeParticles;
	int _N_polyspheres;
    double _sphereoffset;
    double _lastcheck[3];
    double _cutoffMMsq;
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
    std::vector<Eigen::Vector3d> _kvec_arr;
    std::vector<Eigen::Matrix3d> _Mreciprocal_arr;
	//Lubrication parameters
	int _mmax;
	double _gX[3];
	double _gY[3];
	double _g[3];
    std::vector<double> _fXm;
    std::vector<double> _fYm;
    double _cutofflubSq;
    double _Vinv;
    Eigen::Matrix<double, 6, 6> _RP2p;
    std::vector<double> _fitpIs;
    std::vector<double> _fitprrs;
    std::array<double, 20> _minv;
    std::array<double, 20> _m2inv;
    
    
    // parameters for precomputing the long-range resistance matrices
    /*TODO list:
     * - check memory usage for n_cells > 2, 
     * - if neccessay, delete mobility matrix after initializing the lookup table, so that less memory is used by the program
     * - i could potentially make testoverlap much faster for _n_cellsAlongb > 1, if I use the fmod function to reduce the problem to one cell.
     * - set frac to higher value. size of table for frac = 0.2  is 8*9*25**3 * 10^-6 = 1.125 MB (megabytes), so not so large
    */
    double _frac = 0.025; //fraction of the cell width which is used for creating precalculated resistances default should be 0.025
    unsigned int _nxbins;
    std::array<std::array<std::array<Eigen::Matrix3d, 2>, 2>, 2> _T_arr;
    std::vector<std::vector<std::vector<Eigen::Matrix3d>>> _resM_precomp;
    
    void initTrafoMatrixes(){
        _T_arr[0][0][0] = Eigen::Matrix3d::Identity();
        _T_arr[1][0][0] = Eigen::Matrix3d::Identity();
        _T_arr[1][0][0](0,0) = -1;
        _T_arr[0][1][0] = Eigen::Matrix3d::Identity();
        _T_arr[0][1][0](1,1) = -1;
        _T_arr[0][0][1] = Eigen::Matrix3d::Identity();
        _T_arr[0][0][1](2,2) = -1;
        _T_arr[1][1][0] = - _T_arr[0][0][1];
        _T_arr[0][1][1] = - _T_arr[1][0][0];
        _T_arr[1][0][1] = - _T_arr[0][1][0];
        _T_arr[1][1][1] = - Eigen::Matrix3d::Identity(); 
    }

    Eigen::Matrix3d getTableResM(){
        // for my tranformations, T is self adjoint, i.e. T^-1 = T
        int ni[3] = {0,0,0};
        double bcell = _boxsize/_n_cellsAlongb; //width of the cell (should be 10)
        Eigen::Vector3d cellpos; //stores the position r of the particle relative to the cell
        Eigen::Vector3d subcellpos; //stores the position of the tracer in the subcell, to obtain the right mobilty matrix
        for (int i=0; i<3; i++){
            cellpos(i) = fmod(_ppos(i),bcell); //position of the particle inside the cell
            subcellpos(i) = cellpos(i); 
            if (cellpos(i) >= bcell/2){
                ni[i] = 1; //check, in which eigth of the cell the particle resides. This is needed for the transformation
                if (cellpos(i) == bcell/2) subcellpos(i) = abs(subcellpos(i) - bcell) -0.00001; // this case could otherwise lead to a segmentation fault
                else subcellpos(i) = abs(subcellpos(i) - bcell);
            }
        }
        ifdebugPreComp(cout << "cellpos\n" << cellpos << "\nsubcellpos\n" << subcellpos << endl;)
        //find the appropriate index for the stored matrices for the particle position cellpos
        int x = (int)(subcellpos(0) / (bcell * _frac));
        int y = (int)(subcellpos(1) / (bcell * _frac));
        int z = (int)(subcellpos(2) / (bcell * _frac));
        ifdebugPreComp(cout << "x y z:\n" << x << " " << y  << " " << z << "\n ni[]: \n" << ni[0] << " " << ni[1]  << " " << ni[2] << endl;)
        Eigen::Matrix3d resMT = _resM_precomp[x][y][z]; //get the resistance matrix in that spot
        Eigen::Matrix3d T = _T_arr[ni[0]][ni[1]][ni[2]];
        ifdebugPreComp(cout << "T * resMT * T\n" << T * resMT * T << endl;)
        return T * resMT * T;
    }
    
    void compareTableToEwald(){
        _noLub = true;
        Eigen::Matrix3d tmp = getTableResM();
        calcTracerMobilityMatrix(true);
        cout << "~~~~~~\nTableResM\n" << tmp << "\nEwaldResM\n" << _resMNoLub << endl;
        cout << "difference\n" << tmp - _resMNoLub << endl;
        cout << "relative difference\n" << (tmp - _resMNoLub).cwiseQuotient(_resMNoLub) << endl;
        cout << "relative to diagonal\n" << (tmp - _resMNoLub)/_resMNoLub(0,0) << endl;
    }
    
    double getlubfrac(){
        // return the fraction of the lubrication effect on the mobilitymatrix
        // of course, this only works if lubrication is used to calculate _tracerMM !
        double sumMnoLub  = invert3x3(_resMNoLub).sum();
        double sumMLub    = _tracerMM.sum(); //this is with lub
        return  1. - sumMLub/sumMnoLub;
    }

    void getLubDiffM(double &sumMdiff, double &sumMLub){
        sumMLub  = _tracerMM.sum();
        double sumMnoLub  = invert3x3(_resMNoLub).sum();
        sumMdiff = sumMnoLub-sumMLub;
    }
    
    void storeResistanceTable(){
        cout << "storing precomputed Ewald table to file.." << endl;
        
        ofstream resEwaldTable;
        resEwaldTable.open((_files.ewaldTable ).c_str());

        for (int n1 = 0; n1 < _nxbins; n1++){
            for (int n2 = 0; n2 < _nxbins; n2++){
                for (int n3 = 0; n3 < _nxbins; n3++){
                    for (int i=0; i<3; i++){
                        resEwaldTable << _resM_precomp[n1][n2][n3](i,0) << " " << _resM_precomp[n1][n2][n3](i,1) << " " << _resM_precomp[n1][n2][n3](i,2) << " ";
                    }
                    resEwaldTable << endl;
                }
            }
        }
        resEwaldTable.close();
    }
    
    void readResistanceMatrixTable(){
        // source http://stackoverflow.com/questions/3946558/c-read-from-text-file-and-separate-into-variable
        cout << "ewaldTablefile: " << _files.ewaldTable << endl;
        
        ifstream fin(_files.ewaldTable.c_str());
        _resM_precomp = vector<vector<vector<Eigen::Matrix3d>>> (_nxbins, vector<vector<Eigen::Matrix3d>>(_nxbins, vector<Eigen::Matrix3d>(_nxbins,Eigen::Matrix3d::Zero())));
        for (int n1 = 0; n1 < _nxbins; n1++){
            for (int n2 = 0; n2 < _nxbins; n2++){
                for (int n3 = 0; n3 < _nxbins; n3++){
                    //getline(inputFile, line);
                    //istringstream ss(line);
                    double m11, m12, m13, m21, m22, m23, m31, m32, m33;

                    fin >> m11 >> m12 >> m13 >> m21 >> m22 >> m23 >> m31 >> m32 >> m33;
                    _resM_precomp[n1][n2][n3] << m11, m12, m13, m21, m22, m23, m31, m32, m33;
                }
            }
        }
    }



    boost::mt19937 *m_igen;                      //generate instance of random number generator "twister".

    double zerotoone(){
        boost::uniform_01<boost::mt19937&> dist(*m_igen);
        return dist();
    }

    int ran_sign(){
        // the variate generator uses _igen (int rand number generator),
        // samples from uniform integer distribution 0, 1
        boost::variate_generator<boost::mt19937&, boost::uniform_int<>> zeroone(*m_igen, boost::uniform_int<>(0, 1));
    	return (zeroone() * 2) - 1; //this calculation makes value either -1 or 1
    }




private:
    void setRanNumberGen(double seed);
    void countWallCrossing(int crossaxis, int exitmarker);
    void calculateExpHPI(const double& r, double& U, double& Fr);
    void calculateExpPotential(const double& r, double& U, double& Fr);
    void modifyPot(double& U, double& Fr, const double dist);
    void calcLJPot(const double &rSq, double &U, double &dU);
    void initConstMobilityMatrix();
    Eigen::Matrix3d CholInvertPart (const Eigen::MatrixXd &A);
    Eigen::Matrix3d ConjGradInvert(const Eigen::MatrixXd &A);
    Eigen::Matrix3d Cholesky3x3(const Eigen::Matrix3d &mat);
    Eigen::Matrix3d invert3x3 (Eigen::Matrix3d m);
	Eigen::Matrix3d realSpcSm( const Eigen::Vector3d & rij, const bool &self, const double &asq );
    void initMreciprocalTracer();
	Eigen::Matrix3d reciprocalSpcSmTracer( const Eigen::Vector3d & rij );
	Eigen::Matrix3d reciprocalSpcSm( const Eigen::Vector3d & rij, const double &asq );
	Eigen::Matrix3d realSpcM(const double & rsq, const Eigen::Vector3d & rij, const double &asq);
	Eigen::Matrix3d reciprocalSpcM(const double &ksq, const Eigen::Vector3d & kij,  const double &asq);
    Eigen::Matrix3d RotnePrager( Eigen::Vector3d rij, double asq);
    Eigen::Matrix3d RPYamakawaPS(const Eigen::Vector3d & rij, const double asq);
    //Eigen::Matrix3d calcHI2MobMat( Eigen::Vector3d rij, double asq );
    void computeMobilityMatrixHere(Eigen::Vector3d rpos, double xinterval);
    void precomputeResistanceMatrix();


	Eigen::Matrix3d lub2p(const Eigen::Vector3d &rij, const double &rsq );
	Eigen::Matrix3d lubricate( const Eigen::Vector3d & rij );
    Eigen::Vector3d midpointScheme(const Eigen::Vector3d & V0dt,const  Eigen::Vector3d & F);
    void updateMobilityMatrix();
    void initPolySpheres();


    template<typename T>
    string toString(const T& value){
        ostringstream oss;
        oss << value;
        return oss.str();
    }
    
    unsigned int linecount(string file){
        std::ifstream f(file.c_str());
        std::string line;
        unsigned int i;
        for (i = 0; std::getline(f, line); ++i)
            ;
        return i;
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
    void calcMobForcesBeads();
    void saveXYZTraj(string name, const int& move, string a_w);
    void positionHistogram(double block, double possq, double pposprev, int l, int posHisto[]);

    double getPosVariance();
	std::vector<double> getppos();
    double get1DPosVariance(int dim);
    double getUpot(){ return _upot; }
  //  double getDisplacement();
    unsigned int getwallcrossings(int i){ return _wallcrossings[i]; }
    bool checkFirstPassage(double mfpPos, int dim);
    bool testOverlap(double stericrSq_preE=0);
    void moveBack();
    void initPosHisto();
    void addHistoValue();
    void printHistoMatrix(string folder);
    void checkDisplacementforMM();
    string getTestCue(){ return _testcue; };
    void calcTracerMobilityMatrix(bool full);

    void report(std::string);
    

    double getStepDisplacement(){
        return _vdisp.squaredNorm();
    }


    double _binv;
    Eigen::Vector3d minImage(Eigen::Vector3d rij){
        // returns disctance vector with minimal image convention.
        int abc[3];
        Eigen::Vector3d minvec = rij;
        for (int p=0;p<3;p++){
            abc[p]= minvec(p)*(_binv);
            minvec(p) -= abc[p] * _boxsize;
            abc[p]= minvec(p)*(_binv);
            minvec(p) -= abc[p] * _boxsize;
        }
        return minvec;
    }

private:
    void initLubStuff(double particlesize, double polymersize){
        const double lam = _polyrad/_pradius; //Jeffrey1984: lambda = a2/a1. Here, a1 is always the tracer
        const double c1 = pow(1+lam, -3);
        
        double lampow[11];
        for (int i=0;i<11;i++){
            lampow[i] = pow(lam,i);
        }

        _gX[0] = 2 * pow(lam, 2) * c1;
        _gX[1] = lam/5 * ( 1 + 7*lam + lam*lam ) * c1;
        _gX[2] = 1./42 * ( 1 + lam*18 - lampow[2]*29 + lampow[3]*18 + lampow[4]) * c1;

        _gY[1] = 4./15 * lam * (2 + lam + 2*lam*lam) * c1;
        _gY[2] = 2./375 * ( 16 - lam*45 + lampow[2]*58 - lampow[3]*45 + lampow[4]*16) * c1;
        
//         for (int i=0;i<3;i++){
//             cout << _gX[i] << "  --  " << _gY[i] << endl;
//         }

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
        double fy2 = 9./4*lam;
        double fy4 = 6*lam+81./16*lampow[2]+ 18*lampow[3];
        double fy6 = 4*lam + 54*lampow[2]+ 1241./64 *lampow[3 ]+ 81*lampow[4] + 72*lampow[5];
        double fy8 = 279*lampow[2] + 4261./8*lampow[3] + 126369./256*lampow[4] - 117./8*lampow[5] + 648*lampow[6] + 288*lampow[7];
        double fy10 = 1152*lampow[2] +7857./4*lampow[3] +98487./16*lampow[4] + 10548393./1024*lampow[5] +67617./8*lampow[6] - 351./2*lampow[7 ]+ 3888*lampow[8] + 1152*lampow[9];
        const double tmpArrY[] =  {fy0,fy2,fy4,fy6,fy8,fy10};
        _fYm.insert(_fYm.end(), &tmpArrY[0], &tmpArrY[ArrSize]);
        for (int i=0; i<ArrSize; i++){
             //cout << "fYm[i] =" << _fYm[i] << endl;
             //cout << "fXm[i] =" << _fXm[i] << endl;
            _fYm[i] /= pow(2*(1+lam),2*i);
            _fXm[i] /= pow(2*(1+lam),2*i);
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
        
        //m sum precalculations
        _minv[0] = 0; //zero, since this should not be used in any case
        _m2inv[0] = 0;
        for (int m = 1; m < 20; m++){
            _minv[m] = 1./m;
            _m2inv[m] =  1./(m*(m-1));
        }
        _m2inv[1]=-1; //IMPORTANT distinction, due to m_1! Check Jeffrey1984 if unclear
    }

    // ************ RANROD ************
    /*TODO :
     * Implement random number rods depending on n_rods.
     * Test if I can use C++11 on sheldon. If so, use:
     * for(auto& s: _rodvec[plane]){
     *    do smth with s;
     *}
     * instead of
     * (int i=0;i<_rodvec[plane].size();i++){
     *     do smth with _rodvec[plane][i]
     * }
     * for looping over arrys/vectors
     */
    std::array<vector<CRod> , 3> _rodvec; // vector to store polymer rods in cell, one vector stores polymers that are parallel to the same axis
    double _n_rods = 1;
    double _n_rods_max = 1;  // multiplier for maximum number of rods in one cell _n_max = _n_rods_max * _n_rods
    int _n_tries = (int) (ceil(_n_rods) * _n_rods_max + 0.001);
    double _reln = _n_rods / _n_tries;// _n_rods / n_tries is probability of placing one of n_tries new
    unsigned int _N_rods = 0;
    unsigned int avrods =0;
    unsigned int avcount =0;




    void copySphereRods(){
        _polySpheres.resize(0);
        for (int k = 0; k < 3; k++){
            for (int r=0;r<_rodvec[k].size();r++){
                for (int s=0;s<_rodvec[k][r].spheres.size(); s++){
                    //TODO Maybe do this with an &, or something, to copy reference to _rodvec[k].spheres[s]
                    _polySpheres.push_back( _rodvec[k][r].spheres[s] );
                }
            }
        }
        // Update the mobility matrix for new polySpheres
        updateMobilityMatrix();

        //cout << "_polySpheres.size() = " << _polySpheres.size() << endl;
        ifdebug(if (_polySpheres.size() != _rodvec[0][0].spheres.size() * (_rodvec[0].size() + _rodvec[1].size() + _rodvec[2].size())){cout << "Error: number of polyspheres doesnt correspond to number of rods\n";};)
    }

public:
    void initRodsVec(){
        Eigen::Vector3d initVec;
        int ortho1, ortho2, retry =0;

        double _spheredist = _boxsize/(2*_polyrad);
        int Nrods[3]; // number of rods in certain plane, i.e. parallel to a certain axis.

        for (int i=0;i<3;i++){
            // allocate enough space for vector resize.
            _rodvec[i].reserve(_n_rods_max * 9. * _n_rods ); //9 cells per plane
            //TODO if zerotoone()> ... MAKE Nrods[i] random
            Nrods[i] =  9. * _n_rods;// 3x3 cells in one plane -> 9*_n_rods
        }
        bool overlaps = true;
        while (overlaps==true){
            ++retry;
            cout << "RETRY for ranRods Init" << endl;
            for (int axis=0;axis<3;axis++){//axis 0 is x axis.
                ortho1 = axis+1; // No specific order of ortho1 and ortho2
                if (ortho1 == 3) ortho1 = 0;
                ortho2 = 3-(axis+ortho1);
                _rodvec[axis].clear();
                for (int i=0; i<Nrods[axis];i++){
                    for (int count=0; count < 200; count++){
                        overlaps = false;
                        // SO FAR SET INITAL RODS; SUCH THAT CENTRAL CELL * 1.7 IS EMPTY!
                        //TODO ran_sign()*.15*zerotoone() + 0.5 + 1.35 * ran_sign() is between .15+1.35+0.5 = 2 and -.15+1.35+0.5 = 1.7 or between .15-1.35+0.5 = -0.7 and -.15-1.35+0.5 = -1
                        const double centerspace = (_pradius+_polyrad) / _boxsize;
                        const double sidespace = (1.5 - centerspace);
                        //double halfsidespace = (1.5 - centerspace)/2;
                        //double middleofsidespace = halfsidespace + centerspace;
                        initVec(ortho1) = (0.5 + ran_sign()*(  centerspace + zerotoone()*sidespace)) *_boxsize;
                        initVec(ortho2) = (0.5 + ran_sign()*(  centerspace + zerotoone()*sidespace)) *_boxsize;
                        //initVec(ortho1) = (0.5 + ran_sign()*zerotoone()*halfsidespace + middleofsidespace * ran_sign() ) *_boxsize;
                        //initVec(ortho2) = (0.5 + ran_sign()*zerotoone()*halfsidespace + middleofsidespace * ran_sign() ) *_boxsize;
                        initVec(axis) = 0;
                        overlaps = testRodOverlap(initVec,axis,ortho1,ortho2,0);
                        if (overlaps==false) break;
                    }
                    if (overlaps==true){//if there is still overlap after the 'count' for loop -- break out of i loop, and do a retry.
                        if (retry == 500){
                            cout << "could not find suitable rod positions at init after 500 retries." << endl;
                            cout << "Rod numbers were:\n" << _rodvec[0].size() << " ++ " << _rodvec[1].size() << " ++ " << _rodvec[2].size() << endl;
                            throw 18;
                        }
                        else break;
                    }
                    CRod newRod = CRod(axis, initVec, 3*_edgeParticles, _sphereoffset, _boxsize );
                    _rodvec[axis].push_back(newRod);
                }
                if (overlaps==true) break;//if there is still overlap after the 'i' for loop -- break out of axis loop, and do a retry.
            }
        }
        _N_rods = _rodvec[0].size() + _rodvec[1].size() + _rodvec[2].size();
    }

    struct outOfRange{
        const int crossaxis;
        const double b;
        outOfRange(int ca, double boxsize) : crossaxis(ca), b(boxsize) {}
        bool operator()(const CRod& rod) const {
            return abs(rod.coord(crossaxis) - b/2.  ) > 1.5*b;
        }
    };

    void updateRodsVec(const int& crossaxis, const int& exitmarker){//exitmarker is -1 for negative direction, or 1 for positive
        //delete all polymers orthogonal to crossaxis, that are outside the box now
        //update other polymer positions
        ifdebug(cout << "=======================\ncrossaxis and exitmarker: "<< crossaxis << " " << exitmarker << endl;  cout << "_N_rods = " << _N_rods << endl;)
        int ortho[2];
        ortho[0] = crossaxis + 1;
        if (ortho[0] == 3) ortho[0] = 0;
        ortho[1] = 3-(crossaxis + ortho[0]);
        // shift positions of rods
        int plane;
        Eigen::Vector3d tmpvec;
        //First delete
        for (int oa=0;oa<2;oa++){
            plane = ortho[oa];
            //cout << "plane " << plane << endl;
            // More efficient way --- NOT WORKING
            for (int i=0;i<_rodvec[plane].size();i++){
                //shift rod positions parallel to crossaxis. plane is direction that the shifted rods are parallel to.
                _rodvec[plane].at(i).coord(crossaxis) += - exitmarker * _boxsize;
            }
            //efficient erase function: https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
            _rodvec.at(plane).erase( std::remove_if(_rodvec.at(plane).begin(), _rodvec.at(plane).end(), outOfRange(crossaxis,_boxsize)), _rodvec.at(plane).end() );
            // int nrods = _rodvec[plane].size();
            // for (int i=nrods-1;i>=0;i--){//need to count beckwards, due to erase function!
            //     //shift rod positions parallel to crossaxis. plane is direction that the shifted rods are parallel to.
            //     _rodvec[plane].at(i).coord(crossaxis) += - exitmarker * _boxsize;
            //     if (abs(_rodvec[plane].at(i).coord(crossaxis) - _boxsize/2.  ) > 1.5*_boxsize){
            //         // erase rods that have left the simulation box.
            //         _rodvec[plane].erase(_rodvec[plane].begin() + i);
            //     }
            // }
            for (int i=0;i<_rodvec[plane].size();i++){
                //shift rod positions parallel to crossaxis. plane is direction that the shifted rods are parallel to.
                _rodvec[plane].at(i).shiftspheres( crossaxis, - exitmarker * _boxsize);
            }
        }
        //Then create new
        for (int oa=0;oa<2;oa++){
            bool overlaps = true;
            plane = ortho[oa];
            int retry=0, sizeb4new = _rodvec[plane].size();
            while (overlaps==true){
                _rodvec[plane].resize(sizeb4new);//resize vector, so that in case of retry, the newly created rods are deleted
                ++retry;
                //for (int j=0;j<3*_n_tries;j++){// factor 3, since I reassign 3 cells per plane
                //TODO newrods -  Only create as many newrods as were deleted
                //TODO fixed rod number!!
                for (int j=sizeb4new;j< (9*_n_rods);j++){
                    if (zerotoone() < _reln ){
                        for (int count=0; count < 200; count++){
                            tmpvec = Eigen::Vector3d::Zero();//Reset
                            //in direction parallel to crossaxis, choose new position in side cell
                            tmpvec(crossaxis) = (zerotoone()  + exitmarker) * _boxsize;
                            int ortho2 = 3 - (plane + crossaxis);
                            // in direction orthogonal to both plane and crossaxis
                            tmpvec(ortho2) = (zerotoone() * 3 -1) * _boxsize;
                            //cout << _rodvec[plane].size() << endl;
                            overlaps = testRodOverlap(tmpvec,plane,crossaxis,ortho2,0);
                            if (overlaps==false) break;
                        }
                        if (overlaps==true){//if there is still overlap after the 'count' for loop -- break out of j loop, and do a retry.
                            if (retry == 300){
                                cout << "could not find suitable rod positions at update after 300 retries." << endl;
                                cout << "Rod numbers were:\n" << _rodvec[0].size() << " ++ " << _rodvec[1].size() << " ++ " << _rodvec[2].size() << endl;
                                throw 18;
                            }
                            else break;
                        }
                        _rodvec[plane].push_back(  CRod( plane, tmpvec, 3*_edgeParticles, _sphereoffset, _boxsize )  );
                        ifdebug(if (testRodTracerOverlap(tmpvec,plane) == true){  cout << "Overlap between new rod and tracer \nrodaxis = " << plane << " -- rodpos\n" << tmpvec << endl; };)
                    }
                    // NO AXIS LOOP HERE! if (overlaps==true) break;//if there is still overlap after the 'i' for loop -- break out of axis loop, and do a retry.
                }
            }
        }
        copySphereRods();

        _N_rods = _rodvec[0].size() + _rodvec[1].size() + _rodvec[2].size();
        avrods += _N_rods;
        avcount += 1;
    }

    bool testRodTracerOverlap(Eigen::Vector3d rodpos, int axis){
        Eigen::Vector3d testpos = _ppos;
        testpos(axis) = 0.;
        if ((testpos - rodpos).squaredNorm() < _stericrSq + 0.000001){
            //cout << "overlaps! ppos = " << _ppos(0) << ", " << _ppos(1) << ", " << _ppos(2) << endl;
            return true;
        }
        return false;
    }


    bool testRodOverlap(const Eigen::Vector3d& testpos, const int& rodaxis, const int& ortho1, const int& ortho2, const int& sizeb4new){
        // function to test, whether the new rod overlaps with existing rods
        // i_parallel takes care, that only  newly created rods parallel to the testpos rod get tested, since the other lie in different cells.
        double polydiam = _polyrad + _polyrad + 0.0001;
        for (int l = sizeb4new; l < _rodvec.at(rodaxis).size(); l++){
            if ((testpos - _rodvec[rodaxis].at(l).coord).squaredNorm() < _polydiamSq + 0.0001){
                return true;
            }
        }
        // for (int l = 0; l < _rodvec[ortho1].size(); l++){
        //     if (abs(testpos(ortho2) - _rodvec[ortho1].at(l).coord(ortho2)) < polydiam ){
        //         return true;
        //     }
        // }
        // for (int l = 0; l < _rodvec[ortho2].size(); l++){
        //     if (abs(testpos(ortho1) - _rodvec[ortho2].at(l).coord(ortho1)) < polydiam ){
        //         return true;
        //     }
        // }
        return false;
    }

    void testIfSpheresAreOnRods(){
        Eigen::Vector3d spherevec;
        double dot;
        for (int axis=0;axis<3;axis++){
            for (int l = 0; l < _rodvec[axis].size(); l++){
                spherevec=_rodvec[axis][l].spheres[0].pos - _rodvec[axis][l].spheres[1].pos;
                dot = spherevec.dot(_rodvec[axis][l].coord);
                //cout << "dot product: " << dot << endl;
                if (dot != 0){
                    cout << "Bad Sphere placement in rod: axis =" << axis << "  -- index = " << l << endl;
                }
                if (_rodvec[axis][l].coord(axis) != 0){
                    cout << "BAD ROD coord! " << endl;
                }
                Eigen::Vector3d tmpvec = Eigen::Vector3d::Zero();
                tmpvec(axis) = 1.;
                for (int s = 0; s < _rodvec[axis][l].spheres.size(); s++){
                    if ((_rodvec[axis][l].spheres[s].pos -_rodvec[axis][l].coord).cross(tmpvec) != Eigen::Vector3d::Zero()){
                        cout << "sphere does not lie on rod! " << endl;
                    }
                }
            }
        }
        if (_polySpheres.size() == _rodvec[0][0].spheres.size() * (_rodvec[0].size() + _rodvec[1].size() + _rodvec[2].size())  ){
            cout << "Rod and sphere number agree! " << endl;
        }
        else cout << "Rod and sphere number DO NOT agree! " << endl;
    }



    void printAvRods(){
        cout << "nrods in yz plane mean: " << avrods/(3*avcount) << endl;
    }

    void prinRodPos(int axis){
        for (int irod=0;irod<_rodvec[axis].size();irod++){
            double rx = _rodvec[axis][irod].coord(0);
            double ry =_rodvec[axis][irod].coord(1);
            double rz =_rodvec[axis][irod].coord(2);
            cout << ",[" << rx << "," << ry << "," << rz << "]";
        }
        cout << "]," << endl;
    }

    void overlapreport(){
        Eigen::Vector3d vrij;
        for (unsigned int j = 0; j < _polySpheres.size(); j++){
            vrij = (_ppos - _polySpheres[j].pos);
            if (vrij.squaredNorm() <= _stericrSq ){
                cout << vrij.norm() << " !!!!!!!!!!!!! OVErLAP!!  -- TRACER - POLYSPHERE\nsphereindex " << j << endl;
            }
        }
        for (unsigned int i = 0; i < _polySpheres.size(); i++){
            for (unsigned int j = i+1; j < _polySpheres.size(); j++){
                vrij = (_polySpheres[i].pos - _polySpheres[j].pos);
                if (vrij.squaredNorm() < pow(2*_polyrad,2)){
                    cout << "!!!!!!!!!!!!! OVErLAP!!  -- POLYSPHERE - POLYSPHERE" << endl;
                }
            }
        }
        Eigen::Vector3d testpos;
        for (int axis=0;axis<3;axis++){
            testpos = _ppos;
            testpos(axis) = 0.;
            for (int l = 0; l < _rodvec[axis].size(); l++){
                if ((testpos - _rodvec[axis][l].coord).squaredNorm() < _stericrSq + 0.000001){
                    cout << (testpos - _rodvec[axis][l].coord).norm() << " !!!!!!!!!!!!! OVErLAP!!  -- TRACER - ROD\nIn plane " << axis << "\nrod number " << l << endl;
                }
            }
        }
        for (int rodaxis = 0;rodaxis<3;rodaxis++){
            for (int t = 0; t < _rodvec[rodaxis].size(); t++){
                for (int l = t+1; l < _rodvec.at(rodaxis).size(); l++){
                    if ((_rodvec[rodaxis].at(t).coord - _rodvec[rodaxis].at(l).coord).squaredNorm() < _polydiamSq + 0.0001){
                        cout << "!!!!!!!!!!!!! OVErLAP!!  -- ROD - ROD" << endl;
                    }
                }
            }
        }
    }

    void moveBackReport(){
        cout << "_ppos: " << _ppos(0) << ", " << _ppos(1) << ", " << _ppos(2) << endl;
        cout << "_prevpos: " << _prevpos(0) << ", " << _prevpos(1) << ", " << _prevpos(2) << endl;
        cout << "Cholesky3x3(_RMLub)\n" << Cholesky3x3(_RMLub) << endl;
        cout << "_f_mob\n" << _f_mob << endl << "_f_sto\n" << _f_sto << endl;
    }


    void testLub2p();

    //TODO testEwald
    void testEwald(){
        cout << "b4 _tracerMM \n"<< _tracerMM << endl;
        // test Ewald summation for single monomer system
        _noLub = true;
        _EwaldTest=1;
        _edgeParticles = _EwaldTest;
        //---- Ewald summation --------
        _n_cellsAlongb = 1;
        _boxsize=10*_n_cellsAlongb;
        _binv=2/_boxsize;
        _Vinv = 1./pow( _boxsize, 3 );
        _sphereoffset = (_boxsize/_n_cellsAlongb) / _edgeParticles;
        for (int i = 0; i < 3; i++){
            _ppos(i) = (_boxsize/_n_cellsAlongb)/2;
        }
        initPolySpheres();
        initConstMobilityMatrix();
        calcTracerMobilityMatrix(true);
        cout << "_ppos\n" << _ppos << endl;
        cout << "spherepos\n" << _polySpheres[0].pos << endl;
        cout << "Mobility Matrix \n" << _mobilityMatrix << endl;
        cout << "Ewald _tracerMM \n"<< _tracerMM << endl;
        //----- No Ewald ----------
        _noEwald = true;
        _n_cellsAlongb = 1;
        _boxsize=10*_n_cellsAlongb;
        _binv=2./_boxsize;
        _Vinv = 1./pow( _boxsize, 3 );
        _sphereoffset = (_boxsize/_n_cellsAlongb) / _edgeParticles;
        for (int i = 0; i < 3; i++){
            _ppos(i) = _boxsize/2.;
        }
        initPolySpheres();
        initConstMobilityMatrix();
        calcTracerMobilityMatrix(true);
        cout << "No Ewald Mobility Matrix \n" << _mobilityMatrix << endl;
        cout << "No Ewald _tracerMM \n"<< _tracerMM << endl;
        
        //============================================
        // -------------- Normal System ---------------
        //============================================
        
        // TODO ABORT () !!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        abort();
        // TODO ABORT () !!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        _noEwald = false;
        _EwaldTest=0;
        _edgeParticles = (int) ( ( _boxsize/_n_cellsAlongb )/(2*_polyrad) + 0.0001);// round down
        //---- Ewald summation --------
        _n_cellsAlongb = 1;
        _boxsize=10*_n_cellsAlongb;
        _binv=2/_boxsize;
        _Vinv = 1./pow( _boxsize, 3 );
        _sphereoffset = (_boxsize/_n_cellsAlongb) / _edgeParticles;
        _ppos = Eigen::Vector3d(_boxsize/2.,_boxsize/2.,_boxsize/2.);
        initPolySpheres();
        initConstMobilityMatrix();
        calcTracerMobilityMatrix(true);
        cout << "Normal System Ewald _tracerMM \n"<< _tracerMM << endl;
        //----- No Ewald ----------
        _noEwald = true;
        for (int i = 0; i < 3; i++){
            _ppos(i) = _boxsize/2.;
        }
        initPolySpheres();
        initConstMobilityMatrix();
        calcTracerMobilityMatrix(true);
        cout << "Normal System No Ewald _tracerMM \n"<< _tracerMM << endl;
    }
    
    void testRealSpcM(){
        cout <<"\n==========\n" << endl;
        double asq = (_pradius*_pradius+_polyrad*_polyrad)/2;
        Eigen::Vector3d rij(0,0,(_polyrad + _pradius)+1);
        double rsq = rij.squaredNorm();
        cout << "rsq " << rsq << " - rij " << rij << " - asq " << asq << endl;
        cout << "Mreal \n" <<realSpcM( rsq, rij, asq);
        
        cout <<"\n==========\n\n" << endl;
        //abort();
    }
    
    void writePolySpheres(string directory){
        // this function writes the polyspheres positions to a file. This is only used for ranSpheres
        ofstream polyPos;
        polyPos.open((directory + "/Coordinates/polySpherepos.txt" ).c_str());

        for (unsigned int j = 0; j < _N_polyspheres; j++){
            polyPos << _polySpheres[j].pos(0) << " " << _polySpheres[j].pos(1) << " " << _polySpheres[j].pos(2) << " " << endl;
        }
        polyPos.close();
    }


};



#endif /* CCONFIGURATION_H_ */
