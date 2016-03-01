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
#include <Eigen/IterativeLinearSolvers>

#include "CPolymers.h"
#include "CPolySphere.h"
#include "CRod.h"
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
    bool _ranRod;

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
    Eigen::MatrixXd _prevCG;
	bool _HI;
	double _polyrad;
        double _polydiamSq;
	int _edgeParticles;
        double _sphereoffset;
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
    
    double zerotoone(){
        boost::uniform_01<boost::mt19937&> dist(*m_igen);
        return dist();
    }



private:
    void setRanNumberGen(double seed);
    void countWallCrossing(int crossaxis, int exitmarker);
    void calculateExpHPI(const double& r, double& U, double& Fr);
    void calculateExpPotential(const double& r, double& U, double& Fr);
    void modifyPot(double& U, double& Fr, const double dist);
    void calcLJPot(const double &r, double &U, double &dU);
    void initConstMobilityMatrix();
    Eigen::Matrix3d CholInvertPart (const Eigen::MatrixXd &A);
    Eigen::Matrix3d ConjGradInvert(const Eigen::MatrixXd &A);
    Eigen::Matrix3d Cholesky3x3(const Eigen::Matrix3d &mat);
    Eigen::Matrix3d invert3x3 (const Eigen::Matrix3d &A);
	Eigen::Matrix3d realSpcSm( const Eigen::Vector3d & rij, const bool &self, const double &asq );
	Eigen::Matrix3d reciprocalSpcSm( const Eigen::Vector3d & rij, const double &asq );
	Eigen::Matrix3d realSpcM(const double & rsq, const Eigen::Vector3d & rij, const double &asq);
	Eigen::Matrix3d reciprocalSpcM(const double &ksq, const Eigen::Vector3d & kij,  const double &asq);
    Eigen::Matrix3d RotnePrager( const Eigen::Vector3d & rij, const double & asq);

	Eigen::Matrix3d lub2p( Eigen::Vector3d &rij, double &rsq );
	Eigen::Matrix3d lubricate( const Eigen::Vector3d & rij );
    Eigen::Vector3d midpointScheme(Eigen::Vector3d & V0dt, Eigen::Vector3d & F);
    void calcTracerMobilityMatrix(const bool& full);
    void updateMobilityMatrix();
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
        Eigen::Vector3d rij_rem;
        for (int i=0;i<3;i++){
            rij_rem(i) = remainder(rij(i),_boxsize);
        }
        return rij_rem;
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
    double _ntry = 1.;  // multiplier for maximum number of rods in one cell _n_max = _ntry * _n_rods
    unsigned int avrods =0;
    unsigned int avcount =0;
    
        
    void copySphereRods(){
        _polySpheres.clear();
        for (int k = 0; k < 3; k++){
            for (int r=0;r<_rodvec[k].size();r++){
                for (int s=0;s<_rodvec[k][r].spheres.size(); s++){
                    //TODO Maybe do this with an &, or something, to copy reference to _rodvec[k].spheres[s]
                    _polySpheres.push_back( _rodvec[k][r].spheres[s] );
                }
            }
        }
        
        //cout << "_polySpheres.size() = " << _polySpheres.size() << endl;
        assert(_polySpheres.size() == _rodvec[0][0].spheres.size() * (_rodvec[0].size() + _rodvec[1].size() + _rodvec[2].size())    && "Error: number of polyspheres doesnt correspond to number of rods");
    }
    
public:
    int ran_sign(){
    // the variate generator uses _igen (int rand number generator),
    // samples from uniform integer distribution 0, 1
        boost::variate_generator<boost::mt19937&, boost::uniform_int<>> zeroone(*m_igen, boost::uniform_int<>(0, 1));
	return (zeroone() * 2) - 1; //this calculation makes value either -1 or 1 
}

    void initRodsVec(){
        Eigen::Vector3d initVec;
        int ortho1, ortho2, retry =0;
        
        double _spheredist = _boxsize/(2*_polyrad);
        int Nrods[3]; // number of rods in certain plane, i.e. parallel to a certain axis.
        
        for (int i=0;i<3;i++){
            // allocate enough space for vector resize.
            _rodvec[i].reserve(_ntry * 9. * _n_rods ); //9 cells per plane
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
                        double centerspace = (_pradius+_polyrad) / _boxsize;
                        double halfsidespace = (1.5 - centerspace)/2;
                        double middleofsidespace = halfsidespace + centerspace;
                        initVec(ortho1) = (ran_sign()*zerotoone()*halfsidespace + 0.5 + middleofsidespace * ran_sign() ) *_boxsize;
                        initVec(ortho2) = (ran_sign()*zerotoone()*halfsidespace + 0.5 + middleofsidespace * ran_sign() ) *_boxsize;
                        initVec(axis) = 0;
                        overlaps = testRodOverlap(initVec,axis,ortho1,ortho2,_rodvec[axis].size());
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
    }
    
    void updateRodsVec(const int& crossaxis, const int& exitmarker){//exitmarker is -1 for negative direction, or 1 for positive
        //delete all polymers orthogonal to crossaxis, that are outside the box now
        //update other polymer positions
        int n_tries = (int) (ceil(_n_rods) * _ntry + 0.001); 
        int ortho[2] = {1,2};
        if (crossaxis == 1){
            ortho[0]=2; 
            ortho[1]=0;
        }
        else if (crossaxis == 2){
            ortho[0]=0; 
            ortho[1]=1;
        }
        // shift positions of rods
        int plane;
        bool overlaps = true;
        int retry=0, newrods = 0;
        Eigen::Vector3d tmpvec;
        for (int oa=0;oa<2;oa++){
            plane = ortho[oa];
            //cout << "plane " << plane << endl;
            int nrods = _rodvec[plane].size();
            for (int i=nrods-1;i>=0;i--){//need to count beckwards, due to erase function!
                //shift rod positions parallel to crossaxis. plane is direction that the shifted rods are parallel to.
                _rodvec[plane][i].coord[crossaxis] -= exitmarker * _boxsize;
                if (abs(_rodvec[plane][i].coord[crossaxis] - _boxsize/2.  ) > 1.5*_boxsize){
                    // erase rods that have left the simulation box.
                    //cout << _rodvec[plane].size() << " XXX ";
                    _rodvec[plane].erase(_rodvec[plane].begin() + i);
                    //cout << _rodvec[plane].size() << endl;
                }
            }
            double reln = _n_rods / n_tries;// _n_rods / n_tries is probability of placing one of n_tries new
            assert((reln < 1.) && "Error: reln in updateRodsVec must be smaller than 1!");
            // erase the first 3 elements:
            while (overlaps==true){
                _rodvec[plane].erase (_rodvec[plane].end()-newrods,_rodvec[plane].end());
                newrods=0;
                ++retry;
                //cout << "RETRY for ranRods update" << endl;
                for (int j=0;j<3*n_tries;j++){// factor 3, since I reassign 3 cells per plane
                    if (zerotoone() < reln ){  
                        retry = 0;
                        for (int count=0; count < 200; count++){
                            overlaps = false;
                            tmpvec = Eigen::Vector3d::Zero();//Reset
                            //in direction parallel to crossaxis, choose new position in side cell 
                            tmpvec(crossaxis) = (zerotoone()  + exitmarker) * _boxsize;
                            int ortho2 = 3 - (plane + crossaxis);
                            // in direction orthogonal to both plane and crossaxis
                            tmpvec(ortho2) = (zerotoone() * 3 -1) * _boxsize;
                            //cout << _rodvec[plane].size() << endl;
                            overlaps = testRodOverlap(tmpvec,plane,crossaxis,ortho2,newrods);
                            if (overlaps==false) break;
                        }
                        if (overlaps==true){//if there is still overlap after the 'count' for loop -- break out of j loop, and do a retry.
                            if (retry == 100){
                                cout << "could not find suitable rod positions at update after 100 retries." << endl;
                                cout << "Rod numbers were:\n" << _rodvec[0].size() << " ++ " << _rodvec[1].size() << " ++ " << _rodvec[2].size() << endl;
                                throw 18;
                            }
                            else break;
                        }
                        _rodvec[plane].push_back(  CRod( plane, tmpvec, 3*_edgeParticles, _sphereoffset, _boxsize )  );
                        ++newrods;
                    }
                    // NO AXIS LOOP HERE! if (overlaps==true) break;//if there is still overlap after the 'i' for loop -- break out of axis loop, and do a retry.
                }
            }
        }
        copySphereRods();
        
        avrods += _rodvec[0].size() + _rodvec[1].size() + _rodvec[2].size();
        avcount += 1;
    }
    
    
    bool testRodOverlap(Eigen::Vector3d& testpos, const int& rodaxis, const int& ortho1, const int& ortho2, const int i_parallel){
        // function to test, whether the new rod overlaps with existing rods
        // i_parallel takes care, that only  newly created rods parallel to the testpos rod get tested, since the other lie in different cells.
        double polydiam = _polyrad + _polyrad + 0.000001;
        for (int l = _rodvec[rodaxis].size() - i_parallel; l < _rodvec[rodaxis].size(); l++){
            if ((testpos - _rodvec[rodaxis][l].coord).squaredNorm() < _polydiamSq + 0.000001){
                return true;
            }
        }
        for (int l = 0; l < _rodvec[ortho1].size(); l++){
            if (abs(testpos(ortho2) - _rodvec[ortho1][l].coord(ortho2)) < polydiam ){
                return true;
            }
        }
        for (int l = 0; l < _rodvec[ortho2].size(); l++){
            if (abs(testpos(ortho1) - _rodvec[ortho2][l].coord(ortho1)) < polydiam ){
                return true;
            }
        }
        return false;
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
    

};



#endif /* CCONFIGURATION_H_ */
