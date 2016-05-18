#include "headers/CConfiguration.h"


using namespace Eigen;

using namespace std;



const double _6root2 = 1.122462;
const double _srqtPi = sqrt (M_PI);


CConfiguration::CConfiguration(){
}

CConfiguration::CConfiguration( double timestep, model_param_desc modelpar, sim_triggers triggers, file_desc files)
        {
    // seed = 0:  use time, else any integer
    // init random number generator
    setRanNumberGen(0);
    _potRange = modelpar.urange;
    _potStrength = modelpar.ustrength;
    _pradius = modelpar.particlesize/2;   //_pradius is now the actual radius of the particle. hence, I need to change the definition of the LJ potential to include (_pradius +_polyrad)   -  or maybe leave LJ pot out
    _polyrad = modelpar.polymersize / 2;   //This is needed for testOverlap for steric and HI stuff !!
    _polydiamSq = pow(modelpar.polymersize,2);
    _Invpolyrad = 1./_polyrad;
    _fitRPinv = triggers.fitRPinv;
    _boxsize = modelpar.boxsize;
    _n_cellsAlongb = modelpar.n_cells;
    _resetpos = (_boxsize/_n_cellsAlongb)/2; // puts the particle in the center of the cell at the origin
    //_resetpos = (_boxsize)/2;
    _timestep = timestep;
    _rodDistance = modelpar.rodDist;
    _ranSpheres = triggers.ranSpheres;
    _trueRan = triggers.trueRan;
    _ranRod = triggers.ranRod;
    _noLub = triggers.noLub;
    _LJPot = (triggers.includeSteric == false) && (modelpar.particlesize != 0);
    _ranU = triggers.ranPot;
    _poly = CPolymers();
    _upot = 0;
    _f_mob = Vector3d::Zero();
    _f_sto = Vector3d::Zero();
    _tracerMM = Matrix3d::Identity();
    _mu_sto = sqrt( 2 * _timestep );                 //timestep for stochastic force
    _hpi = triggers.hpi;
    _hpi_u = modelpar.hpi_u;
    _hpi_k = modelpar.hpi_k;
    _binv = 2./_boxsize;//this needs to be 2/b apparently
    for (int i = 0; i < 3; i++){
        _ppos(i) = _resetpos;
        _startpos(i) = _resetpos;
        _entryside[i] = 0;
        _wallcrossings[i] = 0;
        _boxnumberXYZ[i] = 0;
        _prevpos(i) = _resetpos;
        _lastcheck[i] = _startpos(i);
    }


    // This is for inclusion of 2nd Order rods if k is 0.2b or larger
    _min = -1, _max = 3;
    if (_ranU || _hpi || (_potRange < 2)){
        _min = 0;
        _max = 2;
    }

    //Ewald sum stuff
    _nmax = modelpar.nmax; // This corresponds to the r_cutoff = _nmax * _boxsize
    _alpha = 1 * sqrt(M_PI) / _boxsize; // This value for alpha corresponds to the suggestion in Beenakker1986
    _k_cutoff = 2. * _alpha * _alpha * _nmax * _boxsize;   /* This corresponds to suggestion by Jain2012 ( equation 15 and 16 ). */
    _nkmax = (int) (_k_cutoff * _boxsize / (2. * M_PI) + 0.5);   /* The last bit (+0.5) may be needed to ensure that the next higher integer
                                                           * k value is taken for kmax */
    _r_cutoffsq = pow(_nmax * _boxsize, 2);   // Like Jain2012, r_cutoff is chosen such that exp(-(r_cutoff * alpha)^2) is small

    // lubrication stuff
    initLubStuff(modelpar.particlesize,modelpar.polymersize);
    _cutofflubSq = pow(modelpar.lubcutint*(_polyrad + _pradius),2);
    _stericrSq = pow(_pradius + _polyrad, 2);


    // init HI vectors matrices, etc
    // Configurations
    _EwaldTest = modelpar.EwaldTest; // Predefine _edgeParticles. Ewaldtest = 0 runs normal. Ewaldtest = 1 runs the program with only spheres in the corners of the cells, i.e. _edgeParticles = 1, EwaldTest = 2 with 2 edgeparticles, and so on
    _noEwald = false;       // noEwald to use normal Rotne Prager instead of Ewald summed one
    _2DLattice = false;     //This creates a 2D lattice like in Phillips1990. I will be able to compare my data to settle whether my HI is correct.

    _Vinv = 1./pow( _boxsize, 3 );
    _cutoffMMsq = pow(0.05*_boxsize/_n_cellsAlongb,2);
    if (_ranSpheres && _boxsize > 10) _cutoffMMsq = pow(0.025*_boxsize,2);
    if (_polyrad != 0) _HI = true;
    if (_EwaldTest != 0) _edgeParticles = _EwaldTest;
    else _edgeParticles = (int) ( ( _boxsize/_n_cellsAlongb )/modelpar.polymersize + 0.0001);// round down
    _sphereoffset = (_boxsize/_n_cellsAlongb) / _edgeParticles;
    if (_ranRod) initRodsVec();//Needs to be before initPolySpheres!
    initPolySpheres();
    if (_HI) {
        _LJPot = false;
    // THIS NEEDS TO COME LAST !!!!!!!
        initConstMobilityMatrix();
        calcTracerMobilityMatrix(true);
    //_pbc_corr = 1 - 1/(_polySpheres.size()+1); // assign system size dependent value to pbc correction for tracer displacement in porous medium
    }

    // TEST CUE to modify the directory the output data is written to!!
    _testcue = "";
    //TODO newrods -  Only create as many newrods as were deleted
    if (_ranRod) cout << "NOTE: Fixed number of rods in plane as 9*nrods!" << endl;
    if (_ranRod) _testcue += "/RPYOverlap/fixnrods1";
    if ( _noEwald ) _testcue += "/noEwald";
    if (_2DLattice) _testcue += "/2DLattice";
    if ( _EwaldTest > 0 ) _testcue += "/EwaldTest" + toString(_EwaldTest);
    if ( _n_cellsAlongb != 1 ){
        _testcue += "/n" + toString(_n_cellsAlongb);
        cout << "Set _edgeParticles = " << _edgeParticles << endl; //<< " --- and _mobilityMatrix.rows() = " <<_mobilityMatrix.rows() << endl;
    }
    if (!_testcue.empty()) cout << "***********************************\n****  WARNING: String '_testcue' is not empty   ****\n***********************************" << endl;
    if ( _boxsize/_n_cellsAlongb != 10 ) cout << "***********************************\n****  WARNING: boxsize b != 10 * n_cell !!!  ****\n***********************************" << endl;
    if (_noEwald) cout << "~~~~~~~~~~~~~~~~~~~~~~~~~!!!!!!!!!! Warning, noEwald is activated ! !!!!!!!!!~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

    // for (int l=0;l<10;l++){
 //        Vector3d rn_vec = Vector3d::Zero();
 //        rn_vec(0) = 3 + 0.001 * l;
 //        cout << rn_vec(0) << "    \n" << lub2p(rn_vec, rn_vec.squaredNorm(), 8) <<"\n -----" << endl ;
 //    }
    iftestEwald(testRealSpcM(); testEwald(); )
    iftestLub2p(testLub2p();)

}


void CConfiguration::updateStartpos(){
    //This function is used if the particle should keep moving after each run, and not start at _resetpos again, like in the first run
    //This (hopefully) will give better averages without having to spend a lot of steps in the beginning of each run to get away from _resetpos
    for (int i = 0; i < 3; i++){
    _startpos(i) = _ppos(i) + _boxsize * _boxnumberXYZ[i];
    }
}

void CConfiguration::checkDisplacementforMM(){
    double movedsq = 0;
    for (int i=0; i<3; i++){
        movedsq += pow(_ppos(i) + _boxsize *  _boxnumberXYZ[i] - _lastcheck[i] , 2);
    }
    if ( movedsq > _cutoffMMsq ){
        calcTracerMobilityMatrix(true);
        // new start point for displacement check
        _lastcheck[0] = _ppos(0) + _boxsize *  _boxnumberXYZ[0];
        _lastcheck[1] = _ppos(1) + _boxsize *  _boxnumberXYZ[1];
        _lastcheck[2] = _ppos(2) + _boxsize *  _boxnumberXYZ[2];
    }
    else calcTracerMobilityMatrix(false); //TODO fullMM: set to true

}

Vector3d CConfiguration::midpointScheme(const Vector3d & V0dt, const Vector3d & F){
    // Implementation of midpoint scheme according to Banchio2003
    Vector3d ppos_prime;
    const double n = 500;
    ppos_prime = _ppos + V0dt / n;

    Vector3d vec_rij;
    Matrix3d lubM = Matrix3d::Zero();

// loop over tracer - particle mobility matrix elements
    for (unsigned int j = 0; j < _polySpheres.size(); j++){
        vec_rij = ppos_prime - _polySpheres[j].pos;
        lubM += lubricate(vec_rij);
    }
    //cout << "############# _mobilityMatrix #########\n" << _mobilityMatrix << endl;
    //cout << "############# lubM #########\n" << lubM << endl;


    // Add lubrication Part to tracer Resistance Matrix and invert
    Matrix3d tracerMM_prime = invert3x3(_resMNoLub + lubM);

    // Note that _f_sto has variance sqrt( _tracerMM ), as it is advised in the paper of Banchio2003
    Vector3d V_primedt = tracerMM_prime * (F * _timestep + _f_sto * _mu_sto);

    // return V_drift * dt
    return n/2. * (V_primedt - V0dt);

}


int CConfiguration::makeStep(){
    //move the particle according to the forces and record trajectory like watched by outsider
//  if (_HI){
    _Vdriftdt = Vector3d::Zero();

    _prevpos = _ppos;

    bool v_nan = true;
    bool lrgDrift = true;
    int tmp  = 0;
    while (v_nan || lrgDrift){  // This loop repeats until the the midpoint scheme does not result in a NaN anymore!
        _V0dt = _tracerMM * (_f_mob * _timestep + _f_sto * _mu_sto);

        if (_HI && !_noLub) _Vdriftdt = midpointScheme(_V0dt, _f_mob);   //TODO  Enable mid-point-scheme / backflow

        v_nan = std::isnan(_Vdriftdt(0));
        lrgDrift = (_Vdriftdt.squaredNorm() > 0.3);

        if (v_nan || lrgDrift) {
            cout << "is Nan = " << v_nan << "    -   largeDrift? = " << lrgDrift << endl;
            cout << "_Vdriftdt.squaredNorm() = " << _Vdriftdt.squaredNorm() << endl;
            calcStochasticForces();
            //cout << "drifCorrection: Vdrift^2 = " << _Vdriftdt.squaredNorm() <<  endl;
            tmp++;
            cout << tmp << endl;
        }
    }

    // Update particle position
    _ppos += (_V0dt + _Vdriftdt) ;

    //cout << (_prevpos - _ppos).norm() << endl << "--------" << endl;
    if (std::isnan(_ppos(0))){
        report("NaN found!");
        return 1;
    }
    //report("TEST");

    //if (_V0dt(0) > 0.1 ) cout << "############ " << _V0dt << endl;
    return 0;
}

void CConfiguration::report(string reason){
    cout << "~*~*~*~*~*~     ALERT ALERT !  ~*~*~*~*~*~" << endl;
    cout << "----------------------------------------------------\n-------------------  "<<reason<<"  -------------\n----------------------------------------------------\n";
    cout << "_edgeParticles " << _edgeParticles << endl;
    cout << "_ppos: " << _ppos(0) << ", " << _ppos(1) << ", " << _ppos(2) << endl;
    cout << "_prevpos: " << _prevpos(0) << ", " << _prevpos(1) << ", " << _prevpos(2) << endl;
    cout << "_tracerMM: " << endl << _tracerMM << endl << " - " << endl;
    cout << "Vdriftdt:\n" << _Vdriftdt << endl << "V0dt:\n" << _V0dt << endl << " **************** " << endl;
    cout << "_f_mob\n" << _f_mob << endl << "_f_sto\n" << _f_sto << endl;
    //cout << "-----mobility Matrix---\n" << _mobilityMatrix << endl;
    cout << "_resMNoLub:\n" << _resMNoLub << endl;
    cout << "_RMLub\n" << _RMLub << endl;
    cout << "Cholesky3x3(_RMLub)\n" << Cholesky3x3(_RMLub) << endl;
    overlapreport();
    testIfSpheresAreOnRods();
}

int CConfiguration::checkBoxCrossing(){
    //should the particle cross the confinement of the cube, let it appear on the opposite side of the box
    for (int i = 0; i < 3; i++){
        int exitmarker = 0;
        if (_ppos(i) < 0.){
            _ppos(i) += _boxsize;
            _boxnumberXYZ[i] -= 1;
            exitmarker = -1;

        }
        else if (_ppos(i) >= _boxsize){
            _ppos(i) -= _boxsize;
            _boxnumberXYZ[i] += 1;
            exitmarker = 1;
        }
        if (exitmarker!=0){
            if (_ranU) _poly.shiftPolySign(i, exitmarker);
            if (_ranRod){
                updateRodsVec(i, exitmarker);
                //cout << "[["<<i<<"," << exitmarker <<"] ";
                //prinRodPos(0); // cout print rod pos!
            }
            if (_ppos(i) > _boxsize){
                cout << "Error: _ppos is outside of allowed range 0 < _ppos < _boxsize!";
                report("Invalid _ppos range!");
                return 1;
            }
        }
    }
    ifdebug(overlapreport();)
    return 0;
}




void CConfiguration::countWallCrossing(int crossaxis, int exitmarker){

    if (_entryside[crossaxis] == exitmarker){                // case: entryside same as exitside
        _wallcrossings[0] += 1;
    }
    else if (_entryside[crossaxis] == (-1.0 * exitmarker)){  // case: entryside opposite exitside
        _wallcrossings[1] += 1;
    }
    else {                                                  // case: exiting "sideways", i.e. through one of the four other sides
        _wallcrossings[2] += 1;
    }
    for (int i = 0; i < 3; i++){                            // mark new entryside
        _entryside[i] = 0;
    }
    _entryside[crossaxis] = -1 * exitmarker;
}



void CConfiguration::calcStochasticForces(){
    // the variate generator uses m_igen (int rand number generator),
    // samples from normal distribution with variance 1 (later sqrt(2) is multiplied)
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > ran_gen(
            *m_igen, boost::normal_distribution<double>(0, 1));

    Vector3d ran_v = Vector3d::Zero();

    ran_v(0) = ran_gen();
    ran_v(1) = ran_gen();
    ran_v(2) = ran_gen();

    if (_HI){
        // return correlated random vector, which is scaled later by sqrt(2 dt)
        //_f_sto = _RMLub.llt().matrixL() * ran_v;
        _f_sto = Cholesky3x3(_RMLub)  * ran_v;
        //cout << "\n" << _tracerMM * Cholesky3x3(_RMLub) << "\n===========\n" << Cholesky3x3(_tracerMM) << endl;
    }

    else { // no HI
        _f_sto = ran_v;
    }
}


void CConfiguration::calcMobilityForces(){
    //calculate mobility forces from potential Epot - Unified version that includes 2ndOrder if k is larger than or equal 0.2 b , except if ranPot is activated.
    double r_abs = 0;
    double r_i = 0, r_k = 0;
    double utmp = 0, frtmp = 0;     //temporary "hilfsvariables"
    double Epot = 0;
    double z1, z2;
    if (_ranU){
        z1 = 1/4 * _boxsize;
        z2 = _boxsize - z1;   //z is in cylindrical coordinates. This indicates above/below which value the exp potential is modifed for random signs.
    }
    //reset mobility forces to zero
    _f_mob = Vector3d::Zero();
    const double r_c = _6root2 * _pradius;    //cutoff for Lennard-Jones calculation (at minimum)

    for (int i = 0; i < 3; i++){
        int k = i + 1;   //k always one direction "further", i.e. if i = 0 = x-direction, then k = 1 = y-direction
        if ( k == 3 ) k = 0;
        int plane = 3 - (i+k); //this is the current plane of the cylindrical coordinates
        int n = 0;     // reset counter for index of next rod in plane  n = 0, 1, 2, 3 -> only needed for ranPot
        for (int nk = _min; nk < _max; nk++){
            for (int ni = _min; ni < _max; ni++){
            r_i = _ppos(i) - ni*_boxsize;
            r_k = _ppos(k) - nk*_boxsize;
            //this is needed if we dont want the rods to cross each other to create a strong potential well
            if (plane == 0){
                r_i -= _rodDistance;
            }
            else if (plane == 1){
                r_k -= _rodDistance;
                r_i -= _rodDistance;
            }

            r_abs = sqrt(r_i * r_i + r_k * r_k); //distance to the rods

            calculateExpPotential(r_abs, utmp, frtmp);

            if (_hpi) calculateExpHPI(r_abs, utmp, frtmp);

            if (_ranU){
                utmp = utmp * _poly.get_sign(plane, n);
                frtmp = frtmp * _poly.get_sign(plane, n);
                if (_ppos(plane) > z2){
                    if (! _poly.samesign(1, plane, n)){
                        _f_mob(plane) += utmp * 4 / _boxsize;              //this takes care of the derivative of the potential modification and resulting force
                        modifyPot(utmp, frtmp, _boxsize - _ppos(plane));
                    }
                }
                else if (_ppos(plane) < z1){
                    if (! _poly.samesign(-1, plane, n)){
                        _f_mob(plane) -= utmp * 4 / _boxsize;              //this takes care of the derivative of the potential modification and resulting force
                        modifyPot(utmp, frtmp, _ppos(plane));
                    }
                }
                n++;  //index of next rod in curent plane
            }


            if (_LJPot && ( r_abs < r_c || _hpi )) calcLJPot(r_abs, utmp, frtmp);


            Epot += utmp;
            _f_mob(i) += frtmp * r_i;
            _f_mob(k) += frtmp * r_k;
            }
        }
    }
    _upot = Epot;
}

void CConfiguration::saveXYZTraj(string name, const int& move, string a_w){
    Vector3d boxCoordinates;
    boxCoordinates << _boxsize *_boxnumberXYZ[0], _boxsize *_boxnumberXYZ[1], _boxsize *_boxnumberXYZ[2];
    Vector3d rtmp;
    FILE *f = fopen(name.c_str(), a_w.c_str());
    bool writeSpheres= ( _ranRod || _2DLattice || _EwaldTest != 0);

    if (!writeSpheres){
        fprintf(f, "%d\n%s (%8.3f %8.3f %8.3f) t=%d \n", 1, "sim_name", _boxsize, _boxsize, _boxsize, move);
    }
    else fprintf(f, "%d\n%s (%8.3f %8.3f %8.3f) t=%d \n", _polySpheres.size() + 1, "sim_name", _boxsize, _boxsize, _boxsize, move);

    // tracer particle
    rtmp = _ppos;//+boxCoordinates;
    fprintf(f, "%3s%9.3f%9.3f%9.3f \n","O", rtmp(0), rtmp(1),  rtmp(2));   // relative position in box

    if (writeSpheres){
        for (unsigned int i = 0; i < _polySpheres.size(); i++) {
            rtmp = _polySpheres[i].pos;//+boxCoordinates;
            fprintf(f, "%3s%9.3f%9.3f%9.3f \n","H", rtmp(0), rtmp(1),  rtmp(2));
        }
    }
    fclose(f);
}


void CConfiguration::setRanNumberGen(double seed){
    if (seed == 0) {
        m_igen = new boost::mt19937(static_cast<unsigned int>(time(NULL)));
        cout << "random seed is time!" << endl;
    } else {
        m_igen = new boost::mt19937(static_cast<unsigned int>(seed));
        cout << "random seed is " << seed << endl;
    }
}




double CConfiguration::getPosVariance(){
    //function to return the variance, to save it as an instant value
    double var = 0;
    for (int m = 0; m < 3; m++){                    //calculate (r_m - r_0)^2|_i
        var += pow((_ppos(m) + _boxsize *_boxnumberXYZ[m] - _startpos(m)) , 2);
    }
    return var;
}

double CConfiguration::get1DPosVariance(int dim){
    //function to return the variance of just one dimension x = 0, y  = 1 or z = 2
    return pow((_ppos(dim) + _boxsize *_boxnumberXYZ[dim] - _startpos(dim)) , 2);
}




bool CConfiguration::checkFirstPassage(double mfpPos, int dim){    //TODO
    //function returns true if in the last step the momentary "way point" (_mfpIndex * xinterval) has been crossed in dim-direction
    // where dim = 0, 1, 2 -> x, y, z. Otherwise it returns false
    if ((_ppos(dim) + _boxsize*_boxnumberXYZ[dim] - _startpos(dim)) > (mfpPos) ) return true;
    else return false;
}



void CConfiguration::moveBack(){
    //moves particle back to previous position
    _ppos = _prevpos;
}
//****************************POTENTIALS**********************************************************



void CConfiguration::calculateExpPotential(const double &r, double& U, double& Fr){
    //function to calculate an exponential Potential U = U_0 * exp(-1 * r * k)
    // k is the interaction range. U_0 is the strength of the potential
    //which is attractive if direction = -1, and repulsive if direction = 1
    //The potential is weighted with kT!

    U = _potStrength * exp(-1 * r / _potRange);
    Fr = U / (_potRange * r);  //This is the force divided by the distance to the rod!
}


void CConfiguration::calculateExpHPI(const double &r, double& U, double& Fr){
    double u = _hpi_u * exp( - r / _hpi_k);
    U += u;
    Fr += u / (_hpi_k * r);
}

/*void CConfiguration::calculateDHPotential(const double r, double& U, double& Fr){
    //function to calculate an exponential Potential U = U_0 * exp(-1 * r * k)
    // k is the interaction range. U_0 is the strength of the potential
    //which is attractive if direction = -1, and repulsive if direction = 1
    //The potential is weighted with kT!

    U = _potStrength * exp(-1 * r / _potRange) / r;
    Fr = U * (1/_potRange + 1/r) / r;  //This is the force divided by the distance to the rod!
}
*/


void CConfiguration::modifyPot(double& U, double& Fr,const double dist){
    //function to modify the potential according to the distance along the polymer axis to the next neighbor,
    //in case the next neighboring polymer part is of opposite sign
    U = U * 4 * dist/_boxsize;
    Fr = Fr * 4 * dist/_boxsize;
}

//****************************STERIC HINDRANCE****************************************************//

bool CConfiguration::testOverlap(){
    //Function to check, whether the diffusing particle of size psize is overlapping with any one of the rods (edges of the box)
    //mostly borrowed from moveParticleAndWatch()


    //"PROPER" METHOD FOR EwaldTest, where the overlap is calculated for spheres in the corners, not rods.
    if ((_EwaldTest != 0 || _2DLattice)){
        Vector3d vrij;
        for (unsigned int j = 0; j < _polySpheres.size(); j++){
            vrij = minImage(_ppos - _polySpheres[j].pos);
            if (vrij.squaredNorm() <= _stericrSq + 0.00001){
                return true;
            }
        }
        return false;
    }

    if (_ranRod){
        Vector3d testpos;
        for (int axis=0;axis<3;axis++){
            testpos = _ppos;
            testpos(axis) = 0.;
            for (int l = 0; l < _rodvec[axis].size(); l++){
                if ((testpos - _rodvec[axis][l].coord).squaredNorm() < _stericrSq + 0.0001){
                    //cout << "overlaps! ppos = " << _ppos(0) << ", " << _ppos(1) << ", " << _ppos(2) << endl;
                    return true;
                }
            }
        }
        return false;
    }


    double r_i = 0, r_k = 0;
    double r_sq = 0;
    double cellwidth = _boxsize/_n_cellsAlongb;

    for (int i = 0; i < 2; i++){
        for (int k = i+1; k < 3; k++){
            for (int n_i = 0; n_i <= _n_cellsAlongb; n_i++ ){
                r_i = _ppos(i) - cellwidth * n_i;

                for (int n_k = 0; n_k <= _n_cellsAlongb; n_k++ ){
                    r_k = _ppos(k) - cellwidth * n_k;
                    //this is needed if we dont want the rods to cross each other to create a strong potential well
                    if (k == 2){
                        r_i -= _rodDistance;
                        if (i == 0) r_k -= _rodDistance;
                    }
                    r_sq = r_i * r_i + r_k * r_k; //distance to the rods
                    if (r_sq <= _stericrSq + 0.00001){
                        return true;
                    }
                }
            }
        }
    }
    return false;
}



void CConfiguration::calcLJPot(const double& r, double& U, double& Fr){
    //Function to calculate the Lennard-Jones Potential
    double  por6 = pow((_pradius / r ), 6);      //por6 stands for "p over r to the power of 6" . The 2 comes from the fact, that I need the particle radius, not the particle size
    U += 4 * ( por6*por6 - por6 + 0.25 );
    Fr +=  24 / ( r * r ) * ( 2 * por6*por6 - por6 );
}


void CConfiguration::initPolySpheres(){
    _polySpheres.clear();
    // store the edgeParticle positions, so that I can simply loop through them later
    std::vector<Vector3d> zeroPos( 3 * _edgeParticles - 2 , Vector3d::Zero() );
    Vector3d nvec;
    // store the edgeParticle positions in first cell in zeroPos
    for (int i = 1; i < _edgeParticles; i++){
        double tmp = i * _sphereoffset;
        zeroPos[i](0) = tmp;
        zeroPos[i + (_edgeParticles - 1)](1) = tmp;
        zeroPos[i + 2 * (_edgeParticles - 1)](2) = tmp;
    }
    for (int nx=0; nx < _n_cellsAlongb; nx++){
        for (int ny=0; ny < _n_cellsAlongb; ny++){
            for (int nz=0; nz < _n_cellsAlongb; nz++){
                nvec << nx, ny, nz; // Position of 0 corner of the simulation box
                for (unsigned int i = 0; i < zeroPos.size(); i++){
                    _polySpheres.push_back( CPolySphere( zeroPos[i] + nvec * _boxsize/_n_cellsAlongb ) );
                    //cout << "----" << _polySpheres[i].pos << endl;
                }
            }
        }
    }
    
    // ******** Phillips1990 *********
    // rods are arranged on a 2D square lattice. The rods point into the z-direction. The lattice is in the x-y-plane.
    // I ONLY CHANGE THE zeroPos Vector3d Array!
    if (_2DLattice){
        _polySpheres.clear();
        nvec = Vector3d::Zero();
        zeroPos.resize(_edgeParticles);
        for (int i = 0; i < _edgeParticles; i++){
            double tmp = i * _sphereoffset;
            zeroPos[i](0) = 0;
            zeroPos[i](1) = 0;
            zeroPos[i](2) = tmp;
        }
        for (int nx=0; nx < _n_cellsAlongb; nx++){
            for (int ny=0; ny < _n_cellsAlongb; ny++){
                for (int nz=0; nz < _n_cellsAlongb; nz++){
                    nvec << nx, ny, nz; // Position of 0 corner of the simulation box
                    for (unsigned int i = 0; i < zeroPos.size(); i++){
                        _polySpheres.push_back( CPolySphere( zeroPos[i] + nvec * _boxsize/_n_cellsAlongb ) );
                        //cout << "----" << _polySpheres[i].pos << endl;
                    }
                }
            }
        }
    }

    // ****** RANDOM RODS *******
    //_edgeParticles stays the same
    if (_ranRod){
        copySphereRods();
    }



    // *******   RANSPHERE  ********
    bool overlap = true;
    Vector3d vrij;
    Vector3d spos;

    int Nspheres = _n_cellsAlongb * _n_cellsAlongb * _n_cellsAlongb;

    if (_ranSpheres){ // This serves to assign a random position to the edgeparticles. Note: Only use it with _EwaldTest = 1!
        int retry = 0;
        while (overlap==true && retry < 100){
            _polySpheres.clear();
            //create sphere in corner at origin
            _polySpheres.push_back( CPolySphere( Vector3d::Zero() ) );
            ++retry;
            //cout << "retry " << retry << endl;
            for (int nx=0; nx < _n_cellsAlongb; nx++){
                for (int ny=0; ny < _n_cellsAlongb; ny++){
                    for (int nz=0; nz < _n_cellsAlongb; nz++){
                        nvec << nx, ny, nz; // Position of 0 corner of the simulation box
                        if (nvec == Vector3d::Zero()) continue;
                        else{
                            Vector3d cellcoordinate = nvec * _boxsize/_n_cellsAlongb;
                            for (int count=0; count < 200; count++){
                                overlap = false;
                                for (int k = 0; k<3; k++){
                                    //offset takes care that there is less chance of overlap for large polyrads
                                    //it increases from 0 to 2*polyrad, for the cell at the origin and the one moste far away from it, respectively
                                    double offset = 2*_polyrad*(nvec.norm() / (sqrt(3)*(_n_cellsAlongb-1)));
                                    //cout << "offset = " << offset << endl;
                                    // Calculate random position of sphere spos
                                    spos(k) = (_boxsize/_n_cellsAlongb - offset) *zerotoone()  + cellcoordinate(k);
                                    if (_trueRan) spos(k) = _boxsize * zerotoone();
                                    //cout << "spos\n" <<spos << endl;
                                }
                                for (int l = 1; l < _polySpheres.size(); l++){
                                    vrij = minImage(spos - _polySpheres[l].pos);
                                    if (vrij.norm() < 2*_polyrad + 0.000001){
                                        overlap = true;
                                        //cout << "found overlap! rij = " << vrij.norm() << endl;
                                        break;
                                    }
                                }
                                // if there is no overlap after the previous for loop to test overlap, then stop this for loop and accept the sphere pos
                                if (overlap==false) break;
                            }
                            if (overlap==true){
                                //cout << "* Overlap between polyspheres. Retry!" << endl;
                                break;//break out of count loop
                            }
                            _polySpheres.push_back( CPolySphere( spos ) );
                        }
                        if (overlap==true) break;//break out of nz loop
                    }
                    if (overlap==true) break;
                }
                if (overlap==true) break;
            }
        }
        cout << "ranSphere initiation retries: " << retry << endl;
        if (overlap==true){
            cout << "ERROR: Overlap between polyspheres could not be avoided." << endl;
            throw 2;
        }
    }

}

//**************************** HYDRODYNAMICS ****************************************************//

void CConfiguration::initConstMobilityMatrix(){
    double rxi = 0.0, ryi = 0.0, rzi = 0.0;
    double rij = 0.0, rijsq = 0.0;
    const double asq = _polyrad * _polyrad;

    // create mobility matrix - Some elements remain constant throughout the simulation. Those are stored here.
    _mobilityMatrix = MatrixXd::Identity( 3 * (_polySpheres.size() + 1) , 3 * (_polySpheres.size() + 1) );

    // Create matrix that stores previous result of conjugate gradient method, for faster conversion.
    _prevCG = Eigen::MatrixXd::Identity(3 * (_polySpheres.size() + 1) , 3 );


    // Eigen-vector for interparticle distance. Needs to initialized as zero vector for self mobilities.
    Vector3d vec_rij = Vector3d::Zero();

    // on-diagonal elements of mobility matrix: TRACER PARTICLE
    Matrix3d selfmob = realSpcSm( vec_rij, true, _pradius * _pradius ) + reciprocalSpcSm( vec_rij, _pradius * _pradius );
    double self_plus = 1. + _pradius / sqrt (M_PI) * ( - 6. * _alpha + 40. * pow(_alpha,3) * _pradius * _pradius / 3. );
    if (!_ranRod && !_noEwald) _mobilityMatrix.block<3,3>(0,0) = selfmob + Matrix3d::Identity() * self_plus; // ewaldCorr

    // now on diagonals for POLYMER SPHERES
    selfmob = realSpcSm( vec_rij, true, asq ) + reciprocalSpcSm( vec_rij, asq );
    // double self_plus = _pradius/_polyrad + _pradius / sqrt (M_PI) * ( - 6. * _alpha );  //TODO  alt
    self_plus = _pradius/_polyrad + _pradius / sqrt (M_PI) * ( - 6. * _alpha + 40. * pow(_alpha,3) * _polyrad * _polyrad / 3. );
    selfmob(0,0) += self_plus;
    selfmob(1,1) += self_plus;
    selfmob(2,2) += self_plus;


    for (unsigned int i = 1; i < _polySpheres.size() + 1 ; i++) {
        unsigned int i_count = 3 * i;
        /*_mobilityMatrix(0,0) = 1;  already taken care of due to identity matrix. The extra F_i0 part in the Ewald sum is only correction to reciprocal sum.
         * I leave out the reciprocal summation for the tracer particle self diff, because it is in only in the simulation box unit cell. */

        if (_ranRod || _noEwald) _mobilityMatrix.block<3,3>(i_count,i_count) = _pradius/_polyrad * Matrix3d::Identity();
        else _mobilityMatrix.block<3,3>(i_count,i_count) = selfmob;
    }


    // The mobility matrix elements for the interacting edgeParticles are stored here, since they do not change position
    Matrix3d muij = Matrix3d::Zero();
    for (unsigned int i = 0; i < _polySpheres.size(); i++){
        unsigned int i_count = 3 * (i + 1);       // plus 1 is necessary due to omitted tracer particle

        for (unsigned int j = i + 1; j < _polySpheres.size(); j++) {
            vec_rij = _polySpheres[i].pos - _polySpheres[j].pos;
            unsigned int  j_count = 3 * (j + 1);    // plus 1 is necessary due to omitted tracer particle

        /*
         * Calculation of RP Ewald sum
         */
            if (_ranRod) muij = RPYamakawaPS(vec_rij, asq);
            else if (_noEwald) muij = RotnePrager( minImage(vec_rij), asq );
            else muij = realSpcSm( vec_rij, false, asq ) + reciprocalSpcSm( vec_rij, asq );

            // both lower and upper triangular of symmetric matrix need be filled
            _mobilityMatrix.block<3,3>(j_count,i_count) = muij;
            _mobilityMatrix.block<3,3>(i_count,j_count) = muij;
            //cout << "----------------\n" << selfmob << endl;
        }
    }
    initMreciprocalTracer();
}


void CConfiguration::calcTracerMobilityMatrix(const bool& full){
    const double asq = 0.5*(_polyrad * _polyrad + _pradius * _pradius);

    // vector for outer product in muij
    Vector3d vec_rij;
    double rij_sq;
    Matrix3d lubM = Matrix3d::Zero();
    Matrix3d muij;

// loop over tracer - particle mobility matrix elements
    for (unsigned int j = 0; j < _polySpheres.size(); j++){
        vec_rij = _ppos - _polySpheres[j].pos;
        if (_noEwald) vec_rij = minImage(vec_rij);
        //vec_rij = minImage(vec_rij); // tmpminim
        //rij_sq = vec_rij.squaredNorm();
        // if (rij_sq <= (_stericrSq + 0.00001)){ // If there is overlap between particles, the distance is set to a very small value, according to Brady and Bossis in Phung1996
        //     // set distance to 0.00000001 + _stericr but preserve direction
        //     double corr = (0.00001 + sqrt(_stericrSq)) / sqrt(rij_sq);
        //     vec_rij *= corr;
        // }

        if (full){
            // Calculation of different particle width Rotne-Prager Ewald sum
            unsigned int j_count = 3 * (j + 1); //  plus 1 is necessary due to omitted tracer particle
            if (_ranRod || _noEwald) muij = RotnePrager( vec_rij, asq );
            // TODO Mreci Tracer
            // else muij = realSpcSm( vec_rij, false, asq ) + reciprocalSpcSm( vec_rij, asq );
            else muij = realSpcSm( vec_rij, false, asq ) + reciprocalSpcSmTracer( vec_rij );
            //TODO del
            //cout << "realSm\n" << realSpcSm( vec_rij, false, asq )<< "\nreciSm\n" << reciprocalSpcSmTracer( vec_rij ) << endl;

            _mobilityMatrix.block<3,3>(j_count,0) = muij;
            _mobilityMatrix.block<3,3>(0,j_count) = muij;
        }
        if (!_noLub ) lubM += lubricate(vec_rij);
    }
    //cout << "############# _mobilityMatrix #########\n" << _mobilityMatrix << endl;
        //cout << "############# lubM #########\n" << lubM << endl;

    // create resistance matrix - Some elements remain constant throughout the simulation. Those are stored here.
    if (full){
         _resMNoLub = ConjGradInvert(_mobilityMatrix);
         //_resMNoLub = CholInvertPart(_mobilityMatrix);
    }
    // Add lubrication Part to tracer Resistance Matrix and invert
    _RMLub = _resMNoLub + lubM;
    _tracerMM = invert3x3(_RMLub);
    //cout << _tracerMM << endl << "........................." << endl;

    // cout << "tracerMM\n" << _tracerMM << endl;
//     cout << "inv(tracerMM)\n" << CholInvertPart(_tracerMM) << endl;
//     cout << "_RMLub\n" << _RMLub << endl;
//     cout << "---XXXX---\n";
//     Matrix3d tmp = (_RMLub.llt().matrixL() * Matrix3d::Identity()).transpose();
//     cout <<  tmp * _RMLub.llt().matrixL() << endl;
//     cout << "------\n";
//     cout << _RMLub.llt().matrixL() * Matrix3d::Identity() << endl;
//     tmp = (_tracerMM.llt().matrixL() * Matrix3d::Identity()).transpose();
//     cout <<  tmp  << endl;
//     cout << "******************************" << endl;
    //cout << "$$$$$$$ lubM $$$$$$$$$\n" << lubM << endl;
    //cout << "############# _tracerMM #########\n" << _tracerMM << endl;
    //cout << _ppos(0) << " " << _ppos(1) << " "<< _ppos(2) << endl;
}



void CConfiguration::updateMobilityMatrix(){// ONLY for ranRod, since I dont have pbc.
    const double asqPoly = _polyrad * _polyrad;
    const double asqTracer = 0.5*(_polyrad * _polyrad + _pradius * _pradius);
    Matrix3d muij = Matrix3d::Zero();
    unsigned int i_count, j_count;

    // create mobility matrix - Some elements remain constant throughout the simulation. Those are stored here.
    _mobilityMatrix = MatrixXd::Identity( 3 * (_polySpheres.size() + 1) , 3 * (_polySpheres.size() + 1) );

    // resize matrix that stores previous result of conjugate gradient method, for faster conversion.
    _prevCG = Eigen::MatrixXd::Identity(3 * (_polySpheres.size() + 1) , 3 );

    // Eigen-vector for interparticle distance. Needs to initialized as zero vector for self mobilities.
    Vector3d vec_rij = Vector3d::Zero();

    //TODO maybe store i_count and j_count first.
    //const int Nspheres = _polySpheres.size()
    //std::array<int, _polySpheres.size()> icount_arr;

    //Update the polymer Spheres mobility matrix elements
    for (unsigned int i = 0; i < _polySpheres.size(); i++){
        i_count = 3 * (i + 1);       // plus 1 is necessary due to omitted tracer particle
        // PolySphere Self interaction
        _mobilityMatrix.block<3,3>(i_count,i_count) *= _pradius/_polyrad;

        // Tracer-polySphere RP
        vec_rij = _ppos - _polySpheres[i].pos;
        muij = RotnePrager( vec_rij, asqTracer );
        _mobilityMatrix.block<3,3>(i_count,0) = muij;
        _mobilityMatrix.block<3,3>(0,i_count) = muij;


        // PolySphere - PolySphere RP
        for (unsigned int j = i + 1; j < _polySpheres.size(); j++) {
            vec_rij = _polySpheres[i].pos - _polySpheres[j].pos;
            j_count = 3 * (j + 1);    // plus 1 is necessary due to omitted tracer particle

            muij = RPYamakawaPS(vec_rij, asqPoly);

            // both lower and upper triangular of symmetric matrix need be filled
            _mobilityMatrix.block<3,3>(j_count,i_count) = muij;
            _mobilityMatrix.block<3,3>(i_count,j_count) = muij;
        }
    }
}


//----------------------- Ewald sum ------------------------

Matrix3d  CConfiguration::realSpcSm( const Vector3d & rij, const bool & self, const double& asq ){  // This should be distance of particles
    const int nmax = _nmax;
    int maxIter = 2*nmax;
    const double r_cutoffsq = _r_cutoffsq;
    double rsq;
    // TODO nmax !!!
    Vector3d rn_vec(3);
    Matrix3d  Mreal = Matrix3d::Zero();
    double v1[maxIter+1] , v2[maxIter+1], v3[maxIter+1];
    for (int n = 0; n <= maxIter; n++){
        v1[n] = (rij(0) + _boxsize * (n-nmax));
        v2[n] = (rij(1) + _boxsize * (n-nmax));
        v3[n] = (rij(2) + _boxsize * (n-nmax));
    }
    for (int n1 = 0; n1 <= maxIter; n1++){

        rn_vec(0) = v1[n1];
        for (int n2 = 0; n2 <= maxIter; n2++){

            rn_vec(1) = v2[n2];
            for (int n3 = 0; n3 <= maxIter; n3++){

                // CASE n ==(0, 0, 0) and nu == eta
                if ((n1==nmax) && (n2==nmax) && (n3==nmax) && self) continue;  // (n1 == nmax) means (n1-nmax == 0) but faster !!
                else{
                    rn_vec(2) = v3[n3];
                    rsq = rn_vec.squaredNorm();
                    if ( rsq <= r_cutoffsq ){
                        if (rsq <= _stericrSq + 0.000001){ // If there is overlap between particles, the distance is set to a very small value, according to Brady and Bossis in Phung1996
                            rsq = 0.000001 + _stericrSq;
                            //cout << "corrected rsq realSpcSm" << endl;
                        }
                        Mreal += realSpcM(rsq, rn_vec, asq);

                    }
                }
            }
        }
    }
    return Mreal;
}

void CConfiguration::initMreciprocalTracer(){
    _Mreciprocal_arr.clear();
    _kvec_arr.clear();
    const double asq = 0.5*( _polyrad*_polyrad + _pradius*_pradius );
    const double k_cutoffsq = pow(_k_cutoff, 2);
    const int nkmax = _nkmax;
    const double ntok = 2 * M_PI / _boxsize;
    Vector3d kvec(3);
    for (int n1 = -nkmax; n1 <= nkmax; n1++){
        for (int n2 = -nkmax; n2 <= nkmax; n2++){
            for (int n3 = -nkmax; n3 <= nkmax; n3++){
                if (n1 == 0 && n2 == 0 && n3 == 0)  continue;
                else{
                    kvec(0) = n1 * ntok, kvec(1) = n2 * ntok, kvec(2) = n3 * ntok;
                    const double ksq = kvec.squaredNorm();
                    if ( ksq <= k_cutoffsq ){
                        _kvec_arr.push_back( kvec );
                        _Mreciprocal_arr.push_back( reciprocalSpcM(ksq, kvec, asq) );
                    }
                }
            }
        }
    }
}

Matrix3d  CConfiguration::reciprocalSpcSmTracer( const Vector3d & rij ){  // This should be distance of particles
    const double k_cutoffsq = pow(_k_cutoff, 2);
    const int nkmax = _nkmax;
    const double ntok = 2 * M_PI / _boxsize;
    Vector3d kvec(3);
    Matrix3d  Mreciprocal = Matrix3d::Zero();
    const int imax = _kvec_arr.size();
    for (int i = 0; i < imax ; i++){
        Mreciprocal += _Mreciprocal_arr[i] * cos((_kvec_arr[i].dot(rij)));
        //TODO del
        //cout << "!!!!! _Mreciprocal_arr[i]\n" << _Mreciprocal_arr[i] << "\ncos =" << cos((_kvec_arr[i].dot(rij))) <<  endl;
    }

    return Mreciprocal * _Vinv;
}

Matrix3d  CConfiguration::reciprocalSpcSm( const Vector3d & rij, const double& asq ){  // This should be distance of particles
    const double k_cutoffsq = pow(_k_cutoff, 2);
    const int nkmax = _nkmax;
    const double ntok = 2 * M_PI / _boxsize;
    Vector3d kvec(3);
    Matrix3d  Mreciprocal = Matrix3d::Zero();
    for (int n1 = -nkmax; n1 <= nkmax; n1++){
        for (int n2 = -nkmax; n2 <= nkmax; n2++){
            for (int n3 = -nkmax; n3 <= nkmax; n3++){
                if (n1 == 0 && n2 == 0 && n3 == 0)  continue;
                else{
                    kvec(0) = n1 * ntok, kvec(1) = n2 * ntok, kvec(2) = n3 * ntok;
                    const double ksq = kvec.squaredNorm();
                    if ( ksq <= k_cutoffsq ){
                        Mreciprocal += reciprocalSpcM(ksq, kvec, asq) * cos((kvec.dot(rij)));
                    }
                }
            }
        }
    }
    return Mreciprocal * _Vinv;
}

Matrix3d  CConfiguration::realSpcM(const double & rsq, const Vector3d & rij, const double& asq) {
    // Idea: To only account for tracer in one box, leave out eta = eta_tracer in sums entirely or something...
    Matrix3d  I = Matrix3d::Identity();
    const double r = sqrt(rsq);
    const double alpha = _alpha;
    const double alphasq = alpha * alpha;
    const double c1 = alphasq * rsq;
    const double c2 = 4. * alphasq * alphasq * alphasq * rsq * rsq;
    const double c3 = alphasq * alphasq * rsq;
    const double c4 = 1./rsq;
    const double c5 = asq*c4;
    const double expc = exp(-c1)/_srqtPi * alpha;
    const double erfcc = erfc(alpha * r) / r;
    

    return _pradius * ( I * ( erfcc *  ( 0.75 + 0.5 * c5 )
        + expc * (3. * c1 - 4.5 + asq * ( c2 - 20. * c3 + 14. * alphasq + c4 )))
            + (rij*rij.transpose()) * c4 * (
                erfcc * ( 0.75 - 1.5 * c5 ) +
                    expc * ( 1.5 - 3. * c1 + asq * (- c2 + 16. * c3 - 2. * alphasq - 3. * c4))));
}


Matrix3d  CConfiguration::reciprocalSpcM(const double& ksq, const Vector3d & kij, const double& asq) {
    Matrix3d  I = Matrix3d::Identity();
    const double alphasq = _alpha * _alpha;
    const double c1 = ksq / ( 4. * alphasq );
    const double ksqInv = 1./ksq;

    return _pradius * ( 1. - 0.333333 * asq * ksq ) * ( 1. + c1 + 2. * c1 * c1 ) * ( 6. * M_PI * ksqInv ) * exp( -c1 )
        * ( I - (kij*kij.transpose()) * ksqInv);
    // return _pradius * ( 1. - 0.333333 * asq * ksq ) * ( 1. + c1 + 2. * c1 * c1 ) * ( 6. * M_PI / (ksq * _V)) * exp( -c1 )
    //     * ( I - (kij*kij.transpose())/ksq);
        //alternative way from Brady1988 TODO alt
    // return _pradius  * ( 1. + c1 + 2. * c1 * c1 ) * ( 6. * M_PI / (ksq * V)) * exp( - c1 ) * ( I - (kij*kij.transpose())/ksq);
}

//------------------------------------- Rotne-Prager ---------------------------------------

Matrix3d CConfiguration::RotnePrager(const Vector3d & rij, const double & asq) {
    // source http://www.fuw.edu.pl/~piotrek/publications/different_sized.pdf.
    Matrix3d  I = Matrix3d::Identity();
    const double rsq = rij.squaredNorm();
    const double r = sqrt(rsq);
    const double rInv = 1./r;

    ifdebug(if (r<2*_polyrad){
        cout << "Too small poly distance in RP!\nasq = "<< asq << "\nr = "<< r <<  endl;
        testIfSpheresAreOnRods();
        overlapreport();};)

    const double c1 = 2. * asq * rInv * rInv;
    return _pradius * .75 * rInv  * ( ( 1. + 0.3333 * c1 ) * I + ( 1. - c1 ) * rij*rij.transpose() * rInv* rInv );
    // return _pradius * 3. / ( 4 * r ) * ( ( 1. + 2. * asq / (3. * rsq )) * I + ( 1. - 2.* asq / rsq ) * rij*rij.transpose() / rsq );
}


Matrix3d CConfiguration::RPYamakawaPS(const Vector3d & rij, const double asq) {
    // ONLY FOR SAME SIZE PARTICLES, i.e. polyspheres
    // source http://www.fuw.edu.pl/~piotrek/publications/different_sized.pdf.
    // overlap testing if activated..
    const double rSq = rij.squaredNorm();
    if (rSq < _polydiamSq) {
        const double r = sqrt(rSq);
        Matrix3d  I = Matrix3d::Identity();
        // overlap detected               //(1. - (9. / 32.) * _r / ( m_a));
        // _pradius * _Invpolyrad  -- term is needed to arrive at self-diffusion of a polySphere
        return _pradius * _Invpolyrad * ( I * (  1. - 0.28125 * r * _Invpolyrad  )  + 0.09375 * _Invpolyrad / r * rij*rij.transpose());    //c1*I + (3./(32.*m_a*_r)) * outer_prod(_rij, _rij);
    }
    else {
        // no overlap
        return RotnePrager(rij, asq);
    }
}

//------------------------------------- Lubrication Resistance Matrix ------------------------------

Matrix3d  CConfiguration::lubricate( const Vector3d & rij ){
    // Addition to lubrication Matrix for nearest edgeparticles
    Matrix3d  lubPart = Matrix3d::Zero();
    // TODO nmax !!!
    Vector3d rn_vec(3);
    double rsq;
//     for (int n1 = -1; n1 <= 1; n1++){
    if (_ranRod){
        rsq = rij.squaredNorm();
        if (rsq <= 0.00001 + _stericrSq){
            rsq = 0.00001 + _stericrSq;
        }
        lubPart = lub2p(rij, rsq); // Only this, to avoid problems with different mmax for Long-Range lubrication in lub2p (i.e. Sum3 and Sum4)
    }
    else{
        double values_1[3] = {rij(0) - _boxsize, rij(0), rij(0) + _boxsize};
        double values_2[3] = {rij(1) - _boxsize, rij(1), rij(1) + _boxsize};
        double values_3[3] = {rij(2) - _boxsize, rij(2), rij(2) + _boxsize};
        for (int n1 = 0; n1 < 3; n1++) {
            rn_vec(0) = values_1[n1];
            for (int n2 = 0; n2 < 3; n2++) {
                rn_vec(1) = values_2[n2];
                for (int n3 = 0; n3 < 3; n3++) {
                    rn_vec(2) = values_3[n3];
                    rsq = rn_vec.squaredNorm();
                    if ( rsq <= _cutofflubSq ){
                        // if (rsq < _cutofflubSq/2){

                        //         //double corr = sqrt((0.0000001 + _stericrSq)/rsq);
                        //
                        //         //cout << "######################## correction! ################" << endl;
                        //         //rn_vec *= corr;
                        //     }
                        // If there is overlap between particles, the distance is set to a very small value, according to Brady and Bossis in Phung1996
                        if (rsq <= 0.00001 + _stericrSq){
                            rsq = 0.00001 + _stericrSq;
                        }
                        lubPart += lub2p(rn_vec, rsq); // Only this, to avoid problems with different mmax for Long-Range lubrication in lub2p (i.e. Sum3 and Sum4)
                    }
                }
            }
        }
    }
    return lubPart;
}



void CConfiguration::testLub2p(){
    Vector3d rij(_polyrad+_pradius,0,0.1);
    double rsq = rij.squaredNorm();

    const unsigned int mmax = _fXm.size();
    
    const double sSq = 4.*rsq/pow(_pradius + _polyrad,2.);
    const Matrix3d rrT = rij * rij.transpose() / rsq;
    const Matrix3d I = Matrix3d::Identity();
    cout << "\ns " << sqrt(sSq) << endl;
    const double sinvSq = 1./sSq;
    const double c1 = 4.*sinvSq;//= pow(2/s,2)
    double c1pows[mmax];
    c1pows[0] = 1;
    for (int m = 1; m < mmax; m++){
        c1pows[m] = c1 * c1pows[m-1];
    }
    const double c2 = (1-c1);
    double c3 = 0;//this needs to be zero, if s>3 and log not computed
    double SumX = 0;
    double SumY = 0;
    //TODO BIG QUESTION; Since the sum is infinite it should actually completely vanish in combination with the log terms. As julian about this!
    if (sSq<9.){
        c3 = log(c2);
        for (int m = 1; m < mmax; m++){
            SumX += ( - _minv[m] * _gX[1] + _m2inv[m] * _gX[2] ) * c1pows[m];
            //SumX += (_fXm[m] - _gX[0] - 1./m * _gX[1] + 1./(m*(m-1)) * _gX[2] ) * c1pows[m];

            SumY += ( - _minv[m] * _gY[1] +  _m2inv[m] * _gY[2] ) * c1pows[m];
        }
    }
    for (int m = 1; m < mmax; m++){
        //TODO precalc 1/m and 1/(m*(m-1)) and put into table, it's always the same at each iteration here and below in SumY
        SumX += (_fXm[m] - _gX[0] ) * c1pows[m];

        SumY += (_fYm[m] ) * c1pows[m];
    }

    for (int m = 1; m < mmax; m++){
        //TODO precalc 1/m and 1/(m*(m-1)) and put into table, it's always the same at each iteration here and below in SumY
        cout << ", " << _fYm[m];
    }
    cout << endl;
    for (int m = 1; m < mmax; m++){
        //TODO precalc 1/m and 1/(m*(m-1)) and put into table, it's always the same at each iteration here and below in SumY
        cout << ", " << _fXm[m];
    }

    double X = _gX[0]/(c2)  -  _gX[1]*c3  -  _gX[2]*(c2)*c3  +  1.  -  _gX[0] + SumX;

    double Y = -_gY[1]*c3 - _gY[2]*(c2)*c3  +  1.  + SumY;

    cout << "X = "<< X << "\nY = " << Y<< endl;
    abort();

}

Matrix3d CConfiguration::lub2p(const Vector3d &rij, const double &rsq){
    // This function returns the 3x3 SELF-lubrication part of the resistance matrix of the tracer particle, i.e. A_{11} in Jeffrey1984
    const unsigned int mmax = _fXm.size();

    const double sSq = 4*rsq/pow(_pradius + _polyrad,2.);
    const Matrix3d rrT = rij * rij.transpose() / rsq;
    const Matrix3d I = Matrix3d::Identity();
//     cout << "\ns " << sqrt(sSq) << endl;
    const double sinvSq = 1./sSq;
    const double sinv = sqrt(sinvSq);
    const double c1 = 4.*sinvSq;//= pow(2/s,2)
    double c1pows[mmax];
    c1pows[0] = 1;
    for (int m = 1; m < mmax; m++){
        c1pows[m] = c1 * c1pows[m-1];
    }
    const double c2 = (1-c1);
    double c3 = 0;//this needs to be zero, if s>3 and log not computed
    double SumX = 0;
    double SumY = 0;
    /*[SOLVED] BIG QUESTION; Since the sum is infinite it should actually completely vanish in combination with the log terms.
    But! The log terms in combination with a limited m is what makes this form correct and converge for small s, as stated in Jeffrey1984
    */
    if (sSq<10.){
        c3 = log(c2);
        for (int m = 1; m < mmax; m++){
            SumX += ( - _minv[m] * _gX[1] + _m2inv[m] * _gX[2] ) * c1pows[m];
            //SumX += (_fXm[m] - _gX[0] - 1./m * _gX[1] + 1./(m*(m-1)) * _gX[2] ) * c1pows[m];

            SumY += ( - _minv[m] * _gY[1] +  _m2inv[m] * _gY[2] ) * c1pows[m];
        }
    }
    for (int m = 1; m < mmax; m++){
        //TODO precalc 1/m and 1/(m*(m-1)) and put into table, it's always the same at each iteration here and below in SumY
        SumX += (_fXm[m] - _gX[0] ) * c1pows[m];

        SumY += (_fYm[m] ) * c1pows[m];
    }

    double X = _gX[0]/(c2)  -  _gX[1]*c3  -  _gX[2]*(c2)*c3  +  1.  -  _gX[0] + SumX;

    double Y = -_gY[1]*c3 - _gY[2]*(c2)*c3  +  1.  + SumY;

    // End Long-Range part
    const Matrix3d lubR = (X - Y) * rrT + Y * I;
    //
    // ------------------ OLD ----------------------


    //Here, i am subtracting the 2paricle RP part
    Matrix3d RPinv;
    // invRP fit function. Calculate fit Polymer
    if (_fitRPinv){
        double c5 = sinv;
        double pI = _fitpIs[0]; double prr = _fitprrs[0];
        int end = _fitpIs.size();
        for (int m = 1; m < end; m++){
            pI += _fitpIs[m] * c5;
            prr += _fitprrs[m] * c5;
            c5 *= sinv;
        }
        RPinv = I * pI + rrT * prr;
    }
    else {  // Matrix inversion
        _RP2p.block<3,3>(0,3) = RotnePrager(rij, (_polyrad * _polyrad + _pradius * _pradius)/2 );
        _RP2p.block<3,3>(3,0) = _RP2p.block<3,3>(0,3);
        RPinv = CholInvertPart( _RP2p );
    }

    return lubR - RPinv;
}


Matrix3d CConfiguration::ConjGradInvert(const MatrixXd &A){
    MatrixXd I = MatrixXd::Identity(A.rows(),3);
    ConjugateGradient<MatrixXd, Lower|Upper > cg;
    cg.setTolerance( 0.00001 );
    cg.compute(A);
    _prevCG = cg.solveWithGuess(I,_prevCG);
    //_prevCG = x;//temp. make this _prevCG = cg.solveWithGuess(I,_prevCG);

    return _prevCG.block<3,3>(0,0);
}


Matrix3d CConfiguration::CholInvertPart (const MatrixXd &A) {
    MatrixXd I = MatrixXd::Identity(A.rows(),3);

    // make sure partInv is 3x3 matrix and A is NxN with N larger 2
    assert( A.rows() > 2 );
    assert( A.rows() == A.cols() );

    // perform inversion by cholesky decompoposition and return upper left 3x3 block
    // The following expression means: solve the system of linear equations A * x = I with the LLT method. x is then the inverse of A!
    // Check http://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html for comparison!
    return A.llt().solve(I).block<3,3>(0,0);
}

Matrix3d CConfiguration::Cholesky3x3(const Matrix3d & mat){
    // source http://rosettacode.org/wiki/Cholesky_decomposition
    Matrix3d L = Matrix3d::Zero();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < (i+1); j++) {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += L(i, k) * L(j, k);
            L(i, j) = (i == j) ?
                           sqrt(mat(i, i) - s) :
                           (1.0 / L(j, j) * (mat(i, j) - s));
        }

    return L;
}

Matrix3d CConfiguration::invert3x3 (const Matrix3d& A) {
    //analytical 3x3 matrix inversion - source wikipedia / stackoverflow
    Matrix3d result;

    double determinant =    +A(0,0)*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
                            -A(0,1)*(A(1,0)*A(2,2)-A(1,2)*A(2,0))
                            +A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
    double invdet = 1/determinant;
    result(0,0) =  (A(1,1)*A(2,2)-A(2,1)*A(1,2))*invdet;
    result(0,1) = -(A(0,1)*A(2,2)-A(0,2)*A(2,1))*invdet;
    result(0,2) =  (A(0,1)*A(1,2)-A(0,2)*A(1,1))*invdet;
    result(1,0) = -(A(1,0)*A(2,2)-A(1,2)*A(2,0))*invdet;
    result(1,1) =  (A(0,0)*A(2,2)-A(0,2)*A(2,0))*invdet;
    result(1,2) = -(A(0,0)*A(1,2)-A(1,0)*A(0,2))*invdet;
    result(2,0) =  (A(1,0)*A(2,1)-A(2,0)*A(1,1))*invdet;
    result(2,1) = -(A(0,0)*A(2,1)-A(2,0)*A(0,1))*invdet;
    result(2,2) =  (A(0,0)*A(1,1)-A(1,0)*A(0,1))*invdet;

    return result;
}




//****************************POS HISTOGRAM****************************************************//

void CConfiguration::initPosHisto(){
    for (int i = 0; i < 100; i++){
        for (int j = 0; j < 100; j++){
            for (int k = 0; k < 100; k++){
                _posHistoM[i][j][k] = 0;  //initialize all the 100*100*100 matrix elements to zero!
            }
        }
    }
}

void CConfiguration::addHistoValue(){
    //adds a value to the position histogram
    int x = (int)(_ppos(0) / _boxsize * 100);     //CAREFUL: THIS CAN'T BE DONE AT A POINT WHERE X MIGHT BE ZERO!!!
    int y = (int)(_ppos(1) / _boxsize * 100);
    int z = (int)(_ppos(2) / _boxsize * 100);
    // if ((x < 0) || (y < 0) || (z < 0) || (x > 99) || (y > 99) || (z > 99)){
    //     cout << "The Position Histogram function 'conf.addHisto()' is in a bad place, since there is a negative position _ppos()" << endl;
    //     cout << "The current position is: " << _ppos(0) << " " << _ppos(1) << " " << _ppos(2) << endl;
    // }
    _posHistoM[x][y][z] += 1;
}

void CConfiguration::printHistoMatrix(string folder){
    //function to print the positionHistogram to a file called posHistoMatrix.txt
    //The elements of the matrix are M[x][y][z]. First the z is counted from 1 to 100 in one row, then the y follows, then, after the 'X', the 100 x elements.

    ofstream matrixfile;
    matrixfile.open((folder + "/InstantValues/posHistoMatrix.txt").c_str());
    int maxval = 0;

    for (int i = 0; i < 100; i++){
        for (int j = 0; j < 100; j++){
            for (int k = 0; k < 100; k++){
                matrixfile << _posHistoM[i][j][k] << " ";
                if (maxval < _posHistoM[i][j][k] ) maxval = _posHistoM[i][j][k];
            }
            matrixfile << endl;
        }
        matrixfile << "X" << endl;
    }
    matrixfile << "MaxVal " << maxval;   //THis does not affect the grid files, since data is only copied to them when a "X" line comes!


    matrixfile.close();
}


std::vector<double> CConfiguration::getppos(){ // returns pointer to current particle position array
    std::vector<double> pos (3);
    for (int i = 0; i < 3; i++){
        pos[i] = _ppos(i) + _boxsize * _boxnumberXYZ[i];
    }
    return pos;
}


//****************************OLD CODE****************************************************



/*
double CConfiguration::getDisplacement(){
   FAULTY, _f_mob is now a ublas::vector !!!!  double d = 0;
    for (int i = 0; i < 3; i++){
        d += pow((_timestep * _f_mob(i) + _mu_sto * _f_sto[i]), 2);
    }
    return sqrt(d);
} */

int CConfiguration::resetposition(){
    //Reset the position to random (allowed) position in box.
    bool overlap = true;
    for (int s =0; s<2000; s++){
        for (int i = 0; i < 3; i++){
            double ranPos = zerotoone() * _boxsize;
            _startpos(i) = ranPos;
            _ppos(i) = ranPos;
            _boxnumberXYZ[i] = 0;
            _prevpos(i) = ranPos;
            _lastcheck[i] = ranPos;
        }
        if (!testOverlap()){
                overlap = false;
                break;
            }
    }
        if (overlap){
            cout << "ERROR !!!!!!!!!!!!!!!!\nCould not find random start position without overlapt between the tracer and monomers.";
            return 1;
        }
    return 0;
}



//----------------------- Rotne Prager ------------------------


//
// matrix<double> CConfiguration::RotnePrager(const double & _r, const double & _rsq,
//                                            const ublas::vector<double> & _rij) {
//     matrix<double> I = identity_matrix<double>(3);
//
//
//     double c1 = 0.75 / _r;     //3. * m_a / (4. * _r);
//     double c3 = 2. / _rsq;        //2. * m_a * m_a / _rsq;
//     double c2 = c3 / 3.;
//     return c1*( (1.+c2)*I + prod( (1.-c3)*I , matrix<double>(outer_prod(_rij,_rij)) / _rsq));
// }


/*
void CConfiguration::calcMobilityForces(){
    //calculate mobility forces from potential Epot
    double r_abs = 0;
    double r_i = 0, r_k = 0;
    double L1[] = {0, _boxsize, 0, _boxsize};  //these two arrays are only needed for the iteration.
    double L2[] = {0, 0, _boxsize, _boxsize};
    double utmp = 0, frtmp = 0;    //temporary "hilfsvariables"
    double Epot = 0;
    double z1, z2;
    if (_ranU){
        z1 = 1/4 * _boxsize;
        z2 = _boxsize - z1;   //z is in cylindrical coordinates. This indicates above/below which value the exp potential is modifed for random signs.
    }
    //reset mobility forces to zero
    for (int l = 0; l < 3; l++) {
        _f_mob(l) = 0;
    }
    const double r_c = _6root2 * _pradius;    //cutoff for Lennard-Jones calculation (at minimum)

    for (int i = 0; i < 3; i++){
            int k = i + 1;   //k always one direction "further", i.e. if i = 0 = x-direction, then k = 1 = y-direction
            if ( k == 3 ) k = 0;
            int plane = 3 - (i+k); //this is the current plane of the cylindrical coordinates
            for (int n = 0; n < 4; n++){
                r_i = _ppos(i) - L1[n];
                r_k = _ppos(k) - L2[n];
                //this is needed if we dont want the rods to cross each other to create a strong potential well
                if (plane == 0){
                    r_i -= _rodDistance;
                }
                else if (plane == 1){
                    r_k -= _rodDistance;
                    r_i -= _rodDistance;
                }
                r_abs = sqrt(r_i * r_i + r_k * r_k); //distance to the rods


                if (_potMod) calculateExpPotentialMOD(r_abs, utmp, frtmp, plane);
                else calculateExpPotential(r_abs, utmp, frtmp);


                if (_ranU){
                    utmp = utmp * _poly.get_sign(plane, n);
                    frtmp = frtmp * _poly.get_sign(plane, n);
                    if (_ppos(plane) > z2){
                        if (! _poly.samesign(1, plane, n)){
                            modifyPot(utmp, frtmp, _boxsize - _ppos(plane));
                            _f_mob(plane) += utmp * 4 / _boxsize;              //this takes care of the derivative of the potential modification and resulting force
                        }
                    }
                    else if (_ppos(plane) < z1){
                        if (! _poly.samesign(-1, plane, n)){
                            modifyPot(utmp, frtmp, _ppos(plane));
                            _f_mob(plane) -= utmp * 4 / _boxsize;              //this takes care of the derivative of the potential modification and resulting force
                        }
                    }
                }


                if (_LJPot && ( r_abs < r_c )) calcLJPot(r_abs, utmp, frtmp);


                Epot += utmp;
                _f_mob(i) += frtmp * r_i;
                _f_mob(k) += frtmp * r_k;
            }

    }
    _upot = Epot;
}


void CConfiguration::calcMobilityForces(){
    //calculate mobility forces from potential Epot
    double r_abs = 0;
    double r_i = 0, r_k = 0;
    double utmp = 0, frtmp = 0;     //temporary "hilfsvariables"
    double Epot = 0;
    //reset mobility forces to zero
    for (int l = 0; l < 3; l++) {
        _f_mob(l) = 0;
    }
    const double r_c = _6root2 * _pradius;    //cutoff for Lennard-Jones calculation (at minimum)

    for (int i = 0; i < 3; i++){
            int k = i + 1;   //k always one direction "further", i.e. if i = 0 = x-direction, then k = 1 = y-direction
            if ( k == 3 ) k = 0;
            int plane = 3 - (i+k); //this is the current plane of the cylindrical coordinates
            for (int ni = -1; ni < 3; ni++){
                for (int nk = -1; nk < 3; nk++){
                r_i = _ppos(i) - ni*_boxsize;
                r_k = _ppos(k) - nk*_boxsize;

                r_abs = sqrt(r_i * r_i + r_k * r_k); //distance to the rods


                if (_potMod) calculateExpPotentialMOD(r_abs, utmp, frtmp, plane);
                else calculateExpPotential(r_abs, utmp, frtmp);


                if (_LJPot && ( r_abs < r_c )) calcLJPot(r_abs, utmp, frtmp);


                Epot += utmp;
                _f_mob(i) += frtmp * r_i;
                _f_mob(k) += frtmp * r_k;
                }
            }

    }
    _upot = Epot;
}


template<class T>
bool CConfiguration::CholInvertPart (const matrix<T>& input, matrix<T>& partInv) {
    using namespace boost::numeric::ublas;
    if (partInv.size1() != 3) { cout << "Error Bad partInv size" << endl; exit(EXIT_FAILURE); }
    // create a working copy of the input
    matrix<T> A(input);

    triangular_matrix<double, lower> L(A.size1(), A.size1());
    L.clear();
    // perform cholesky decomposition
    size_t res = cholesky_decompose(A, L);

    // create 3x3 identity matrix of "partInv"
 //   partInv.assign(identity_matrix<T>(3));

    ublas::vector<double> b (A.size1());

    for (int i = 0; i < 3; i++){
        std::fill(b.begin(), b.end(), 0.0);  //make b zero vector
        b(i)=1;                              //make b unit vector in direction i
        inplace_solve(L, b, lower_tag() );
        inplace_solve(trans(L), b, upper_tag() );
        partInv(i, 0) = b(0);
        partInv(i, 1) = b(1);
        partInv(i, 2) = b(2);
    }
    return true;
}


*/
