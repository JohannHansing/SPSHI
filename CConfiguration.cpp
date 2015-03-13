#include "headers/CConfiguration.h"
#include <boost/timer.hpp>


using namespace Eigen;

using namespace std;

const double _6root2 = 1.122462;


CConfiguration::CConfiguration(){
}

CConfiguration::CConfiguration(
        double timestep,  double potRange,  double potStrength,  double boxsize, double rodDistance, const bool ewaldCorr,
        double psize, const bool noLub, const bool steric, const bool ranU, bool hpi, double hpi_u, double hpi_k,
		double polymersize)
		{
    _potRange = potRange;
    _potStrength = potStrength;
    _pradius = psize/2;   //_pradius is now the actual radius of the particle. hence, I need to change the definition of the LJ potential to include (_pradius + _polyrad)   -  or maybe leave LJ pot out
    _polyrad = polymersize / 2;   //This is needed for testOverlap for steric and HI stuff !!
	_boxsize = boxsize;
    _resetpos = _boxsize/2;
    _timestep = timestep;
    _rodDistance = rodDistance;
    _ewaldCorr = ewaldCorr;
    _noLub = noLub;
    _LJPot = (steric == false) && (psize != 0);
    _ranU = ranU;
    _poly = CPolymers();
    _hpi = hpi;
    _upot = 0;
    _f_mob = Vector3d::Zero();
    _f_sto = Vector3d::Zero();
    _mu_sto = sqrt( 2 * _timestep );                 //timestep for stochastic force
	_hpi = hpi; 
	_hpi_u = hpi_u;
	_hpi_k = hpi_k;
    for (int i = 0; i < 3; i++){
        _ppos(i) = _resetpos;
        _startpos(i) = _resetpos;
        _entryside[i] = 0;
        _wallcrossings[i] = 0;
        _boxnumberXYZ[i] = 0;
        _prevpos(i) = _resetpos;
        _lastcheck[i] = _startpos(i);
    }
    
    
    bool posHisto = false;
    if (posHisto) initPosHisto();
	

    // This is for inclusion of 2nd Order rods if k is 0.2b or larger
    _min = -1, _max = 3;
    if (_ranU || _hpi || (_potRange < 2)){
        _min = 0;
        _max = 2;
    }

    // seed = 0:  use time, else any integer
    // init random number generator
    setRanNumberGen(0);
    
    //Ewald sum stuff
    _nmax = 2; // This corresponds to a cutoff of r_cutoff = 1 * _boxsize
	_alpha = 1 * sqrt(M_PI) / _boxsize; // This value for alpha corresponds to the suggestion in Beenakker1986
	_k_cutoff = 2. * _alpha * _alpha * _nmax * _boxsize;   /* This corresponds to suggestion by Jain2012 ( equation 15 and 16 ). */
	_nkmax = (int) (_k_cutoff * _boxsize / (2. * M_PI) + 0.5);   /* The last bit (+0.5) may be needed to ensure that the next higher integer 
		                                                   * k value is taken for kmax */
	_r_cutoffsq = pow(_nmax * _boxsize, 2);   // Like Jain2012, r_cutoff is chosen such that exp(-(r_cutoff * alpha)^2) is small
	
	// lubrication stuff
	const double lam = polymersize/psize;
	const double c1 = pow(1+lam, -3);
	_g[0] = 2 * pow(lam, 2) * c1;
	_g[1] = lam/5 * ( 1 + 7*lam + lam*lam ) * c1;
	_g[2] = 1/42 * ( 1 + lam*(18 - lam*(29 + lam*(18 + lam)))) * c1;
    _cutofflubSq = pow(7,2)*(pow(_polyrad,2) + pow(_pradius,2));
    _stericrSq = pow(_pradius + _polyrad, 2);
	
	
	// init HI vectors matrices, etc
    _V = pow( _boxsize, 3);
    _cutoffMMsq = pow(0.05*_boxsize,2);
    _n_cellsAlongb = 5;
    if (polymersize != 0) _HI = true;
	if (_HI) {
		_edgeParticles = (int) ( ( _boxsize/_n_cellsAlongb )/polymersize + 0.001);
		_LJPot = false;
        // THIS NEEDS TO COME LAST !!!!!!!
		initConstMobilityMatrix();
		calcTracerMobilityMatrix(true);
	}


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
	else calcTracerMobilityMatrix(false);

}

Vector3d CConfiguration::midpointScheme(Vector3d V0dt, Vector3d F){
    // Implementation of midpoint scheme according to Banchio2003
    Vector3d ppos_prime;
    int n = 200;
    ppos_prime = _ppos + V0dt / n;

	Vector3d vec_rij;
	Matrix3d lubM = Matrix3d::Zero(); 
	
// loop over tracer - particle mobility matrix elements
	for (unsigned int j = 0; j < 3 * _edgeParticles - 2; j++){
		vec_rij = ppos_prime - _polySpheres[j].getPosition();
		lubM += lubricate(vec_rij);
	}
	//cout << "############# _mobilityMatrix #########\n" << _mobilityMatrix << endl;
    //cout << "############# lubM #########\n" << lubM << endl;


	// Add lubrication Part to tracer Resistance Matrix and invert	
    Matrix3d tracerMM_prime = invert3x3(_resMNoLub + lubM);
    
    // Note that _f_sto has variance sqrt( _tracerMM ), as it is advised in the paper of Banchio2003
	Vector3d V_primedt = tracerMM_prime * (F * _timestep + _f_sto * _mu_sto);
    
    // return V_drift * dt
    return n/2 * (V_primedt - V0dt);
    
}


void CConfiguration::makeStep(){
    //move the particle according to the forces and record trajectory like watched by outsider
//	if (_HI){
    Vector3d Vdriftdt = Vector3d::Zero();
    Vector3d V0dt = Vector3d::Zero();

    _prevpos = _ppos;

    bool v_nan = true;
    while (v_nan){  // This loop repeats until the the midpoint scheme does not result in a NaN anymore!
        V0dt = _tracerMM * (_f_mob * _timestep + _f_sto * _mu_sto);

        if (_HI && !_noLub) Vdriftdt = midpointScheme(V0dt, _f_mob);   //TODO  Enable mid-point-scheme / backflow
        v_nan = std::isnan(Vdriftdt(0));
        if (v_nan) {
            calcStochasticForces();
            //cout << "NaN found" << endl;
        }
    }

	// Update particle position
	_ppos += V0dt + Vdriftdt;

/*	}
	else{
	    for (int i = 0; i < 3; i++){
	        _prevpos(i) = _ppos(i);
	        _ppos(i) += _timestep * _f_mob(i) + _mu_sto * _f_sto[i];
	    }
	}
  */  
}

void CConfiguration::checkBoxCrossing(){
    //should the particle cross the confinement of the cube, let it appear on the opposite side of the box
    for (int i = 0; i < 3; i++){
        if (_ppos(i) < 0){
            _ppos(i) += _boxsize;
            _boxnumberXYZ[i] -= 1;
            countWallCrossing(i, -1);
            if (_ranU) _poly.shiftPolySign(i, -1);

        }
        else if (_ppos(i) > _boxsize){
            _ppos(i) -= _boxsize;
            _boxnumberXYZ[i] += 1;
            countWallCrossing(i, 1);
            if (_ranU) _poly.shiftPolySign(i, 1);
        }
    }
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
	    _f_sto = _RMLub.llt().matrixL() * ran_v;
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
	
	// calc HI mobility matrix here, since it needs to be defined for random force normalisation
	//if (_HI){ calcTracerMobilityMatrix(); }
}


void CConfiguration::saveXYZTraj(string name, const int& move, string a_w){
    FILE *f = fopen(name.c_str(), a_w.c_str());

    fprintf(f, "%d\n%s (%8.3f %8.3f %8.3f) t=%d \n", 1, "sim_name", _boxsize, _boxsize, _boxsize, move);

    // polymer particles
    fprintf(f, "%3s%9.3f%9.3f%9.3f \n","H", _ppos(0)+ _boxsize *_boxnumberXYZ[0], _ppos(1)+ _boxsize *_boxnumberXYZ[1], _ppos(2)+ _boxsize *_boxnumberXYZ[2]);   // absolute position
	// fprintf(f, "%3s%9.3f%9.3f%9.3f \n","H", _ppos(0), _ppos(1), _ppos(2));  // relative position in box
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



void CConfiguration::calculateExpPotential(const double r, double& U, double& Fr){
    //function to calculate an exponential Potential U = U_0 * exp(-1 * r * k)
    // k is the interaction range. U_0 is the strength of the potential
    //which is attractive if direction = -1, and repulsive if direction = 1
    //The potential is weighted with kT!

    U = _potStrength * exp(-1 * r / _potRange);
    Fr = U / (_potRange * r);  //This is the force divided by the distance to the rod!
}


void CConfiguration::calculateExpHPI(const double r, double& U, double& Fr){
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


void CConfiguration::modifyPot(double& U, double& Fr, double dist){
    //function to modify the potential according to the distance along the polymer axis to the next neighbor,
    //in case the next neighboring polymer part is of opposite sign
    U = U * 4 * dist/_boxsize;
    Fr = Fr * 4 * dist/_boxsize;
}

//****************************STERIC HINDRANCE****************************************************//

bool CConfiguration::testOverlap(){
    //Function to check, whether the diffusing particle of size psize is overlapping with any one of the rods (edges of the box)
    //mostly borrowed from moveParticleAndWatch()
    bool overlaps = false;

    // for (unsigned int i=0; i < _polySpheres.size(); i++ ){
//         // First update tracer position in Polysphere class instances to
//         _polySpheres[i].updateTracerVec( _ppos );
//         if ( _polySpheres[i].getTracerdistsq() <= _stericrSq ){
//             overlaps=true;
//             continue;
//         }
//     }
//     return overlaps;

    // KEEP THIS TO MAYBE IMPLEMENT MORE EFFICIENT testOVERLAP
    // --> store cell of tracerparticle and adjust to L1 and L2

    double r_i = 0, r_k = 0;
    double r_sq = 0;
    double L1[] = {0, _boxsize/_n_cellsAlongb, 0, _boxsize/_n_cellsAlongb};  //these two arrays are only needed for the iteration.
    double L2[] = {0, 0, _boxsize/_n_cellsAlongb, _boxsize/_n_cellsAlongb};
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
                    if (r_sq <= _stericrSq){
                        overlaps = true;
                        return overlaps;
                    }
                }
            }
            if (overlaps == true) continue;
        }
        if (overlaps == true) continue;
    }
    return overlaps;
}


void CConfiguration::calcLJPot(const double r, double& U, double& Fr){
    //Function to calculate the Lennard-Jones Potential
    double  por6 = pow((_pradius / r ), 6);      //por6 stands for "p over r to the power of 6" . The 2 comes from the fact, that I need the particle radius, not the particle size
    U += 4 * ( por6*por6 - por6 + 0.25 );
    Fr +=  24 / ( r * r ) * ( 2 * por6*por6 - por6 );
}

//**************************** HYDRODYNAMICS ****************************************************//

void CConfiguration::initConstMobilityMatrix(){
	double rxi = 0.0, ryi = 0.0, rzi = 0.0;
	double rij = 0.0, rijsq = 0.0;
	const double asq = _polyrad * _polyrad;
	
	
	// store the edgeParticle positions, so that I can simply loop through them later
    std::vector<Vector3d> zeroPos( 3 * _edgeParticles - 2 , Vector3d::Zero() );
    Vector3d nvec;
    double offset = (_boxsize/_n_cellsAlongb - 2 * _polyrad * _edgeParticles)/_edgeParticles;
	// store the edgeParticle positions in first cell in zeroPos
	for (int i = 1; i < _edgeParticles; i++){
		double tmp = i * (_boxsize/_n_cellsAlongb) / _edgeParticles;
		zeroPos[i](0) = tmp;
		zeroPos[i + (_edgeParticles - 1)](1) = tmp;  
		zeroPos[i + 2 * (_edgeParticles - 1)](2) = tmp;
	}
	for (int nx=0; nx < _n_cellsAlongb; nx++){
	    for (int ny=0; ny < _n_cellsAlongb; ny++){
            for (int nz=0; nz < _n_cellsAlongb; nz++){
                nvec << nx, ny, nz; // Position of 0 corner of the simulation box
                for (unsigned int i = 0; i < zeroPos.size(); i++){
            		double tmp = i * 2 * _polyrad;
                    _polySpheres.push_back( CPolySphere( zeroPos[i] + nvec * _boxsize/_n_cellsAlongb, _startpos ) );
            	}
            }
        }
	}
    

	// create mobility matrix - Some elements remain constant throughout the simulation. Those are stored here.
	_mobilityMatrix = MatrixXd::Identity( 3 * (_polySpheres.size() + 1) , 3 * (_polySpheres.size() + 1) );
     

    // Eigen-vector for interparticle distance. Needs to initialized as zero vector for self mobilities.
	Vector3d vec_rij = Vector3d::Zero();
	
	// on-diagonal elements of mobility matrix: Tracer first
    Matrix3d selfmob = realSpcSm( vec_rij, true, _pradius * _pradius ) + reciprocalSpcSm( vec_rij, _pradius * _pradius );
    double self_plus = 1. + _pradius / sqrt (M_PI) * ( - 6. * _alpha + 40. * pow(_alpha,3) * _pradius * _pradius / 3. );

	if (_ewaldCorr) _mobilityMatrix.block<3,3>(0,0) = selfmob + Matrix3d::Identity() * self_plus; // ewaldCorr
	
	// now on diagonals for edgeparticles
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

		_mobilityMatrix.block<3,3>(i_count,i_count) = selfmob;
	}


	// The mobility matrix elements for the interacting edgeParticles are stored here, since they do not change position
	for (unsigned int i = 0; i < _polySpheres.size(); i++){
		unsigned int i_count = 3 * (i + 1);       // plus 1 is necessary due to omitted tracer particle
		
		for (unsigned int j = i + 1; j < _polySpheres.size(); j++) {
			vec_rij = _polySpheres[i].getPosition() - _polySpheres[j].getPosition();
			unsigned int  j_count = 3 * (j + 1);    // plus 1 is necessary due to omitted tracer particle
		
		/*
		 * Calculation of RP Ewald sum
		 */
	        Matrix3d muij = Matrix3d::Zero();
			muij = realSpcSm( vec_rij, false, asq ) + reciprocalSpcSm( vec_rij, asq );
	//		cout << vec_rij << endl;
	//		cout << muij << endl; // TODO del
	//		cout << "real: " << realSpcSm( vec_rij, false ) << " ---- reciprocal: " << reciprocalSpcSm( vec_rij, false ) << endl;
			
			// both lower and upper triangular of symmetric matrix need be filled
			_mobilityMatrix.block<3,3>(j_count,i_count) = muij;
			_mobilityMatrix.block<3,3>(i_count,j_count) = muij;
            //cout << "----------------\n" << selfmob << endl;
		}
	}	
}



void CConfiguration::calcTracerMobilityMatrix(bool full){
	const double asq = (_polyrad * _polyrad + _pradius * _pradius)/2;
		
	// vector for outer product in muij
	Vector3d vec_rij(3);
	Matrix3d lubM = Matrix3d::Zero(); 
	
// loop over tracer - particle mobility matrix elements
	for (unsigned int j = 0; j < 3 * _edgeParticles - 2; j++){
        vec_rij = _ppos - _polySpheres[j].getPosition();
		
		if (full){
		    // Calculation of different particle width Rotne-Prager Ewald sum 
			unsigned int j_count = 3 * (j + 1); //  plus 1 is necessary due to omitted tracer particle
			Matrix3d muij  = realSpcSm( vec_rij, false, asq ) + reciprocalSpcSm( vec_rij, asq );
	
			// only lower triangular of symmetric matrix is filled and used
			_mobilityMatrix.block<3,3>(j_count,0) = muij;
			_mobilityMatrix.block<3,3>(0,j_count) = muij;
			//noalias(subrange(_mobilityMatrix, j_count, j_count + 3, 0, 3)) = muij;
		}
		if (!_noLub) lubM += lubricate(vec_rij);
	}
	//cout << "############# _mobilityMatrix #########\n" << _mobilityMatrix << endl;
        //cout << "############# lubM #########\n" << lubM << endl;
		
	// create resistance matrix - Some elements remain constant throughout the simulation. Those are stored here.
	if (full){ 
		 _resMNoLub = CholInvertPart(_mobilityMatrix); 
         //cout << "$$$$$$$ _resMNoLub $$$$$$$$$\n" << _resMNoLub << endl;
	}
	// Add lubrication Part to tracer Resistance Matrix and invert	
    _RMLub = _resMNoLub + lubM;
	_tracerMM = invert3x3(_RMLub);
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



//----------------------- Ewald sum ------------------------

Matrix3d  CConfiguration::realSpcSm( const Vector3d & rij, const bool self, const double asq ){  // This should be distance of particles
    const int nmax = _nmax;
	int maxIter = 2*nmax;
	const double r_cutoffsq = _r_cutoffsq;
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
				if ((n1-nmax) == 0 && (n2-nmax) == 0 && (n3-nmax) == 0 && self) continue;
				else{  
					rn_vec(2) = v3[n3];
					const double rsq = (rn_vec.transpose() * rn_vec);
	                if ( rsq <= r_cutoffsq ){ 
						Mreal += realSpcM(rsq, rn_vec, asq);
						
					}
				}
			}
		}
	}
	return Mreal;
}



Matrix3d  CConfiguration::reciprocalSpcSm( const Vector3d & rij, const double asq ){  // This should be distance of particles
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
					const double ksq = (kvec.transpose()*kvec);
                    if ( ksq <= k_cutoffsq ){  
						Mreciprocal += reciprocalSpcM(ksq, kvec, asq) * cos((kvec.transpose()* rij));
					}
				}
            }
        }
    }

    return Mreciprocal;   // V is taken care of in reciprocalSpcM 
        
}

Matrix3d  CConfiguration::realSpcM(const double & rsq, const Vector3d & rij, const double asq) {
    // Idea: To only account for tracer in one box, leave out eta = eta_tracer in sums entirely or something...
	Matrix3d  I = Matrix3d::Identity();   
    const double r = sqrt(rsq);
	const double alpha = _alpha;
	const double alphasq = alpha * alpha;
	const double c1 = alphasq * rsq;
	const double c2 = 4. * alphasq * alphasq * alphasq * rsq * rsq;
	const double c3 = alphasq * alphasq * rsq;
	const double c4 = 1./rsq;
	const double c5 = asq/rsq;
	const double expc = exp(-c1)/sqrt (M_PI) * alpha;
	const double erfcc = erfc(alpha * r) / r;
	
	return _pradius * ( I * ( erfcc *  ( 0.75 + 0.5 * c5 )
		+ expc * (3. * c1 - 4.5 + asq * ( c2 - 20. * c3 + 14. * alphasq + c4 )))
			+ (rij*rij.transpose()) / rsq * (
				erfcc * ( 0.75 - 3. * c5 ) +
					expc * ( 1.5 - 3. * c1 + asq * (- c2 + 16. * c3 - 2. * alphasq - 3. * c4))));
}


Matrix3d  CConfiguration::reciprocalSpcM(const double ksq, const Vector3d & kij, const double asq) {
    Matrix3d  I = Matrix3d::Identity();
	const double alphasq = _alpha * _alpha;
    const double c1 = ksq / ( 4. * alphasq );

	return _pradius * ( 1. - 0.333333 * asq * ksq ) * ( 1. + c1 + 2. * c1 * c1 ) * ( 6. * M_PI / (ksq * _V)) * exp( -c1 ) 
        * ( I - (kij*kij.transpose())/ksq);
		//alternative way from Brady1988 TODO alt
    // return _pradius  * ( 1. + c1 + 2. * c1 * c1 ) * ( 6. * M_PI / (ksq * V)) * exp( - c1 ) * ( I - (kij*kij.transpose())/ksq);
}




//------------------------------------- Lubrication Resistance Matrix ------------------------------

Matrix3d  CConfiguration::lubricate( const Vector3d & rij ){
	// Addition to lubrication Matrix for nearest edgeparticles
	Matrix3d  lubPart = Matrix3d::Zero();
    // TODO nmax !!!
	Vector3d rn_vec(3);
	Matrix3d  Mreal = Matrix3d::Zero();
//     for (int n1 = -1; n1 <= 1; n1++){
    double values_1[3] = {rij(0) - _boxsize, rij(0), rij(0) + _boxsize};
    double values_2[3] = {rij(1) - _boxsize, rij(1), rij(1) + _boxsize};
    double values_3[3] = {rij(2) - _boxsize, rij(2), rij(2) + _boxsize};
    for (int n1 = 0; n1 < 3; n1++) {
		rn_vec(0) = values_1[n1];
		for (int n2 = 0; n2 < 3; n2++) {
			rn_vec(1) = values_2[n2];
		    for (int n3 = 0; n3 < 3; n3++) {
				rn_vec(2) = values_3[n3];
				const double rsq = (rn_vec.transpose() * rn_vec);
                if ( rsq <= _cutofflubSq ){ 
					if (rsq < _cutofflubSq/2) lubPart += lub2p(rn_vec, rsq, 8);
					else lubPart += lub2p(rn_vec, rsq, 5);
				}
			}
		}
	}
	return lubPart;
}

Matrix3d CConfiguration::lub2p( Vector3d rij, double rsq, unsigned int mmax ){
	// This function returns the 3x3 lubrication part of the resistance matrix of the tracer particle
    
    //If there is particle overlap, set particle surface distance for lubrication to be very small, according to Phung1996
    if (rsq <= _stericrSq) rsq = _stericrSq + 0.0000000001;
    
	double s = 2*sqrt(rsq)/(_pradius + _polyrad);
	// cout << "\ns " << s << endl;
	double c1 = pow(2/s, 2);
	// cout << "c1 " << c1 << endl;
	double Sum1 = - c1 * ( _g[2] + _g[1] );
	double Sum2 = c1;
	// cout << "Sum1: " << Sum1 << " $$$$ Sum2: " << Sum2 << endl;
	for (int m = 2; m < mmax; m++){
		double c2 = pow(c1, m);
		Sum1 += c2/m * ( _g[2]/(m-1) - _g[1]);
		Sum2 += c2;
	}
	Sum2 = Sum2 * _g[0];
	// cout << "Sum1: " << Sum1 << " $$$$ Sum2: " << Sum2 << endl;
	double c3 = - ( _g[1] + _g[2] * ( 1 - c1 ) ) * log( 1 - c1 );
	double c4 = ( _g[0]/(1-c1) - _g[0] +  2*c3  +  2*Sum1  +  Sum2 ) / (rsq);
	Matrix3d lubR = Matrix3d::Identity() * (c3 + Sum1) + rij * rij.transpose() * c4;

	//cout << "\n------- lubR -----" << endl;
	//cout << lubR << endl;
	return lubR; // * (_pradius + _polyrad)/(2*_pradius);  // TODO this correction for lubrication normalization is wrong I think
}







Matrix3d CConfiguration::CholInvertPart (const MatrixXd A) {
	MatrixXd I = MatrixXd::Identity(A.rows(),A.rows());
	
	// make sure partInv is 3x3 matrix and A is NxN with N larger 2
	assert( A.rows() > 2 );
	assert( A.rows() == A.cols() );

    // perform inversion by cholesky decompoposition and return upper left 3x3 block
    return A.ldlt().solve(I).block<3,3>(0,0);
}

Matrix3d CConfiguration::invert3x3 (const Matrix3d A) { 
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
    _posHistoM.resize(100);
    for (int i = 0; i < 100; i++){
        _posHistoM[i].resize(100);
        for (int j = 0; j < 100; j++){
            _posHistoM[i][j].resize(100, 0);  //initialize all the 100*100*100 matrix elements to zero!
        }
    }
}

void CConfiguration::addHistoValue(){
    //adds a value to the position histogram
    int x = _ppos(0) / _boxsize * 100;     //CAREFUL: THIS CAN'T BE DONE AT A POINT WHERE X MIGHT BE ZERO!!!
    int y = _ppos(1) / _boxsize * 100;
    int z = _ppos(2) / _boxsize * 100;
    if ((x < 0) || (y < 0) || (z < 0) || (x > 99) || (y > 99) || (z > 99)){
        cout << "The Position Histogram function 'conf.addHisto()' is in a bad place, since there is a negative position _ppos()" << endl;
        cout << "The current position is: " << _ppos(0) << " " << _ppos(1) << " " << _ppos(2) << endl;
    }
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

void CConfiguration::resetposition(){
    //Reset the position to random (allowed) position in cell.
    boost::mt19937 rng;    
	boost::uniform_01<boost::mt19937&> zeroone(*m_igen);
	bool overlap = true;
	while (overlap == true){
	    for (int i = 0; i < 3; i++){
	        _entryside[i] = 0;
			double ranPos = zeroone();
	        _startpos(i) = ranPos;
	        _ppos(i) = ranPos;
	        _boxnumberXYZ[i] = 0;
	    }
		overlap = testOverlap();
	}
}



//----------------------- Rotne Prager ------------------------


 Matrix3d  CConfiguration::RotnePrager(const double & r, const double & rsq,
 										   const Vector3d & rij) {
// 	Matrix3d  I = Matrix3d::Identity();
//
// 	double c1 = 0.75 * _polyrad / r; 	//3. * a / (4. * r);
// 	double c2 = 2. * _polyrad * _polyrad / rsq;
// 	return c1*( (1.+c2/3)*I + prod( (1.-c2)*I , matrix<double>((rij*rij.transpose())) / rsq));
// 	                         // = outer_prod(v1,v2) * (1-c2)  works also for the second part of this eq.
}

Matrix3d CConfiguration::RotnePragerDiffRad(const double & r, const double & rsq,
										          const Vector3d & rij) {
	// Matrix3d  I = Matrix3d::Identity();
	//
	// double c1 = 0.75 * _pradius / r; 	//3. * _p/2 / (4. * _r);
	// double c2 = (_polyrad * _polyrad + _pradius * _pradius) / rsq;
	// return c1*( (1.+c2/3)*I + prod( (1.-c2)*I , matrix<double>((rij*rij.transpose())) / rsq));
}



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



