#include "headers/CConfiguration.h"
#include "headers/CholDecomp.h"


using namespace ublas;

using namespace std;

const double _6root2 = 1.122462;


CConfiguration::CConfiguration(){
}

CConfiguration::CConfiguration(
        double timestep,  double potRange,  double potStrength,  double boxsize, double rodDistance, const bool potMod,
        double psize, const bool posHisto, const bool steric, const bool ranU, bool hpi, double hpi_u, double hpi_k,
		double polymersize)
		{
    _potRange = potRange;
    _potStrength = potStrength;
    _pradius = psize/2;   //_pradius is now the actual radius of the particle. hence, I need to change the definition of the LJ potential to include (_pradius + _polyrad)   -  or maybe leave LJ pot out
    _boxsize = boxsize;
    _resetpos = _boxsize/2;
    _timestep = timestep;
    _rodDistance = rodDistance;
    _potMod = potMod;
    _LJPot = (steric == false) && (psize != 0);
    _ranU = ranU;
    _poly = CPolymers();
    _hpi = hpi;
    _upot = 0;
    _mu_sto = sqrt( 2 * _timestep );                 //timestep for stochastic force
	_hpi = hpi; 
	_hpi_u = hpi_u;
	_hpi_k = hpi_k;
    for (int i = 0; i < 3; i++){
        _ppos[i] = _resetpos;
        _startpos[i] = _resetpos;
        _entryside[i] = 0;
        _wallcrossings[i] = 0;
        _boxnumberXYZ[i] = 0;
        _prevpos[i] = _resetpos;
    }

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
	
	// init HI vectors matrices, etc
	_tracerMM = identity_matrix<double> (3);   
	_polyrad = polymersize / 2;   //This is needed for testOverlap for steric
    if (polymersize != 0) _HI = true;
	if (_HI) {
		_edgeParticles = (int) (10/polymersize);
		_epos.resize(3 * _edgeParticles - 2);
		initConstMobilityMatrix();
		_LJPot = false;
	}

}

void CConfiguration::updateStartpos(){
    //This function is used if the particle should keep moving after each run, and not start at _resetpos again, like in the first run
    //This (hopefully) will give better averages without having to spend a lot of steps in the beginning of each run to get away from _resetpos
    for (int i = 0; i < 3; i++){
    _startpos[i] = _ppos[i] + _boxsize * _boxnumberXYZ[i];
    }
}



void CConfiguration::makeStep(){
    //move the particle according to the forces and record trajectory like watched by outsider
//	if (_HI){
	cout << _tracerMM << endl;  //TODO del
	
		ublas::vector<double> F(3,0.0);
	    for (int i = 0; i < 3; i++){
			F(i) = _f_mob[i];
	        _prevpos[i] = _ppos[i];
	    }
		ublas::vector<double> MMdotF (3, 0.0);
		
		noalias(MMdotF) = prod(_tracerMM,F);
		
		// Update particle position
		for (int i = 0; i < 3; i++){
		    _ppos[i] += MMdotF(i) * _timestep + _f_sto(i) * _mu_sto;
		}
		
/*	}
	else{
	    for (int i = 0; i < 3; i++){
	        _prevpos[i] = _ppos[i];
	        _ppos[i] += _timestep * _f_mob[i] + _mu_sto * _f_sto[i];
	    }
	}
  */  
}

void CConfiguration::checkBoxCrossing(){
    //should the particle cross the confinement of the cube, let it appear on the opposite side of the box
    for (int i = 0; i < 3; i++){
        if (_ppos[i] < 0){
            _ppos[i] += _boxsize;
            _boxnumberXYZ[i] -= 1;
            countWallCrossing(i, -1);
            if (_ranU) _poly.shiftPolySign(i, -1);

        }
        else if (_ppos[i] > _boxsize){
            _ppos[i] -= _boxsize;
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
	
	ublas::vector<double> ran_v(3);
	ran_v.clear();
	
    ran_v(0) = ran_gen();
	ran_v(1) = ran_gen();
	ran_v(2) = ran_gen();
    
	if (_HI){  		
		// taken from matthias:
		// init lower triangular matrix L(3Nx3N) and stand norm dist rand vector ran_v(3N)
		triangular_matrix<double, lower> L(3, 3);
		L.clear();

		size_t res = cholesky_decompose(_tracerMM, L);
		// test wheter decomposition was correct...
		if(res != 0) {
			cout << "Decomposition failed in row "<< res - 1 << " !"<< endl;
		}

		// return correlated random vector, which is scaled later by sqrt(2 dt)
		noalias(_f_sto) = prod(L, ran_v);		
	}

    else { // no HI
        noalias(_f_sto) = ran_v;
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
    for (int l = 0; l < 3; l++) {
        _f_mob[l] = 0;
    }
    const double r_c = _6root2 * _pradius;    //cutoff for Lennard-Jones calculation (at minimum)

    for (int i = 0; i < 3; i++){
        int k = i + 1;   //k always one direction "further", i.e. if i = 0 = x-direction, then k = 1 = y-direction
        if ( k == 3 ) k = 0;
        int plane = 3 - (i+k); //this is the current plane of the cylindrical coordinates
        int n = 0;     // reset counter for index of next rod in plane  n = 0, 1, 2, 3 -> only needed for ranPot
        for (int nk = _min; nk < _max; nk++){
            for (int ni = _min; ni < _max; ni++){
            r_i = _ppos[i] - ni*_boxsize;
            r_k = _ppos[k] - nk*_boxsize;
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
			
			if (_hpi) calculateExpHPI(r_abs, utmp, frtmp);

            if (_ranU){
                utmp = utmp * _poly.get_sign(plane, n);
                frtmp = frtmp * _poly.get_sign(plane, n);
                if (_ppos[plane] > z2){
                    if (! _poly.samesign(1, plane, n)){
                        _f_mob[plane] += utmp * 4 / _boxsize;              //this takes care of the derivative of the potential modification and resulting force
                        modifyPot(utmp, frtmp, _boxsize - _ppos[plane]);
                    }
                }
                else if (_ppos[plane] < z1){
                    if (! _poly.samesign(-1, plane, n)){
                        _f_mob[plane] -= utmp * 4 / _boxsize;              //this takes care of the derivative of the potential modification and resulting force
                        modifyPot(utmp, frtmp, _ppos[plane]);
                    }
                }
                n++;  //index of next rod in curent plane
            }


            if (_LJPot && ( r_abs < r_c || _hpi )) calcLJPot(r_abs, utmp, frtmp);


            Epot += utmp;
            _f_mob[i] += frtmp * r_i;
            _f_mob[k] += frtmp * r_k;
            }
        }
    }
    _upot = Epot;
	
	// calc HI mobility matrix here, since it needs to be defined for random force normalisation
	if (_HI){ calcTracerMobilityMatrix(); }
}


void CConfiguration::saveXYZTraj(string name, const int& move, string a_w){
    FILE *f = fopen(name.c_str(), a_w.c_str());

    fprintf(f, "%d\n%s (%8.3f %8.3f %8.3f) t=%d \n", 1, "sim_name", _boxsize, _boxsize, _boxsize, move);

    // polymer particles
    fprintf(f, "%3s%9.3f%9.3f%9.3f \n","H", _ppos[0]+ _boxsize *_boxnumberXYZ[0], _ppos[1]+ _boxsize *_boxnumberXYZ[1], _ppos[2]+ _boxsize *_boxnumberXYZ[2]);   // absolute position
	// fprintf(f, "%3s%9.3f%9.3f%9.3f \n","H", _ppos[0], _ppos[1], _ppos[2]);  // relative position in box
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
        var += pow((_ppos[m] + _boxsize *_boxnumberXYZ[m] - _startpos[m]) , 2);
    }
    return var;
}

double CConfiguration::get1DPosVariance(int dim){
    //function to return the variance of just one dimension x = 0, y  = 1 or z = 2
    return pow((_ppos[dim] + _boxsize *_boxnumberXYZ[dim] - _startpos[dim]) , 2);
}




bool CConfiguration::checkFirstPassage(double mfpPos, int dim){    //TODO
    //function returns true if in the last step the momentary "way point" (_mfpIndex * xinterval) has been crossed in dim-direction
    // where dim = 0, 1, 2 -> x, y, z. Otherwise it returns false
    if ((_ppos[dim] + _boxsize*_boxnumberXYZ[dim] - _startpos[dim]) > (mfpPos) ) return true;
    else return false;
}



void CConfiguration::moveBack(){
    //moves particle back to previous position
    for (int i = 0; i < 3; i++) {_ppos[i] = _prevpos[i];}
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

void CConfiguration::calculateExpPotentialMOD(const double r, double& U, double& Fr, int plane){
    //function to calculate an exponential Potential U = U_0 * exp(-1 * r * k)
    // k is the interaction range. U_0 is the strength of the potential
    //which is attractive if direction = -1, and repulsive if direction = 1
    //The potential is weighted with kT!
    //index, is the dimension (x, y or z) which is normal to the plane that the potential is being calculated in.

	 // U = _potStrength * exp(-1 * r / _potRange) * ( 1 - 2 / 3 * pow((2*_ppos[3-index]/_boxsize - 1), 2));
	 U = _potStrength * exp( -r / _potRange) * (1 - abs(2 * _ppos[plane]/_boxsize - 1) * 0.66667);
	 Fr = U / (_potRange * r);  //This is the force divided by the distance to the rod!
	 
    //DEBYE!!!
   // U = _potStrength * exp(-1 * r / _potRange) / r;
    //    Fr = U * (1/_potRange + 1/r) / r;  //This is the force divided by the distance to the rod!
	 
    // BESSEL
  //  U = 2 * _potStrength * boost::math::cyl_bessel_k(0, r / _potRange);
  //  Fr = 2 * _potStrength * boost::math::cyl_bessel_k(1, r / _potRange) / (_potRange * r);
}

void CConfiguration::modifyPot(double& U, double& Fr, double dist){
    //function to modify the potential according to the distance along the polymer axis to the next neighbor,
    //in case the next neighboring polymer part is of opposite sign
    U = U * 4 * dist/_boxsize;
    Fr = Fr * 4 * dist/_boxsize;
}

//****************************STERIC HINDRANCE****************************************************//

bool CConfiguration::testOverlap(){
    //Function to check, whether the diffusing particle of size psize is overlapping with any one of the rods (edges of the box)
    //most if borrowed from moveParticleAndWatch()
    bool overlaps = false;
    double r_i = 0, r_k = 0;
    double r_abs = 0;
    double L1[] = {0, _boxsize, 0, _boxsize};  //these two arrays are only needed for the iteration.
    double L2[] = {0, 0, _boxsize, _boxsize};

    for (int i = 0; i < 2; i++){
        for (int k = i+1; k < 3; k++){
            for (int n = 0; n < 4; n++){
                r_i = _ppos[i] - L1[n];
                r_k = _ppos[k] - L2[n];
                //this is needed if we dont want the rods to cross each other to create a strong potential well
                if (k == 2){
                    r_i -= _rodDistance;
                    if (i == 0) r_k -= _rodDistance;
                }
                r_abs = sqrt(r_i * r_i + r_k * r_k); //distance to the rods
                if (r_abs < (_pradius + _polyrad)) overlaps = true;
            }
        }
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
	double rxij = 0.0, ryij = 0.0, rzij = 0.0;
	double rij = 0.0, rijsq = 0.0;
	
	// ublas-vector for outer product in muij
	ublas::vector<double> vec_rij(3, 0.0);
	
	// create mobility matrix - Some elements remain constant throughout the simulation. Those are stored here.
	_mobilityMatrix = identity_matrix<double>(3 * (3 * _edgeParticles - 1));
	
	// set _edgeParticles positions _pos to zero
	for (int i = 0; i < 3 * _edgeParticles - 2; i++){
	        _epos[i].fill(0);
	}
	
	// store the edgeParticle positions, so that I can simply loop through them later
	for (int i = 1; i < _edgeParticles; i++){
		double tmp = i * 2 * _polyrad;
		_epos[i][0] = tmp;
		_epos[i + (_edgeParticles - 1)][1] = tmp;  
		_epos[i + 2 * (_edgeParticles - 1)][2] = tmp;
	}
	
	// on-diagonal elements are self mobilities. Thus = 1 for tracer and lambda = p/a for edgeparticles
	//_mobilityMatrix(0,0) = 1;  //already taken care of due to identity matrix?
	//_mobilityMatrix(1,1) = 1;
	//_mobilityMatrix(2,2) = 1;
		
	for (unsigned int i = 1; i < 3 * _edgeParticles - 1 ; i++) {
		unsigned int i_count = 3 * i;

		double povera = _pradius/_polyrad;
		_mobilityMatrix(i_count,i_count) = povera;
		_mobilityMatrix(i_count+1,i_count+1) = povera;
		_mobilityMatrix(i_count+2,i_count+2) = povera;
	}
	
	// The mobility matrix elements for the interacting edgeParticles are stored here, since they do not change
	for (unsigned int i = 0; i < 3 * _edgeParticles - 2; i++){
		rxi = _epos[i][0];
		ryi = _epos[i][1];
		rzi = _epos[i][2];
		unsigned int i_count = 3 * (i + 1);       // plus 1 is necessary due to omitted tracer particle
		
		for (unsigned int j = i + 1; j < 3 * _edgeParticles - 2; j++) {
			rxij = rxi - _epos[j][0];
			ryij = ryi - _epos[j][1];
			rzij = rzi - _epos[j][2];
			unsigned int  j_count = 3 * (j + 1);    // plus 1 is necessary due to omitted tracer particle

			rijsq = rxij * rxij + ryij * ryij + rzij * rzij;

			rij = sqrt(rijsq);
		
		/*
		 * Calculation of Rotne-Prager Tensor
		 */
			matrix<double> muij (3,3,0.0);
			vec_rij(0) = rxij;
			vec_rij(1) = ryij;
			vec_rij(2) = rzij;
			
			noalias(muij) = RotnePrager(rij,rijsq,vec_rij);
			
			// only lower triangular of symmetric matrix is filled and used
			noalias(subrange(_mobilityMatrix, j_count, j_count + 3, i_count, i_count + 3)) = muij;
		}
	}		
}


void CConfiguration::calcTracerMobilityMatrix(){
	double rxij = 0.0, ryij = 0.0, rzij = 0.0;
	double rij = 0.0, rijsq = 0.0;
		
	// ublas-vector for outer product in muij
	ublas::vector<double> vec_rij(3, 0.0);
	
// loop over tracer - particle mobility matrix elements
	for (unsigned int j = 0; j < 3 * _edgeParticles - 2; j++){
		rxij = _ppos[0] - _epos[j][0];
		ryij = _ppos[1] - _epos[j][1];
		rzij = _ppos[2] - _epos[j][2];
		unsigned int j_count = 3 * (j + 1); //  plus 1 is necessary due to omitted tracer particle

		rijsq = rxij * rxij + ryij * ryij + rzij * rzij;

		rij = sqrt(rijsq);
		
	/*
	 * Calculation of different particle width Rotne-Prager Tensor
	 */
		matrix<double> muij (3,3,0.0);
		vec_rij(0) = rxij;
		vec_rij(1) = ryij;
		vec_rij(2) = rzij;
		
		noalias(muij) = RotnePragerDiffRad(rij,rijsq,vec_rij);
		
		// only lower triangular of symmetric matrix is filled and used
		noalias(subrange(_mobilityMatrix, j_count, j_count + 3, 0, 3)) = muij;
	}		
	// create resistance matrix - Some elements remain constant throughout the simulation. Those are stored here.
	matrix<double> resistanceMatrix = identity_matrix<double> (3 * (3 * _edgeParticles - 1));
	matrix<double> mobM = _mobilityMatrix;

	//cout << _mobilityMatrix << endl;  //TODO del
	
	
	bool inverted;
	// invert full mobility matrix to obtain resistance matrix
	cout << "invert mobMat" << endl;
	inverted = InvertMatrix(mobM, resistanceMatrix);
	cout << "inverted" << endl;
	if (!inverted){
		cout << "ERROR: Could not invert mobility matrix!" << endl;
		exit (EXIT_FAILURE);
	}
	// invert 3x3 submatrix of resistance matrix to obtain single particle mobility matrix
	matrix<double> submobil (3,3);
	matrix<double> submatrix = subrange(resistanceMatrix,0,3,0,3);
	inverted = InvertMatrix(submatrix,submobil);
	if (!inverted){
		cout << "ERROR: Could not invert resistance submatrix!" << endl;
		exit (EXIT_FAILURE);
	}
	noalias(_tracerMM) = submobil;

	cout << _tracerMM << endl;  //TODO del
	
	
}


matrix<double> CConfiguration::RotnePrager(const double & r, const double & rsq,
										   const ublas::vector<double> & rij) {
	matrix<double> I = identity_matrix<double>(3);

	double c1 = 0.75 * _polyrad / r; 	//3. * a / (4. * r);
	double c2 = 2. * _polyrad * _polyrad / rsq;	
	return c1*( (1.+c2/3)*I + prod( (1.-c2)*I , matrix<double>(outer_prod(rij,rij)) / rsq));
}

matrix<double> CConfiguration::RotnePragerDiffRad(const double & r, const double & rsq,
										          const ublas::vector<double> & rij) {
	matrix<double> I = identity_matrix<double>(3);

	double c1 = 0.75 * _pradius / r; 	//3. * _p/2 / (4. * _r);
	double c2 = (_polyrad * _polyrad + _pradius * _pradius) / rsq;  
	return c1*( (1.+c2/3)*I + prod( (1.-c2)*I , matrix<double>(outer_prod(rij,rij)) / rsq));
}

template<class T> 
bool CConfiguration::InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) { 
    typedef ublas::permutation_matrix<std::size_t> pmatrix; 
    // create a working copy of the input 
    matrix<T> A(input); 
    // create a permutation matrix for the LU-factorization 
    pmatrix pm(A.size1()); 
    // perform LU-factorization 
    int res = lu_factorize(A,pm); 
    if( res != 0 )
        return false; 
    // create identity matrix of "inverse" 
    inverse.assign(ublas::identity_matrix<T>(A.size1())); 
    // backsubstitute to get the inverse 
    lu_substitute(A, pm, inverse); 
    return true; 
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
    int x = _ppos[0] / _boxsize * 100;     //CAREFUL: THIS CAN'T BE DONE AT A POINT WHERE X MIGHT BE ZERO!!!
    int y = _ppos[1] / _boxsize * 100;
    int z = _ppos[2] / _boxsize * 100;
    if ((x < 0) || (y < 0) || (z < 0) || (x > 99) || (y > 99) || (z > 99)){
        cout << "The Position Histogram function 'conf.addHisto()' is in a bad place, since there is a negative position _ppos[]" << endl;
        cout << "The current position is: " << _ppos[0] << " " << _ppos[1] << " " << _ppos[2] << endl;
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


//****************************OLD CODE****************************************************



/*
double CConfiguration::getDisplacement(){
   FAULTY, _f_mob is now a ublas::vector !!!!  double d = 0;
    for (int i = 0; i < 3; i++){
        d += pow((_timestep * _f_mob[i] + _mu_sto * _f_sto[i]), 2);
    }
    return sqrt(d); 
} */

void CConfiguration::resetposition(){
    //Reset the position after every run.
    for (int i = 0; i < 3; i++){
        _entryside[i] = 0;
        _startpos[i] = _resetpos;
        _ppos[i] = _resetpos;
        _boxnumberXYZ[i] = 0;
    }
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
        _f_mob[l] = 0;
    }
    const double r_c = _6root2 * _pradius;    //cutoff for Lennard-Jones calculation (at minimum)

    for (int i = 0; i < 3; i++){
            int k = i + 1;   //k always one direction "further", i.e. if i = 0 = x-direction, then k = 1 = y-direction
            if ( k == 3 ) k = 0;
            int plane = 3 - (i+k); //this is the current plane of the cylindrical coordinates
            for (int n = 0; n < 4; n++){
                r_i = _ppos[i] - L1[n];
                r_k = _ppos[k] - L2[n];
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
                    if (_ppos[plane] > z2){
                        if (! _poly.samesign(1, plane, n)){
                            modifyPot(utmp, frtmp, _boxsize - _ppos[plane]);
                            _f_mob[plane] += utmp * 4 / _boxsize;              //this takes care of the derivative of the potential modification and resulting force
                        }
                    }
                    else if (_ppos[plane] < z1){
                        if (! _poly.samesign(-1, plane, n)){
                            modifyPot(utmp, frtmp, _ppos[plane]);
                            _f_mob[plane] -= utmp * 4 / _boxsize;              //this takes care of the derivative of the potential modification and resulting force
                        }
                    }
                }


                if (_LJPot && ( r_abs < r_c )) calcLJPot(r_abs, utmp, frtmp);


                Epot += utmp;
                _f_mob[i] += frtmp * r_i;
                _f_mob[k] += frtmp * r_k;
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
        _f_mob[l] = 0;
    }
    const double r_c = _6root2 * _pradius;    //cutoff for Lennard-Jones calculation (at minimum)

    for (int i = 0; i < 3; i++){
            int k = i + 1;   //k always one direction "further", i.e. if i = 0 = x-direction, then k = 1 = y-direction
            if ( k == 3 ) k = 0;
            int plane = 3 - (i+k); //this is the current plane of the cylindrical coordinates
            for (int ni = -1; ni < 3; ni++){
                for (int nk = -1; nk < 3; nk++){
                r_i = _ppos[i] - ni*_boxsize;
                r_k = _ppos[k] - nk*_boxsize;

                r_abs = sqrt(r_i * r_i + r_k * r_k); //distance to the rods


                if (_potMod) calculateExpPotentialMOD(r_abs, utmp, frtmp, plane);
                else calculateExpPotential(r_abs, utmp, frtmp);


                if (_LJPot && ( r_abs < r_c )) calcLJPot(r_abs, utmp, frtmp);


                Epot += utmp;
                _f_mob[i] += frtmp * r_i;
                _f_mob[k] += frtmp * r_k;
                }
            }

    }
    _upot = Epot;
}
*/


