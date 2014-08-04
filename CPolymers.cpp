#include "headers/CPolymers.h"

using namespace std;


CPolymers::CPolymers(){
    //initiate Polymers by giving all 3*3*4 = 36 a random sign.
    setRanNumberGen(0);
    _sign.resize(3);
    for (int i = 0; i < 3; i++){
        _sign[i].resize(3);
        for (int j = 0; j < 3; j++){
            _sign[i][j].resize(4);  //initialize all the 100*100*100 matrix elements to zero!
            for (int r = 0; r < 4; r++){
                _sign[i][j][r] = ran_sign();
            }
        }
    }
	//printSignFile("w");
}


void CPolymers::shiftPolySign(int crossaxis, int exitmarker){
    int dir = exitmarker;   //is either -1 or 1
    int k = crossaxis;      //is either 0, 1 or 2
    int m = 1+dir;
    int l = 1-dir;
    _sign[l][k] = _sign[1][k];
    _sign[1][k] = _sign[m][k];
    for (int i=0; i < 4; i++){
        _sign[m][k][i] = ran_sign();
    }
    //two integers that stand for the direction above and below the current, with periodic order 0,1,2,0 = x, y, z, x
    int above = k + 1;
    if (above==3) above = 0;

    int below = k - 1;
    if ( below == -1) below = 2;

    for (int j = 0; j < 3; j++){
    _sign[j][above][l] = _sign[j][above][m];
    _sign[j][above][2-dir] = _sign[j][above][2+dir];
    _sign[j][below][(int)((l)/2)] = _sign[j][below][(int)((m)/2)];
    _sign[j][below][(int)(2+(l)/2)] = _sign[j][below][(int)(2+(m)/2)];

    //reassign new random signs
    _sign[j][above][m] = ran_sign();
    _sign[j][above][2+dir] = ran_sign();
    _sign[j][below][(int)((m)/2)] = ran_sign();
    _sign[j][below][(int)(2+(m)/2)] = ran_sign();
    }
	//printSignFile("a");
}


/*int CPolymers::ran_sign(){
    boost::random::uniform_int_distribution<> onetwo(-1,1);  // distribution that maps to -1,1
    return onetwo(_rng);
}
*/

int CPolymers::ran_sign(){
    // the variate generator uses _igen (int rand number generator),
    // samples from uniform integer distribution 0, 1
    boost::variate_generator<boost::mt19937&, 
	                         boost::uniform_int<> 
							 > zeroone(
                             *_igen, boost::uniform_int<>(0, 1));
	
	int plusminus = (zeroone() * 2) - 1;
	return plusminus; //this calculation makes value either -1 or 1 
}

void CPolymers::setRanNumberGen(double seed){
    if (seed == 0) {
        _igen = new boost::mt19937(static_cast<unsigned int>(time(NULL)));
    } else {
        _igen = new boost::mt19937(static_cast<unsigned int>(seed));
    }
}

bool CPolymers::samesign(const int direction, const int plane_axis, const int edgenumber){
    //return true, if two adjacent polymer "pieces" are of the same sign (charge)
    return (_sign[1][plane_axis][edgenumber] == _sign[1+direction][plane_axis][edgenumber]);
}

int CPolymers::get_sign(const int plane_axis, const int edgenumber){
    return _sign[1][plane_axis][edgenumber];
}



void CPolymers::printSignFile(const string a_w){
    string name = "polymercharge.xyz";
    FILE *f = fopen(name.c_str(), a_w.c_str());

    fprintf(f, "%d\n%s (%8.3f %8.3f %8.3f)\n", 36, "sim_name", 10, 10, 10);
	
    double L1[] = {0, 10, 0, 10};  //these two arrays are only needed for the iteration.
    double L2[] = {0, 0, 10, 10};
    double coordinate[3];
    string color = "P";
	
    // polymer chains
	
	for (int i = 0; i < 3; i++){
		for (int axis=0; axis < 3; axis++){
			for (int pos=0; pos<4; pos++){
			    int above = axis + 1;
			    if (above==3) above = 0;

			    int below = axis - 1;
			    if ( below == -1) below = 2;
				coordinate[axis] = 10*i - 5;
				coordinate[above] = L1[pos];
				coordinate[below] = L2[pos];
				if (_sign[i][axis][pos] == -1) color = "H";
				else if (_sign[i][axis][pos] == 1) color = "O";
				else color = "P";
				fprintf(f, "%3s%9.3f%9.3f%9.3f \n",color.c_str(), coordinate[0], coordinate[1], coordinate[2]);
			}
		}
	}   
    fclose(f);
}



