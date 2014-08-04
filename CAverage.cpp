/*
 * CAverage.cpp
 *
 *  Created on: May 4, 2010
 *      Author: radtke
 */
/*
 * Version modified by Johann Hansing. Original Version by Matthias Radtke
 */

#include "headers/CAverage.h"


using namespace std;

CAverage::CAverage() {
}

CAverage::CAverage(string name, string folder, int datapoints, int runs)
    : _instantValues (datapoints), _errors (datapoints)      //vectors will automatically be initialized to zero
{
	_points = datapoints;
	_name = name;
	_runs = runs;

	clear();

	// create empty save file..
	ofstream myfile;
	string s;
	ostringstream outStream;
	outStream << folder << "/InstantValues/" << _name << ".dat";
	s = outStream.str();
	_file_instV = s;
	myfile.open(_file_instV.c_str());
	myfile.close();
}


/* adds Value v as instant Value to vector at datapoint n*/
void CAverage::addInstantValue(double v, int n) {
	_instantValues[n] += v;
	//Formula from Wikipedia: Standartabweichung: Bearbeitung fŸr auflaufende Messwerte + Standardfehler
	//continues saveAverageInstantValues()!
	_errors[n] += v*v;
}

/*
void CAverage::addInstantValue(double v, int n) {                                   //TODO
    _instantValues.push_back(v);
}
*/

/* adds Value v as instant Value to vector at datapoint n*/
void CAverage::addFPValue(double fpt, int n) {    //fpt is first passage time
    if (_instantValues.size() < n){
        _instantValues.push_back(fpt);
        _errors.push_back(fpt*fpt);
        _instances.push_back(1);
    }
    else {
        _instantValues[n-1] += fpt;
        _errors[n-1] += fpt*fpt;
        _instances[n-1] += 1;
    }
    //Formula from Wikipedia: Standardabweichung: Bearbeitung fŸr auflaufende Messwerte + Standardfehler
    //continues saveAverageInstantValues()!

}


/* calculates an average value of all saved instant values
double CAverage::getAverage() {
	double sum = 0.0;
	for (unsigned int i = 0; i < _instantValues.size(); i++) {
		sum += _instantValues[i];
	}

	return (sum / double(_instantValues.size()));
}
*/

void CAverage::clear() {
	_instantValues.clear();
//	_averageValues.clear();
}

/* save instant values to file */

void CAverage::saveAverageInstantValues(double timeInt) {
    //timeInt is the interval between the times that an average instant value was saved
	if (_points > 0) {
		ofstream myfile;
		myfile.open(_file_instV.c_str());


		for (unsigned int i = 0; i < _points; i++) {
			myfile << ((i+1)*timeInt)<< " " << (_instantValues[i]/_runs) << " " <<
			//the next formula is to calculate the error from the standard deviation
			        //CAREFUL No error for _runs = 1, due to divide by zero!
			sqrt((_errors[i] - _instantValues[i] * _instantValues[i] / _runs)/(_runs*(_runs - 1))) <<
			endl;
		}
		myfile.close();

		_instantValues.clear();
	}
}

/* This works, but I would have to have 3 of these for _ppos[i] -> i=0,1,2
void CAverage::saveAverageInstantValues(double timeInt) {                                    //TODO
    //timeInt is the interval between the times that an average instant value was saved
    if (_points > 0) {                                                                       //TODO
        ofstream myfile;
        myfile.open(_file_instV.c_str());


        for (unsigned int i = 1; i < (_instantValues.size()/1000); i++){

            double value=0;
            double error=0;
            int n = 0;
            while ((n+i) < _instantValues.size()){
                value +=  pow((_instantValues[n] - _instantValues[n+i]), 2);
                error += value*value;
                n++;
            }

            myfile << (i*timeInt)<< " " << value/n << " " <<
                        //the next formula is to calculate the error from the standard deviation
                                //CAREFUL No error for _runs = 1, due to divide by zero!
                        sqrt((error - value*value / n)/(n*(n - 1))) <<   endl;
        }

        myfile.close();

        _instantValues.clear();
    }
}
/*


/* save Pushed Back instant values to file */
void CAverage::saveAverageFPValue(double posInt) {
    //timeInt is the interval between the times that an average instant value was saved
    int size = _instantValues.size();
    if (size > 0) {
        ofstream myfile;
        myfile.open(_file_instV.c_str());


        for (unsigned int i = 0; i < size; i++) {
            myfile << ((i+1)*posInt)<< " " << (_instantValues[i]/_instances[i]) << " " <<
            //the next formula is to calculate the error from the standard deviation
                    //CAREFUL No error for _runs = 1, due to divide by zero!
            sqrt((_errors[i] - _instantValues[i] * _instantValues[i] / _instances[i])/(_instances[i]*(_instances[i] - 1))) << " " <<
            _instances[i] <<
            endl;
        }
        myfile.close();

        _instantValues.clear();
    }
}
