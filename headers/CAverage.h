/*
 * CAverage.h
 *
 *  Created on: May 4, 2010
 *      Author: Johann Hansing
 */

#ifndef CAVERAGE_H_
#define CAVERAGE_H_

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>


class CAverage {
private:
	int _points;
	int _runs;
	std::string _name;
	std::string _file_instV;
	std::vector<double> _instantValues;
	std::vector<double> _errors;
	std::vector<int> _instances;      //instances records how many particle runs have contributed to take the average of the mfp time
//	std::vector<double> _averageValues;
	bool saveFileExists;

public:
	CAverage();
	CAverage(std::string name, std::string folder, int datapoints, int runs);

	void addInstantValue(double, int);
	void addFPValue(double fpt, int n);

//	double getAverage();
	void clear();

	void saveAverageInstantValues(double timeInt);
	void saveAverageFPValue(double posInt);
	std::string getName() {	return _name; }
};

#endif /* CAVERAGE_H_ */
