
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <string.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sstream>
#include <math.h>
#include <boost/filesystem.hpp>
#include "CAverage.h"
#include "CConfiguration.h"
#include "parameter_structs.h"

#include <stdlib.h>     //for using the function sleep

// Create global instances of structs
struct sim_param_desc _simpar;
struct model_param_desc _modelpar;
struct file_desc _files;
struct sim_triggers _triggers;


//Function declarations
template<typename T>
string toString(const T& value){
    ostringstream oss;
    oss << value;
    return oss.str();
}

template <typename T, size_t N>
inline
size_t sizeOfArray( const T(&)[ N ] )
{
  return N;
}

void createDataFolder(string testcue){
    //NOTE: Maybe I can leave out dt, as soon as I settled on a timestep
    //NOTE: As soon as I create input-list with variables, I must change this function
    char range[5];
    sprintf(range, "%.3f", _modelpar.urange);
    //In the definition of folder, the addition has to START WITH A STRING! for the compiler to know what to do (left to right).
    _files.folder = "sim_data/noreset";
    if (_triggers.fitRPinv) _files.folder = _files.folder + "/fitRPinv";
    if (!testcue.empty()) _files.folder = _files.folder + "/test" + testcue;
    _files.folder += "/ewaldCorr";
    if (_triggers.ranSpheres) _files.folder = _files.folder +  "/ranSpheres";
    if (_triggers.noLub) _files.folder = _files.folder +  "/noLub";
    if (_triggers.ranPot) _files.folder = _files.folder + "/ranPot";
    if (_triggers.includeSteric) _files.folder = _files.folder + "/steric";    //TODO steric2
    if (_triggers.hpi) _files.folder = _files.folder + "/HPI/hpiu" + toString(_modelpar.hpi_u) + "/hpik" + toString(_modelpar.hpi_k);
    _files.folder = _files.folder
            + "/dt" + toString(_simpar.timestep)
            + "/t" + toString(_simpar.simtime)
            + "/a" + toString(_modelpar.polymersize)
            + "/d" + toString(_modelpar.rodDist)
            + "/b" + toString(_modelpar.boxsize)
            + "/p" + toString(_modelpar.particlesize)
            + "/k" + range
            + "/u" + toString(_modelpar.ustrength);
    boost::filesystem::create_directories(_files.folder);
    boost::filesystem::create_directory(_files.folder + "/InstantValues");
    boost::filesystem::create_directory(_files.folder + "/Coordinates");
    cout << "Writing data to folder:\n" << _files.folder << endl;
}


void parameterFile(string testcue){
    //Creates a file where the simulation settings are stored
    ofstream parameterFile;
    parameterFile.open((_files.folder + "/parameters.txt").c_str());

    // Print time and date
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    parameterFile << "date " << (now->tm_year + 1900) << '-'
         << (now->tm_mon + 1) << '-'
         <<  now->tm_mday
         << endl;
    parameterFile << "starttime " << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << endl;
    parameterFile << "Sim_dir " << _files.folder << endl;
    if (!testcue.empty()) parameterFile << "Test cue " << testcue << endl;
    parameterFile << "fitRPinv " << _triggers.fitRPinv << endl;
    parameterFile << "ewaldCorr " << true << endl;
    parameterFile << "ranSpheres " << _triggers.ranSpheres << endl;
    parameterFile << "noLub " << _triggers.noLub << endl;
    parameterFile << "recordPosHisto " << _triggers.recordPosHisto << endl;
    parameterFile << "steric " << _triggers.includeSteric << endl;
    parameterFile << "ranPot " << _triggers.ranPot  << endl;
    parameterFile << "HPI " << _triggers.hpi  << endl;
    if (_triggers.hpi == true){
        parameterFile << "hpi_u " << _modelpar.hpi_u  << endl;
        parameterFile << "hpi_k " << _modelpar.hpi_k  << endl;
    }
    parameterFile << "d " << _modelpar.rodDist << endl;
    parameterFile << "p " << _modelpar.particlesize << endl;
    parameterFile << "b " << _modelpar.boxsize << endl;
    parameterFile << "dt " << _simpar.timestep << endl;
    parameterFile << "runs " << _simpar.runs << endl;
    parameterFile << "steps " << _simpar.steps << endl;
    parameterFile << "time " << _simpar.timestep*_simpar.steps << endl;
    parameterFile << "k " << _modelpar.urange << endl;
    parameterFile << "U_0 " << _modelpar.ustrength << endl;
    parameterFile << "a " << _modelpar.polymersize << endl;

    parameterFile.close();
}

void parameterFileAppend(double executiontime){
    // Appends parameters to parameters file
    ofstream parameterFile;
    parameterFile.open((_files.folder + "/parameters.txt").c_str(), std::ios_base::app);
    parameterFile << "ExecutionTime " << executiontime << " s" << endl;
    parameterFile.close();
}


