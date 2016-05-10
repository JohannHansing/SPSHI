
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
#include "misc/SimpleIni.h"

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
    _files.folder = "sim_data/newlub/noreset";
    if (_triggers.fitRPinv) _files.folder = _files.folder + "/fitRPinv";
    if (_triggers.ranRod) _files.folder = _files.folder +  "/ranRod";
    if (_triggers.trueRan) _files.folder = _files.folder +  "/trueRan";
    else if (_triggers.ranSpheres) _files.folder = _files.folder +  "/ranSpheres";
    if (!testcue.empty()) _files.folder = _files.folder + "/test" + testcue;
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
    parameterFile << "ranRod " << _triggers.ranRod << endl;
    parameterFile << "ranSpheres " << _triggers.ranSpheres << endl;
    parameterFile << "trueRan " << _triggers.trueRan << endl;
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
    parameterFile << "n_cells " << _modelpar.n_cells << endl;
    parameterFile << "EwaldTest " << _modelpar.EwaldTest << endl;
    parameterFile << "nmax " << _modelpar.nmax << endl;
    parameterFile << "lubcutint " << _modelpar.lubcutint << endl;

    parameterFile.close();
}

void parameterFileAppend(double executiontime){
    // Appends parameters to parameters file
    ofstream parameterFile;
    parameterFile.open((_files.folder + "/parameters.txt").c_str(), std::ios_base::app);
    parameterFile << "ExecutionTime " << executiontime << " s" << endl;
    parameterFile.close();
}


void read_run_params(const std::string& filename) {// Inspired by function from Richard
    CSimpleIniA ini;
    ini.SetUnicode();
    if (ini.LoadFile(filename.c_str()) < 0) {
      std::cout << "Could not open run parameters file " << filename << std::endl;
      abort();
    }
    
    //TRIGGERS:
    _triggers.writeTrajectory = (strcmp(ini.GetValue("", "writeTraj", "") , "true") == 0 ) ;    // relative position TODO
    _triggers.fitRPinv = (strcmp(ini.GetValue("", "fitRPinv", "") , "true") == 0 ) ;
    _triggers.ranSpheres = (strcmp(ini.GetValue("", "ranSpheres", "") , "true") == 0 ) ;
    _triggers.recordPosHisto = (strcmp(ini.GetValue("", "recordPosHisto", "") , "true") == 0 ) ;
    _triggers.noLub = (strcmp(ini.GetValue("", "noLub", "") , "true") == 0 ) ;
    _triggers.includeSteric = (strcmp(ini.GetValue("", "includeSteric", "") , "true") == 0 ) ;  // steric 2
	_triggers.ranPot = (strcmp(ini.GetValue("", "ranPot", "") , "true") == 0 ) ;
	_triggers.hpi = (strcmp(ini.GetValue("", "hpi", "") , "true") == 0 ) ;          // hpi exp
	int boolpar = 8;
    
	// Checking for correct structure of input arguments
	// for (int k= 0; k < argc; k++ ) cout << "parameter " << k << " " << argv[k] << endl;
//     for (int b_i=1; b_i<=boolpar; b_i++){
//         if ((strcmp(argv[b_i] , "true") == 1 )  && (strcmp(argv[b_i] , "false") == 1 )){
//             cerr << "Error; Bool parameter " << b_i << " is not either 'true' or 'false'!" << endl;
//             exit(1);
//         }
//     }
    std::cout << ini.GetSection("") << endl;
    
    
    // SIMULATION PARAMS
    _simpar.runs = atoi(ini.GetValue("", "runs", ""));                       // Number of Simulation runs to get mean values from
    _simpar.timestep = atof(ini.GetValue("", "dt", ""));
    _simpar.simtime = atoi(ini.GetValue("", "t", ""));                   // simulation time
    _simpar.instantvalues = 200;


    // MODEL PARAMS
    _modelpar.rodDist = atof(ini.GetValue("", "d", ""));                //Distance of the rods depending on boxsize. for zero do rodDist[]={0.0}
    _modelpar.boxsize = atof(ini.GetValue("", "b", ""));
    _modelpar.particlesize = atof(ini.GetValue("", "p", ""));
    _modelpar.urange = atof(ini.GetValue("", "k", ""));
    _modelpar.ustrength = atof(ini.GetValue("", "u", ""));
	_modelpar.hpi_u = atof(ini.GetValue("", "hpi_u", ""));
	_modelpar.hpi_k = atof(ini.GetValue("", "hpi_k", ""));
	_modelpar.polymersize = atof(ini.GetValue("", "a", ""));   // diameter of polymer chains, i.e. edgeparticles
}
