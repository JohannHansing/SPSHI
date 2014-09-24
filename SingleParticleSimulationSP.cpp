
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <string.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sstream>
#include <math.h>
#include <boost/filesystem.hpp>
#include "headers/CAverage.h"
#include "headers/CConfiguration.h"


using namespace std;

//ISSUES:
//Size of int large enough for large numbers such as steps???
//fopen() fclose() each time that write XYZtraj is called? -> slow
//rounding errors when calculating int steps = int simtime / double timestep



//Function declarations
string createDataFolder(
        bool respos, double dt, double simtime, double potRange, double potStrength, double boxsize,
        double particlesize, double rDist, bool potMod, bool steric, bool randomPot, bool hpi, double hpi_u, double hpi_k, double polymersize);
void settingsFile(
        string folder, bool respos, double particlesize, double boxsize, double timestep, double runs, double steps,
        double potStrength, double potRange, double rDist, bool potMod, bool recordMFP, bool steric, bool randomPot, bool hpi, double hpi_u,
		double hpi_k, double polymersize, double executiontime);
void printReport(bool resetPos, int entry, int opposite, int sides, const double timestep[], const double urange[], const double ustrength[], 
        const double rodDist[], const double particlesize[], unsigned int runs, int tsize, int rsize, int ssize, int dsize, int psize, bool potMod, 
        bool steric, bool randomPot, bool hpi, double hpi_u, double hpi_k, double polymersize);

template<typename T>
string toString(const T& value);

template <typename T, size_t N>
inline
size_t sizeOfArray( const T(&)[ N ] );





int main(int argc, const char* argv[]){
	// measure runtime of the simulation
	clock_t start, end;
	
    //Main includes the iteration loops for the simulation

    //NOTE: so far wallcrossings is added for all runs!!! makes it kind of incorrect, since after each run, ppos is reset.
    //NOTE: so far saving Instant Values for each tenth step!

    //TRIGGERS:
    bool writeTrajectory = (strcmp(argv[1] , "true") == 0 ) ;    // relative position TODO
    bool resetPos = (strcmp(argv[2] , "true") == 0 ) ;
    bool potentialMod = (strcmp(argv[3] , "true") == 0 ) ;       //BESSEL TODO
    bool recordMFP = (strcmp(argv[4] , "true") == 0 ) ;
    bool recordPosHisto = (strcmp(argv[5] , "true") == 0 ) ;
    bool includeSteric = (strcmp(argv[6] , "true") == 0 ) ;  // steric 2
	bool ranPot = (strcmp(argv[7] , "true") == 0 ) ;
	bool hpi = (strcmp(argv[8] , "true") == 0 ) ;          // hpi exp
	int boolpar = 8;
	
	// Checking for correct structure of input arguments
	for (int k= 0; k < argc; k++ ) cout << "parameter " << k << " " << argv[k] << endl;
	for (int b_i=1; b_i<=boolpar; b_i++){
		if ((strcmp(argv[b_i] , "true") == 1 )  && (strcmp(argv[b_i] , "false") == 1 )){
			cerr << "Error; Bool parameter " << b_i << " is not either 'true' or 'false'!" << endl;
			exit(1);
		}
	}


    int runs = atoi( argv[boolpar+1] );                       // Number of Simulation runs to get mean values from
    double timestep = atof( argv[boolpar+2] );
    int simtime = atoi( argv[boolpar+3] );                   // simulation time
    int instantvalues = 200;
    unsigned int steps;
    
    double rodDist = atof( argv[boolpar+4] );                //Distance of the rods depending on boxsize. for zero do rodDist[]={0.0}
    double boxsize = atof( argv[boolpar+5] );
    double particlesize = atof( argv[boolpar+6] );
    double urange = atof( argv[boolpar+7] );
    double ustrength = atof( argv[boolpar+8] );
	double hpi_u = atof( argv[boolpar+9] );
	double hpi_k = atof( argv[boolpar+10] );
	double polymersize = atof( argv[boolpar+11] );   // diameter of polymer chains, i.e. edgeparticles
    unsigned int saveInt;
    int instValIndex;                             //Counter for addInstantValue
	double HI = false;

	//HI
	if (polymersize != 0){
		if (fmod(10, polymersize) != 0) {
			cerr << "Error; bad polymersize! (Nonzero modulus when dividing 10)" << endl;
			exit(1);
		}
		HI = true;
		includeSteric = true;
	}
	
	
	
	
    //MFP
    double fpInt = boxsize/10;

    //printReport(resetPos, conf.getwallcrossings(0), conf.getwallcrossings(1), conf.getwallcrossings(2), timestep, urange, ustrength, rodDist, particlesize, runs,
    //            sizeOfArray(timestep), sizeOfArray(urange), sizeOfArray(ustrength), sizeOfArray(rodDist), sizeOfArray(particlesize), potentialMod, includeSteric);


    
    steps = simtime/timestep;
    saveInt = steps/instantvalues;
	const int trajout = (int)(1/timestep);
	const int MMcalcStep = (int)(0.05/timestep);
        
    //Create data folders and print location as string to string "folder"
    string folder = createDataFolder(resetPos, timestep, simtime, urange, ustrength, boxsize, particlesize, rodDist, potentialMod, 
	                                 includeSteric, ranPot, hpi, hpi_u, hpi_k, polymersize);


    //initialize averages
    CAverage energyU = CAverage("Upot", folder, instantvalues, runs);
    CAverage squareDisp = CAverage("squaredisp", folder, instantvalues, runs);
    //CAverage displacem = CAverage("displacement", folder, instantvalues, runs);
    //CAverage squareDisp_x = CAverage("squaredisp_x", folder, instantvalues, runs);
    //CAverage squareDisp_y = CAverage("squaredisp_y", folder, instantvalues, runs);
    //CAverage squareDisp_z = CAverage("squaredisp_z", folder, instantvalues, runs);
    //CAverage mfp_x = CAverage("mfp_x", folder, 1, 1);
    CAverage mfp_xyz;
    if ( recordMFP ) mfp_xyz = CAverage("mfp_xyz", folder, 1, 1);

    
    //initialize instance of configuration
    CConfiguration conf = CConfiguration(timestep, urange, ustrength, boxsize, rodDist, potentialMod, particlesize, recordPosHisto, includeSteric,
    ranPot, hpi , hpi_u, hpi_k, polymersize);
		

    //create file to save the trajectory
    string traj_file = folder + "/Coordinates/single_traj.xyz";
    if (writeTrajectory) conf.saveXYZTraj(traj_file,0,"w");


    //cout << "Starting Run Number: " << simcounter << " out of " << totalsims << endl;
	start = clock();
    cout << "Starting Simulation!" << endl;
    
	unsigned int stepcount = 0;
    ofstream trajectoryfile;
    trajectoryfile.open((folder + "/Coordinates/trajectory.txt").c_str());
    
// **************START OF RUNS-LOOP
    for (int l = 0; l<runs; l++){

        if (resetPos) conf.resetposition();         
        else conf.updateStartpos();

        instValIndex = 0;
        int fpCounter[3] = {0};                  //counter for first passage times (next point to pass first: fpCounter*fpInt


        if (l%100 == 0) cout << "run " << l << endl;

        for (int i = 0; i < steps; i++){  //calculate stochastic force first, then mobility force!!						
		    // calc HI mobility matrix here, since it needs to be defined for random force normalisation
		    if (HI){
				if ( i%MMcalcStep == 0 ){ conf.calcTracerMobilityMatrix(true); }
				else { conf.calcTracerMobilityMatrix(false); }
			}

            conf.calcStochasticForces();


            conf.calcMobilityForces();


            if (((i+1)%100 == 0) && (l == 0) && writeTrajectory){       //Save the first trajectory to file
                conf.saveXYZTraj(traj_file, i, "a");                    // TODO change back ((i+1)%XXX == 0) to 100
            }




            if (((i+1)%saveInt) == 0){       //saving Instant Values for each saveInt'th step!
                energyU.addInstantValue(conf.getUpot(), instValIndex);
                squareDisp.addInstantValue(conf.getPosVariance(), instValIndex);
            //    displacem.addInstantValue(conf.getDisplacement(), instValIndex);
            //    squareDisp_x.addInstantValue(conf.get1DPosVariance(0), instValIndex);
            //    squareDisp_y.addInstantValue(conf.get1DPosVariance(1), instValIndex);
            //    squareDisp_z.addInstantValue(conf.get1DPosVariance(2), instValIndex);

                instValIndex += 1;
            }

            /*if ((conf.checkFirstPassage(fpInt*(fpCounter+1))) && recordMFP) {
                fpCounter+=1;
                mfp_x.addFPValue(i*timestep, fpCounter);
            }
            */
            if (recordMFP){
                for (int a=0; a < 3; a++){
                    if (conf.checkFirstPassage(fpInt*(fpCounter[a]+1), a)) {
                        fpCounter[a]+=1;
                        mfp_xyz.addFPValue(i*timestep, fpCounter[a]);
                    }
                }
            }
			
			stepcount++;
            conf.makeStep();    //move particle at the end of iteration

            
            //if (includeSteric && conf.testOverlap()) conf.moveBack();   //TODO steric2
            //else conf.checkBoxCrossing();
            

                //TODO steric
            while (includeSteric && conf.testOverlap()){
                conf.moveBack();
                conf.calcStochasticForces();
                conf.makeStep();
            }
            conf.checkBoxCrossing(); //check if particle has crossed the confinement of the box
            


            if (stepcount%trajout == 0) {
				std::vector<double> ppos = conf.getppos();
				trajectoryfile << stepcount * timestep << "\t" << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;
			}
            if (((i % 5) == 0) && recordPosHisto) conf.addHistoValue();
        }
    }//----------END OF RUNS-LOOP




    //watch out: Save average instant values at timeInterval: timestep * saveinterval saveInt!!!
    energyU.saveAverageInstantValues(saveInt*timestep);
    squareDisp.saveAverageInstantValues(saveInt*timestep);
    //displacem.saveAverageInstantValues(saveInt*timestep);
    //squareDisp_x.saveAverageInstantValues(saveInt*timestep);
    //squareDisp_y.saveAverageInstantValues(saveInt*timestep);
    //squareDisp_z.saveAverageInstantValues(saveInt*timestep);
    //if (recordMFP) mfp_x.saveAverageFPValue(fpInt);
    if (recordMFP) mfp_xyz.saveAverageFPValue(fpInt);

    if ( recordPosHisto ) conf.printHistoMatrix(folder);



    //printReport(resetPos, conf.getwallcrossings(0), conf.getwallcrossings(1), conf.getwallcrossings(2), timestep, urange, ustrength, rodDist, particlesize, runs,
    //        sizeOfArray(timestep), sizeOfArray(urange), sizeOfArray(ustrength), sizeOfArray(rodDist), sizeOfArray(particlesize), potentialMod, includeSteric);

    
	cout << "Simulation Finished" << endl;
	end = clock();
	double runtime = (double)((end-start)/(CLOCKS_PER_SEC));
	cout << runtime << " seconds runtime." << endl;
	
	//If settingsFile is saved, then the simulation was successfull
    settingsFile(folder, resetPos, particlesize, boxsize, timestep, runs, steps, ustrength, urange, rodDist, 
	potentialMod, recordMFP, includeSteric, ranPot ,hpi, hpi_u, hpi_k, polymersize, (runtime/3600));
	
	trajectoryfile.close();
	
	
    return 0;
}



//--------------------------------------------------------------------------
//**************************************************************************
//--------------------------------------------------------------------------






string createDataFolder(bool resetpos, double timestep, double simtime, double potRange, double potStrength,
        double boxsize, double particlesize, double rDist, bool potMod, bool steric, bool randomPot, bool hpi, double hpi_u, double hpi_k,
		double polymersize){
    //NOTE: Maybe I can leave out dt, as soon as I settled on a timestep
    //NOTE: As soon as I create input-list with variables, I must change this function
    char range[5];
    sprintf(range, "%.3f", potRange);
    //In the definition of folder, the addition has to START WITH A STRING! for the compiler to know what to do (left to right).
    string folder = "sim_data";
    if (!resetpos) folder = folder + "/noreset";
    if (randomPot) folder = folder + "/ranPot";
    if (steric) folder = folder + "/steric";    //TODO steric2
    if (potMod) folder = folder +  "/potMod";   //"/potMod";  TODO!!! Bessel
	if (hpi) folder = folder + "/HPI/hpiu" + toString(hpi_u) + "/hpik" + toString(hpi_k);
    folder = folder
            + "/dt" + toString(timestep)
            + "/t" + toString(simtime)
			+ "/a" + toString(polymersize)
            + "/d" + toString(rDist)
            + "/b" + toString(boxsize)
            + "/p" + toString(particlesize)
            + "/k" + range
            + "/u" + toString(potStrength);
    boost::filesystem::create_directories(folder);
    boost::filesystem::create_directory(folder + "/InstantValues");
    boost::filesystem::create_directory(folder + "/Coordinates");
    return folder;
}


void settingsFile(string folder, bool resetpos, double particlesize, double boxsize, double timestep, double runs, double steps, 
    double potStrength, double potRange, double rDist, bool potMod, bool recordMFP, bool steric, bool randomPot, bool hpi, double hpi_u, double hpi_k, 
	double polymersize, double executiontime){
    //Creates a file where the simulation settings are stored
    //MAYBE ALSO INCLUDE TIME AND DATE!!
    ofstream settingsfile;
    settingsfile.open((folder + "/sim_Settings.txt").c_str());
    settingsfile << "Time required for execution: "
    << executiontime 
    << " hours." << "\n\n";
    settingsfile << "Sim dir: " << folder << endl;
    settingsfile << "ResetPos " << resetpos << endl;
    settingsfile << "potMod " << potMod << endl;//" (Bessel)" << endl;  //TODO Bessel!
    settingsfile << "recordMFP " << recordMFP << endl;
    settingsfile << "includesteric " << steric << endl;
	settingsfile << "ranPot " << randomPot  << endl;
	settingsfile << "HPI " << hpi  << endl;
	if (hpi == true){
		settingsfile << "hpi_u " << randomPot  << endl;
		settingsfile << "hpi_k " << randomPot  << endl;
    }
    settingsfile << "d " << rDist << endl;
    settingsfile << "p " << particlesize << endl;
	settingsfile << "b " << boxsize << endl;
    settingsfile << "dt " << timestep  << endl << "runs " << runs << endl << "steps " << steps << endl << "time: " << timestep*steps << endl;
    settingsfile << "k " << potRange << endl << "U_0 " << potStrength << endl;
	settingsfile << "a" << polymersize << endl;
	    
    settingsfile.close();
}




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
