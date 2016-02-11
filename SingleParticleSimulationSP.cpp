
#include "headers/SingleParticleSimulationSP.h"

using namespace std;

//ISSUES:
//Size of int large enough for large numbers such as steps???
//fopen() fclose() each time that write XYZtraj is called? -> slow
//rounding errors when calculating int steps = int simtime / double timestep

int main(int argc, const char* argv[]){
	// measure runtime of the simulation
	clock_t start, end;

    //Main includes the iteration loops for the simulation

    //NOTE: so far wallcrossings is added for all runs!!! makes it kind of incorrect, since after each run, ppos is reset.
    //NOTE: so far saving Instant Values for each tenth step!

    //TODO struct delete 'bool' add struct name, e.g. _triggers.writeTrajectory
    //TRIGGERS:
    _triggers.writeTrajectory = (strcmp(argv[1] , "true") == 0 ) ;    // relative position TODO
    _triggers.fitRPinv = (strcmp(argv[2] , "true") == 0 ) ;
    _triggers.ranSpheres = (strcmp(argv[3] , "true") == 0 ) ;
    _triggers.recordPosHisto = (strcmp(argv[4] , "true") == 0 ) ;
    _triggers.noLub = (strcmp(argv[5] , "true") == 0 ) ;
    _triggers.includeSteric = (strcmp(argv[6] , "true") == 0 ) ;  // steric 2
	_triggers.ranPot = (strcmp(argv[7] , "true") == 0 ) ;
	_triggers.hpi = (strcmp(argv[8] , "true") == 0 ) ;          // hpi exp
	int boolpar = 8;

	// Checking for correct structure of input arguments
	for (int k= 0; k < argc; k++ ) cout << "parameter " << k << " " << argv[k] << endl;
	for (int b_i=1; b_i<=boolpar; b_i++){
		if ((strcmp(argv[b_i] , "true") == 1 )  && (strcmp(argv[b_i] , "false") == 1 )){
			cerr << "Error; Bool parameter " << b_i << " is not either 'true' or 'false'!" << endl;
			exit(1);
		}
	}

    //TODO struct delete int, double, etc, add _simpar
    _simpar.runs = atoi( argv[boolpar+1] );                       // Number of Simulation runs to get mean values from
    _simpar.timestep = atof( argv[boolpar+2] );
    _simpar.simtime = atoi( argv[boolpar+3] );                   // simulation time
    _simpar.instantvalues = 200;

    _modelpar.rodDist = atof( argv[boolpar+4] );                //Distance of the rods depending on boxsize. for zero do rodDist[]={0.0}
    _modelpar.boxsize = atof( argv[boolpar+5] );
    _modelpar.particlesize = atof( argv[boolpar+6] );
    _modelpar.urange = atof( argv[boolpar+7] );
    _modelpar.ustrength = atof( argv[boolpar+8] );
	_modelpar.hpi_u = atof( argv[boolpar+9] );
	_modelpar.hpi_k = atof( argv[boolpar+10] );
	_modelpar.polymersize = atof( argv[boolpar+11] );   // diameter of polymer chains, i.e. edgeparticles
    int instValIndex;                             //Counter for addInstantValue
    
	double HI = false;

	//HI
	if (_modelpar.polymersize != 0){
		if (fmod(10, _modelpar.polymersize) != 0 && (1 == fmod(10, _modelpar.polymersize)/_modelpar.polymersize)) {  // The second comparison is needed cause fmod()  sometimes does not work properly and gives fmod(a,b) = b, which of course is sensless
			cerr << "Error; bad polymersize! (Nonzero modulus when dividing 10)" << endl;
			exit(1);
		}
		HI = true;
        
        std::string filename = "../rfitInvRP.py"; // The python script should lie in the folder "above" Release, i.e. where the C++ files are, such that it is tracked by git
        std::string command = "python ";
        std::string parameters = ' ' + toString(_modelpar.particlesize) + ' ' + toString(_modelpar.polymersize);
        command += filename + parameters;
        cout << "running: $" + command << endl;
        system(command.c_str());
        sleep(5);         //make the programme wait for 5 seconds, such that fit script can finish
	}

    if (!_triggers.includeSteric) {
        _triggers.noLub = false;
        cout << "!!! WARNING: Steric is disabled.\nActivating Lubrication, since it needs to be activated to replace steric interaction!!";
    }




    //MFP
    double fpInt = _modelpar.boxsize/10;



    _simpar.steps = _simpar.simtime/_simpar.timestep;
    _simpar.saveInt = _simpar.steps/_simpar.instantvalues;
    const int trajout = (int)(10/_simpar.timestep);
    const int MMcalcStep = (int)(0.05/_simpar.timestep);


    //initialize instance of configuration
    CConfiguration conf = CConfiguration(_simpar.timestep, _modelpar, _triggers, _files);
    if (_triggers.recordPosHisto) conf.initPosHisto();

    //Create data folders and print location as string to string "folder"
    //TODO struct createDataFolder does not return a folder anymore
    createDataFolder(conf.getTestCue());
    

    //initialize averages
    CAverage energyU = CAverage("Upot", _files.folder, _simpar.instantvalues, _simpar.runs);
    CAverage squareDisp = CAverage("squaredisp", _files.folder, _simpar.instantvalues, _simpar.runs);
    //CAverage displacem = CAverage("displacement", _files.folder, _simpar.instantvalues, _simpar.runs);
    //CAverage squareDisp_x = CAverage("squaredisp_x", _files.folder, _simpar.instantvalues, _simpar.runs);
    //CAverage squareDisp_y = CAverage("squaredisp_y", _files.folder, _simpar.instantvalues, _simpar.runs);
    //CAverage squareDisp_z = CAverage("squaredisp_z", _files.folder, _simpar.instantvalues, _simpar.runs);
    //CAverage mfp_x = CAverage("mfp_x", _files.folder, 1, 1);
    //CAverage mfp_xyz;
    //if ( recordMFP ) mfp_xyz = CAverage("mfp_xyz", _files.folder, 1, 1);



    //create file to save the trajectory
    string traj_file = _files.folder + "/Coordinates/single_traj.xyz";
    if (_triggers.writeTrajectory) conf.saveXYZTraj(traj_file,0,"w");


    //cout << "Starting Run Number: " << simcounter << " out of " << totalsims << endl;
    start = clock();
    cout << "Starting Simulation!" << endl;

    unsigned int stepcount = 0;
    ofstream trajectoryfile;
    trajectoryfile.open((_files.folder + "/Coordinates/trajectory.txt").c_str());

    // Write parameter file parameters.txt
    parameterFile(conf.getTestCue());

    if (_triggers.includeSteric && conf.testOverlap()){
        cout << "ERROR !!!!!!!!!!!!!!!!\nThere is an OVERLAP between the polymer network and the particle start position!" << endl;
        return 1;
    }

// ************** START OF RUNS-LOOP *****************
    for (int l = 0; l<_simpar.runs; l++){

        conf.updateStartpos();

        instValIndex = 0;
        int fpCounter[3] = {0};                  //counter for first passage times (next point to pass first: fpCounter*fpInt
        int boxcheck = 0; // stores 0 or 1 value for conf.checkBoxCrossing().
        int stepcheck = 0; // " for conf.makeStep()."

        //if (l%100==0) cout << "run " << l << endl;

        for (int i = 0; i < _simpar.steps; i++){  //calculate stochastic force first, then mobility force!!
            // calc HI mobility matrix here, since it needs to be defined for random force normalisation

            if (HI){ // Calculate displacement and update mobilityMatrix if it is larger than 0.1*tracer Radius
                conf.checkDisplacementforMM();
            }
            //    if ( i%MMcalcStep == 0 ){ conf.calcTracerMobilityMatrix(true); }
            //    else { conf.calcTracerMobilityMatrix(false); }
            //}

            conf.calcStochasticForces();


            if (_modelpar.ustrength != 0) conf.calcMobilityForces();


            if (((i+1)%100 == 0) && (l == 0) && _triggers.writeTrajectory){       //Save the first trajectory to file
                conf.saveXYZTraj(traj_file, i, "a");                    // TODO change back ((i+1)%XXX == 0) to 100
            }




            if (((i+1)%_simpar.saveInt) == 0){       //saving Instant Values for each saveInt'th step!
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
            // if (recordMFP){
            //     for (int a=0; a < 3; a++){
            //         if (conf.checkFirstPassage(fpInt*(fpCounter[a]+1), a)) {
            //             fpCounter[a]+=1;
            //             mfp_xyz.addFPValue(i*timestep, fpCounter[a]);
            //         }
            //     }
            // }

            stepcheck = conf.makeStep();    //move particle at the end of iteration


            /* //TODO steric2
            if (includeSteric && conf.testOverlap()) conf.moveBack();
            else boxcheck = conf.checkBoxCrossing();
            */ // end steric2

                //TODO steric
            while (_triggers.includeSteric && conf.testOverlap()){
                conf.moveBack();
                conf.calcStochasticForces();
                stepcheck = conf.makeStep();
            }
            boxcheck = conf.checkBoxCrossing(); //check if particle has crossed the confinement of the box
            // end steric


            if (boxcheck==1 || stepcheck==1) return 1; // If boxcrossing is out of range stop program execution!

            stepcount++;
            if (stepcount%trajout == 0) {
                std::vector<double> ppos = conf.getppos();
                trajectoryfile << fixed << stepcount * _simpar.timestep << "\t" << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;
            }
            if (_triggers.recordPosHisto && ((i % 5) == 0)) conf.addHistoValue();

        }
        if ( _triggers.recordPosHisto ) conf.printHistoMatrix(_files.folder);
        
        if (l%100 == 0)  cout << "run " << toString(l) << endl;  
    }//----------END OF RUNS-LOOP ----------------




    //watch out: Save average instant values at timeInterval: timestep * saveinterval saveInt!!!
    energyU.saveAverageInstantValues(_simpar.saveInt*_simpar.timestep);
    squareDisp.saveAverageInstantValues(_simpar.saveInt*_simpar.timestep);
    //displacem.saveAverageInstantValues(saveInt*timestep);
    //squareDisp_x.saveAverageInstantValues(saveInt*timestep);
    //squareDisp_y.saveAverageInstantValues(saveInt*timestep);
    //squareDisp_z.saveAverageInstantValues(saveInt*timestep);
    //if (recordMFP) mfp_x.saveAverageFPValue(fpInt);
    //if (recordMFP) mfp_xyz.saveAverageFPValue(fpInt);

    //if ( _triggers.recordPosHisto ) conf.printHistoMatrix(_files.folder);




	cout << "Simulation Finished" << endl;
	end = clock();
	double runtime = (double)((end-start)/(CLOCKS_PER_SEC));
	cout << runtime << " seconds runtime." << endl;
    
    //TODO struct parametereFileAppend does not take folder as argument anymore
    parameterFileAppend(runtime);

	trajectoryfile.close();


    return 0;
}