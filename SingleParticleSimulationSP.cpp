
#include "headers/SingleParticleSimulationSP.h"

using namespace std;
#define ifdebug(x)

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

    //TRIGGERS:
    _triggers.writeTrajectory = (strcmp(argv[1] , "true") == 0 ) ;    // relative position TODO
    _triggers.fitRPinv = (strcmp(argv[2] , "true") == 0 ) ;
    _triggers.ranSpheres = (strcmp(argv[3] , "ranSpheres") == 0  || strcmp(argv[3] , "trueRan") == 0) ;
    _triggers.trueRan = (strcmp(argv[3] , "trueRan") == 0 ) ;
    _triggers.ranRod = (strcmp(argv[3] , "ranRod") == 0 ) ;
    _triggers.recordPosHisto = (strcmp(argv[4] , "true") == 0 ) ;
    _triggers.noLub = (strcmp(argv[5] , "true") == 0 ) ;
    _triggers.includeSteric = (strcmp(argv[6] , "steric") == 0 || strcmp(argv[6] , "steric2") == 0 || strcmp(argv[6] , "LJ") == 0 || strcmp(argv[6] , "LJ025") == 0) ;
    _triggers.stericType = argv[6];
    _triggers.ranPot = (strcmp(argv[7] , "true") == 0 ) ;
    _triggers.HI2 = (strcmp(argv[8] , "true") == 0 ) ;          // hpi exp
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

    _modelpar.rodDist = 0;                //Deprecated -- can use argv[boolpar+4] for sth else 
    _modelpar.boxsize = atof( argv[boolpar+5] );
    _modelpar.particlesize = atof( argv[boolpar+6] );
    _modelpar.urange = atof( argv[boolpar+7] );
    _modelpar.ustrength = atof( argv[boolpar+8] );
    _modelpar.hpi_u = atof( argv[boolpar+9] );
    _modelpar.hpi_k = atof( argv[boolpar+10] );
    _modelpar.polymersize = atof( argv[boolpar+11] );   // diameter of polymer chains, i.e. edgeparticles
    _modelpar.n_cells = atoi( argv[boolpar+12] );   // number of cells along one axis of the simulation box
    _modelpar.EwaldTest = atoi( argv[boolpar+13] );  // index to set the number of polymer particles. Default is EwaldTest = 0!
    _modelpar.nmax = atoi( argv[boolpar+14] ); // This corresponds to the r_cutoff = _nmax * _boxsize
    _modelpar.lubcutint = atoi( argv[boolpar+15] );

    //default values! atoi returns zero, if no valid number is given (which is good).
    if (_modelpar.lubcutint == 0) _modelpar.lubcutint = 9;
    if (_modelpar.nmax == 0) _modelpar.nmax = 3;
    if (_modelpar.n_cells == 0) _modelpar.n_cells = 1;
    if (_modelpar.boxsize == 0) _modelpar.boxsize = 10*_modelpar.n_cells;

    // first argument must be the name of the parameter file
    // read_run_params(("../jobs/input/" + toString(argv[1]) + ".txt").c_str());


    int instValIndex;                             //Counter for addInstantValue
    const string sterictype = _triggers.stericType;

    bool HI = false;

    //HI
    if (_modelpar.polymersize != 0){
        HI = true;

        if (_triggers.fitRPinv == true){
            std::string filename = "../rfitInvRP.py"; // The python script should lie in the folder "above" Release, i.e. where the C++ files are, such that it is tracked by git
            std::string command = "python ";
            std::string parameters = ' ' + toString(_modelpar.particlesize) + ' ' + toString(_modelpar.polymersize);
            command += filename + parameters;
            cout << "running: $" + command << endl;
            system(command.c_str());
            sleep(3);         //make the programme wait for 5 seconds, such that fit script can finish
        }
    }
    //TODO LJ vs steric
    if (_triggers.stericType!="steric" && _triggers.stericType!="steric2" && _triggers.noLub) {
        cout << "!!! Error: Lubrication is disabled.\nActivating steric interaction is vital!!";
        return 4;
    }
    
    
    if (_triggers.ranRod && _modelpar.n_cells !=1){
        cout << "Error: _triggers.ranRod && _modelpar.n_cells != 1" << endl;
        return 4;
    }
    if (_triggers.ranRod && _modelpar.ustrength !=0){
            cout << "Error: _triggers.ranRod && _modelpar.ustrength !=0.\ncalcmobility forces not yet implemented" << endl;
            return 4;
    }
    





    //MFP
    double fpInt = _modelpar.boxsize/10;



    _simpar.steps = _simpar.simtime/_simpar.timestep;
    _simpar.saveInt = _simpar.steps/_simpar.instantvalues;
    const int trajout = (int)(10/_simpar.timestep);
    const int MMcalcStep = (int)(0.05/_simpar.timestep);


    //initialize instance of configuration
    CConfiguration conf;
    try{
        conf = CConfiguration(_simpar.timestep, _modelpar, _triggers, _files);
    }
    catch(int e){
        cout << "An exception occurred. Exception Nr. " << e << '\n';
        return 1;
    }
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


    if (_triggers.ranSpheres && conf.testOverlap()){
        cout << "Resetting particle position to random allowed place in box." << endl;
        int oltest = conf.resetposition(); //assigns a random allowed position to tracer particle in box.
        if (oltest == 1){
            return 1;
        }
    }

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
            ifdebug(if (i%50==0){cout << i;};)
            // calc HI mobility matrix here, since it needs to be defined for random force normalisation

            if (HI){ // Calculate displacement and update mobilityMatrix if it is larger than 0.1*tracer Radius
                conf.checkDisplacementforMM();
            }

            conf.calcStochasticForces();

            if (_modelpar.EwaldTest==0) conf.calcMobilityForces();
            else conf.calcMobForcesBeads();
            
            

            if (((i+1)%100 == 0) && (l == 0) && _triggers.writeTrajectory){       //Save the first trajectory to file
                conf.saveXYZTraj(traj_file, i, "a");                    // TODO change back ((i+1)%XXX == 0) to 100
            }


            if (((i+1)%_simpar.saveInt) == 0){       //saving Instant Values for each saveInt'th step!
                energyU.addInstantValue(conf.getUpot(), instValIndex);
                squareDisp.addInstantValue(conf.getPosVariance(), instValIndex);

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


            //TODO steric2
            if (sterictype=="steric2"){
                if (_triggers.includeSteric && conf.testOverlap()) conf.moveBack();
                else boxcheck = conf.checkBoxCrossing();
            }// end steric2
            

            //TODO steric
            else if (sterictype=="steric"){
                int cnt=0;
                while (_triggers.includeSteric && conf.testOverlap()){
                    conf.moveBack();
                    conf.calcStochasticForces();
                    stepcheck = conf.makeStep();
                    ifdebug ((cout << "moveBack!");) //conf.moveBackReport();)
                    cnt++;
                    if (cnt==3000){
                    cout << "Bad particle position. Cannot avoid overlap with moveBack." << endl;
                    conf.saveXYZTraj(("tmptrajMakeStepErr.xyz"), 0, "w");
                    return 4;
                    }
                }
                boxcheck = conf.checkBoxCrossing(); //check if particle has crossed the confinement of the box
            }// end steric
            else{// case of LJ potential for steric
                boxcheck = conf.checkBoxCrossing();
            }


            if (boxcheck==1 || stepcheck==1){
                cout << "Error: boxcheck==1 or stepcheck==1." << endl;
                return 1; // If boxcrossing is out of range stop program execution!
            }

            stepcount++;
            if (stepcount%trajout == 0) {
                std::vector<double> ppos = conf.getppos();
                trajectoryfile << fixed << stepcount * _simpar.timestep << "\t" << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;
            }
            if (_triggers.recordPosHisto && ((i % 5) == 0)) conf.addHistoValue();
            
            //TODO testing
//             conf._ppos(0) = (i+1)*0.3*_modelpar.boxsize + _modelpar.boxsize/2;
//             conf._ppos(1) = _modelpar.boxsize/2; conf._ppos(2) = _modelpar.boxsize/2;
//             if (i==3) return 0;

        }
        if ( (l%20 == 0) && _triggers.recordPosHisto ) conf.printHistoMatrix(_files.folder);

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
