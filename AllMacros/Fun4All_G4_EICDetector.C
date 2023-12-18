#ifndef MACRO_FUN4ALLG4EICDETECTOR_C
#define MACRO_FUN4ALLG4EICDETECTOR_C

#include <GlobalVariables.C>

#include <DisplayOn.C>
#include <G4Setup_EICDetector.C>
#include <G4_DSTReader_EICDetector.C>
#include <G4_EventEvaluator.C>
#include <G4_FwdJets.C>
#include <G4_Global.C>
#include <G4_Input.C>
#include <G4_Production.C>
#include <G4_User.C>

#include <TROOT.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <RooUnblindPrecision.h>

R__LOAD_LIBRARY(libfun4all.so)

int Fun4All_G4_EICDetector(
    const string &inputParamsFile = "params.config",
    const string &settingsFile = "settings.config",
    const int nEvents = -1,
    const double etamin = -3.4,
    const double etamax = 3.4,
    const double mommin = 1.0,
    const double mommax = 20.0,
    const string &outputFile = "G4EICDetector.root",
    bool debug = true)
{
  //---------------
  // Fun4All server
  //---------------
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);
  Enable::DEV_DEBUG = debug;
  //Opt to print all random seed used for debugging reproducibility. Comment out to reduce stdout prints.
  //PHRandomSeed::Verbosity(1);

  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();
  // By default every random number generator uses
  // PHRandomSeed() which reads /dev/urandom to get its seed
  // if the RANDOMSEED flag is set its value is taken as initial seed
  // which will produce identical results so you can debug your code
  // rc->set_IntFlag("RANDOMSEED", 12345);

  // User initialization
  Enable::USER = true;
  UserInit(settingsFile, inputParamsFile);

  //===============
  // Input options
  //===============

  // switching IPs by comment/uncommenting the following lines
  // used for both beamline setting and for the event generator crossing boost
  Enable::IP6 = true;
  //Enable::IP8 = true;

   
  //===============
  // The following Ion energy and electron energy setting needs to be speficied
  // The setting options for e-p high divergence setting (p energy x e energy):
  // Option: 275x18, 275x10, 100x10, 100x5, 41x5
  //
  // The setting options for e-p high divergence setting (p energy x e energy):
  // Option: 275x18, 275x10, 100x10, 100x5, 41x5
  //
  // The setting options for e-p high divergence setting (A energy x e energy):
  // Option: 110x18, 110x10, 110x5, 41x5

  // Setting proton beam pipe energy. If you don't know what to set here, leave it at 275
  Enable::HFARFWD_ION_ENERGY = 275;

  // Setting electron beam pipe energy. If you don't know what to set here, leave it at 18
  Enable::HFARBWD_E_ENERGY = 18;

  // Beam Scattering configuration setting specified by CDR
  //
  // Option 1: ep-high-acceptance
  // Option 2: ep-high-divergence
  // Option 3: eA
  //
  // Enable::BEAM_COLLISION_SETTING = "ep-high-divergence";
  // If you don't know what to put here, set it to ep-high-divergence   
  //
  // Enable::BEAM_COLLISION_SETTING = "eA";
  Enable::BEAM_COLLISION_SETTING = "ep-high-divergence";

  // Simple multi particle generator in eta/phi/pt ranges
  Input::SIMPLE = true;
  // Input::SIMPLE_NUMBER = 2; // if you need 2 of them
  // Input::SIMPLE_VERBOSITY = 1;

  // Particle gun (same particles in always the same direction)
  // Input::GUN = true;
  // Input::GUN_NUMBER = 3; // if you need 3 of them
  // Input::GUN_VERBOSITY = 0;

  // Particle ion gun
  // Input::IONGUN = true; 

  // Upsilon generator
  // Input::UPSILON = true;
  // Input::UPSILON_NUMBER = 3; // if you need 3 of them
  // Input::UPSILON_VERBOSITY = 0;

  //-----------------
  // Initialize the selected Input/Event generation
  //-----------------
  InputInit();
  //--------------
  // Set generator specific options
  //--------------
  // can only be set after InputInit() is called

  // Simple Input generator:
  // if you run more than one of these Input::SIMPLE_NUMBER > 1
  // add the settings for other with [1], next with [2]...
  if (Input::SIMPLE)
  {
    int PDG = (int)(*G4USER::settingParams)["PDG"];
    int N = (int)(*G4USER::settingParams)["N_PER_EVENT"];
    INPUTGENERATOR::SimpleEventGenerator[0]->add_particles(G4USER::PDGDictionary[PDG], N);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                                                              PHG4SimpleEventGenerator::Uniform,
                                                                              PHG4SimpleEventGenerator::Uniform);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., (int)(*G4USER::settingParams)["VERTEX_DISTRIBUTION_WIDTH"]);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(etamin, etamax);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_p_range(mommin, mommax);
  }
  
  // particle gun
  // if you run more than one of these Input::GUN_NUMBER > 1
  // add the settings for other with [1], next with [2]...
  if (Input::GUN)
  {
    INPUTGENERATOR::Gun[0]->AddParticle("pi-", 0, 1, 0);
    INPUTGENERATOR::Gun[0]->set_vtx(0, 0, 0);
  }

  // register all input generators with Fun4All
  InputRegister();

  // set up production relatedstuff
  //   Enable::PRODUCTION = true;

  //======================
  // Write the DST
  //======================

  // Enable::DSTOUT = true;
  DstOut::OutputDir = ".";
  DstOut::OutputFile = outputFile;
  Enable::DSTOUT_COMPRESS = true;  // Compress DST files

  //Option to convert DST to human command readable TTree for quick poke around the outputs
  // Enable::DSTREADER = true;

  // turn the display on (default off)
  Enable::DISPLAY = false;

  //======================
  // What to run
  //======================
  // Global options (enabled for all subsystems - if implemented)
  //  Enable::ABSORBER = true;
  //  Enable::OVERLAPCHECK = true;
  //  Enable::VERBOSITY = 1;

  // whether to simulate the Be section of the beam pipe
  Enable::PIPE = true;
  // If need to disable EIC beam pipe extension beyond the Be-section:
  //G4PIPE::use_forward_pipes = true;
  
  // barrel tracker
  Enable::TrackingService = true;
  // Enable::TrackingService_VERBOSITY = INT_MAX - 10;
  Enable::BARREL = true;
  // fst
  Enable::FST = true;

  Enable::TRACKING = true;
  Enable::TRACKING_EVAL = Enable::TRACKING && true;
  G4TRACKING::DISPLACED_VERTEX = true;  // this option exclude vertex in the track fitting and use RAVE to reconstruct primary and 2ndary vertexes
                                        // projections to calorimeters
  
  Enable::PLUGDOOR = false;

  // Other options
  Enable::GLOBAL_RECO = G4TRACKING::DISPLACED_VERTEX;  // use reco vertex for global event vertex
  Enable::GLOBAL_FASTSIM = true;

  // new settings using Enable namespace in GlobalVariables.C
  Enable::BLACKHOLE = true;
  //Enable::BLACKHOLE_SAVEHITS = false; // turn off saving of bh hits
  //BlackHoleGeometry::visible = true;

  //

  //---------------
  // World Settings
  //---------------
  //  G4WORLD::PhysicsList = "FTFP_BERT"; //FTFP_BERT_HP best for calo
  //  G4WORLD::WorldMaterial = "G4_AIR"; // set to G4_GALACTIC for material scans
  //  G4WORLD::WorldMaterial = "G4_Galactic"; // set to G4_GALACTIC for material scans

  //---------------
  // Magnet Settings
  //---------------

  //  const string magfield = "1.5"; // alternatively to specify a constant magnetic field, give a float number, which will be translated to solenoidal field in T, if string use as fieldmap name (including path)
  //  G4MAGNET::magfield = string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root");  // default map from the calibration database
  double mag_field = (double)(*G4USER::designParams)["B_FIELD"];
  G4MAGNET::magfield_rescale = -mag_field / 1.5;  // make consistent with expected Babar field strength of 1.4T

  //---------------
  // Pythia Decayer
  //---------------
  // list of decay types in
  // $OFFLINE_MAIN/include/g4decayer/EDecayType.hh
  // default is All:
  // G4P6DECAYER::decayType = EDecayType::kAll;

  // Initialize the selected subsystems
  G4Init();

  //---------------------
  // GEANT4 Detector description
  //---------------------
  G4Setup();
  //--------------
  // Tracking and PID
  //--------------

  if (Enable::TRACKING) Tracking_Reco();

  //-----------------
  // Global Vertexing
  //-----------------

  if (Enable::GLOBAL_RECO)
  {
    Global_Reco();
  }
  else if (Enable::GLOBAL_FASTSIM)
  {
    Global_FastSim();
  }

  string outputroot = outputFile;
  string remove_this = ".root";
  size_t pos = outputroot.find(remove_this);
  if (pos != string::npos)
  {
    outputroot.erase(pos, remove_this.length());
  }

  if (Enable::DSTREADER) G4DSTreader_EICDetector(outputroot + "_DSTReader.root");

  //----------------------
  // Simulation evaluation
  //----------------------

  if (Enable::TRACKING_EVAL) Tracking_Eval(outputroot + "_g4tracking_eval.root");

  if (Enable::USER) UserAnalysisInit();

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();

  //--------------
  // Set up Output Manager
  //--------------
  if (Enable::PRODUCTION)
  {
    Production_CreateOutputDir();
  }

  if (Enable::DSTOUT)
  {
    string FullOutFile = DstOut::OutputDir + "/" + DstOut::OutputFile;
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
    if (Enable::DSTOUT_COMPRESS) DstCompress(out);
    se->registerOutputManager(out);
  }

  //-----------------
  // Event processing
  //-----------------
  if (Enable::DISPLAY)
  {
    DisplayOn();

    gROOT->ProcessLine("Fun4AllServer *se = Fun4AllServer::instance();");
    gROOT->ProcessLine("PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");");

    cout << "-------------------------------------------------" << endl;
    cout << "You are in event display mode. Run one event with" << endl;
    cout << "se->run(1)" << endl;
    cout << "Run Geant4 command with following examples" << endl;
    gROOT->ProcessLine("displaycmd()");

    return 0;
  }
  // if we use a negative number of events we go back to the command line here
  if (nEvents < 0)
  {
    //DisplayOn();
    //gROOT->ProcessLine("PHG4Reco *g4 = QTGui();");
    QTGui();
    gSystem->Exit(0);
    return 0;
  }
  // if we run any of the particle generators and use 0 it'll run forever
  if (nEvents == 0 && !Input::READHITS && !Input::HEPMC && !Input::READEIC)
  {
    cout << "using 0 for number of events is a bad idea when using particle generators" << endl;
    cout << "it will run forever, so I just return without running anything" << endl;
    return 0;
  }

  se->run(nEvents);

  //-----
  // Exit
  //-----

  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  if (Enable::PRODUCTION)
  {
    Production_MoveOutput();
  }
  gSystem->Exit(0);
  return 0;
}
#endif
