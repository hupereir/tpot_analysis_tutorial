#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstOutputManager.h>

#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <g4main/PHG4SimpleEventGenerator.h>
#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <micromegas/MicromegasRawDataDecoder.h>
#include <micromegas/MicromegasRawDataCalibration.h>
#include <micromegas/MicromegasRawDataEvaluation.h>

// local macros
#include <G4Setup_sPHENIX.C>
#include <G4_Global.C>

#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_Reco.C>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libqa_modules.so)

R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)

//____________________________________________________________________
int Fun4All_ReadRawData(
  int nEvents = 500,
  const char* inputFile = "/sphenix/lustre01/sphnxpro/commissioning/TPOT/junk/TPOT_ebdc39_junk-00029863-0000.evt",
  const char* evaluationFile = "MicromegasRawDataEvaluation-00029863-0000.root"
  )
{
  // print inputs
  std::cout << "Fun4All_ReadRawData - nEvents: " << nEvents << std::endl;
  std::cout << "Fun4All_ReadRawData - inputFile: " << inputFile << std::endl;
  std::cout << "Fun4All_ReadRawData - evaluationFile: " << evaluationFile << std::endl;

  // options
  Enable::PIPE = true;
  Enable::MBD = true;
  Enable::MAGNET = true;
  Enable::PLUGDOOR = false;

  // enable all absorbers
  // this is equivalent to the old "absorberactive" flag
  Enable::ABSORBER = true;

  // central tracking
  Enable::MVTX = true;
  Enable::INTT = true;
  Enable::TPC = true;
  Enable::MICROMEGAS = true;
  Enable::BLACKHOLE = true;

  // server
  auto se = Fun4AllServer::instance();
  se->Verbosity(1);

  // make sure to printout random seeds for reproducibility
  PHRandomSeed::Verbosity(1);

  // reco const
  auto rc = recoConsts::instance();
  // rc->set_IntFlag("RANDOMSEED",PHRandomSeed());
  // rc->set_IntFlag("RANDOMSEED",1);

  // Geant4 initialization
  G4Init();

  if( true )
  {
    // raw data evaluation
    auto micromegasRawDataEvaluation = new MicromegasRawDataEvaluation;
    micromegasRawDataEvaluation->Verbosity(1);
    micromegasRawDataEvaluation->set_evaluation_outputfile(evaluationFile);
    se->registerSubsystem( micromegasRawDataEvaluation );
  }

  // input manager
  auto in = new Fun4AllPrdfInputManager;
  in->fileopen(inputFile);
  se->registerInputManager(in);

  // process events
  se->run(nEvents);

  // terminate
  se->End();
  se->PrintTimer();

  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}

// This function is only used to test if we can load this as root6 macro
// without running into unresolved libraries and include files
void RunFFALoadTest() {}
