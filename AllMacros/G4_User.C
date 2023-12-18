#ifndef MACRO_G4USER_C
#define MACRO_G4USER_C

#include <fun4all/Fun4AllServer.h>
#include <G4_ReadParams.C>

R__LOAD_LIBRARY(libfun4all.so)

class PHG4Reco;

namespace Enable
{
// if you want this to run by default, initialize this to true
// Otherwise you have to use Enable::USER = true; in your macro
  bool USER = false;
  int USER_VERBOSITY = 0;
}

namespace G4USER
{
// here you can set parameters in your macro via
// G4USER::myparam = 1;
// add as many as you need
  int myparam = 0;
  std::string settingName = "HOLY_GRAIL";
  std::map<signed int, std::string> PDGDictionary;
  map<string, double>* settingParams = new map<string, double>();
  map<string, double>* designParams = new map<string, double>();
}

// This initializes the G4 part if you have a detector implemented
// You need to tell its dimensions to the surrounding black hole
void UserInit(string settingsFile, string designFile)
{
  // set the black hole dimensions surrounding the detector
  // XXX: maximum radius of your detector
  // YYY: maximum extension in z
  // ZZZ: maximum extension in -z (use -ZZZ)
  //BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, XXX);
  //BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, YYY);
  //BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, ZZZ);
  G4USER::PDGDictionary.insert({11, "e-"});
  G4USER::PDGDictionary.insert({-11, "e+"});
  G4USER::PDGDictionary.insert({22, "gamma"});
  G4USER::PDGDictionary.insert({211, "pi+"});
  G4USER::PDGDictionary.insert({-211, "pi-"});
  ReadParams(settingsFile, G4USER::settingParams);
  ReadParams(designFile, G4USER::designParams);
  
  std::cout << "USER INTITALIZED DESIGN" << std::endl;
  std::cout << "Set Magnetic Field scale factor is = " << (*G4USER::designParams)[Form("B_FIELD")] << std::endl; 
  std::cout << "Number of BARREL layers = " << (*G4USER::settingParams)[Form("NLAYERS_SI_BAR")] << std::endl;
  std::cout << "Number of Electron Side Disks = " << (*G4USER::settingParams)[Form("NLAYERS_SI_EDISK")] << std::endl;
  std::cout << "Number of Hadron Side Disks = " << (*G4USER::settingParams)[Form("NLAYERS_SI_HDISK")] << std::endl;
  std::cout << "Maximum radius of your detector = " << (*G4USER::settingParams)[Form("MAX_RADIUS")] << std::endl;
  std::cout << "Maximum extension in z = " << (*G4USER::settingParams)[Form("MAX_HZ")] << std::endl;
  std::cout << "Maximum extension in -z = " << (*G4USER::settingParams)[Form("MAX_EZ")] << std::endl;
  std::cout << "Average support thickness = " << (*G4USER::settingParams)[Form("AVG_SUPRT_THICKNESS")] << std::endl;
  std::cout << "Maximum Vertex support radius = " << (*G4USER::settingParams)[Form("MAX_VTX_SUPRT_RADIUS")] << std::endl;

  
}

// If you have a detector - here goes the setup
void UserDetector(PHG4Reco *g4Reco)
{
  return;
}

// Any analysis goes here (registering your module with Fun4All)
void UserAnalysisInit()
{
  Fun4AllServer* se = Fun4AllServer::instance();

  return;
}

#endif
