#ifndef MACRO_G4BARREL_C
#define MACRO_G4BARREL_C

#include <g4lblvtx/AllSiliconTrackerSubsystem.h>
#include <g4lblvtx/SimpleNtuple.h>
#include <g4main/PHG4Reco.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>
#include "GlobalVariables.C"
#include <G4_User.C>

#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

using namespace std;

namespace Enable
{
  bool BARREL = false;
  bool BARREL_OVERLAPCHECK = false;
}  // namespace Enable
//-----------------------------------------------------------------------------------//
void BarrelInit()
{
}

void BarrelFastKalmanFilterConfig(PHG4TrackFastSim * kalman_filter, int ilay, double radius, double pitch, bool addproj)
{

  // import Kalman filter config (lines 226 to 246 here: https://github.com/eic/g4lblvtx/blob/master/macros/auxiliary_studies/simplified_geometry/Fun4All_G4_simplified_v2.C):

  // add Vertexing Layers
  kalman_filter->add_phg4hits(
      Form("G4HIT_BARR_%d",ilay),  // const std::string& phg4hitsNames,
      PHG4TrackFastSim::Cylinder,
      999.,                      // radial-resolution [cm]
      pitch / 10000. / sqrt(12.),  // azimuthal-resolution [cm]
      pitch / 10000. / sqrt(12.),  // z-resolution [cm]
      0.9,                         // efficiency,
      0                          // noise hits
  );
  kalman_filter->add_cylinder_state(Form("CYLINDER_%d",ilay), radius);
  if(addproj)TRACKING::ProjectionNames.insert(Form("CYLINDER_%d",ilay));
}

void BarrelSetup(PHG4Reco* g4Reco) //
{
  cout << "##-------------BarrelSetup--------------##" << endl;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::BARREL_OVERLAPCHECK;
  PHG4CylinderSubsystem *cyl(nullptr);
  double radius = 0.;
  const int nLayers = (*G4USER::settingParams)[Form("NLAYERS_SI_BAR")];
  cout << "nLayers is : " << nLayers << endl;
  for(unsigned int i = 1; i < nLayers + 1; i++)
  {
    radius = (*G4USER::designParams)[Form("SI_L%i_RADIUS",i)];
    double pitch = (*G4USER::designParams)[Form("SI_L%i_PITCH",i)]*Units::um;
    double e_length = (*G4USER::designParams)[Form("SI_L%i_E_LENGTH",i)];
    double h_length = (*G4USER::designParams)[Form("SI_L%i_H_LENGTH",i)];
    double thickness = (*G4USER::designParams)[Form("SI_L%i_THICKNESS",i)] * 9.37 / 100.;
    double z_length = h_length + e_length;
    double z_pos = h_length - e_length;
    cout << "Layer " << i << "Radius = " << radius 
    << " cms, Pitch = " << pitch << 
    " cms, Thickness = " << thickness << " cms, Z-length = " 
    << z_length << " cms, Z-position = " << z_pos << endl;
    cyl = new PHG4CylinderSubsystem("BARR", i);
    cyl->set_string_param("material", "G4_Si");
    cyl->set_double_param("radius", radius);
    cyl->set_double_param("thickness", thickness); // It XX0 = thickess * 100 / 9.37
    cyl->set_double_param("place_z", z_pos/2.);
    cyl->set_double_param("length", z_length);
    cyl->SetActive();
    cyl->OverlapCheck(OverlapCheck);
    g4Reco->registerSubsystem(cyl);

    BarrelFastKalmanFilterConfig(TRACKING::FastKalmanFilter, i, radius, pitch, false);
    //BarrelFastKalmanFilterConfig(TRACKING::FastKalmanFilterInnerTrack, i, radius, pitch, false);
    //BarrelFastKalmanFilterConfig(TRACKING::FastKalmanFilterSiliconTrack, i, radius, pitch, false);
  }

  return radius;
}

#endif

