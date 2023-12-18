#ifndef MACRO_G4FSTEIC_C
#define MACRO_G4FSTEIC_C

#include "GlobalVariables.C"
#include "G4_User.C"
#include "G4_Pipe_EIC.C"

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4detectors/PHG4SectorSubsystem.h>
#include <g4trackfastsim/PHG4TrackFastSim.h>

#include <g4main/PHG4Reco.h>
#include <TMath.h>
#include <string>

R__LOAD_LIBRARY(libg4detectors.so)

int make_LANL_FST_station(string name, PHG4Reco *g4Reco, double zpos, double Rmin,
                          double Rmax, double tSilicon, double pitch);
int make_supportCyl(string name, PHG4Reco *g4Reco,
                    double r, double t, double length);
//-----------------------------------------------------------------------------------//
namespace Enable
{
  static bool FST = false;
  bool FST_OVERLAPCHECK = false;
}  // namespace Enable

namespace G4FST
{
  namespace SETTING
  {
    bool SUPPORTCYL = false;
  }  // namespace SETTING
  namespace SERVICE
  {
  double G4_Al = 15.0; //um
  double G4_KAPTON = 20.0; //um
  double G4_WATER = 100.0; //um
  double G4_GRAPHITE = 50.0; //um
  double G4_AIR = 1.0; //cm
  } //namespace SERVICE
}  // namespace G4FST

//-----------------------------------------------------------------------------------//
void FST_Init()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 48.);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 127.);
  if (G4FST::SETTING::SUPPORTCYL)
  {
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -127.);
  }
}
//-----------------------------------------------------------------------------------//

/* 

The idea is to calculate all the rmax rmin based on z and theta1 and theta2 in a python script (obj_fun.py)
Once calculated just pass the disk dimensions into this script

*/

double computeRMax(double zpos, char Side)
{
  double rMax = 0;
  const int nLayers = (int)(*G4USER::settingParams)[Form("NLAYERS_SI_BAR")];
  if (fabs(zpos) > fabs((*G4USER::designParams)[Form("SI_L%i_%c_LENGTH", nLayers, Side)])) 
  {
    return (*G4USER::settingParams)[Form("MAX_RADIUS")];
  }
  cout << "Average support thickness is " << (*G4USER::settingParams)["AVG_SUPRT_THICKNESS"] << endl;
  for (unsigned int j = 1; j < nLayers; j++)
  {
    double length1 = (*G4USER::designParams)[Form("SI_L%i_%c_LENGTH", j, Side)];
    double length2 = (*G4USER::designParams)[Form("SI_L%i_%c_LENGTH", j + 1, Side)];
    if (fabs(zpos) >= length1 && fabs(zpos) <= length2)
    {
      // compute the angle
      double R2 = (*G4USER::designParams)[Form("SI_L%i_RADIUS", j + 1)] - 
                  (*G4USER::designParams)[Form("SI_L%i_THICKNESS", j + 1)]*Units::um -
                  (*G4USER::settingParams)["AVG_SUPRT_THICKNESS"];
      cout << "Actual radius is " << (*G4USER::designParams)[Form("SI_L%i_RADIUS", j)] << endl;
      double R1 = (*G4USER::designParams)[Form("SI_L%i_RADIUS", j)] - 
                  (*G4USER::designParams)[Form("SI_L%i_THICKNESS", j)]*Units::um -
                  (*G4USER::settingParams)["AVG_SUPRT_THICKNESS"];
      double Z2 = (*G4USER::designParams)[Form("SI_L%i_%c_LENGTH", j + 1, Side)];
      double Z1 = (*G4USER::designParams)[Form("SI_L%i_%c_LENGTH", j, Side)];
      double theta = atan2(fabs(fabs(R2) - fabs(R1)), fabs(fabs(Z2) - fabs(Z1)));
      double r0 = R1 - tan(theta)*Z1; // c = y1 - tan(Theta)*x1
      rMax = fabs(zpos - Z1) * tan(theta) + R1 - (*G4USER::settingParams)["AVG_SUPRT_THICKNESS"];
      cout << "j is : " << j << ", pos is : " << zpos << ", rMax is : " << rMax << ", R2 is : " << R2 << ", R1 is : " << R1 << ", Z1 is : " << Z1 << ", Z2 is : " << Z2 << ", theta is : " << theta << endl;
      break;
    }
  }
  return rMax;
}

void FSTSetup(PHG4Reco *g4Reco)
{
  const int nEDisks = (int)(*G4USER::settingParams)[Form("NLAYERS_SI_EDISK")];
  const int nHDisks = (int)(*G4USER::settingParams)[Form("NLAYERS_SI_HDISK")];
  const int nLayers = (int)(*G4USER::settingParams)[Form("NLAYERS_SI_BAR")];
  const double min_Radius = (*G4USER::settingParams)[Form("MIN_RADIUS")];
  const double max_Radius = (*G4USER::settingParams)[Form("MAX_RADIUS")];
  const double min_Hz = (*G4USER::settingParams)[Form("MIN_HZ")];
  const double min_Ez = (*G4USER::settingParams)[Form("MIN_EZ")]; // note this is negative
  const double max_Hz = (*G4USER::settingParams)[Form("MAX_HZ")];
  const double max_Ez = (*G4USER::settingParams)[Form("MAX_EZ")]; // note this is negative
  
  cout << "#------------- Silicon Disks -------------#" << endl;
  
  for (unsigned int i = 1; i < nEDisks + 1; i++){
    double z = (*G4USER::designParams)[Form("SI_ED%i_Z", i)]; //SI_ED1_Z
    double rMin = G4PIPE::be_pipe_radius + G4PIPE::be_pipe_thickness;
    if(fabs(z) > fabs(G4PIPE::be_pipe_length_neg)) rMin = 0.0521*fabs(z) + G4PIPE::be_pipe_thickness;
    double rMax = computeRMax(fabs(z), 'E');
    double thickness = (*G4USER::designParams)[Form("SI_ED%i_THICKNESS", i)];
    double pitch = (*G4USER::designParams)[Form("SI_ED%i_PITCH", i)]*Units::um;
    cout << "Design Params for EST_" << i << " : ZPos = " << z << 
    ", RMin = " << rMin << ", RMax = " << rMax << ", Thickness = " << 
    thickness << ", Pitch = " << pitch << endl;
    make_LANL_FST_station(Form("EST_%i", i), g4Reco, z, rMin, rMax, thickness, pitch);
  }
  

  for (unsigned int i = 1; i < nHDisks + 1; i++){
    double z = (*G4USER::designParams)[Form("SI_HD%i_Z", i)];
    double rMin = G4PIPE::be_pipe_radius + G4PIPE::be_pipe_thickness;
    if(fabs(z) > G4PIPE::be_pipe_length_plus) rMin = 0.0521*fabs(z) + G4PIPE::be_pipe_thickness;
    double rMax = computeRMax(fabs(z), 'H');
    double thickness = (*G4USER::designParams)[Form("SI_HD%i_THICKNESS", i)];
    double pitch = (*G4USER::designParams)[Form("SI_HD%i_PITCH", i)]*Units::um;
    cout << "Design Params for FST" << i << " : ZPos = " << z << 
    ", RMin = " << rMin << ", RMax = " << rMax << ", Thickness = " << 
    thickness << ", Pitch = " << pitch << endl;
    make_LANL_FST_station(Form("FST_%i", i), g4Reco, z, rMin, rMax, thickness, pitch);
  }

  cout << "\n<><><><><><><> End of EST/FST Information <><><><><><><>\n" << endl;
}

/*---------------------------------------------------------------------*
 * Barrel tracker designed by LANL EIC team                            *
 * See technical notes for details: arXiv:2009.02888                   *
 * Contact Ping and Xuan @LANL for questions:                          *
 *   Xuan: xuanli@lanl.gov                                             *
 *   Ping: cpwong@lanl.gov                                             *
 *---------------------------------------------------------------------*/

//-----------------------------------------------------------------------------------//
int make_LANL_FST_station(string name, PHG4Reco *g4Reco,
                          double zpos, double Rmin, double Rmax, double tSilicon, double pitch)  //silicon thickness
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::FST_OVERLAPCHECK;
  Rmin = Rmin * Units::cm;
  Rmax = Rmax * Units::cm;
  zpos = zpos * Units::cm;
  double overall_thickness = tSilicon * Units::um  + G4FST::SERVICE::G4_Al * Units::um 
                             + G4FST::SERVICE::G4_KAPTON * Units::um  + G4FST::SERVICE::G4_WATER * Units::um 
                             + G4FST::SERVICE::G4_GRAPHITE * Units::um  + G4FST::SERVICE::G4_AIR * Units::cm
                             + G4FST::SERVICE::G4_GRAPHITE * Units::um ;
  //cout << Units::um  << endl;
  //cout << tSilicon << " " << G4FST::SERVICE::G4_Al << " " << G4FST::SERVICE::G4_KAPTON << " " << G4FST::SERVICE::G4_WATER << " " << G4FST::SERVICE::G4_GRAPHITE << " " << G4FST::SERVICE::G4_AIR << " " << G4FST::SERVICE::G4_GRAPHITE << endl;
  //cout << tSilicon * Units::um  << " " << G4FST::SERVICE::G4_Al * Units::um  << " " << G4FST::SERVICE::G4_KAPTON * Units::um  << " " << G4FST::SERVICE::G4_WATER * Units::um  << " " << G4FST::SERVICE::G4_GRAPHITE * Units::um  << " " << G4FST::SERVICE::G4_AIR * Units::cm << " " << G4FST::SERVICE::G4_GRAPHITE * Units::um  << endl;
  //cout << "THickness of Disk is : " << overall_thickness << endl;
  if (zpos < 0) zpos = zpos + overall_thickness;
  else zpos = zpos - overall_thickness;
  //cout << "ZPos of Disk is : " << zpos << endl;
  double min_polar_angle = atan2(Rmin, zpos);
  double max_polar_angle = atan2(Rmax, zpos);

  // always facing the interaction point
  double polar_angle = 0;
  if (zpos < 0)
  {
    zpos = -zpos;
    polar_angle = M_PI;
  }
  if (max_polar_angle < min_polar_angle)
  {
    double t = max_polar_angle;
    max_polar_angle = min_polar_angle;
    min_polar_angle = t;
  }
  PHG4SectorSubsystem *fst;
  fst = new PHG4SectorSubsystem(name);

  fst->SuperDetector(name);

  fst->get_geometry().set_normal_polar_angle(polar_angle);
  fst->get_geometry().set_normal_start(zpos);
  fst->get_geometry().set_min_polar_angle(min_polar_angle);
  fst->get_geometry().set_max_polar_angle(max_polar_angle);
  fst->get_geometry().set_max_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  fst->get_geometry().set_min_polar_edge(PHG4Sector::Sector_Geometry::ConeEdge());
  fst->get_geometry().set_N_Sector(1);
  fst->get_geometry().set_material("G4_AIR");
  fst->OverlapCheck(OverlapCheck);  //true);//overlapcheck);

  // build up layers

  fst->get_geometry().AddLayer("SiliconSensor", "G4_Si", tSilicon * Units::um, true, 100);
  fst->get_geometry().AddLayer("Metalconnection", "G4_Al", G4FST::SERVICE::G4_Al * Units::um, false, 100);
  fst->get_geometry().AddLayer("HDI", "G4_KAPTON", G4FST::SERVICE::G4_KAPTON * Units::um, false, 100);
  fst->get_geometry().AddLayer("Cooling", "G4_WATER", G4FST::SERVICE::G4_WATER * Units::um, false, 100);
  fst->get_geometry().AddLayer("Support", "G4_GRAPHITE", G4FST::SERVICE::G4_GRAPHITE * Units::um, false, 100);
  fst->get_geometry().AddLayer("Support_Gap", "G4_AIR", G4FST::SERVICE::G4_AIR * Units::cm, false, 100);
  fst->get_geometry().AddLayer("Support2", "G4_GRAPHITE", G4FST::SERVICE::G4_GRAPHITE * Units::um, false, 100);

  g4Reco->registerSubsystem(fst);

  if (TRACKING::FastKalmanFilter)
  {
    TRACKING::FastKalmanFilter->add_phg4hits(string("G4HIT_") + name,           //      const std::string& phg4hitsNames,
                                             PHG4TrackFastSim::Vertical_Plane,  //      const DETECTOR_TYPE phg4dettype,
                                             pitch / sqrt(12.),                 //      const float radres,
                                             pitch / sqrt(12.),                 //      const float phires,
                                             50e-4 / sqrt(12.),                 //      const float lonres, *ignored in plane detector*
                                             0.9,                                 //      const float eff,
                                             0);                                //      const float noise
  }
  zpos = 0.;
  return 0;
}

#endif

