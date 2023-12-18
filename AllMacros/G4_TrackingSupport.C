#ifndef MACRO_G4TrackingService_C
#define MACRO_G4TrackingService_C

#include <GlobalVariables.C>

#include <g4detectors/PHG4ConeSubsystem.h>
#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4main/PHG4Reco.h>

#include <fun4all/Fun4AllServer.h>

#include <cmath>
#include <vector>

//ECCE Tracking Services
//Should be 5 barrels and 4 cones

using namespace std;

class ServiceProperties
{
 public:
  ServiceProperties();

  explicit ServiceProperties(const string &name,
                             const double &rad_len_copper,
                             const double &rad_len_aluminum,
                             const double &rad_len_water,
                             const double &rad_len_plastic,
                             const double &rad_len_carbon,
                             const double &rad_len_iron,
                             const double &z_south,
                             const double &z_north,
                             const double &r_south,
                             const double &r_north);

  virtual ~ServiceProperties(){};

  const string get_name();
  const double get_rad_len_copper();
  const double get_rad_len_aluminum();
  const double get_rad_len_water();
  const double get_rad_len_plastic();
  const double get_rad_len_carbon();
  const double get_rad_len_iron();
  const double get_z_south();
  const double get_z_north();
  const double get_r_south();
  const double get_r_north();

 private:
  const string m_name = "service";
  const double m_rad_len_copper = 0.0;
  const double m_rad_len_aluminum = 0.0;
  const double m_rad_len_water = 0.0;
  const double m_rad_len_plastic = 0.0;
  const double m_rad_len_carbon = 0.0;
  const double m_rad_len_iron = 0.0;
  const double m_z_south = 0.0;
  const double m_z_north = 0.0;
  const double m_r_south = 0.0;
  const double m_r_north = 0.0;
};
ServiceProperties::ServiceProperties(const string &name,
                                     const double &rad_len_copper,
                                     const double &rad_len_aluminum,
                                     const double &rad_len_water,
                                     const double &rad_len_plastic,
                                     const double &rad_len_carbon,
                                     const double &rad_len_iron,
                                     const double &z_south,
                                     const double &z_north,
                                     const double &r_south,
                                     const double &r_north)
  : m_name(name)
  , m_rad_len_copper(rad_len_copper)
  , m_rad_len_aluminum(rad_len_aluminum)
  , m_rad_len_water(rad_len_water)
  , m_rad_len_plastic(rad_len_plastic)
  , m_rad_len_carbon(rad_len_carbon)
  , m_rad_len_iron(rad_len_iron)
  , m_z_south(z_south)
  , m_z_north(z_north)
  , m_r_south(r_south)
  , m_r_north(r_north)
{
}

const string ServiceProperties::get_name() { return m_name; }
const double ServiceProperties::get_rad_len_copper() { return m_rad_len_copper; }
const double ServiceProperties::get_rad_len_aluminum() { return m_rad_len_aluminum; }
const double ServiceProperties::get_rad_len_water() { return m_rad_len_water; }
const double ServiceProperties::get_rad_len_plastic() { return m_rad_len_plastic; }
const double ServiceProperties::get_rad_len_carbon() { return m_rad_len_carbon; }
const double ServiceProperties::get_rad_len_iron() { return m_rad_len_iron; }
const double ServiceProperties::get_z_south() { return m_z_south; }
const double ServiceProperties::get_z_north() { return m_z_north; }
const double ServiceProperties::get_r_south() { return m_r_south; }
const double ServiceProperties::get_r_north() { return m_r_north; }

namespace Enable
{
  bool TrackingService = false;
  bool TrackingService_ABSORBER = false;
  bool TrackingService_OVERLAPCHECK = false;
  int TrackingService_VERBOSITY = 0;

}  // namespace Enable

vector<double> get_thickness(ServiceProperties *object)
{
  vector<double> thickness = {(object->get_rad_len_copper() / 100) * G4TrackingService::materials[0].second
                             ,(object->get_rad_len_aluminum() / 100) * G4TrackingService::materials[1].second
                             ,(object->get_rad_len_water() / 100) * G4TrackingService::materials[2].second
                             ,(object->get_rad_len_plastic() / 100) * G4TrackingService::materials[3].second
                             ,(object->get_rad_len_carbon() / 100) * G4TrackingService::materials[4].second
                             ,(object->get_rad_len_iron() / 100) * G4TrackingService::materials[5].second};
  return thickness;
}

void TrackingServiceInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, 280.);
  // extends only to -z
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -450.);

}

double TrackingServiceCone(ServiceProperties *object, PHG4Reco *g4Reco, double radius)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::TrackingService_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TrackingService_OVERLAPCHECK;
  int verbosity = max(Enable::VERBOSITY, Enable::TrackingService_VERBOSITY);

  PHG4ConeSubsystem *cone;

  double innerRadiusSouth = object->get_r_south();
  double innerRadiusNorth = object->get_r_north();
  double length = abs(object->get_z_north() - object->get_z_south());
  vector<double> thickness = get_thickness(object);

  for (int i = 0; i < G4TrackingService::nMaterials; ++i)
  {
    if (thickness[i] == 0) continue;
    cone = new PHG4ConeSubsystem(object->get_name(), G4TrackingService::subsysID);
    cone->Verbosity(verbosity);
    cone->SetR1(innerRadiusSouth, innerRadiusSouth + thickness[i]);
    cone->SetR2(innerRadiusNorth, innerRadiusNorth + thickness[i]);
    cone->SetPlaceZ(object->get_z_south() + length / 2 + G4TrackingService::GlobalOffset);
    cone->SetZlength(length / 2);
    cone->SetMaterial(G4TrackingService::materials[i].first);
    cone->SuperDetector("TrackingService");
    if (AbsorberActive) cone->SetActive();
    cone->OverlapCheck(OverlapCheck);
    g4Reco->registerSubsystem(cone);
    ++G4TrackingService::subsysID;
    innerRadiusSouth += thickness[i];
    innerRadiusNorth += thickness[i];
  }
  radius = max(innerRadiusSouth, innerRadiusNorth);

  return radius;
}

double TrackingServiceCylinder(ServiceProperties *object, PHG4Reco *g4Reco, double radius)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::TrackingService_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TrackingService_OVERLAPCHECK;
  int verbosity = max(Enable::VERBOSITY, Enable::TrackingService_VERBOSITY);

  PHG4CylinderSubsystem *cyl;

  double innerRadius = object->get_r_south();
  double length = abs(object->get_z_north() - object->get_z_south());
  vector<double> thickness = get_thickness(object);

  for (int i = 0; i < G4TrackingService::nMaterials; ++i)
  {
    if (thickness[i] == 0) continue;
    cyl = new PHG4CylinderSubsystem(object->get_name(), G4TrackingService::subsysID);
    cyl->Verbosity(verbosity);
    cyl->set_double_param("place_z", object->get_z_south() + length / 2 + G4TrackingService::GlobalOffset);
    cyl->set_double_param("radius", innerRadius);
    cyl->set_double_param("length", length);
    cyl->set_string_param("material", G4TrackingService::materials[i].first);
    cyl->set_double_param("thickness", thickness[i]);
    cyl->SuperDetector("TrackingService");
    if (AbsorberActive) cyl->SetActive();
    cyl->OverlapCheck(OverlapCheck);
    g4Reco->registerSubsystem(cyl);
    ++G4TrackingService::subsysID;
    innerRadius += thickness[i];
  }
  radius = innerRadius;

  return radius;
}

double TrackingService(PHG4Reco *g4Reco, double radius)
{
  /* Support structure thickness*/
  double shellX0 = 100 * G4TrackingService::ShellThickness / G4TrackingService::materials[4].second;
  double CuThickness = 13.;          // 0.18668 cms
  double AlThickness = 0.;           // 0. cms
  double WaterThickness = 0.70;      // 0.25256 cms
  double PlasticThickness = 0.48;    // 0.241488 cms
  double CarbonThickness = shellX0;  // 0.3 cms

  double X_X0_percent = 0.1; // % X0
  
  vector<ServiceProperties *> cylinders, cones;

  ServiceProperties *tmp = new ServiceProperties("tmp", 
                                                  CuThickness, AlThickness, 
                                                  WaterThickness, PlasticThickness, 
                                                  CarbonThickness, 0, 
                                                  0., 0., 
                                                  0., 0.
                                                );
  double avg_thickness_inner = get_thickness(tmp)[0] + 
                                get_thickness(tmp)[1] + 
                                get_thickness(tmp)[2] + 
                                get_thickness(tmp)[3] + 
                                get_thickness(tmp)[4];
  avg_thickness_inner *= Units::mm;
  
  (*G4USER::settingParams)[Form("AVG_SUPRT_THICKNESS")] = avg_thickness_inner;

  double vtx_sprt_rad = (*G4USER::settingParams)["MAX_VTX_SPRT_RADIUS"];
  double vtx_e_length = 0.;
  double vtx_h_length = 0.;

  int tmpLyr = 0;
  double close = 9999.;
  for (int ilyr = 1; ilyr < (*G4USER::settingParams)["NLAYERS_SI_BAR"] + 1; ilyr++)
  {
    if(fabs((*G4USER::designParams)[Form("SI_L%i_RADIUS", ilyr)] - (*G4USER::settingParams)["MAX_VTX_SPRT_RADIUS"]) < close)
    {
      vtx_sprt_rad = (*G4USER::designParams)[Form("SI_L%i_RADIUS", ilyr)] + (*G4USER::designParams)[Form("SI_L%i_THICKNESS", ilyr)] * Units::um ;
      vtx_e_length = -1*(*G4USER::designParams)[Form("SI_L%i_E_LENGTH", ilyr)];
      vtx_h_length = (*G4USER::designParams)[Form("SI_L%i_H_LENGTH", ilyr)];
      tmpLyr = ilyr;
      close = fabs((*G4USER::designParams)[Form("SI_L%i_RADIUS", ilyr)] - 
                    (*G4USER::settingParams)["MAX_VTX_SPRT_RADIUS"]
                  );
    }
  }
  cout << "Vertex Support Radius = " << vtx_sprt_rad << 
  ", Vertex E Length = " << vtx_e_length << 
  ", Vertex H Length = " << vtx_h_length << endl;

  // Vertex Cylindrical Support Structure 
  cylinders.push_back(new ServiceProperties("VtxSuprt_Cyl", 
                                            0, 0,
                                            0, 0,
                                            X_X0_percent, 0, 
                                            vtx_e_length, vtx_h_length, 
                                            vtx_sprt_rad, 0
                                            )
                      );
 
  // vtx radius to tmpLyr + 1 support cone
  double cone_north_rad = vtx_sprt_rad - avg_thickness_inner;
  double cone_south_rad = (*G4USER::designParams)[Form("SI_L%i_RADIUS", tmpLyr + 1)] - avg_thickness_inner;
  double cone_north_ez = vtx_e_length;
  double cone_south_ez = -1.*(*G4USER::designParams)[Form("SI_L%i_E_LENGTH", tmpLyr + 1)];
  cones.push_back(new ServiceProperties(Form("Support_ECone%i", tmpLyr), 
                                        CuThickness, AlThickness, 
                                        WaterThickness, PlasticThickness, 
                                        CarbonThickness, 0, 
                                        cone_south_ez, cone_north_ez, 
                                        cone_south_rad, cone_north_rad
                                        )
                  );
  cout << "Cone support radius from : " << 
  cone_south_rad << " to " << cone_north_rad << ", e length from : " <<
  cone_south_ez << " to " << cone_north_ez << endl;
  
  for (int ilyr = tmpLyr + 1; ilyr < (*G4USER::settingParams)["NLAYERS_SI_BAR"]; ilyr++)
  {
    cone_north_rad = (*G4USER::designParams)[Form("SI_L%i_RADIUS", ilyr)] - avg_thickness_inner;
    cone_south_rad = (*G4USER::designParams)[Form("SI_L%i_RADIUS", ilyr + 1)] - avg_thickness_inner;
    cone_north_ez = -1.*(*G4USER::designParams)[Form("SI_L%i_E_LENGTH", ilyr)];
    cone_south_ez = -1.*(*G4USER::designParams)[Form("SI_L%i_E_LENGTH", ilyr + 1)];
    cones.push_back(new ServiceProperties(Form("Support_ECone%i", ilyr), 
                                          CuThickness, AlThickness, 
                                          WaterThickness, PlasticThickness, 
                                          CarbonThickness, 0, 
                                          cone_south_ez, cone_north_ez, 
                                          cone_south_rad, cone_north_rad
                                          )
                    );
    cout << "Cone support radius from : " <<
    cone_south_rad << " to " << cone_north_rad << ", e length from : " <<
    cone_south_ez << " to " << cone_north_ez << endl;
  }
  // Lets move to the hadron side

  // vtx to tmpLyr + 1 support cone for hadron
  cone_south_rad = vtx_sprt_rad - avg_thickness_inner;
  cone_north_rad = (*G4USER::designParams)[Form("SI_L%i_RADIUS", tmpLyr + 1)] - avg_thickness_inner;
  cone_south_ez = vtx_h_length;
  cone_north_ez = (*G4USER::designParams)[Form("SI_L%i_H_LENGTH", tmpLyr + 1)];
  cones.push_back(new ServiceProperties(Form("Support_HCone%i", tmpLyr), 
                                        CuThickness, AlThickness, 
                                        WaterThickness, PlasticThickness, 
                                        CarbonThickness, 0, 
                                        cone_south_ez, cone_north_ez, 
                                        cone_south_rad, cone_north_rad
                                        )
                  );
  
  cout << "Cone support radius Hadron side from : " <<
  cone_south_rad << " to " << cone_north_rad << ", h length from : " <<
  cone_south_ez << " to " << cone_north_ez << endl;

  for (int ilyr = tmpLyr + 1; ilyr < (*G4USER::settingParams)["NLAYERS_SI_BAR"]; ilyr++)
  {
    cone_south_rad = (*G4USER::designParams)[Form("SI_L%i_RADIUS", ilyr)] - avg_thickness_inner;
    cone_north_rad = (*G4USER::designParams)[Form("SI_L%i_RADIUS", ilyr + 1)] - avg_thickness_inner;
    cone_south_ez = (*G4USER::designParams)[Form("SI_L%i_H_LENGTH", ilyr)];
    cone_north_ez = (*G4USER::designParams)[Form("SI_L%i_H_LENGTH", ilyr + 1)];
    cones.push_back(new ServiceProperties(Form("Support_HCone%i", ilyr), 
                                          CuThickness, AlThickness, 
                                          WaterThickness, PlasticThickness, 
                                          CarbonThickness, 0, 
                                          cone_south_ez, cone_north_ez, 
                                          cone_south_rad, cone_north_rad
                                          )
                    );
    cout << "Cone support radius Hadron side from : " <<
    cone_south_rad << " to " << cone_north_rad << ", h length from : " <<
    cone_south_ez << " to " << cone_north_ez << endl;
  }

  // Finally in a cylinder covering the last cylinder and the entire Z length
  double cyl_north_ez = (*G4USER::settingParams)["MAX_HZ"];
  double cyl_south_ez = (*G4USER::settingParams)["MAX_EZ"];
  int lastLyr = (int)(*G4USER::settingParams)["NLAYERS_SI_BAR"];
  double cyl_rad = (*G4USER::designParams)[Form("SI_L%i_RADIUS", lastLyr)] +
                    (*G4USER::designParams)[Form("SI_L%i_THICKNESS", lastLyr)]* 9.37 / 100.;
  (*G4USER::settingParams)[Form("MAX_RADIUS")] = cyl_rad;
  cylinders.push_back(new ServiceProperties("Support_Cyl_Final", 
                                            CuThickness, AlThickness, 
                                            WaterThickness, PlasticThickness, 
                                            CarbonThickness, 0, 
                                            cyl_south_ez, cyl_north_ez,
                                            cyl_rad, 0
                                            )
                      );
  cout << "Cylinder support radius = " << cyl_rad << ", e length from : " <<
  cyl_south_ez << " to " << cyl_north_ez << endl;

  for (ServiceProperties *cylinder : cylinders) radius += TrackingServiceCylinder(cylinder, g4Reco, radius);
  for (ServiceProperties *cone : cones) radius += TrackingServiceCone(cone, g4Reco, radius);

  return radius;
}

#endif
