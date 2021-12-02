#include "PHG4TpcDetector.h"
#include "PHG4TpcDefs.h"
#include "PHG4TpcDisplayAction.h"

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phparameter/PHParameters.h>

#include <phool/recoConsts.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <cmath>
#include <iostream>  // for basic_ostream::operator<<
#include <map>       // for map
#include <sstream>

class G4VSolid;
class PHCompositeNode;

//_______________________________________________________________
PHG4TpcDetector::PHG4TpcDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4TpcDisplayAction *>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_InnerCageRadius(m_Params->get_double_param("gas_inner_radius") * cm - m_Params->get_double_param("cage_layer_9_thickness") * cm - m_Params->get_double_param("cage_layer_8_thickness") * cm - m_Params->get_double_param("cage_layer_7_thickness") * cm - m_Params->get_double_param("cage_layer_6_thickness") * cm - m_Params->get_double_param("cage_layer_5_thickness") * cm - m_Params->get_double_param("cage_layer_4_thickness") * cm - m_Params->get_double_param("cage_layer_3_thickness") * cm - m_Params->get_double_param("cage_layer_2_thickness") * cm - m_Params->get_double_param("cage_layer_1_thickness") * cm)
  , m_OuterCageRadius(m_Params->get_double_param("gas_outer_radius") * cm + m_Params->get_double_param("cage_layer_9_thickness") * cm + m_Params->get_double_param("cage_layer_8_thickness") * cm + m_Params->get_double_param("cage_layer_7_thickness") * cm + m_Params->get_double_param("cage_layer_6_thickness") * cm + m_Params->get_double_param("cage_layer_5_thickness") * cm + m_Params->get_double_param("cage_layer_4_thickness") * cm + m_Params->get_double_param("cage_layer_3_thickness") * cm + m_Params->get_double_param("cage_layer_2_thickness") * cm + m_Params->get_double_param("cage_layer_1_thickness") * cm)
{
}

//_______________________________________________________________
int PHG4TpcDetector::IsInTpc(G4VPhysicalVolume *volume) const
{
  if (m_ActiveFlag)
  {
    if (m_ActiveVolumeSet.find(volume) != m_ActiveVolumeSet.end())
    {
      return 1;
    }
  }
  if (m_AbsorberActiveFlag)
  {
    if (m_AbsorberVolumeSet.find(volume) != m_AbsorberVolumeSet.end())
    {
      return -1;
    }
  }
  return 0;
}

//_______________________________________________________________
void PHG4TpcDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  // create Tpc envelope
  // tpc consists of (from inside to gas volume, outside is reversed up to now)
  // 1st layer cu
  // 2nd layer FR4
  // 3rd layer HoneyComb
  // 4th layer cu
  // 5th layer FR4
  // 6th layer Kapton
  // 7th layer cu
  // 8th layer Kapton
  // 9th layer cu

  double steplimits = m_Params->get_double_param("steplimits") * cm;
  if (std::isfinite(steplimits))
  {
    m_G4UserLimits = new G4UserLimits(steplimits);
  }

  G4VSolid *tpc_envelope = new G4Tubs("tpc_envelope", m_InnerCageRadius, m_OuterCageRadius, m_Params->get_double_param("tpc_length") * cm / 2., 0., 2 * M_PI);

  recoConsts *rc = recoConsts::instance();
  G4LogicalVolume *tpc_envelope_logic = new G4LogicalVolume(tpc_envelope,
                                                            GetDetectorMaterial(rc->get_StringFlag("WorldMaterial")),
                                                            "tpc_envelope");
  m_DisplayAction->AddVolume(tpc_envelope_logic, "TpcEnvelope");

  ConstructTpcCageVolume(tpc_envelope_logic);
  ConstructTpcGasVolume(tpc_envelope_logic);

  new G4PVPlacement(0, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm),
                    tpc_envelope_logic, "tpc_envelope",
                    logicWorld, 0, false, OverlapCheck());
}

int PHG4TpcDetector::ConstructTpcGasVolume(G4LogicalVolume *tpc_envelope)
{
  static std::map<int, std::string> tpcgasvolname =
      {{PHG4TpcDefs::North, "tpc_gas_north"},
       {PHG4TpcDefs::South, "tpc_gas_south"}};

  // Window / central membrane
  double tpc_window_thickness = m_Params->get_double_param("window_thickness") * cm;
  double tpc_half_length = (m_Params->get_double_param("tpc_length") * cm - tpc_window_thickness) / 2.;

  G4VSolid *tpc_window = new G4Tubs("tpc_window", m_Params->get_double_param("gas_inner_radius") * cm, m_Params->get_double_param("gas_outer_radius") * cm, tpc_window_thickness / 2., 0., 2 * M_PI);
  G4LogicalVolume *tpc_window_logic = new G4LogicalVolume(tpc_window,
                                                          GetDetectorMaterial(m_Params->get_string_param("window_surface1_material")),
                                                          "tpc_window");

  //  G4VisAttributes *visatt = new G4VisAttributes();
  //  visatt->SetVisibility(true);
  //  visatt->SetForceSolid(true);
  //  visatt->SetColor(PHG4TPCColorDefs::tpc_cu_color);
  //  tpc_window_logic->SetVisAttributes(visatt);

  m_DisplayAction->AddVolume(tpc_window_logic, "TpcWindow");
  G4VPhysicalVolume *tpc_window_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                         tpc_window_logic, "tpc_window",
                                                         tpc_envelope, false, PHG4TpcDefs::Window, OverlapCheck());

  m_AbsorberVolumeSet.insert(tpc_window_phys);

  // Window / central membrane core
  double tpc_window_surface1_thickness = m_Params->get_double_param("window_surface1_thickness") * cm;
  double tpc_window_surface2_thickness = m_Params->get_double_param("window_surface2_thickness") * cm;
  double tpc_window_surface2_core_thickness = tpc_window_thickness - 2 * tpc_window_surface1_thickness;
  double tpc_window_core_thickness = tpc_window_surface2_core_thickness - 2 * (tpc_window_surface2_thickness);

  G4VSolid *tpc_window_surface2_core =
      new G4Tubs("tpc_window_surface2_core", m_Params->get_double_param("gas_inner_radius") * cm, m_Params->get_double_param("gas_outer_radius") * cm,
                 tpc_window_surface2_core_thickness / 2., 0., 2 * M_PI);
  G4LogicalVolume *tpc_window_surface2_core_logic = new G4LogicalVolume(tpc_window_surface2_core,
                                                                        GetDetectorMaterial(m_Params->get_string_param("window_surface2_material")),
                                                                        "tpc_window_surface2_core");

  m_DisplayAction->AddVolume(tpc_window_surface2_core_logic, "TpcWindow");
  //  visatt = new G4VisAttributes();
  //  visatt->SetVisibility(true);
  //  visatt->SetForceSolid(true);
  //  visatt->SetColor(PHG4TPCColorDefs::tpc_pcb_color);
  //  tpc_window_surface2_core_logic->SetVisAttributes(visatt);
  G4VPhysicalVolume *tpc_window_surface2_core_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                                       tpc_window_surface2_core_logic, "tpc_window_surface2_core",
                                                                       tpc_window_logic, false, PHG4TpcDefs::WindowCore, OverlapCheck());

  m_AbsorberVolumeSet.insert(tpc_window_surface2_core_phys);

  G4VSolid *tpc_window_core =
      new G4Tubs("tpc_window", m_Params->get_double_param("gas_inner_radius") * cm, m_Params->get_double_param("gas_outer_radius") * cm,
                 tpc_window_core_thickness / 2., 0., 2 * M_PI);
  G4LogicalVolume *tpc_window_core_logic = new G4LogicalVolume(tpc_window_core,
                                                               GetDetectorMaterial(m_Params->get_string_param("window_core_material")),
                                                               "tpc_window_core");

  m_DisplayAction->AddVolume(tpc_window_core_logic, "TpcHoneyComb");
  G4VPhysicalVolume *tpc_window_core_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                              tpc_window_core_logic, "tpc_window_core",
                                                              tpc_window_surface2_core_logic, false, PHG4TpcDefs::WindowCore, OverlapCheck());

  m_AbsorberVolumeSet.insert(tpc_window_core_phys);

  // Gas
  G4VSolid *tpc_gas = new G4Tubs("tpc_gas", m_Params->get_double_param("gas_inner_radius") * cm, m_Params->get_double_param("gas_outer_radius") * cm, tpc_half_length / 2., 0., 2 * M_PI);

  G4LogicalVolume *tpc_gas_logic = new G4LogicalVolume(tpc_gas,
                                                       GetDetectorMaterial(m_Params->get_string_param("tpc_gas")),
                                                       "tpc_gas");

  tpc_gas_logic->SetUserLimits(m_G4UserLimits);
  m_DisplayAction->AddVolume(tpc_gas_logic, "TpcGas");
  G4VPhysicalVolume *tpc_gas_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, (tpc_half_length + tpc_window_thickness) / 2.),
                                                      tpc_gas_logic, tpcgasvolname[PHG4TpcDefs::North],
                                                      tpc_envelope, false, PHG4TpcDefs::North, OverlapCheck());
  std::cout << "north copy no: " << tpc_gas_phys->GetCopyNo() << std::endl;

  m_ActiveVolumeSet.insert(tpc_gas_phys);
  tpc_gas_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -(tpc_half_length + tpc_window_thickness) / 2.),
                                   tpc_gas_logic, tpcgasvolname[PHG4TpcDefs::South],
                                   tpc_envelope, false, PHG4TpcDefs::South, OverlapCheck());

  std::cout << "south copy no: " << tpc_gas_phys->GetCopyNo() << std::endl;
  m_ActiveVolumeSet.insert(tpc_gas_phys);

#if G4VERSION_NUMBER >= 1033
  const G4RegionStore *theRegionStore = G4RegionStore::GetInstance();
  G4Region *tpcregion = theRegionStore->GetRegion("REGION_TPCGAS");
  tpc_gas_logic->SetRegion(tpcregion);
  tpcregion->AddRootLogicalVolume(tpc_gas_logic);
#endif
  return 0;
}

int PHG4TpcDetector::ConstructTpcCageVolume(G4LogicalVolume *tpc_envelope)
{
  // 8th layer cu
  // 9th layer FR4
  // 10th layer HoneyComb
  // 11th layer FR4
  // 12th layer Kapton
  // 13th layer FR4
  static const int nlayers = 9;
  static const double thickness[nlayers] = {m_Params->get_double_param("cage_layer_1_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_2_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_3_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_4_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_5_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_6_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_7_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_8_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_9_thickness") * cm};

  static const std::string material[nlayers] = {m_Params->get_string_param("cage_layer_1_material"),
                                                m_Params->get_string_param("cage_layer_2_material"),
                                                m_Params->get_string_param("cage_layer_3_material"),
                                                m_Params->get_string_param("cage_layer_4_material"),
                                                m_Params->get_string_param("cage_layer_5_material"),
                                                m_Params->get_string_param("cage_layer_6_material"),
                                                m_Params->get_string_param("cage_layer_7_material"),
                                                m_Params->get_string_param("cage_layer_8_material"),
                                                m_Params->get_string_param("cage_layer_9_material")};

  double tpc_cage_radius = m_InnerCageRadius;
  std::ostringstream name;
  for (int i = 0; i < nlayers; i++)
  {
    name.str("");
    int layerno = i + 1;
    name << "tpc_cage_layer_" << layerno;
    G4VSolid *tpc_cage_layer = new G4Tubs(name.str(), tpc_cage_radius, tpc_cage_radius + thickness[i], m_Params->get_double_param("tpc_length") * cm / 2., 0., 2 * M_PI);
    G4LogicalVolume *tpc_cage_layer_logic = new G4LogicalVolume(tpc_cage_layer,
                                                                GetDetectorMaterial(material[i]),
                                                                name.str());
    m_DisplayAction->AddTpcInnerLayer(tpc_cage_layer_logic);
    G4VPhysicalVolume *tpc_cage_layer_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                               tpc_cage_layer_logic, name.str(),
                                                               tpc_envelope, false, layerno, OverlapCheck());
    m_AbsorberVolumeSet.insert(tpc_cage_layer_phys);
    tpc_cage_radius += thickness[i];
  }
  // outer cage
  tpc_cage_radius = m_OuterCageRadius;
  for (int i = 0; i < nlayers; i++)
  {
    tpc_cage_radius -= thickness[i];
    name.str("");
    int layerno = 10 + 1 + i;  // so the accompanying inner layer is layer - 10
    name << "tpc_cage_layer_" << layerno;
    G4VSolid *tpc_cage_layer = new G4Tubs(name.str(), tpc_cage_radius, tpc_cage_radius + thickness[i], m_Params->get_double_param("tpc_length") * cm / 2., 0., 2 * M_PI);
    G4LogicalVolume *tpc_cage_layer_logic = new G4LogicalVolume(tpc_cage_layer,
                                                                GetDetectorMaterial(material[i]),
                                                                name.str());
    m_DisplayAction->AddTpcOuterLayer(tpc_cage_layer_logic);
    G4VPhysicalVolume *tpc_cage_layer_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                               tpc_cage_layer_logic, name.str(),
                                                               tpc_envelope, false, layerno, OverlapCheck());
    m_AbsorberVolumeSet.insert(tpc_cage_layer_phys);
  }

  return 0;
}
