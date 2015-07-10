#include "PHG4InnerHcalSubsystem.h"
#include "PHG4InnerHcalDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4InnerHcalSteppingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <fun4all/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4InnerHcalSubsystem::PHG4InnerHcalSubsystem( const std::string &name, const int lyr ):
  PHG4Subsystem( name ),
  detector_( 0 ),
  steppingAction_( NULL ),
  eventAction_(NULL),
  place_in_x(0),
  place_in_y(0),
  place_in_z(0),
  rot_in_x(0),
  rot_in_y(0),
  rot_in_z(0),
  material("SS310"),  // default - stainless steel ss310
  active(0),
  absorberactive(0),
  layer(lyr),
  blackhole(0),
  tilt_angle(NAN),
  detector_type(name),
  superdetector("NONE")
{

  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str().c_str());
}

//_______________________________________________________________________
int PHG4InnerHcalSubsystem::Init( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  detector_ = new PHG4InnerHcalDetector(topNode, Name(), layer);
  detector_->SetPlace(place_in_x, place_in_y, place_in_z);
  detector_->SetXRot(rot_in_x);
  detector_->SetYRot(rot_in_y);
  detector_->SetZRot(rot_in_z);
  detector_->SetActive(active);
  detector_->SetAbsorberActive(absorberactive);
  detector_->BlackHole(blackhole);
  detector_->SetMaterial(material);
  detector_->SuperDetector(superdetector);
  detector_->OverlapCheck(overlapcheck);
  if (isfinite(tilt_angle)) detector_->SetTilt(tilt_angle);
  if (active)
    {
      ostringstream nodename;
      if (superdetector != "NONE")
	{
	  nodename <<  "G4HIT_" << superdetector;
	}
      else
	{
	  nodename <<  "G4HIT_" << detector_type << "_" << layer;
	}
      // create hit list
      PHG4HitContainer* block_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str());
      if ( !block_hits )
	{

	  dstNode->addNode( new PHIODataNode<PHObject>( block_hits = new PHG4HitContainer(), nodename.str().c_str(), "PHObject" ));

	}
      if (absorberactive)
	{
	  nodename.str("");
	  if (superdetector != "NONE")
	    {
	      nodename <<  "G4HIT_ABSORBER_" << superdetector;
	    }
	  else
	    {
	      nodename <<  "G4HIT_ABSORBER_" << detector_type << "_" << layer;
	    }
	}
      block_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str());
      if ( !block_hits )
	{

	  dstNode->addNode( new PHIODataNode<PHObject>( new PHG4HitContainer(), nodename.str().c_str(), "PHObject" ));

	}
      // create stepping action
      steppingAction_ = new PHG4InnerHcalSteppingAction(detector_);

      eventAction_ = new PHG4EventActionClearZeroEdep(topNode, nodename.str());
    }
  if (blackhole && !active)
    {
      steppingAction_ = new PHG4InnerHcalSteppingAction(detector_);
    }
  return 0;

}

//_______________________________________________________________________
int
PHG4InnerHcalSubsystem::process_event( PHCompositeNode * topNode )
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
    if (steppingAction_)
        {
            steppingAction_->SetInterfacePointers( topNode );
        }
    return 0;
}


//_______________________________________________________________________
PHG4Detector* PHG4InnerHcalSubsystem::GetDetector( void ) const
{
    return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4InnerHcalSubsystem::GetSteppingAction( void ) const
{
    return steppingAction_;
}

void
PHG4InnerHcalSubsystem::SetTiltAngle(const double tilt)
{
  tilt_angle = tilt*deg;
}
