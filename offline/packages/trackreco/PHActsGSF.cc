#include "PHActsGSF.h"
#include "MakeSourceLinks.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsGsfTrackFittingAlgorithm.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/TrackFitting/BetheHeitlerApprox.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

#include "ActsEvaluator.h"

#include <TDatabasePDG.h>

//____________________________________________________________________________..
PHActsGSF::PHActsGSF(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int PHActsGSF::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "PHActsGSF::InitRun begin" << std::endl;
  }

  if (m_actsEvaluator)
  {
    PHNodeIterator iter(topNode);

    PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

    if (!dstNode)
    {
      std::cerr << "DST node is missing, quitting" << std::endl;
      throw std::runtime_error("Failed to find DST node in PHActsTrkFitter::createNodes");
    }

    PHNodeIterator dstIter(topNode);
    PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));

    if (!svtxNode)
    {
      svtxNode = new PHCompositeNode("SVTX");
      dstNode->addNode(svtxNode);
    }
    m_seedTracks = findNode::getClass<SvtxTrackMap>(topNode, _seed_track_map_name);

    if (!m_seedTracks)
    {
      m_seedTracks = new SvtxTrackMap_v2;

      PHIODataNode<PHObject>* seedNode =
          new PHIODataNode<PHObject>(m_seedTracks, _seed_track_map_name, "PHObject");
      svtxNode->addNode(seedNode);
    }
  }

  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto bha = Acts::makeDefaultBetheHeitlerApprox();
  ActsGsfTrackFittingAlgorithm gsf;
  m_fitCfg.fit = gsf.makeGsfFitterFunction(
      m_tGeometry->geometry().tGeometry,
      m_tGeometry->geometry().magField,
      bha,
      12, 1e-4, 
      MixtureReductionAlgorithm::KLDistance, false, false);

  if (m_actsEvaluator)
  {
    m_evaluator = std::make_unique<ActsEvaluator>(m_evalname);
    m_evaluator->Init(topNode);
    m_evaluator->verbosity(Verbosity());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHActsGSF::process_event(PHCompositeNode* topNode)
{
  auto logLevel = Acts::Logging::FATAL;

  if (m_actsEvaluator)
  {
    m_evaluator->next_event(topNode);
  }

  if (Verbosity() > 4)
  {
    logLevel = Acts::Logging::VERBOSE;
  }

  /// Fill an additional track map if using the acts evaluator
  /// for proto track comparison to fitted track
  if (m_actsEvaluator)
  {
    /// wipe at the beginning of every new fit pass, so that the seeds
    /// are whatever is currently in SvtxTrackMap
    m_seedTracks->clear();
    for (const auto& [key, track] : *m_trackMap)
    {
      m_seedTracks->insert(track);
    }
  }

  auto logger = Acts::getDefaultLogger("PHActsGSF", logLevel);

  for (const auto& [key, track] : *m_trackMap)
  {
    auto pSurface = makePerigee(track);
    if (!pSurface)
    {
      //! If no vertex was assigned to track, just skip it
      continue;
    }
    const auto seed = makeSeed(track, pSurface);

    auto svtxseed = new TrackSeed_v2();
    std::map<TrkrDefs::cluskey, Acts::Vector3> clusterPositions;
    for (auto& cKey : get_cluster_keys(track))
    {
      auto cluster = m_clusterContainer->findCluster(cKey);
      auto globalPosition = m_tGeometry->getGlobalPosition(cKey, cluster);
      clusterPositions.insert(std::make_pair(cKey, globalPosition));
      svtxseed->insert_cluster_key(cKey);
    }
    svtxseed->set_phi(track->get_phi());
    TrackSeedHelper::circleFitByTaubin(svtxseed,clusterPositions, 0, 57);
    TrackSeedHelper::lineFit(svtxseed, clusterPositions, 7, 57);

    ActsTrackFittingAlgorithm::MeasurementContainer measurements;
    TrackSeed* tpcseed = track->get_tpc_seed();
    TrackSeed* silseed = track->get_silicon_seed();

    /// We only fit full sPHENIX tracks
    if (!silseed or !tpcseed)
    {
      continue;
    }

    auto crossing = silseed->get_crossing();
    if (crossing == SHRT_MAX)
    {
      continue;
    }

    /*
    auto sourceLinks = getSourceLinks(tpcseed, measurements, crossing);
    auto silSourceLinks = getSourceLinks(silseed, measurements, crossing);
    */

    // loop over modifiedTransformSet and replace transient elements modified for the previous track with the default transforms
    MakeSourceLinks makeSourceLinks;
    makeSourceLinks.setVerbosity(Verbosity());
    makeSourceLinks.set_pp_mode(m_pp_mode);

    makeSourceLinks.resetTransientTransformMap(
        m_alignmentTransformationMapTransient,
        m_transient_id_set,
        m_tGeometry);

    // TPC source links
    auto sourceLinks = makeSourceLinks.getSourceLinks(
        tpcseed,
        measurements,
        m_clusterContainer,
        m_tGeometry,
        m_globalPositionWrapper,
        m_alignmentTransformationMapTransient,
        m_transient_id_set,
        crossing);

    // silicon source links
    auto silSourceLinks = makeSourceLinks.getSourceLinks(
        silseed,
        measurements,
        m_clusterContainer,
        m_tGeometry,
        m_globalPositionWrapper,
        m_alignmentTransformationMapTransient,
        m_transient_id_set,
        crossing);

    // copy transient map for this track into transient geoContext
    m_transient_geocontext = m_alignmentTransformationMapTransient;

    for (auto& siSL : silSourceLinks)
    {
      sourceLinks.push_back(siSL);
    }

    auto calibptr = std::make_unique<Calibrator>();
    CalibratorAdapter calibrator(*calibptr, measurements);
    auto magcontext = m_tGeometry->geometry().magFieldContext;
    auto calcontext = m_tGeometry->geometry().calibContext;

    auto ppoptions = Acts::PropagatorPlainOptions();

    ActsTrackFittingAlgorithm::GeneralFitterOptions options{
        m_transient_geocontext,
        magcontext,
        calcontext,
        &(*pSurface),
        ppoptions};
    if (Verbosity() > 2)
    {
      std::cout << "calling gsf with position "
                << seed.position(m_transient_geocontext).transpose()
                << " and momentum " << seed.momentum().transpose()
                << std::endl;
    }
    auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
    auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
    ActsTrackFittingAlgorithm::TrackContainer tracks(trackContainer, trackStateContainer);
    auto result = fitTrack(sourceLinks, seed, options, calibrator, tracks);

    if (result.ok())
    {
      updateTrack(result, track, tracks, svtxseed, measurements);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::shared_ptr<Acts::PerigeeSurface> PHActsGSF::makePerigee(SvtxTrack* track) const
{
  SvtxVertex* vertex = m_vertexMap->get(track->get_vertex_id());
  if (!vertex)
  {
    return nullptr;
  }

  Acts::Vector3 vertexpos(vertex->get_x() * Acts::UnitConstants::cm,
                          vertex->get_y() * Acts::UnitConstants::cm,
                          vertex->get_z() * Acts::UnitConstants::cm);

  return Acts::Surface::makeShared<Acts::PerigeeSurface>(
      vertexpos);
}

ActsTrackFittingAlgorithm::TrackParameters PHActsGSF::makeSeed(SvtxTrack* track,
                                                               const std::shared_ptr<Acts::PerigeeSurface>& psurf) const
{
  Acts::Vector4 fourpos(track->get_x() * Acts::UnitConstants::cm,
                        track->get_y() * Acts::UnitConstants::cm,
                        track->get_z() * Acts::UnitConstants::cm,
                        10 * Acts::UnitConstants::ns);

  int charge = track->get_charge();
  Acts::Vector3 momentum(track->get_px(),
                         track->get_py(),
                         track->get_pz());

  ActsTransformations transformer;
  auto cov = transformer.rotateSvtxTrackCovToActs(track);

  return ActsTrackFittingAlgorithm::TrackParameters::create(psurf,
                                                            m_tGeometry->geometry().getGeoContext(),
                                                            fourpos,
                                                            momentum,
                                                            charge / momentum.norm(),
                                                            cov,
                                                            Acts::ParticleHypothesis::electron(),
                                                            1 * Acts::UnitConstants::cm)
      .value();
}

ActsTrackFittingAlgorithm::TrackFitterResult PHActsGSF::fitTrack(
    const std::vector<Acts::SourceLink>& sourceLinks,
    const ActsTrackFittingAlgorithm::TrackParameters& seed,
    const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
    const CalibratorAdapter& calibrator,
    ActsTrackFittingAlgorithm::TrackContainer& tracks)
{
  return (*m_fitCfg.fit)(sourceLinks, seed, options, calibrator, tracks);
}

void PHActsGSF::updateTrack(FitResult& result, SvtxTrack* track,
                            ActsTrackFittingAlgorithm::TrackContainer& tracks,
                            const TrackSeed* seed,
                            const ActsTrackFittingAlgorithm::MeasurementContainer& measurements)
{
  std::vector<Acts::MultiTrajectoryTraits::IndexType> trackTips;
  trackTips.reserve(1);
  auto& outtrack = result.value();
  trackTips.emplace_back(outtrack.tipIndex());
  ActsExamples::Trajectories::IndexedParameters indexedParams;

  indexedParams.emplace(std::pair{outtrack.tipIndex(),
                                  ActsExamples::TrackParameters{outtrack.referenceSurface().getSharedPtr(),
                                                                outtrack.parameters(), outtrack.covariance(), Acts::ParticleHypothesis::electron()}});

  updateSvtxTrack(trackTips, indexedParams, tracks, track);

  if (m_actsEvaluator)
  {
    m_evaluator->evaluateTrackFit(tracks, trackTips, indexedParams, track,
                                  seed, measurements);
  }
}

void PHActsGSF::updateSvtxTrack(std::vector<Acts::MultiTrajectoryTraits::IndexType>& tips,
                                Trajectory::IndexedParameters& paramsMap,
                                ActsTrackFittingAlgorithm::TrackContainer& tracks,
                                SvtxTrack* track)
{
  const auto& mj = tracks.trackStateContainer();
  const auto& tracktip = tips.front();
  const auto& params = paramsMap.find(tracktip)->second;
  const auto trajState =
      Acts::MultiTrajectoryHelpers::trajectoryState(mj, tracktip);

  if (Verbosity() > 1)
  {
    std::cout << "Old track parameters: " << std::endl
              << "   (" << track->get_x()
              << ", " << track->get_y() << ", " << track->get_z()
              << ")" << std::endl
              << "   (" << track->get_px() << ", " << track->get_py()
              << ", " << track->get_pz() << ")" << std::endl;
    std::cout << "New GSF track parameters: " << std::endl
              << "   " << params.position(m_transient_geocontext).transpose()
              << std::endl
              << "   " << params.momentum().transpose()
              << std::endl;
  }

  /// Will create new states
  track->clear_states();

  // create a state at pathlength = 0.0
  // This state holds the track parameters, which will be updated below
  float pathlength = 0.0;
  SvtxTrackState_v1 out(pathlength);
  out.set_x(0.0);
  out.set_y(0.0);
  out.set_z(0.0);
  track->insert_state(&out);

  track->set_x(params.position(m_transient_geocontext)(0) / Acts::UnitConstants::cm);
  track->set_y(params.position(m_transient_geocontext)(1) / Acts::UnitConstants::cm);
  track->set_z(params.position(m_transient_geocontext)(2) / Acts::UnitConstants::cm);

  track->set_px(params.momentum()(0));
  track->set_py(params.momentum()(1));
  track->set_pz(params.momentum()(2));
  track->set_charge(params.charge());
  track->set_chisq(trajState.chi2Sum);
  track->set_ndf(trajState.NDF);

  ActsTransformations transformer;
  transformer.setVerbosity(Verbosity());

  if (params.covariance())
  {
    auto rotatedCov = transformer.rotateActsCovToSvtxTrack(params);
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        track->set_error(i, j, rotatedCov(i, j));
      }
    }
  }

  transformer.fillSvtxTrackStates(mj, tracktip, track, m_transient_geocontext);
}

//____________________________________________________________________________..
std::vector<TrkrDefs::cluskey> PHActsGSF::get_cluster_keys(SvtxTrack* track)
{
  std::vector<TrkrDefs::cluskey> out;
  for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
  {
    if (seed)
    {
      std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
    }
  }

  return out;
}

//____________________________________________________________________________..
int PHActsGSF::End(PHCompositeNode* /*unused*/)
{
  if (m_actsEvaluator)
  {
    m_evaluator->End();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHActsGSF::getNodes(PHCompositeNode* topNode)
{
  // track map
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!m_trackMap)
  {
    std::cout << PHWHERE << " The input track map is not available. Exiting PHActsGSF" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // cluster map
  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterContainer)
  {
    std::cout << PHWHERE << "The input cluster container is not available. Exiting PHActsGSF" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "The input Acts tracking geometry is not available. Exiting PHActsGSF" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // global position wrapper
  m_globalPositionWrapper.loadNodes(topNode);

  // vertex map
  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!m_vertexMap)
  {
    std::cout << PHWHERE << "Vertex map unavailable, exiting PHActsGSF" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_alignmentTransformationMapTransient = findNode::getClass<alignmentTransformationContainer>(topNode, "alignmentTransformationContainerTransient");
  if (!m_alignmentTransformationMapTransient)
  {
    std::cout << PHWHERE << "alignmentTransformationContainerTransient not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
