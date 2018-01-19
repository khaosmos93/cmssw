//-------------------------------------------------
//
/**  \class L2MuonSeedGeneratorFromL1T
 *
 *   L2 muon seed generator:
 *   Transform the L1 informations in seeds for the
 *   L2 muon reconstruction
 *
 *
 *
 *   \author  A.Everett, R.Bellan, J. Alcaraz
 *
 *    ORCA's author: N. Neumeister
 */
//
//--------------------------------------------------

// Class Header
#include "RecoMuon/L2MuonSeedGenerator/src/L2MuonSeedGeneratorFromL1T.h"


// Framework
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"

using namespace std;
using namespace edm;
using namespace l1t;

bool SortByL1Pt(L2MuonTrajectorySeed &A, L2MuonTrajectorySeed &B) {
  l1t::MuonRef Ref_L1A = A.l1tParticle();
  l1t::MuonRef Ref_L1B = B.l1tParticle();
  return (Ref_L1A->pt() > Ref_L1B->pt());
};

bool SortByL1QandPt(L2MuonTrajectorySeed &A, L2MuonTrajectorySeed &B) {
  l1t::MuonRef Ref_L1A = A.l1tParticle();
  l1t::MuonRef Ref_L1B = B.l1tParticle();

  // Compare quality first
  if (Ref_L1A->hwQual() > Ref_L1B->hwQual())  return true;
  if (Ref_L1A->hwQual() < Ref_L1B->hwQual())  return false;

  // For same quality L1s compare pT
  if (Ref_L1A->pt() > Ref_L1B->pt())  return true;
  if (Ref_L1A->pt() < Ref_L1B->pt())  return false;

  return false;
};

// constructors
L2MuonSeedGeneratorFromL1T::L2MuonSeedGeneratorFromL1T(const edm::ParameterSet& iConfig) :
  theSource(iConfig.getParameter<InputTag>("InputObjects")),
  theL1GMTReadoutCollection(iConfig.getParameter<InputTag>("GMTReadoutCollection")), // to be removed
  thePropagatorName(iConfig.getParameter<string>("Propagator")),
  theL1MinPt(iConfig.getParameter<double>("L1MinPt")),
  theL1MaxEta(iConfig.getParameter<double>("L1MaxEta")),
  theL1MinQuality(iConfig.getParameter<unsigned int>("L1MinQuality")),
  useOfflineSeed(iConfig.getUntrackedParameter<bool>("UseOfflineSeed", false)),
  useUnassociatedL1(iConfig.existsAs<bool>("UseUnassociatedL1") ?
            iConfig.getParameter<bool>("UseUnassociatedL1") : true),
  matchingDR(iConfig.getParameter<std::vector<double> >("MatchDR")),
  etaBins(iConfig.getParameter<std::vector<double> >("EtaMatchingBins")),
  centralBxOnly_( iConfig.getParameter<bool>("CentralBxOnly") ),

  useNewLogic( iConfig.getParameter<bool>("UseNewLogic") ), //Min
  printout( iConfig.getParameter<bool>("Printout") ), //Min
  // MatchType : 0 Min dR(1to1), 1 Higher Q(1to1), 2 All matched L1
  matchType( iConfig.getParameter<int>("MatchType") ), //Min
  // SortType : 0 not sort, 1 Pt, 2 Q and Pt
  sortType( iConfig.getParameter<int>("SortType") )  //Min

  {

    muCollToken_ = consumes<MuonBxCollection>(theSource);

    if(useOfflineSeed) {
      theOfflineSeedLabel = iConfig.getUntrackedParameter<InputTag>("OfflineSeedLabel");
      offlineSeedToken_ = consumes<edm::View<TrajectorySeed> >(theOfflineSeedLabel);

      // check that number of eta bins -1 matches number of dR cones
      if( matchingDR.size()!=etaBins.size() - 1  ) {
        throw cms::Exception("Configuration") << "Size of MatchDR "
                << "does not match number of eta bins." << endl;
      }

  }

  // service parameters
  ParameterSet serviceParameters = iConfig.getParameter<ParameterSet>("ServiceParameters");

  // the services
  theService = new MuonServiceProxy(serviceParameters);

  // the estimator
  theEstimator = new Chi2MeasurementEstimator(10000.);


  produces<L2MuonTrajectorySeedCollection>();
}

// destructor
L2MuonSeedGeneratorFromL1T::~L2MuonSeedGeneratorFromL1T(){
  if (theService)   delete theService;
  if (theEstimator) delete theEstimator;
}

void L2MuonSeedGeneratorFromL1T::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("GMTReadoutCollection", edm::InputTag("")); // to be removed
  desc.add<edm::InputTag>("InputObjects",edm::InputTag("hltGmtStage2Digis"));
  desc.add<string>("Propagator", "");
  desc.add<double>("L1MinPt",-1.);
  desc.add<double>("L1MaxEta",5.0);
  desc.add<unsigned int>("L1MinQuality",0);
  desc.addUntracked<bool>("UseOfflineSeed",false);
  desc.add<bool>("UseUnassociatedL1", true);
  desc.add<std::vector<double>>("MatchDR", {0.3});
  desc.add<std::vector<double>>("EtaMatchingBins", {0., 2.5});
  desc.add<bool>("CentralBxOnly", true);

  desc.add<bool>("UseNewLogic", false); //Min
  desc.add<bool>("Printout", true); //Min
  desc.add<int>("MatchType", 0); //Min
  desc.add<int>("SortType", 0); //Min

  desc.addUntracked<edm::InputTag>("OfflineSeedLabel", edm::InputTag(""));

  edm::ParameterSetDescription psd0;
  psd0.addUntracked<std::vector<std::string>>("Propagators", {
      "SteppingHelixPropagatorAny"
  });
  psd0.add<bool>("RPCLayers", true);
  psd0.addUntracked<bool>("UseMuonNavigation", true);
  desc.add<edm::ParameterSetDescription>("ServiceParameters", psd0);
  descriptions.add("L2MuonSeedGeneratorFromL1T",desc);
}

void L2MuonSeedGeneratorFromL1T::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  const std::string metname = "Muon|RecoMuon|L2MuonSeedGeneratorFromL1T";
  MuonPatternRecoDumper debug;

  if(printout) cout << endl;
  if(printout) cout << endl;
  if(printout) cout << iEvent.id().run() << ":" << iEvent.id().luminosityBlock() << ":" << iEvent.id().event() << endl;
  if(printout) cout << endl;

  auto output = std::make_unique<L2MuonTrajectorySeedCollection>();

  if(useNewLogic) { // New logic
    if(printout) cout << endl;
    if(printout) cout << "            |||||||||||||||||||||||||||" << endl;
    if(printout) cout << "            |||  Testing New Logic  |||" << endl;
    if(printout) cout << "            |||||||||||||||||||||||||||" << endl;

    edm::Handle<MuonBxCollection> muColl;
    iEvent.getByToken(muCollToken_, muColl);
    unsigned int NmuColl = muColl->size();
    /*LogTrace(metname)*/ if(printout) cout << "\n\n[ Number of L1 muons : " << NmuColl << "]" << endl;

    // For New Logic
    vector< vector<double> > dRmtx; //Min
    vector< vector<const TrajectorySeed *> > selOffseeds; //Min
    unsigned int i,j; // for the matrix element

    edm::Handle<edm::View<TrajectorySeed> > offlineSeedHandle;
    //vector<int> offlineSeedMap;
    if(useOfflineSeed) {
      iEvent.getByToken(offlineSeedToken_, offlineSeedHandle);
      unsigned int NofflineSeed = offlineSeedHandle->size();
      if(printout) cout << "[ Number of offline seeds : " << NofflineSeed << "]" << endl;

      //offlineSeedMap = vector<int>(NofflineSeed, 0);

      // Initialize dRmtx and selOffseeds
      dRmtx = vector< vector<double> >(NmuColl, vector<double>(NofflineSeed, 999)); //Min
      selOffseeds = vector< vector <const TrajectorySeed *> >(NmuColl, vector<const TrajectorySeed *>(NofflineSeed, 0)); //Min

      if(printout) cout << endl;   //Min
      if(printout) cout << "|||||||||||||||||" << endl;
      if(printout) cout << "|| init Matrix ||" << endl;
      if(printout) cout << "|||||||||||||||||" << endl;
      for(i=0; i<NmuColl; ++i) {
        for(j=0; j<NofflineSeed; ++j) {
          if(printout) cout << dRmtx[i][j] << "    ";
        }
        if(printout) cout << endl;
      }

      if(printout) cout << endl;
      if(printout) cout << "||||||||||||||||||||||" << endl;
      if(printout) cout << "|| init selOffSeeds ||" << endl;
      if(printout) cout << "||||||||||||||||||||||" << endl;
      for(i=0; i<NmuColl; ++i) {
        for(j=0; j<NofflineSeed; ++j) {
          if(printout) cout << selOffseeds[i][j] << "    ";
        }
        if(printout) cout << endl;
      }
    } // if(useOfflineSeed)

    bool isAsso = false;   // At least one L1mu are associated with offSeed   //Min

    for (int ibx = muColl->getFirstBX(); ibx <= muColl->getLastBX(); ++ibx) {
      if (centralBxOnly_ && (ibx != 0)) continue;
      if(printout) cout << endl;
      if(printout) cout << "[~ loop on L1 : ibx = " << ibx << " ~]" << endl;
      if(printout) cout << endl;

      for (auto it = muColl->begin(ibx); it != muColl->end(ibx); it++){

        unsigned int imu = distance(muColl->begin(muColl->getFirstBX()),it); //Position of given L1mu  //Min

        unsigned int quality = it->hwQual();
        int valid_charge = it->hwChargeValid();

        float pt    =  it->pt();
        float eta   =  it->eta();
        float theta =  2*atan(exp(-eta));
        float phi   =  it->phi();
        int charge  =  it->charge();
        // Set charge=0 for the time being if the valid charge bit is zero
        if (!valid_charge) charge = 0;

        int link = 36 + (int)(it -> tfMuonIndex() / 3.);
        bool barrel = true;
        if ( (link >= 36 && link <= 41) || (link >= 66 && link <= 71)) barrel = false;

        if ( pt < theL1MinPt || fabs(eta) > theL1MaxEta ) continue;

        /*LogTrace(metname)*/ if(printout) cout << "||||||||||||||||||||||||" << endl;
                              if(printout) cout << "||| New L2 Muon Seed |||" << endl;
                              if(printout) cout << "||| ith L1 Muon : " << imu << "  |||" << endl;
                              if(printout) cout << "||||||||||||||||||||||||" << endl;
        /*LogTrace(metname)*/ if(printout) cout << "Pt = "         << pt     << " GeV/c\n";
        /*LogTrace(metname)*/ if(printout) cout << "eta = "        << eta << "\n";
        /*LogTrace(metname)*/ if(printout) cout << "theta = "      << theta  << " rad\n";
        /*LogTrace(metname)*/ if(printout) cout << "phi = "        << phi    << " rad\n";
        /*LogTrace(metname)*/ if(printout) cout << "charge = "     << charge << "\n";
        /*LogTrace(metname)*/ if(printout) cout << "In Barrel? = " << barrel << "\n";
                              if(printout) cout << "quality = "<< quality << "\n";

        if ( quality <= theL1MinQuality ) {
          if(printout) cout << "!( quality <= theL1MinQuality )" << endl;
          continue;
        }

        // Update the services
        theService->update(iSetup);

        const DetLayer *detLayer = 0;
        float radius = 0.;

        CLHEP::Hep3Vector vec(0.,1.,0.);
        vec.setTheta(theta);
        vec.setPhi(phi);

        DetId theid;
        // Get the det layer on which the state should be put
        if ( barrel ){
          /*LogTrace(metname)*/ if(printout) cout << "The seed is in the barrel";

          // MB2
          DetId id = DTChamberId(0,2,0);
          detLayer = theService->detLayerGeometry()->idToLayer(id);
          /*LogTrace(metname)*/ if(printout) cout << "L2 Layer: " << debug.dumpLayer(detLayer);

          const BoundSurface* sur = &(detLayer->surface());
          const BoundCylinder* bc = dynamic_cast<const BoundCylinder*>(sur);

          radius = fabs(bc->radius()/sin(theta));
          theid  = id;

          /*LogTrace(metname)*/ if(printout) cout << "radius "<<radius;

          if ( pt < 3.5 ) pt = 3.5;
        }
        else {
          /*LogTrace(metname)*/ if(printout) cout << "The seed is in the endcap";

          DetId id;
          // ME2
          if ( theta < Geom::pi()/2. )
            id = CSCDetId(1,2,0,0,0);
          else
            id = CSCDetId(2,2,0,0,0);

          detLayer = theService->detLayerGeometry()->idToLayer(id);
          /*LogTrace(metname)*/ if(printout) cout << "L2 Layer: " << debug.dumpLayer(detLayer);

          radius = fabs(detLayer->position().z()/cos(theta));
          theid = id;

          if( pt < 1.0) pt = 1.0;
        }

        // Fallback solution using ME2
        DetId fallback_id;
        theta < Geom::pi()/2. ? fallback_id = CSCDetId(1,2,0,0,0) : fallback_id = CSCDetId(2,2,0,0,0);
        const DetLayer* ME2DetLayer = theService->detLayerGeometry()->idToLayer(fallback_id);

        vec.setMag(radius);

        GlobalPoint pos(vec.x(),vec.y(),vec.z());

        GlobalVector mom(pt*cos(phi), pt*sin(phi), pt*cos(theta)/sin(theta));

        GlobalTrajectoryParameters param(pos,mom,charge,&*theService->magneticField());
        AlgebraicSymMatrix55 mat;

        mat[0][0] = (0.25/pt)*(0.25/pt);  // sigma^2(charge/abs_momentum)
        if ( !barrel ) mat[0][0] = (0.4/pt)*(0.4/pt);

        //Assign q/pt = 0 +- 1/pt if charge has been declared invalid
        if (!valid_charge) mat[0][0] = (1./pt)*(1./pt);

        mat[1][1] = 0.05*0.05;        // sigma^2(lambda)
        mat[2][2] = 0.2*0.2;          // sigma^2(phi)
        mat[3][3] = 20.*20.;          // sigma^2(x_transverse))
        mat[4][4] = 20.*20.;          // sigma^2(y_transverse))

        CurvilinearTrajectoryError error(mat);

        const FreeTrajectoryState state(param,error);

        /*LogTrace(metname)*/ if(printout) cout << "Free trajectory State from the parameters";
        /*LogTrace(metname)*/ if(printout) cout << debug.dumpFTS(state);

        // Propagate the state on the MB2/ME2 surface
        TrajectoryStateOnSurface tsos = theService->propagator(thePropagatorName)->propagate(state, detLayer->surface());

        /*LogTrace(metname)*/ if(printout) cout << "State after the propagation on the layer";
        /*LogTrace(metname)*/ if(printout) cout << debug.dumpLayer(detLayer);
        /*LogTrace(metname)*/ if(printout) cout << debug.dumpTSOS(tsos);

        double dRcone = matchingDR[0];
        if ( fabs(eta) < etaBins.back() ){
          std::vector<double>::iterator lowEdge = std::upper_bound (etaBins.begin(), etaBins.end(), fabs(eta));
          dRcone    =  matchingDR.at( lowEdge - etaBins.begin() - 1);
        }



        if (tsos.isValid()) {

          edm::OwnVector<TrackingRecHit> container;

          if(useOfflineSeed) {

            if( ( !valid_charge || charge == 0 ) ) { //useOfflineSeed && ( !valid_charge || charge == 0)

              // Fill dRmtx and selOffseeds and return true if there exist Associated offSeeds for given L1(imu)
              bool isAssoOffseed = NewAssociateOfflineSeedToL1(offlineSeedHandle, dRmtx, tsos, imu, selOffseeds, dRcone );

              if(isAssoOffseed) {
                if(printout) cout << "\n isAssoOffseed ^_^ \n" << endl;
                isAsso = true;
              }

              else { // Using old way
                if(useUnassociatedL1) {
                  // convert the TSOS into a PTSOD
                  PTrajectoryStateOnDet const & seedTSOS = trajectoryStateTransform::persistentState( tsos, theid.rawId());
                  output->push_back(L2MuonTrajectorySeed(seedTSOS,container,alongMomentum,
                                    MuonRef(muColl,  distance(muColl->begin(muColl->getFirstBX()),it)  )));
                }
              }
            } //end of ( !valid_charge || charge == 0 )
            else if ( valid_charge ) {
            // Get the compatible dets on the layer
              std::vector< pair<const GeomDet*,TrajectoryStateOnSurface> >
                detsWithStates = detLayer->compatibleDets(tsos,
                                  *theService->propagator(thePropagatorName),
                                  *theEstimator);

              if (detsWithStates.size() == 0 && barrel ) {
                // try again to propagate but using ME2 as reference
                tsos = theService->propagator(thePropagatorName)->propagate(state, ME2DetLayer->surface());
                detsWithStates = ME2DetLayer->compatibleDets(tsos,
                                               *theService->propagator(thePropagatorName),
                                               *theEstimator);
              }

              if (detsWithStates.size()){

                TrajectoryStateOnSurface newTSOS = detsWithStates.front().second;
                const GeomDet *newTSOSDet = detsWithStates.front().first;

                if(printout) cout << "Most compatible det" << endl;
                if(printout) cout << debug.dumpMuonId(newTSOSDet->geographicalId()) << endl;

                if(printout) cout << "L1 info: Det and State:" << endl;
                if(printout) cout << debug.dumpMuonId(newTSOSDet->geographicalId()) << endl;

                if (newTSOS.isValid()){

                  //if(printout) cout << "(x, y, z) = (" << newTSOS.globalPosition().x() << ", "
                  //                  << newTSOS.globalPosition().y() << ", " << newTSOS.globalPosition().z() << ")";
                  if(printout) cout << "pos: (r=" << newTSOS.globalPosition().mag() << ", phi="
                                    << newTSOS.globalPosition().phi() << ", eta=" << newTSOS.globalPosition().eta() << ")" << endl;
                  if(printout) cout << "mom: (q*pt=" << newTSOS.charge()*newTSOS.globalMomentum().perp() << ", phi="
                                    << newTSOS.globalMomentum().phi() << ", eta=" << newTSOS.globalMomentum().eta() << ")" << endl;

                  //if(printout) cout << "State on it";
                  //if(printout) cout << debug.dumpTSOS(newTSOS);
                  // Get the compatible dets on the layer : end

                  bool isAssoOffseed = NewAssociateOfflineSeedToL1(offlineSeedHandle, dRmtx, newTSOS, imu, selOffseeds, dRcone );

                  if(isAssoOffseed) {
                    if(printout) cout << "\n isAssoOffseed ^_^ \n" << endl;
                    isAsso = true;
                  }
                  else { // Using old way
                    if(useUnassociatedL1) { // false
                      // convert the TSOS into a PTSOD
                      PTrajectoryStateOnDet const & seedTSOS = trajectoryStateTransform::persistentState( newTSOS,newTSOSDet->geographicalId().rawId());
                      output->push_back(L2MuonTrajectorySeed(seedTSOS,container,alongMomentum,
                                        MuonRef(muColl,  distance(muColl->begin(muColl->getFirstBX()),it)  )));
                    }
                  }
                }
              }
            } //end of ( valid_charge )
          } // end of ( useOfflineSeed )

          else { // Using old way
            // convert the TSOS into a PTSOD
            PTrajectoryStateOnDet const & seedTSOS = trajectoryStateTransform::persistentState( tsos, theid.rawId());
            output->push_back(L2MuonTrajectorySeed(seedTSOS,container,alongMomentum,
                              MuonRef(muColl,  distance(muColl->begin(muColl->getFirstBX()),it)  )));
          }
        } // end of ( tsos.isValid() )

      } // for (auto it = muColl->begin(ibx)
    } // for (int ibx = muColl->getFirstBX()

    // MatchType : 0 Min dR(1to1), 1 Higher Q(1to1), 2 All matched L1
    if( matchType == 0 ) { // Min dR, one-to-one
      if(useNewLogic && useOfflineSeed && isAsso) {
        unsigned int NofflineSeed1 = offlineSeedHandle->size();
        if(printout) cout << endl;
        if(printout) cout << "|||||||||||||||||" << endl;
        if(printout) cout << "|| New  Matrix ||" << endl;
        if(printout) cout << "|||||||||||||||||" << endl;

        for(i=0; i<NmuColl; ++i) {
          for(j=0; j<NofflineSeed1; ++j) {
            if(printout) cout << dRmtx[i][j] << "    ";
          }
          if(printout) cout << endl;
        }

        if(printout) cout << endl;
        if(printout) cout << "|||||||||||||||||||||" << endl;
        if(printout) cout << "|| New SelOffSeeds ||" << endl;
        if(printout) cout << "|||||||||||||||||||||" << endl;
        for(i=0; i<NmuColl; ++i) {
          for(j=0; j<NofflineSeed1; ++j) {
            if(printout) cout << selOffseeds[i][j] << "    ";
          }
          if(printout) cout << endl;
        }

        // vector of matched dR with order (just for printing out)
        vector<double> minOrder;

        unsigned int nL1;
        unsigned int i, j;
        for(nL1=0; nL1 < NmuColl; ++nL1) {
          double tempDR = 999;
          int theL1 = 0;
          int theOffs = 0;

          for(i=0; i<NmuColl; ++i) {
            for(j=0; j<NofflineSeed1; ++j) {
              if(tempDR > dRmtx[i][j]) {
                tempDR = dRmtx[i][j];
                theL1 = i;
                theOffs = j;
              }
            }
          }

          auto Newit = muColl->begin(muColl->getFirstBX()) + theL1;

          double NewdRcone = matchingDR[0];

          if ( fabs(Newit->eta()) < etaBins.back() ){
            std::vector<double>::iterator lowEdge = std::upper_bound (etaBins.begin(), etaBins.end(), fabs(Newit->eta()));
            NewdRcone    =  matchingDR.at( lowEdge - etaBins.begin() - 1);
          }

          if( !(tempDR < NewdRcone) )  break;  // FIXME

          if(printout) cout << "----------------------------" << endl;
          if(printout) cout << endl;
          if(printout) cout << nL1 <<"th minDR : " << tempDR << "     L1Q, Pt : " << Newit -> hwQual() << ", " << Newit -> pt() << "     L1 at " << theL1 << "     Offseed at " << theOffs << endl;

          minOrder.push_back(tempDR);

          // Remove row and column for given (L1mu, offSeed)
          for(i=0; i<NmuColl; ++i) {
            for(j=0; j<NofflineSeed1; ++j) {
              dRmtx[theL1][j] = 999;
              dRmtx[i][theOffs] = 999;
            }
          }

          // Check dRmtx
          if(printout) cout << endl;
          if(printout) cout << "|||||||||||||||||" << endl;
          if(printout) cout << "|| stp  Matrix ||" << endl;
          if(printout) cout << "|||||||||||||||||" << endl;

          for(i=0; i<NmuColl; ++i) {
            for(j=0; j<NofflineSeed1; ++j) {
              if(printout) cout << dRmtx[i][j] << "    ";
            }
            if(printout) cout << "\n";
          }

          //put given L1mu and offSeed to output
          edm::OwnVector<TrackingRecHit> Newcontainer; //Min

          PTrajectoryStateOnDet const & seedTSOS = selOffseeds[theL1][theOffs] ->startingState();
                        TrajectorySeed::const_iterator
                        tsci  = selOffseeds[theL1][theOffs]->recHits().first,
                        tscie = selOffseeds[theL1][theOffs]->recHits().second;
                        for(; tsci!=tscie; ++tsci) {
                          Newcontainer.push_back(*tsci);
                        }
                        output->push_back(L2MuonTrajectorySeed(seedTSOS,Newcontainer,alongMomentum,
                                          MuonRef(muColl,  distance(muColl->begin(muColl->getFirstBX()),Newit)  )));

                        if(printout) cout << theL1 << "th L1mu -> output" << endl;
        }

        // Just for check
        if(printout) cout << "\n minDR( 0th, 1st, 2nd, ...) : ";
        for(unsigned int l = 0; l<minOrder.size(); ++l) {
          if(printout) cout << minOrder[l] << "    ";
        }
        if(printout) cout << "\n\n";
      }

      else {
        if(printout) cout << "!isAssoOffseed -> No new matrix :(\n" << endl;
      }

    } // if( matchType == 0 )

    if( matchType == 1 ) { // Higher Q, one-to-one
      if(useNewLogic && useOfflineSeed && isAsso) {
        unsigned int NofflineSeed1 = offlineSeedHandle->size();
        if(printout) cout << endl;
        if(printout) cout << "|||||||||||||||||" << endl;
        if(printout) cout << "|| New  Matrix ||" << endl;
        if(printout) cout << "|||||||||||||||||" << endl;

        for(i=0; i<NmuColl; ++i) {
          for(j=0; j<NofflineSeed1; ++j) {
            if(printout) cout << dRmtx[i][j] << "    ";
          }
          if(printout) cout << endl;
        }

        if(printout) cout << endl;
        if(printout) cout << "|||||||||||||||||||||" << endl;
        if(printout) cout << "|| New SelOffSeeds ||" << endl;
        if(printout) cout << "|||||||||||||||||||||" << endl;
        for(i=0; i<NmuColl; ++i) {
          for(j=0; j<NofflineSeed1; ++j) {
            if(printout) cout << selOffseeds[i][j] << "    ";
          }
          if(printout) cout << endl;
        }

        // vector of matched dR with order (just for printing out)
        vector<double> minOrder;

        unsigned int nL1;
        unsigned int i, j;
        for(nL1=0; nL1 < NmuColl; ++nL1) {
          double tempDR = 999;
          int theL1 = 0;
          int theOffs = 0;
          auto theit = muColl->begin(muColl->getFirstBX());

          // L1Q > 10 (L1Q = 12)
          for(i=0; i<NmuColl; ++i) {
            theit = muColl->begin(muColl->getFirstBX()) + i;
            if (theit->hwQual() > 10) {
              for(j=0; j<NofflineSeed1; ++j) {
                if(tempDR > dRmtx[i][j]) {
                  tempDR = dRmtx[i][j];
                  theL1 = i;
                  theOffs = j;
                }
              }
            }
          }

          // 6 < L1Q <= 10 (L1Q = 8)
          if (tempDR == 999) {
            for(i=0; i<NmuColl; ++i) {
              theit = muColl->begin(muColl->getFirstBX()) + i;
              if ( (theit->hwQual() <= 10) && (theit->hwQual() > 6) ) {
                for(j=0; j<NofflineSeed1; ++j) {
                  if(tempDR > dRmtx[i][j]) {
                    tempDR = dRmtx[i][j];
                    theL1 = i;
                    theOffs = j;
                  }
                }
              }
            }
          }

          // L1Q <= 6 (L1Q = 4)
          if (tempDR == 999) {
            for(i=0; i<NmuColl; ++i) {
              theit = muColl->begin(muColl->getFirstBX()) + i;
              if (theit->hwQual() <= 6) {
                for(j=0; j<NofflineSeed1; ++j) {
                  if(tempDR > dRmtx[i][j]) {
                    tempDR = dRmtx[i][j];
                    theL1 = i;
                    theOffs = j;
                  }
                }
              }
            }
          }

          auto Newit = muColl->begin(muColl->getFirstBX()) + theL1;

          double NewdRcone = matchingDR[0];

          if ( fabs(Newit->eta()) < etaBins.back() ){
            std::vector<double>::iterator lowEdge = std::upper_bound (etaBins.begin(), etaBins.end(), fabs(Newit->eta()));
            NewdRcone    =  matchingDR.at( lowEdge - etaBins.begin() - 1);
          }

          if( !(tempDR < NewdRcone) )  break;  // FIXME

          if(printout) cout << "----------------------------" << endl;
          if(printout) cout << endl;
          if(printout) cout << nL1 <<"th minDR : " << tempDR << "     L1Q, Pt : " << Newit -> hwQual() << ", " << Newit -> pt() << "     L1 at " << theL1 << "     Offseed at " << theOffs << endl;

          minOrder.push_back(tempDR);

          // Remove row and column for given (L1mu, offSeed)
          for(i=0; i<NmuColl; ++i) {
            for(j=0; j<NofflineSeed1; ++j) {
              dRmtx[theL1][j] = 999;
              dRmtx[i][theOffs] = 999;
            }
          }

          // Check dRmtx
          if(printout) cout << endl;
          if(printout) cout << "|||||||||||||||||" << endl;
          if(printout) cout << "|| stp  Matrix ||" << endl;
          if(printout) cout << "|||||||||||||||||" << endl;

          for(i=0; i<NmuColl; ++i) {
            for(j=0; j<NofflineSeed1; ++j) {
              if(printout) cout << dRmtx[i][j] << "    ";
            }
            if(printout) cout << "\n";
          }

          //put given L1mu and offSeed to output
          edm::OwnVector<TrackingRecHit> Newcontainer; //Min

          PTrajectoryStateOnDet const & seedTSOS = selOffseeds[theL1][theOffs] ->startingState();
                        TrajectorySeed::const_iterator
                        tsci  = selOffseeds[theL1][theOffs]->recHits().first,
                        tscie = selOffseeds[theL1][theOffs]->recHits().second;
                        for(; tsci!=tscie; ++tsci) {
                          Newcontainer.push_back(*tsci);
                        }
                        output->push_back(L2MuonTrajectorySeed(seedTSOS,Newcontainer,alongMomentum,
                                          MuonRef(muColl,  distance(muColl->begin(muColl->getFirstBX()),Newit)  )));

                        if(printout) cout << theL1 << "th L1mu -> output" << endl;
        }

        // Just for check
        if(printout) cout << "\n minDR( 0th, 1st, 2nd, ...) : ";
        for(unsigned int l = 0; l<minOrder.size(); ++l) {
          if(printout) cout << minOrder[l] << "    ";
        }
        if(printout) cout << "\n\n";
      }

      else {
        if(printout) cout << "!isAssoOffseed -> No new matrix :(\n" << endl;
      }

    } // if( matchType == 1 )

    else if( matchType == 2 ) { // Put all the matched L1

      if(useNewLogic && useOfflineSeed && isAsso) {
        unsigned int NofflineSeed1 = offlineSeedHandle->size();
        if(printout) cout << endl;
        if(printout) cout << "|||||||||||||||||" << endl;
        if(printout) cout << "|| New  Matrix ||" << endl;
        if(printout) cout << "|||||||||||||||||" << endl;

        for(i=0; i<NmuColl; ++i) {
          for(j=0; j<NofflineSeed1; ++j) {
            if(printout) cout << dRmtx[i][j] << "    ";
          }
          if(printout) cout << endl;
        }

        if(printout) cout << endl;
        if(printout) cout << "|||||||||||||||||||||" << endl;
        if(printout) cout << "|| New SelOffSeeds ||" << endl;
        if(printout) cout << "|||||||||||||||||||||" << endl;
        for(i=0; i<NmuColl; ++i) {
          for(j=0; j<NofflineSeed1; ++j) {
            if(printout) cout << selOffseeds[i][j] << "    ";
          }
          if(printout) cout << endl;
        }

        vector<double> outputOrder;

        for(i=0; i<NmuColl; ++i) {
          auto Newit = muColl->begin(muColl->getFirstBX()) + i;
          double tempDR = 999;
          int theOffs = 0;

          double NewdRcone = matchingDR[0];
          if ( fabs(Newit->eta()) < etaBins.back() ){
            std::vector<double>::iterator lowEdge = std::upper_bound (etaBins.begin(), etaBins.end(), fabs(Newit->eta()));
            NewdRcone    =  matchingDR.at( lowEdge - etaBins.begin() - 1);
          }

          for(j=0; j<NofflineSeed1; ++j) {
            if( tempDR > dRmtx[i][j] ) {
              tempDR = dRmtx[i][j];
              theOffs = j;
            }
          }

          if( !(tempDR < NewdRcone) )  continue;  // FIXME

          if( selOffseeds[i][theOffs] != 0 ) {
            //put given L1mu and offSeed to output
            edm::OwnVector<TrackingRecHit> Newcontainer; //Min

            PTrajectoryStateOnDet const & seedTSOS = selOffseeds[i][theOffs] ->startingState();
                          TrajectorySeed::const_iterator
                          tsci  = selOffseeds[i][theOffs]->recHits().first,
                          tscie = selOffseeds[i][theOffs]->recHits().second;
                          for(; tsci!=tscie; ++tsci) {
                            Newcontainer.push_back(*tsci);
                          }
                          output->push_back(L2MuonTrajectorySeed(seedTSOS,Newcontainer,alongMomentum,
                                            MuonRef(muColl,  distance(muColl->begin(muColl->getFirstBX()),Newit)  )));

                          if(printout) cout << i << "th L1mu -> output" << endl;

            outputOrder.push_back(dRmtx[i][theOffs]);
          }
        }

        // Just for check
        if(printout) cout << "\n outputDR( 0th, 1st, 2nd, ...) : ";
        for(unsigned int l = 0; l<outputOrder.size(); ++l) {
          if(printout) cout << outputOrder[l] << "    ";
        }
        if(printout) cout << "\n\n";
      }

      else {
        if(printout) cout << "!isAssoOffseed -> No new matrix :(\n" << endl;
      }

    } // else if( matchType == 2 )

    dRmtx.clear();
    selOffseeds.clear();

  } // New Logic

  else if(!useNewLogic) { // Old Logic

    if(printout) cout << endl;
    if(printout) cout << "            |||||||||||||||||||" << endl;
    if(printout) cout << "            |||  Old Logic  |||" << endl;
    if(printout) cout << "            |||||||||||||||||||" << endl;

    edm::Handle<MuonBxCollection> muColl;
    iEvent.getByToken(muCollToken_, muColl);
    /*LogTrace(metname)*/ if(printout) cout << "Number of muons " << muColl->size() << endl;

    edm::Handle<edm::View<TrajectorySeed> > offlineSeedHandle;
    vector<int> offlineSeedMap;
    if(useOfflineSeed) {
      iEvent.getByToken(offlineSeedToken_, offlineSeedHandle);
      /*LogTrace(metname)*/ if(printout) cout << "Number of offline seeds " << offlineSeedHandle->size() << endl;
      offlineSeedMap = vector<int>(offlineSeedHandle->size(), 0);
    }

    for (int ibx = muColl->getFirstBX(); ibx <= muColl->getLastBX(); ++ibx) {
      if (centralBxOnly_ && (ibx != 0)) continue;
      for (auto it = muColl->begin(ibx); it != muColl->end(ibx); it++){

        unsigned int quality = it->hwQual();
        int valid_charge = it->hwChargeValid();

        float pt    =  it->pt();
        float eta   =  it->eta();
        float theta =  2*atan(exp(-eta));
        float phi   =  it->phi();
        int charge  =  it->charge();
        // Set charge=0 for the time being if the valid charge bit is zero
        if (!valid_charge) charge = 0;

        int link = 36 + (int)(it -> tfMuonIndex() / 3.);
        bool barrel = true;
        if ( (link >= 36 && link <= 41) || (link >= 66 && link <= 71)) barrel = false;

        if ( pt < theL1MinPt || fabs(eta) > theL1MaxEta ) continue;

        /*LogTrace(metname)*/ if(printout) cout << "New L2 Muon Seed" ;
        /*LogTrace(metname)*/ if(printout) cout << "Pt = "         << pt     << " GeV/c";
        /*LogTrace(metname)*/ if(printout) cout << "eta = "        << eta;
        /*LogTrace(metname)*/ if(printout) cout << "theta = "      << theta  << " rad";
        /*LogTrace(metname)*/ if(printout) cout << "phi = "        << phi    << " rad";
        /*LogTrace(metname)*/ if(printout) cout << "charge = "     << charge;
        /*LogTrace(metname)*/ if(printout) cout << "In Barrel? = " << barrel;

        if ( quality <= theL1MinQuality ) continue;
        /*LogTrace(metname)*/ if(printout) cout << "quality = "<< quality;

        // Update the services
        theService->update(iSetup);

        const DetLayer *detLayer = 0;
        float radius = 0.;

        CLHEP::Hep3Vector vec(0.,1.,0.);
        vec.setTheta(theta);
        vec.setPhi(phi);

        DetId theid;
        // Get the det layer on which the state should be put
        if ( barrel ){
          /*LogTrace(metname)*/ if(printout) cout << "The seed is in the barrel";

          // MB2
          DetId id = DTChamberId(0,2,0);
          detLayer = theService->detLayerGeometry()->idToLayer(id);
          /*LogTrace(metname)*/ if(printout) cout << "L2 Layer: " << debug.dumpLayer(detLayer);

          const BoundSurface* sur = &(detLayer->surface());
          const BoundCylinder* bc = dynamic_cast<const BoundCylinder*>(sur);

          radius = fabs(bc->radius()/sin(theta));
          theid  = id;

          /*LogTrace(metname)*/ if(printout) cout << "radius "<<radius;

          if ( pt < 3.5 ) pt = 3.5;
        }
        else {
          /*LogTrace(metname)*/ if(printout) cout << "The seed is in the endcap";

          DetId id;
          // ME2
          if ( theta < Geom::pi()/2. )
            id = CSCDetId(1,2,0,0,0);
          else
            id = CSCDetId(2,2,0,0,0);

          detLayer = theService->detLayerGeometry()->idToLayer(id);
          /*LogTrace(metname)*/ if(printout) cout << "L2 Layer: " << debug.dumpLayer(detLayer);

          radius = fabs(detLayer->position().z()/cos(theta));
          theid = id;

          if( pt < 1.0) pt = 1.0;
        }

        // Fallback solution using ME2
        DetId fallback_id;
        theta < Geom::pi()/2. ? fallback_id = CSCDetId(1,2,0,0,0) : fallback_id = CSCDetId(2,2,0,0,0);
        const DetLayer* ME2DetLayer = theService->detLayerGeometry()->idToLayer(fallback_id);

        vec.setMag(radius);

        GlobalPoint pos(vec.x(),vec.y(),vec.z());

        GlobalVector mom(pt*cos(phi), pt*sin(phi), pt*cos(theta)/sin(theta));

        GlobalTrajectoryParameters param(pos,mom,charge,&*theService->magneticField());
        AlgebraicSymMatrix55 mat;

        mat[0][0] = (0.25/pt)*(0.25/pt);  // sigma^2(charge/abs_momentum)
        if ( !barrel ) mat[0][0] = (0.4/pt)*(0.4/pt);

        //Assign q/pt = 0 +- 1/pt if charge has been declared invalid
        if (!valid_charge) mat[0][0] = (1./pt)*(1./pt);

        mat[1][1] = 0.05*0.05;        // sigma^2(lambda)
        mat[2][2] = 0.2*0.2;          // sigma^2(phi)
        mat[3][3] = 20.*20.;          // sigma^2(x_transverse))
        mat[4][4] = 20.*20.;          // sigma^2(y_transverse))

        CurvilinearTrajectoryError error(mat);

        const FreeTrajectoryState state(param,error);

        /*LogTrace(metname)*/ if(printout) cout << "Free trajectory State from the parameters";
        /*LogTrace(metname)*/ if(printout) cout << debug.dumpFTS(state);

        // Propagate the state on the MB2/ME2 surface
        TrajectoryStateOnSurface tsos = theService->propagator(thePropagatorName)->propagate(state, detLayer->surface());

        /*LogTrace(metname)*/ if(printout) cout << "State after the propagation on the layer";
        /*LogTrace(metname)*/ if(printout) cout << debug.dumpLayer(detLayer);
        /*LogTrace(metname)*/ if(printout) cout << debug.dumpTSOS(tsos);

      double dRcone = matchingDR[0];
      if ( fabs(eta) < etaBins.back() ){
      std::vector<double>::iterator lowEdge = std::upper_bound (etaBins.begin(), etaBins.end(), fabs(eta));
      dRcone    =  matchingDR.at( lowEdge - etaBins.begin() - 1);
      }

        if (tsos.isValid()) {

          edm::OwnVector<TrackingRecHit> container;

          if(useOfflineSeed && ( !valid_charge || charge == 0) ) {

            const TrajectorySeed *assoOffseed =
              associateOfflineSeedToL1(offlineSeedHandle, offlineSeedMap, tsos, dRcone );

            if(assoOffseed!=0) {
              PTrajectoryStateOnDet const & seedTSOS = assoOffseed->startingState();
              TrajectorySeed::const_iterator
              tsci  = assoOffseed->recHits().first,
              tscie = assoOffseed->recHits().second;
              for(; tsci!=tscie; ++tsci) {
                container.push_back(*tsci);
              }
              output->push_back(L2MuonTrajectorySeed(seedTSOS,container,alongMomentum,
                                MuonRef(muColl,  distance(muColl->begin(muColl->getFirstBX()),it)  )));
            }
            else {
              if(useUnassociatedL1) {
                // convert the TSOS into a PTSOD
                PTrajectoryStateOnDet const & seedTSOS = trajectoryStateTransform::persistentState( tsos, theid.rawId());
                output->push_back(L2MuonTrajectorySeed(seedTSOS,container,alongMomentum,
                                  MuonRef(muColl,  distance(muColl->begin(muColl->getFirstBX()),it)  )));
              }
            }
          }
          else if (useOfflineSeed && valid_charge){
            // Get the compatible dets on the layer
            std::vector< pair<const GeomDet*,TrajectoryStateOnSurface> >
              detsWithStates = detLayer->compatibleDets(tsos,
                                *theService->propagator(thePropagatorName),
                                *theEstimator);

            if (detsWithStates.size() == 0 && barrel ) {
              // try again to propagate but using ME2 as reference
              tsos = theService->propagator(thePropagatorName)->propagate(state, ME2DetLayer->surface());
              detsWithStates = ME2DetLayer->compatibleDets(tsos,
                                             *theService->propagator(thePropagatorName),
                                             *theEstimator);
            }

            if (detsWithStates.size()){

              TrajectoryStateOnSurface newTSOS = detsWithStates.front().second;
              const GeomDet *newTSOSDet = detsWithStates.front().first;

              /*LogTrace(metname)*/ if(printout) cout << "Most compatible det";
              /*LogTrace(metname)*/ if(printout) cout << debug.dumpMuonId(newTSOSDet->geographicalId());

              /*LogDebug(metname)*/ if(printout) cout << "L1 info: Det and State:";
              /*LogDebug(metname)*/ if(printout) cout << debug.dumpMuonId(newTSOSDet->geographicalId());

              if (newTSOS.isValid()){

                ///*LogDebug(metname)*/ if(printout) cout << "(x, y, z) = (" << newTSOS.globalPosition().x() << ", "
                //                  << newTSOS.globalPosition().y() << ", " << newTSOS.globalPosition().z() << ")";
                /*LogDebug(metname)*/ if(printout) cout << "pos: (r=" << newTSOS.globalPosition().mag() << ", phi="
                                  << newTSOS.globalPosition().phi() << ", eta=" << newTSOS.globalPosition().eta() << ")";
                /*LogDebug(metname)*/ if(printout) cout << "mom: (q*pt=" << newTSOS.charge()*newTSOS.globalMomentum().perp() << ", phi="
                                  << newTSOS.globalMomentum().phi() << ", eta=" << newTSOS.globalMomentum().eta() << ")";

                ///*LogDebug(metname)*/ if(printout) cout << "State on it";
                ///*LogDebug(metname)*/ if(printout) cout << debug.dumpTSOS(newTSOS);

                const TrajectorySeed *assoOffseed =
                  associateOfflineSeedToL1(offlineSeedHandle, offlineSeedMap, newTSOS, dRcone);

                if(assoOffseed!=0) {
                  PTrajectoryStateOnDet const & seedTSOS = assoOffseed->startingState();
                  TrajectorySeed::const_iterator
                    tsci  = assoOffseed->recHits().first,
                    tscie = assoOffseed->recHits().second;
                  for(; tsci!=tscie; ++tsci) {
                    container.push_back(*tsci);
                  }
                  output->push_back(L2MuonTrajectorySeed(seedTSOS,container,alongMomentum,
                                    MuonRef(muColl,  distance(muColl->begin(muColl->getFirstBX()),it)  )));
                }
                else {
                  if(useUnassociatedL1) {
                    // convert the TSOS into a PTSOD
                    PTrajectoryStateOnDet const & seedTSOS = trajectoryStateTransform::persistentState( newTSOS,newTSOSDet->geographicalId().rawId());
                    output->push_back(L2MuonTrajectorySeed(seedTSOS,container,alongMomentum,
                                      MuonRef(muColl,  distance(muColl->begin(muColl->getFirstBX()),it)  )));
                  }
                }
              }
            }
          }
          else {
            // convert the TSOS into a PTSOD
            PTrajectoryStateOnDet const & seedTSOS = trajectoryStateTransform::persistentState( tsos, theid.rawId());
            output->push_back(L2MuonTrajectorySeed(seedTSOS,container,alongMomentum,
                              MuonRef(muColl,  distance(muColl->begin(muColl->getFirstBX()),it)  )));
          }
        }
      }
    }

  } // Old Logic


  if(printout) cout << "Print output" << endl;
  for(unsigned int k=0; k<output->size(); ++k){
    l1t::MuonRef Ref_L1 = output->at(k).l1tParticle();
    if(printout) cout << Ref_L1->pt() << "    ";
  }
  if(printout) cout << "\n\n";

  if(sortType == 1) { // Sort output by L1 pT
    sort(output->begin(),output->end(),SortByL1Pt);
  }

  if(sortType == 2) { // Sort output by L1 pT
    sort(output->begin(),output->end(),SortByL1QandPt);
  }

  for(unsigned int k=0; k<output->size(); ++k){
    l1t::MuonRef Ref_L1 = output->at(k).l1tParticle();
    if(printout) cout << Ref_L1->pt() << "    ";
  }
  if(printout) cout << "\n\n";

  if(printout) cout << endl;
  if(printout) cout << "-----iEvent.put(output)-----" << endl;
  iEvent.put(std::move(output));
  if(printout) cout << endl;
}



bool L2MuonSeedGeneratorFromL1T::NewAssociateOfflineSeedToL1( edm::Handle<edm::View<TrajectorySeed> > & offseeds,
                                     std::vector< std::vector<double> > & dRmtx, //Min
                                     TrajectoryStateOnSurface & newTsos,
                                     unsigned int imu, //Min
                                     std::vector< std::vector<const TrajectorySeed *> > & selOffseeds, //Min
                                     double dRcone ) {   //Min

  const std::string metlabel = "Muon|RecoMuon|L2MuonSeedGeneratorFromL1T";
  bool isAssociated = false;
  MuonPatternRecoDumper debugtmp;

  edm::View<TrajectorySeed>::const_iterator offseed, endOffseed = offseeds->end();
  unsigned int nOffseed(0);

  if(printout) cout << endl;
  if(printout) cout << "[~ loop on L2 off seeds ~]" << endl;
  if(printout) cout << endl;

  for(offseed=offseeds->begin(); offseed!=endOffseed; ++offseed, ++nOffseed) {
    if(printout) cout << "||||||||||||||||||||||||" << endl;
    if(printout) cout << "||| " << nOffseed << "th Offline Seed |||" << endl;
    if(printout) cout << "||||||||||||||||||||||||" << endl;


    GlobalPoint  glbPos = theService->trackingGeometry()->idToDet(offseed->startingState().detId())->surface().toGlobal(offseed->startingState().parameters().position());
    GlobalVector glbMom = theService->trackingGeometry()->idToDet(offseed->startingState().detId())->surface().toGlobal(offseed->startingState().parameters().momentum());

    // Preliminary check
    double preDr = deltaR( newTsos.globalPosition().eta(), newTsos.globalPosition().phi(), glbPos.eta(), glbPos.phi() );
    if(preDr > 1.0) {
      if(printout) cout << "!(preDr > 1.0)" << endl;
      continue;
    }

    const FreeTrajectoryState offseedFTS(glbPos, glbMom, offseed->startingState().parameters().charge(), &*theService->magneticField());
    TrajectoryStateOnSurface offseedTsos = theService->propagator(thePropagatorName)->propagate(offseedFTS, newTsos.surface());
    if(printout) cout << "Offline seed info: Det and State" << std::endl;
    if(printout) cout << debugtmp.dumpMuonId(offseed->startingState().detId()) << std::endl;
    //if(printout) cout << "(x, y, z) = (" << newTSOS.globalPosition().x() << ", "
    //                     << newTSOS.globalPosition().y() << ", " << newTSOS.globalPosition().z() << ")" << std::endl;
    if(printout) cout << "pos: (r=" << offseedFTS.position().mag() <<
                              ", phi=" << offseedFTS.position().phi() <<
                              ", eta=" << offseedFTS.position().eta() << ")" << std::endl;
    if(printout) cout << "mom: (q*pt=" << offseedFTS.charge()*offseedFTS.momentum().perp() <<
                              ", phi=" << offseedFTS.momentum().phi() <<
                              ", eta=" << offseedFTS.momentum().eta() << ")" << std::endl << std::endl;
    //if(printout) cout << debugtmp.dumpFTS(offseedFTS) << std::endl;

    if(offseedTsos.isValid()) {
      if(printout) cout << "Offline seed info after propagation to L1 layer:" << std::endl;
      //if(printout) cout << "(x, y, z) = (" << offseedTsos.globalPosition().x() << ", "
      //                   << offseedTsos.globalPosition().y() << ", " << offseedTsos.globalPosition().z() << ")" << std::endl;
      if(printout) cout << "pos: (r=" << offseedTsos.globalPosition().mag() <<
                                ", phi=" << offseedTsos.globalPosition().phi() <<
                                ", eta=" << offseedTsos.globalPosition().eta() << ")" << std::endl;
      if(printout) cout << "mom: (q*pt=" << offseedTsos.charge()*offseedTsos.globalMomentum().perp() <<
                                ", phi=" << offseedTsos.globalMomentum().phi() <<
                                ", eta=" << offseedTsos.globalMomentum().eta() << ")" << std::endl << std::endl;
      //if(printout) cout << debugtmp.dumpTSOS(offseedTsos) << std::endl;
      double newDr = deltaR( newTsos.globalPosition().eta(),     newTsos.globalPosition().phi(),
                             offseedTsos.globalPosition().eta(), offseedTsos.globalPosition().phi() );


      if(printout) cout << "   -- DR = " << newDr << std::endl;
      if( newDr < dRcone ) {
        if(printout) cout << "          --> OK! " << newDr << std::endl << std::endl;

        dRmtx[imu][nOffseed] = newDr;
        selOffseeds[imu][nOffseed] = &*offseed;

        isAssociated = true;
      }
      else {
        if(printout) cout << "          --> Rejected. " << newDr << std::endl << std::endl;
      }
    }
    else {
      if(printout) cout << "Invalid offline seed TSOS after propagation!" << std::endl << std::endl;
    }
  }

  return isAssociated;
}


// FIXME: does not resolve ambiguities yet!
const TrajectorySeed* L2MuonSeedGeneratorFromL1T::associateOfflineSeedToL1( edm::Handle<edm::View<TrajectorySeed> > & offseeds,
                                     std::vector<int> & offseedMap,
                                     TrajectoryStateOnSurface & newTsos,
                                     double dRcone ) {

  const std::string metlabel = "Muon|RecoMuon|L2MuonSeedGeneratorFromL1T";
  MuonPatternRecoDumper debugtmp;

  edm::View<TrajectorySeed>::const_iterator offseed, endOffseed = offseeds->end();
  const TrajectorySeed *selOffseed = 0;
  double bestDr = 999.;
  unsigned int nOffseed(0);
  int lastOffseed(-1);

  for(offseed=offseeds->begin(); offseed!=endOffseed; ++offseed, ++nOffseed) {
    if(offseedMap[nOffseed]!=0) continue;
    GlobalPoint  glbPos = theService->trackingGeometry()->idToDet(offseed->startingState().detId())->surface().toGlobal(offseed->startingState().parameters().position());
    GlobalVector glbMom = theService->trackingGeometry()->idToDet(offseed->startingState().detId())->surface().toGlobal(offseed->startingState().parameters().momentum());

    // Preliminary check
    double preDr = deltaR( newTsos.globalPosition().eta(), newTsos.globalPosition().phi(), glbPos.eta(), glbPos.phi() );
    if(preDr > 1.0) continue;

    const FreeTrajectoryState offseedFTS(glbPos, glbMom, offseed->startingState().parameters().charge(), &*theService->magneticField());
    TrajectoryStateOnSurface offseedTsos = theService->propagator(thePropagatorName)->propagate(offseedFTS, newTsos.surface());
    /*LogDebug(metlabel)*/ if(printout) cout << "Offline seed info: Det and State" << std::endl;
    /*LogDebug(metlabel)*/ if(printout) cout << debugtmp.dumpMuonId(offseed->startingState().detId()) << std::endl;
    ///*LogDebug(metlabel)*/ if(printout) cout << "(x, y, z) = (" << newTSOS.globalPosition().x() << ", "
    //                     << newTSOS.globalPosition().y() << ", " << newTSOS.globalPosition().z() << ")" << std::endl;
    /*LogDebug(metlabel)*/ if(printout) cout << "pos: (r=" << offseedFTS.position().mag() << ", phi="
               << offseedFTS.position().phi() << ", eta=" << offseedFTS.position().eta() << ")" << std::endl;
    /*LogDebug(metlabel)*/ if(printout) cout << "mom: (q*pt=" << offseedFTS.charge()*offseedFTS.momentum().perp() << ", phi="
               << offseedFTS.momentum().phi() << ", eta=" << offseedFTS.momentum().eta() << ")" << std::endl << std::endl;
    ///*LogDebug(metlabel)*/ if(printout) cout << debugtmp.dumpFTS(offseedFTS) << std::endl;

    if(offseedTsos.isValid()) {
      /*LogDebug(metlabel)*/ if(printout) cout << "Offline seed info after propagation to L1 layer:" << std::endl;
      ///*LogDebug(metlabel)*/ if(printout) cout << "(x, y, z) = (" << offseedTsos.globalPosition().x() << ", "
      //                   << offseedTsos.globalPosition().y() << ", " << offseedTsos.globalPosition().z() << ")" << std::endl;
      /*LogDebug(metlabel)*/ if(printout) cout << "pos: (r=" << offseedTsos.globalPosition().mag() << ", phi="
             << offseedTsos.globalPosition().phi() << ", eta=" << offseedTsos.globalPosition().eta() << ")" << std::endl;
      /*LogDebug(metlabel)*/ if(printout) cout << "mom: (q*pt=" << offseedTsos.charge()*offseedTsos.globalMomentum().perp() << ", phi="
             << offseedTsos.globalMomentum().phi() << ", eta=" << offseedTsos.globalMomentum().eta() << ")" << std::endl << std::endl;
      ///*LogDebug(metlabel)*/ if(printout) cout << debugtmp.dumpTSOS(offseedTsos) << std::endl;
      double newDr = deltaR( newTsos.globalPosition().eta(),     newTsos.globalPosition().phi(),
                 offseedTsos.globalPosition().eta(), offseedTsos.globalPosition().phi() );
      /*LogDebug(metlabel)*/ if(printout) cout << "   -- DR = " << newDr << std::endl;
      if( newDr < dRcone && newDr<bestDr ) {
        /*LogDebug(metlabel)*/ if(printout) cout << "          --> OK! " << newDr << std::endl << std::endl;
        selOffseed = &*offseed;
        bestDr = newDr;
        offseedMap[nOffseed] = 1;
        if(lastOffseed>-1) offseedMap[lastOffseed] = 0;
        lastOffseed = nOffseed;
      }
      else {
        /*LogDebug(metlabel)*/ if(printout) cout << "          --> Rejected. " << newDr << std::endl << std::endl;
      }
    }
    else {
      /*LogDebug(metlabel)*/ if(printout) cout << "Invalid offline seed TSOS after propagation!" << std::endl << std::endl;
    }
  }

  return selOffseed;
}
