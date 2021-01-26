
//for getting detid
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include <map>

#include <TString.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "ggAnalysis/ggNtuplizer/interface/GEDPhoIDTools.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

using namespace std;


//pfpohton for crystal selection from cluster
Int_t nPFPhoton;
Float_t PFPhotonE[100];
Float_t PFPhotonET[100];
Float_t PFPhotonHadEM[100];
Float_t PFPhotonEta[100];
Float_t PFPhotonPhi[100];
Float_t PFPhotonSigmaIetaIeta[100];
Float_t PFPhotonR9[100];
Float_t PFChargedHadronIso[100];
Float_t PFNeutralHadronIso[100];
Float_t PFPhotonIso[100];
Float_t PFMIPTotE[100];
Int_t PFMIPNhit[100];
Float_t PFEcalIso03[100];
Float_t PFHcalIso03[100];
Float_t PFTrkIso03[100];
Int_t PFPixMatch[100];


//detid array   
Int_t   nAllCellsEB;
Int_t   AllCellsIEtaEB[30000];
Int_t   AllCellsIPhiEB[30000];
Float_t AllCellsE_EB[30000];
Float_t AllTimeEB[30000];
Int_t AllClusteredEB[30000][30];
Float_t AllFracEB[30000][30];


vector<float> phoSeedIEta_;
vector<float> phoSeedIPhi_;
vector<float> ophoSeedIEta_;
vector<float> ophoSeedIPhi_;

//ootphoton
Int_t          onPho_;
vector<float>  ophoE_;
vector<float>  ophoSigmaE_;
vector<float>  ophoEt_;
vector<float>  ophoEta_;
vector<float>  ophoPhi_;
vector<float>  ophoCalibE_;
vector<float>  ophoCalibEt_;
vector<float>  ophoSCE_;
vector<float>  ophoSCRawE_;
vector<float>  ophoESEnP1_;
vector<float>  ophoESEnP2_;
vector<float>  ophoSCEta_;
vector<float>  ophoSCPhi_;
vector<float>  ophoSCEtaWidth_;
vector<float>  ophoSCPhiWidth_;
vector<float>  ophoSCBrem_;
vector<int>    ophohasPixelSeed_;
vector<int>    ophoEleVeto_;
vector<float>  ophoR9_;
vector<float>  ophoHoverE_;
vector<float>  ophoESEffSigmaRR_;
vector<float>  ophoSigmaIEtaIEtaFull5x5_;
vector<float>  ophoSigmaIEtaIPhiFull5x5_;
vector<float>  ophoSigmaIPhiIPhiFull5x5_;
vector<float>  ophoE2x2Full5x5_;
vector<float>  ophoE5x5Full5x5_;
vector<float>  ophoR9Full5x5_;
vector<float>  ophoPFChIso_;
vector<float>  ophoPFChPVIso_;
vector<float>  ophoPFPhoIso_;
vector<float>  ophoPFNeuIso_;
vector<float>  ophoPFChWorstVetoIso_;
vector<float>  ophoPFChWorstIso_;
vector<float>  ophoEcalPFClusterIso_;
vector<float>  ophoHcalPFClusterIso_;
vector<float>  ophoPFClusEcalIso_;
vector<float>  ophoPFClusHcalIso_;

//vector<float>  ophoSeedBCE_;
//vector<float>  ophoSeedBCEta_;
vector<float>  ophoIDMVA_;
vector<ULong64_t> ophoFiredSingleTrgs_;
vector<ULong64_t> ophoFiredDoubleTrgs_;
vector<ULong64_t> ophoFiredTripleTrgs_;
vector<ULong64_t> ophoFiredL1Trgs_;
vector<float>  ophoSeedTime_;
vector<float>  ophoSeedEnergy_;
vector<float>  ophoMIPTotEnergy_;
//vector<float>  ophoSeedTimeFull5x5_;
vector<float>  ophoMIPChi2_;
vector<float>  ophoMIPSlope_;
vector<float>  ophoMIPIntercept_;
vector<float>  ophoMIPNhitCone_;
vector<float>  ophoMIPIsHalo_;
vector<UShort_t> ophoxtalBits_;
vector<UShort_t> ophoIDbit_;
vector<float>    ophoScale_stat_up_;
vector<float>    ophoScale_stat_dn_;
vector<float>    ophoScale_syst_up_;
vector<float>    ophoScale_syst_dn_;
vector<float>    ophoScale_gain_up_;
vector<float>    ophoScale_gain_dn_;
vector<float>    ophoResol_rho_up_;
vector<float>    ophoResol_rho_dn_;
vector<float>    ophoResol_phi_up_;
vector<float>    ophoResol_phi_dn_;


Int_t          nPho_;
vector<float>  phoE_;
vector<float>  phoSigmaE_;
vector<float>  phoEt_;
vector<float>  phoEta_;
vector<float>  phoPhi_;
vector<float>  phoCalibE_;
vector<float>  phoCalibEt_;
vector<float>  phoSCE_;
vector<float>  phoSCRawE_;
vector<float>  phoESEnP1_;
vector<float>  phoESEnP2_;
vector<float>  phoSCEta_;
vector<float>  phoSCPhi_;
vector<float>  phoSCEtaWidth_;
vector<float>  phoSCPhiWidth_;
vector<float>  phoSCBrem_;
vector<int>    phohasPixelSeed_;
vector<int>    phoEleVeto_;
vector<float>  phoR9_;
vector<float>  phoHoverE_;
vector<float>  phoESEffSigmaRR_;
vector<float>  phoSigmaIEtaIEtaFull5x5_;
vector<float>  phoSigmaIEtaIPhiFull5x5_;
vector<float>  phoSigmaIPhiIPhiFull5x5_;
vector<float>  phoE2x2Full5x5_;
vector<float>  phoE5x5Full5x5_;
vector<float>  phoR9Full5x5_;
vector<float>  phoPFChIso_;
vector<float>  phoPFChPVIso_;
vector<float>  phoPFPhoIso_;
vector<float>  phoPFNeuIso_;
vector<float>  phoPFChWorstVetoIso_;
vector<float>  phoPFChWorstIso_;
vector<float>  phoEcalPFClusterIso_;
vector<float>  phoHcalPFClusterIso_;
//vector<float>  phoSeedBCE_;
//vector<float>  phoSeedBCEta_;
vector<float>  phoIDMVA_;
vector<ULong64_t> phoFiredSingleTrgs_;
vector<ULong64_t> phoFiredDoubleTrgs_;
vector<ULong64_t> phoFiredTripleTrgs_;
vector<ULong64_t> phoFiredL1Trgs_;
vector<float>  phoSeedTime_;
vector<float>  phoSeedEnergy_;
vector<float>  phoMIPTotEnergy_;
//vector<float>  phoSeedTimeFull5x5_;
vector<float>  phoMIPChi2_;
vector<float>  phoMIPSlope_;
vector<float>  phoMIPIntercept_;
vector<float>  phoMIPNhitCone_;
vector<float>  phoMIPIsHalo_;
vector<UShort_t> phoxtalBits_;
vector<UShort_t> phoIDbit_;
vector<float>    phoScale_stat_up_;
vector<float>    phoScale_stat_dn_;
vector<float>    phoScale_syst_up_;
vector<float>    phoScale_syst_dn_;
vector<float>    phoScale_gain_up_;
vector<float>    phoScale_gain_dn_;
vector<float>    phoResol_rho_up_;
vector<float>    phoResol_rho_dn_;
vector<float>    phoResol_phi_up_;
vector<float>    phoResol_phi_dn_;

//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;
    
  }
  return false;
}


//crystal info
void ggNtuplizer::branchesCrystals(TTree* tree) {

    //pfphoton for crsystal from cluster
    tree->Branch("nPFPhoton", &nPFPhoton, "nPFPhoton/I");
    tree->Branch("PFPhotonE", PFPhotonE, "PFPhotonE[nPFPhoton]/F");
    tree->Branch("PFPhotonHadEM", PFPhotonHadEM, "PFPhotonHadEM[nPFPhoton]/F");
    tree->Branch("PFPhotonEta", PFPhotonEta, "PFPhotonEta[nPFPhoton]/F");
    tree->Branch("PFPhotonPhi",PFPhotonPhi, "PFPhotonPhi[nPFPhoton]/F");
    tree->Branch("PFPhotonSigmaIetaIeta", PFPhotonSigmaIetaIeta,"PFPhotonSigmaIetaIeta[nPFPhoton]/F");
    tree->Branch("PFPhotonR9", PFPhotonR9, "PFPhotonR9[nPFPhoton]/F");
    tree->Branch("PFChargedHadronIso", PFChargedHadronIso, "PFChargedHadronIso[nPFPhoton]/F");
    tree->Branch("PFNeutralHadronIso", PFNeutralHadronIso, "PFNeutralHadronIso[nPFPhoton]/F");
    tree->Branch("PFPhotonIso", PFPhotonIso, "PFPhotonIso[nPFPhoton]/F");

    tree->Branch("PFMIPTotE", PFMIPTotE, "PFMIPTotE[nPFPhoton]/F");
    tree->Branch("PFPhotonET", PFPhotonET, "PFPhotonET[nPFPhoton]/F");
    tree->Branch("PFMIPNhit", PFMIPNhit, "PFMIPNhit[nPFPhoton]/I");
    tree->Branch("PFEcalIso03", PFEcalIso03, "PFEcalIso03[nPFPhoton]/F");
    tree->Branch("PFHcalIso03", PFHcalIso03, "PFHcalIso03[nPFPhoton]/F");
    tree->Branch("PFTrkIso03", PFTrkIso03, "PFTrkIso03[nPFPhoton]/F");
    tree->Branch("PFPixMatch", PFPixMatch, "PFPixMatch[nPFPhoton]/I");


    //detid info
    tree->Branch("nAllCellsEB", &nAllCellsEB, "nAllCellsEB/I");
    tree->Branch("AllCellsIEtaEB", &AllCellsIEtaEB, "AllCellsIEtaEB[nAllCellsEB]/I");
    tree->Branch("AllCellsIPhiEB", &AllCellsIPhiEB, "AllCellsIPhiEB[nAllCellsEB]/I");  
    tree->Branch("AllCellsE_EB", &AllCellsE_EB, "AllCellsE_EB[nAllCellsEB]/F");
    tree->Branch("AllTimeEB", &AllTimeEB, "AllFracEB[nAllCellsEB]/F");
    tree->Branch("AllClusteredEB", AllClusteredEB, "AllClusteredEB[nAllCellsEB][30]/I");
    tree->Branch("AllFracEB", AllFracEB, "AllFracEB[nAllCellsEB][30]/F");

}



void ggNtuplizer::branchesPhotons(TTree* tree) {

    

    //seed geo in detid
    tree->Branch("phoSeedIEta",             &phoSeedIEta_);
    tree->Branch("phoSeedIPhi",             &phoSeedIPhi_);
    tree->Branch("ophoSeedIEta",            &ophoSeedIEta_);
    tree->Branch("ophoSeedIPhi",            &ophoSeedIPhi_);

    //out-of-time
    //in-time
    tree->Branch("onPho",                    &onPho_);
    tree->Branch("ophoE",                    &ophoE_);
    tree->Branch("ophoSigmaE",               &ophoSigmaE_);
    tree->Branch("ophoEt",                   &ophoEt_);
    tree->Branch("ophoEta",                  &ophoEta_);
    tree->Branch("ophoPhi",                  &ophoPhi_);
    tree->Branch("ophoCalibE",               &ophoCalibE_);
    tree->Branch("ophoCalibEt",              &ophoCalibEt_);
    tree->Branch("ophoSCE",                  &ophoSCE_);
    tree->Branch("ophoSCRawE",               &ophoSCRawE_);
    tree->Branch("ophoESEnP1",               &ophoESEnP1_);
    tree->Branch("ophoESEnP2",               &ophoESEnP2_);
    tree->Branch("ophoSCEta",                &ophoSCEta_);
    tree->Branch("ophoSCPhi",                &ophoSCPhi_);
    tree->Branch("ophoSCEtaWidth",           &ophoSCEtaWidth_);
    tree->Branch("ophoSCPhiWidth",           &ophoSCPhiWidth_);
    tree->Branch("ophoSCBrem",               &ophoSCBrem_);
    tree->Branch("ophohasPixelSeed",         &ophohasPixelSeed_);
    tree->Branch("ophoEleVeto",              &ophoEleVeto_);
    tree->Branch("ophoR9",                   &ophoR9_);
    tree->Branch("ophoHoverE",               &ophoHoverE_);
    tree->Branch("ophoESEffSigmaRR",         &ophoESEffSigmaRR_);
    tree->Branch("ophoSigmaIEtaIEtaFull5x5", &ophoSigmaIEtaIEtaFull5x5_);
    tree->Branch("ophoSigmaIEtaIPhiFull5x5", &ophoSigmaIEtaIPhiFull5x5_);
    tree->Branch("ophoSigmaIPhiIPhiFull5x5", &ophoSigmaIPhiIPhiFull5x5_);
    tree->Branch("ophoE2x2Full5x5",          &ophoE2x2Full5x5_);
    tree->Branch("ophoE5x5Full5x5",          &ophoE5x5Full5x5_);
    tree->Branch("ophoR9Full5x5",            &ophoR9Full5x5_);
    //tree->Branch("ophoSeedBCE",              &ophoSeedBCE_);
    //tree->Branch("ophoSeedBCEta",            &ophoSeedBCEta_);
    tree->Branch("ophoPFChIso",              &ophoPFChIso_);
    tree->Branch("ophoPFChPVIso",            &ophoPFChPVIso_);
    tree->Branch("ophoPFPhoIso",             &ophoPFPhoIso_);
    tree->Branch("ophoPFNeuIso",             &ophoPFNeuIso_);

    tree->Branch("ophoPFClusEcalIso",        &ophoPFClusEcalIso_);
    tree->Branch("ophoPFClusHcalIso",        &ophoPFClusHcalIso_);

    tree->Branch("ophoPFChWorstIso",         &ophoPFChWorstIso_);
    tree->Branch("ophoPFChWorstVetoIso",     &ophoPFChWorstVetoIso_);
    tree->Branch("ophoEcalPFClusterIso",     &ophoEcalPFClusterIso_);
    tree->Branch("ophoHcalPFClusterIso",     &ophoHcalPFClusterIso_);
    tree->Branch("ophoIDMVA",                &ophoIDMVA_);
    tree->Branch("ophoFiredSingleTrgs",      &ophoFiredSingleTrgs_);
    tree->Branch("ophoFiredDoubleTrgs",      &ophoFiredDoubleTrgs_);
    tree->Branch("ophoFiredTripleTrgs",      &ophoFiredTripleTrgs_);
    tree->Branch("ophoFiredL1Trgs",          &ophoFiredL1Trgs_);
    tree->Branch("ophoSeedTime",             &ophoSeedTime_);
    tree->Branch("ophoSeedEnergy",           &ophoSeedEnergy_);
    tree->Branch("ophoMIPTotEnergy",         &ophoMIPTotEnergy_);
    //tree->Branch("ophoSeedTimeFull5x5",              &ophoSeedTimeFull5x5_);
    tree->Branch("ophoMIPChi2",                      &ophoMIPChi2_);
    tree->Branch("ophoMIPSlope",                     &ophoMIPSlope_);
    tree->Branch("ophoMIPIntercept",                 &ophoMIPIntercept_);
    tree->Branch("ophoMIPNhitCone",                  &ophoMIPNhitCone_);
    tree->Branch("ophoMIPIsHalo",                    &ophoMIPIsHalo_);

    tree->Branch("ophoxtalBits",      &ophoxtalBits_);
    tree->Branch("ophoIDbit",         &ophoIDbit_);
    tree->Branch("ophoScale_stat_up", &ophoScale_stat_up_);
    tree->Branch("ophoScale_stat_dn", &ophoScale_stat_dn_);
    tree->Branch("ophoScale_syst_up", &ophoScale_syst_up_);
    tree->Branch("ophoScale_syst_dn", &ophoScale_syst_dn_);
    tree->Branch("ophoScale_gain_up", &ophoScale_gain_up_);
    tree->Branch("ophoScale_gain_dn", &ophoScale_gain_dn_);
    tree->Branch("ophoResol_rho_up",  &ophoResol_rho_up_);
    tree->Branch("ophoResol_rho_dn",  &ophoResol_rho_dn_);
    tree->Branch("ophoResol_phi_up",  &ophoResol_phi_up_);
    tree->Branch("ophoResol_phi_dn",  &ophoResol_phi_dn_);



    //in-time
    tree->Branch("nPho",                    &nPho_);
    tree->Branch("phoE",                    &phoE_);
    tree->Branch("phoSigmaE",               &phoSigmaE_);
    tree->Branch("phoEt",                   &phoEt_);
    tree->Branch("phoEta",                  &phoEta_);
    tree->Branch("phoPhi",                  &phoPhi_);
    tree->Branch("phoCalibE",               &phoCalibE_);
    tree->Branch("phoCalibEt",              &phoCalibEt_);
    tree->Branch("phoSCE",                  &phoSCE_);
    tree->Branch("phoSCRawE",               &phoSCRawE_);
    tree->Branch("phoESEnP1",               &phoESEnP1_);
    tree->Branch("phoESEnP2",               &phoESEnP2_);
    tree->Branch("phoSCEta",                &phoSCEta_);
    tree->Branch("phoSCPhi",                &phoSCPhi_);
    tree->Branch("phoSCEtaWidth",           &phoSCEtaWidth_);
    tree->Branch("phoSCPhiWidth",           &phoSCPhiWidth_);
    tree->Branch("phoSCBrem",               &phoSCBrem_);
    tree->Branch("phohasPixelSeed",         &phohasPixelSeed_);
    tree->Branch("phoEleVeto",              &phoEleVeto_);
    tree->Branch("phoR9",                   &phoR9_);
    tree->Branch("phoHoverE",               &phoHoverE_);
    tree->Branch("phoESEffSigmaRR",         &phoESEffSigmaRR_);
    tree->Branch("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5_);
    tree->Branch("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5_);
    tree->Branch("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5_);
    tree->Branch("phoE2x2Full5x5",          &phoE2x2Full5x5_);
    tree->Branch("phoE5x5Full5x5",          &phoE5x5Full5x5_);
    tree->Branch("phoR9Full5x5",            &phoR9Full5x5_);
    //tree->Branch("phoSeedBCE",              &phoSeedBCE_);
    //tree->Branch("phoSeedBCEta",            &phoSeedBCEta_);
    tree->Branch("phoPFChIso",              &phoPFChIso_);
    tree->Branch("phoPFChPVIso",            &phoPFChPVIso_);
    tree->Branch("phoPFPhoIso",             &phoPFPhoIso_);
    tree->Branch("phoPFNeuIso",             &phoPFNeuIso_);
    tree->Branch("phoPFChWorstIso",         &phoPFChWorstIso_);
    tree->Branch("phoPFChWorstVetoIso",     &phoPFChWorstVetoIso_);
    tree->Branch("phoEcalPFClusterIso",     &phoEcalPFClusterIso_);
    tree->Branch("phoHcalPFClusterIso",     &phoHcalPFClusterIso_);
    tree->Branch("phoIDMVA",                &phoIDMVA_);
    tree->Branch("phoFiredSingleTrgs",      &phoFiredSingleTrgs_);
    tree->Branch("phoFiredDoubleTrgs",      &phoFiredDoubleTrgs_);
    tree->Branch("phoFiredTripleTrgs",      &phoFiredTripleTrgs_);
    tree->Branch("phoFiredL1Trgs",          &phoFiredL1Trgs_);
    tree->Branch("phoSeedTime",             &phoSeedTime_);
    tree->Branch("phoSeedEnergy",           &phoSeedEnergy_);
    tree->Branch("phoMIPTotEnergy",         &phoMIPTotEnergy_);
    //tree->Branch("phoSeedTimeFull5x5",              &phoSeedTimeFull5x5_);
    tree->Branch("phoMIPChi2",                      &phoMIPChi2_);
    tree->Branch("phoMIPSlope",                     &phoMIPSlope_);
    tree->Branch("phoMIPIntercept",                 &phoMIPIntercept_);
    tree->Branch("phoMIPNhitCone",                  &phoMIPNhitCone_);
    tree->Branch("phoMIPIsHalo",                    &phoMIPIsHalo_);

    tree->Branch("phoxtalBits",      &phoxtalBits_);
    tree->Branch("phoIDbit",         &phoIDbit_);
    tree->Branch("phoScale_stat_up", &phoScale_stat_up_);
    tree->Branch("phoScale_stat_dn", &phoScale_stat_dn_);
    tree->Branch("phoScale_syst_up", &phoScale_syst_up_);
    tree->Branch("phoScale_syst_dn", &phoScale_syst_dn_);
    tree->Branch("phoScale_gain_up", &phoScale_gain_up_);
    tree->Branch("phoScale_gain_dn", &phoScale_gain_dn_);
    tree->Branch("phoResol_rho_up",  &phoResol_rho_up_);
    tree->Branch("phoResol_rho_dn",  &phoResol_rho_dn_);
    tree->Branch("phoResol_phi_up",  &phoResol_phi_up_);
    tree->Branch("phoResol_phi_dn",  &phoResol_phi_dn_);

}


void ggNtuplizer::fillTracks(const edm::Event& e, const edm::EventSetup& es) {
    /*
    //track info from untime.cc
    edm::Handle<reco::TrackCollection> trackHandle;
    //edm::Handle<std::vector<reco::Track> > trackHandle;
    e.getByToken(tracklabel_, trackHandle);
    const reco::TrackCollection *tracks = trackHandle.product();
    reco::TrackCollection::const_iterator trItr = tracks->begin();
    nAllTracks=0;

    //Fill track bank.
    for(trItr = tracks->begin(); trItr!= tracks->end() && nAllTracks<10000; ++trItr){
        if (trItr->pt()>1.){
            AllTracksPt[nAllTracks]=trItr->pt();
            // if (trItr->pt()>25) n20Trk++;
            AllTracksPhi[nAllTracks]=trItr->phi();
            AllTracksEta[nAllTracks]=trItr->eta();
            AllTracksNhit[nAllTracks]=trItr->numberOfValidHits();
            AllTracksChi2[nAllTracks]=trItr->normalizedChi2();
            AllTracksPtErr[nAllTracks]=trItr->ptError();
            AllTracksQ[nAllTracks]=trItr->charge();
            //AllTracksDxy[nAllTracks]=trItr->dxy(math::XYZPoint(vertexBeamSpot.x0(
            //),vertexBeamSpot.y0(),vertexBeamSpot.z0()));
            //AllTracksZ0[nAllTracks]=trItr->dz(math::XYZPoint(vertexBeamSpot.x0(),
            //                                               vertexBeamSpot.y0(),vertexBeamSpot.z0()));
            AllTracksZ0[nAllTracks]=trItr->vz();
            nAllTracks++;
        }
    }
  
*/

}



void ggNtuplizer::fillCrystals(const edm::Event& e, const edm::EventSetup& es) {

    //cluster crystal from untime.cc
    /*
    edm::Handle<reco::Photon> recoPhotonHandle;
    e.getByToken(recophotonCollection_, recoPhotonHandle);
    */
    edm::Handle<edm::View<pat::Photon> > photonHandle;
    e.getByToken(photonCollection_, photonHandle);
    const edm::View<pat::Photon> *mphotons = photonHandle.product();
    edm::View<pat::Photon>::const_iterator photonit = mphotons->begin();

    std::map<DetId, vector<std::pair<Int_t, Float_t> > > crysclusEB;
    nPFPhoton = 0;
    for(photonit = mphotons->begin(); photonit!= mphotons->end()&&nPFPhoton<100; ++photonit){
        std::vector<std::pair<DetId, float> > clusdet = photonit->superCluster()->hitsAndFractions();
        for (int iii=0;iii<int(clusdet.size());++iii){
        DetId blarg = clusdet[iii].first;
        
        std::map<DetId, vector< pair<Int_t, Float_t> > >::iterator selcrys = crysclusEB.find(blarg); 
        
        if (selcrys == crysclusEB.end()){
        vector< std::pair<Int_t, Float_t> > blerg;
        blerg.push_back( make_pair(nPFPhoton, clusdet[iii].second));
        crysclusEB.insert(make_pair(blarg,blerg));
        }
        else{
        selcrys->second.push_back(make_pair(nPFPhoton, clusdet[iii].second));
        }
    }

    PFPhotonHadEM[nPFPhoton] = photonit->hadTowOverEm();
    PFPhotonE[nPFPhoton]=photonit->energy();
    PFPhotonET[nPFPhoton]=photonit->pt();
    PFMIPTotE[nPFPhoton]=photonit->mipTotEnergy();
    PFMIPNhit[nPFPhoton]=photonit->mipNhitCone();
    //if(nPFPhoton<20)cout << "PhotonE["<< nPFPhoton << "]: " << PFPhotonE[nPFPhoton] << endl; 
    PFPhotonEta[nPFPhoton]=photonit->eta();
    //if(nPFPhoton<20)cout << "PhotonEta["<< nPFPhoton << "]: " << PFPhotonEta[nPFPhoton] << endl; 
    PFPhotonPhi[nPFPhoton]=photonit->phi();
    //if(nPFPhoton<20)cout << "PhotonPhi["<< nPFPhoton << "]: " << PFPhotonPhi[nPFPhoton] << endl; 
    PFPhotonSigmaIetaIeta[nPFPhoton]=photonit->full5x5_sigmaIetaIeta();
    //if(nPFPhoton<20)cout << "PhotonSigmaietaIeta["<< nPFPhoton << "]: " << PFPhotonSigmaIetaIeta[nPFPhoton] << endl; 
    PFPhotonR9[nPFPhoton]=photonit->r9();
    // isolator03_.fGetIsolation(&*photonit,pfH.product(),primVtxRef,vtxH);
    PFChargedHadronIso[nPFPhoton]= photonit->chargedHadronIso();
    //if(nPFPhoton<20)cout << "ChargedIso["<< nPFPhoton << "]: " << PFChargedHadronIso[nPFPhoton] << endl;
    PFNeutralHadronIso[nPFPhoton]= photonit->neutralHadronIso();
    //if(nPFPhoton<20)cout << "NeutralIso["<< nPFPhoton << "]: " << PFNeutralHadronIso[nPFPhoton] << endl;
    PFPhotonIso[nPFPhoton]= photonit->photonIso();
    //if(nPFPhoton<20)cout << "PhotonIso["<< nPFPhoton << "]: " << PFPhotonIso[nPFPhoton] << endl; 
    PFEcalIso03[nPFPhoton] = photonit->ecalRecHitSumEtConeDR03();
    PFHcalIso03[nPFPhoton] = photonit->hcalTowerSumEtConeDR03();
    PFTrkIso03[nPFPhoton] = photonit->trkSumPtHollowConeDR03();
    PFPixMatch[nPFPhoton] = photonit->hasPixelSeed();
    nPFPhoton++;
    }



    //getting the info of detid
    edm::Handle<EcalRecHitCollection> ecalhitsCollHEB;
    e.getByToken(ebReducedRecHitCollection_, ecalhitsCollHEB);
    const EcalRecHitCollection* rechitsCollectionEB_ = ecalhitsCollHEB.product();

    

    nAllCellsEB=0;
    for (EcalRecHitCollection::const_iterator it = rechitsCollectionEB_->begin();it!=rechitsCollectionEB_->end() && nAllCellsEB < 30000;++it)
    {
        DetId blarg = it->detid();
        //Here's where the map is used to put what photon this crystal goes with. 
        std::map<DetId, vector< pair<Int_t, Float_t> > >::const_iterator selcrys = crysclusEB.find(blarg);
        if (selcrys==crysclusEB.end())
        {
            AllClusteredEB[nAllCellsEB][0]=-1;
        }
        else{
            vector< pair<Int_t, Float_t> > blerg = selcrys->second;
            for (int ikk=0;ikk<int(blerg.size()) && ikk<30;ikk++)
            {
                AllClusteredEB[nAllCellsEB][ikk]=blerg[ikk].first;
                AllFracEB[nAllCellsEB][ikk]=blerg[ikk].second;
            }
        }
        
        //array
        EBDetId dit = it->detid();
        AllCellsIEtaEB[nAllCellsEB]=dit.ieta();
        AllCellsIPhiEB[nAllCellsEB]=dit.iphi();
        AllCellsE_EB[nAllCellsEB]=it->energy();
        AllTimeEB[nAllCellsEB]=it->time();
        

        nAllCellsEB++;

    }

}


void ggNtuplizer::fillPhotons(const edm::Event& e, const edm::EventSetup& es) {
  
    ophoSeedIEta_.clear();
    ophoSeedIPhi_.clear();


    // cleanup from previous execution
    ophoE_                   .clear();
    ophoSigmaE_              .clear();
    ophoEt_                  .clear();
    ophoEta_                 .clear();
    ophoPhi_                 .clear();
    ophoCalibE_              .clear();
    ophoCalibEt_             .clear();
    ophoSCE_                 .clear();
    ophoSCRawE_              .clear();
    ophoESEnP1_              .clear();
    ophoESEnP2_              .clear();
    ophoSCEta_               .clear();
    ophoSCPhi_               .clear();
    ophoSCEtaWidth_          .clear();
    ophoSCPhiWidth_          .clear();
    ophoSCBrem_              .clear();
    ophohasPixelSeed_        .clear();
    ophoEleVeto_             .clear();
    ophoR9_                  .clear();
    ophoHoverE_              .clear();
    ophoESEffSigmaRR_        .clear();
    ophoSigmaIEtaIEtaFull5x5_.clear();
    ophoSigmaIEtaIPhiFull5x5_.clear();
    ophoSigmaIPhiIPhiFull5x5_.clear();
    ophoE2x2Full5x5_         .clear();
    ophoE5x5Full5x5_         .clear();
    ophoR9Full5x5_           .clear();
    ophoPFChIso_             .clear();
    ophoPFChPVIso_           .clear();
    ophoPFPhoIso_            .clear();
    ophoPFNeuIso_            .clear();
    ophoPFChWorstVetoIso_    .clear();
    ophoPFChWorstIso_        .clear();
    ophoEcalPFClusterIso_    .clear();
    ophoHcalPFClusterIso_    .clear();
    ophoPFClusEcalIso_       .clear();
    ophoPFClusHcalIso_       .clear();
    //ophoSeedBCE_           .clear();
    //ophoSeedBCEta_         .clear();
    ophoIDMVA_               .clear();
    ophoFiredSingleTrgs_     .clear();
    ophoFiredDoubleTrgs_     .clear();
    ophoFiredTripleTrgs_     .clear();
    ophoFiredL1Trgs_         .clear();
    ophoxtalBits_            .clear();
    ophoSeedTime_            .clear();
    ophoSeedEnergy_          .clear();
    ophoMIPTotEnergy_        .clear();  
    //ophoSeedTimeFull5x5_   .clear();
    ophoMIPChi2_           .clear();
    ophoMIPSlope_          .clear();
    ophoMIPIntercept_      .clear();
    ophoMIPNhitCone_       .clear();
    ophoMIPIsHalo_         .clear();

    ophoIDbit_        .clear();
    ophoScale_stat_up_.clear();
    ophoScale_stat_dn_.clear();
    ophoScale_syst_up_.clear();
    ophoScale_syst_dn_.clear();
    ophoScale_gain_up_.clear();
    ophoScale_gain_dn_.clear();
    ophoResol_rho_up_ .clear();
    ophoResol_rho_dn_ .clear();
    ophoResol_phi_up_ .clear();
    ophoResol_phi_dn_ .clear();  

    onPho_ = 0;

    edm::Handle<EcalRecHitCollection> ecalhitsCollHEB;
    e.getByToken(ebReducedRecHitCollection_, ecalhitsCollHEB);


    //added oPhotonToken in ggNtuplizer.h
    edm::Handle<edm::View<pat::Photon> > OOTPhotonsH;
    e.getByToken(OOTphotonCollection_, OOTPhotonsH);

    if (!OOTPhotonsH.isValid()) {
        edm::LogWarning("ggNtuplizer") << "no pat::Photons in event";
        return;
    }

    /*
    edm::Handle<reco::PhotonCollection> recoPhotonHandle;
    e.getByToken(recophotonCollection_, recoPhotonHandle);
    */

    EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
    noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

    for (edm::View<pat::Photon>::const_iterator iPho = OOTPhotonsH->begin(); iPho != OOTPhotonsH->end(); ++iPho) {

        ophoE_               .push_back(iPho->energy());
        //ophoCalibE_          .push_back(iPho->userFloat("ecalEnergyPostCorr"));
        ophoEt_              .push_back(iPho->et());
        //ophoCalibEt_         .push_back(iPho->et()*iPho->userFloat("ecalEnergyPostCorr")/iPho->energy());
        //ophoSigmaE_          .push_back(iPho->userFloat("ecalEnergyErrPostCorr"));
        ophoEta_             .push_back(iPho->eta());
        ophoPhi_             .push_back(iPho->phi());
        ophoSCE_             .push_back((*iPho).superCluster()->energy());
        ophoSCRawE_          .push_back((*iPho).superCluster()->rawEnergy());
        ophoESEnP1_          .push_back((*iPho).superCluster()->preshowerEnergyPlane1());
        ophoESEnP2_          .push_back((*iPho).superCluster()->preshowerEnergyPlane2());
        ophoSCEta_           .push_back((*iPho).superCluster()->eta());
        ophoSCPhi_           .push_back((*iPho).superCluster()->phi());
        ophoSCEtaWidth_      .push_back((*iPho).superCluster()->etaWidth());
        ophoSCPhiWidth_      .push_back((*iPho).superCluster()->phiWidth());
        ophoSCBrem_          .push_back((*iPho).superCluster()->phiWidth()/(*iPho).superCluster()->etaWidth());
        ophohasPixelSeed_    .push_back((Int_t)iPho->hasPixelSeed());
        ophoEleVeto_         .push_back((Int_t)iPho->passElectronVeto());
        ophoR9_              .push_back(iPho->r9());
        ophoHoverE_          .push_back(iPho->hadTowOverEm());
        ophoESEffSigmaRR_    .push_back(lazyTool.eseffsirir(*((*iPho).superCluster())));
        
        ophoPFClusEcalIso_   .push_back(iPho->ecalPFClusterIso());
        ophoPFClusHcalIso_   .push_back(iPho->hcalPFClusterIso());
        ophoPFChIso_         .push_back(iPho->chargedHadronIso()); 
        ophoPFChPVIso_       .push_back(iPho->chargedHadronPFPVIso()); 
        ophoPFPhoIso_        .push_back(iPho->photonIso());
        ophoPFNeuIso_        .push_back(iPho->neutralHadronIso());
        ophoPFChWorstIso_    .push_back(iPho->chargedHadronWorstVtxIso()); 
        ophoPFChWorstVetoIso_.push_back(iPho->chargedHadronWorstVtxGeomVetoIso());
        ophoEcalPFClusterIso_.push_back(iPho->ecalPFClusterIso());
        ophoHcalPFClusterIso_.push_back(iPho->hcalPFClusterIso());
        //ophoIDMVA_           .push_back(iPho->userFloat("PhotonMVAEstimatorRunIIFall17v2Values"));  

        /*
        // VID decisions     
        UShort_t tmpphoIDbit = 0;        
        bool isPassLoose  = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
        if (isPassLoose)  setbit(tmpphoIDbit, 0);   
        bool isPassMedium = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-medium");
        if (isPassMedium) setbit(tmpphoIDbit, 1);    
        bool isPassTight  = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-tight");
        if (isPassTight)  setbit(tmpphoIDbit, 2);
        
        ophoIDbit_.push_back(tmpphoIDbit);      

        // systematics for energy scale and resolution
        
        ophoScale_stat_up_.push_back(iPho->userFloat("energyScaleStatUp"));
        ophoScale_stat_dn_.push_back(iPho->userFloat("energyScaleStatDown"));
        ophoScale_syst_up_.push_back(iPho->userFloat("energyScaleSystUp"));
        ophoScale_syst_dn_.push_back(iPho->userFloat("energyScaleSystDown"));
        ophoScale_gain_up_.push_back(iPho->userFloat("energyScaleGainUp"));
        ophoScale_gain_dn_.push_back(iPho->userFloat("energyScaleGainDown"));
        ophoResol_rho_up_ .push_back(iPho->userFloat("energySigmaRhoUp"));
        ophoResol_rho_dn_ .push_back(iPho->userFloat("energySigmaRhoDown"));
        ophoResol_phi_up_ .push_back(iPho->userFloat("energySigmaPhiUp"));
        ophoResol_phi_dn_ .push_back(iPho->userFloat("energySigmaPhiDown"));
        */
        ///////////////////////////////SATURATED/UNSATURATED ///from ggFlash////
        DetId seed = (iPho->superCluster()->seed()->hitsAndFractions())[0].first;
        bool isBarrel = seed.subdetId() == EcalBarrel;
        const EcalRecHitCollection * rechits = (isBarrel?lazyTool.getEcalEBRecHitCollection():lazyTool.getEcalEERecHitCollection());
                
        EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
        if (theSeedHit != rechits->end()) {
        //std::cout<<"(*theSeedHit).time()"<<(*theSeedHit).time()<<"seed energy: "<<(*theSeedHit).energy()<<std::endl;  

        ophoSeedTime_  .push_back((*theSeedHit).time());
        ophoSeedEnergy_.push_back((*theSeedHit).energy());

        EBDetId dit   = theSeedHit->detid();
        ophoSeedIEta_  .push_back(dit.ieta());
        ophoSeedIPhi_  .push_back(dit.iphi());

        } else{
        ophoSeedTime_  .push_back(-99.);
        ophoSeedEnergy_.push_back(-99.);

        ophoSeedIEta_  .push_back(-99.);
        ophoSeedIPhi_  .push_back(-99.);
        }

        /*
        unsigned short nSaturated = 0, nLeRecovered = 0, nNeighRecovered = 0, nGain1 = 0, nGain6 = 0, nWeired = 0;
        int isSaturated       = 0;
        int isSaturated_gain6 = 0;
        
        UShort_t tmpxtalbit = 0;

        auto matrix5x5 = lazyTool.matrixDetId(seed,-2,+2,-2,+2);
        for (auto & deId : matrix5x5 ) {
        /// cout << "matrix " << deId.rawId() << endl;
        auto rh = rechits->find(deId);
        if( rh != rechits->end() ) {
        nSaturated += rh->checkFlag( EcalRecHit::kSaturated );
        nLeRecovered += rh->checkFlag( EcalRecHit::kLeadingEdgeRecovered );
        nNeighRecovered += rh->checkFlag( EcalRecHit::kNeighboursRecovered );
        nGain1 += rh->checkFlag( EcalRecHit::kHasSwitchToGain1 );
        nGain6 += rh->checkFlag( EcalRecHit::kHasSwitchToGain6 );
        nWeired += rh->checkFlag( EcalRecHit::kWeird ) || rh->checkFlag( EcalRecHit::kDiWeird );
        
        if( rh->checkFlag( EcalRecHit::kHasSwitchToGain1 ) && rh->checkFlag( EcalRecHit::kSaturated ) && !isSaturated){ //this is to fill only once, i.e. only if xtal has this, no need to check for other xtals

        setbit(tmpxtalbit, 0);
        isSaturated = 1;
        //break;
        }
        
        if( rh->checkFlag( EcalRecHit::kHasSwitchToGain6 ) && rh->checkFlag( EcalRecHit::kSaturated ) && !isSaturated_gain6){ //this is to fill only once, i.e. only if xtal has this, no need to check for other xtals

        setbit(tmpxtalbit, 1);
        isSaturated_gain6 = 1;
        //break;
        }
        
        }//if( rh != rechits->end() ) 
        
        if (nWeired>0) setbit(tmpxtalbit,2);      
        if (nGain6>0) setbit(tmpxtalbit,3); 

        }//for(auto & deId : matrix5x5 )
    
        ophoxtalBits_.push_back(tmpxtalbit);
        */
        ophoFiredSingleTrgs_     .push_back(matchSinglePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
        ophoFiredDoubleTrgs_     .push_back(matchDoublePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
        ophoFiredTripleTrgs_     .push_back(matchTriplePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
        ophoFiredL1Trgs_         .push_back(matchL1TriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));

        ophoSigmaIEtaIEtaFull5x5_ .push_back(iPho->full5x5_sigmaIetaIeta());
        ophoSigmaIEtaIPhiFull5x5_ .push_back(iPho->full5x5_showerShapeVariables().sigmaIetaIphi);
        ophoSigmaIPhiIPhiFull5x5_ .push_back(iPho->full5x5_showerShapeVariables().sigmaIphiIphi);
        ophoE2x2Full5x5_          .push_back(lazyToolnoZS.e2x2(*((*iPho).superCluster()->seed())));
        ophoE5x5Full5x5_          .push_back(iPho->full5x5_e5x5());
        ophoR9Full5x5_            .push_back(iPho->full5x5_r9());
        ophoMIPTotEnergy_         .push_back(iPho->mipTotEnergy());

        //phoSeedBCE_        .push_back((*iPho).superCluster()->seed()->energy());
        //phoSeedBCEta_      .push_back((*iPho).superCluster()->seed()->eta());
        
        //phoSeedTimeFull5x5_.push_back(lazyToolnoZS.SuperClusterSeedTime(*((*iPho).superCluster())));
        ophoMIPChi2_        .push_back(iPho->mipChi2());
        ophoMIPSlope_       .push_back(iPho->mipSlope());
        ophoMIPIntercept_   .push_back(iPho->mipIntercept());
        ophoMIPNhitCone_    .push_back(iPho->mipNhitCone());
        ophoMIPIsHalo_      .push_back(iPho->mipIsHalo());    
        
        onPho_++;
    }


    phoSeedIEta_.clear();
    phoSeedIPhi_.clear();

    // cleanup from previous execution
    phoE_                   .clear();
    phoSigmaE_              .clear();
    phoEt_                  .clear();
    phoEta_                 .clear();
    phoPhi_                 .clear();
    phoCalibE_              .clear();
    phoCalibEt_             .clear();
    phoSCE_                 .clear();
    phoSCRawE_              .clear();
    phoESEnP1_              .clear();
    phoESEnP2_              .clear();
    phoSCEta_               .clear();
    phoSCPhi_               .clear();
    phoSCEtaWidth_          .clear();
    phoSCPhiWidth_          .clear();
    phoSCBrem_              .clear();
    phohasPixelSeed_        .clear();
    phoEleVeto_             .clear();
    phoR9_                  .clear();
    phoHoverE_              .clear();
    phoESEffSigmaRR_        .clear();
    phoSigmaIEtaIEtaFull5x5_.clear();
    phoSigmaIEtaIPhiFull5x5_.clear();
    phoSigmaIPhiIPhiFull5x5_.clear();
    phoE2x2Full5x5_         .clear();
    phoE5x5Full5x5_         .clear();
    phoR9Full5x5_           .clear();
    phoPFChIso_             .clear();
    phoPFChPVIso_           .clear();
    phoPFPhoIso_            .clear();
    phoPFNeuIso_            .clear();
    phoPFChWorstVetoIso_    .clear();
    phoPFChWorstIso_        .clear();
    phoEcalPFClusterIso_    .clear();
    phoHcalPFClusterIso_    .clear();
    //phoSeedBCE_           .clear();
    //phoSeedBCEta_         .clear();
    phoIDMVA_               .clear();
    phoFiredSingleTrgs_     .clear();
    phoFiredDoubleTrgs_     .clear();
    phoFiredTripleTrgs_     .clear();
    phoFiredL1Trgs_         .clear();
    phoxtalBits_            .clear();
    phoSeedTime_            .clear();
    phoSeedEnergy_          .clear();
    phoMIPTotEnergy_        .clear();  
    //phoSeedTimeFull5x5_   .clear();
    phoMIPChi2_           .clear();
    phoMIPSlope_          .clear();
    phoMIPIntercept_      .clear();
    phoMIPNhitCone_       .clear();
    phoMIPIsHalo_         .clear();

    phoIDbit_        .clear();
    phoScale_stat_up_.clear();
    phoScale_stat_dn_.clear();
    phoScale_syst_up_.clear();
    phoScale_syst_dn_.clear();
    phoScale_gain_up_.clear();
    phoScale_gain_dn_.clear();
    phoResol_rho_up_ .clear();
    phoResol_rho_dn_ .clear();
    phoResol_phi_up_ .clear();
    phoResol_phi_dn_ .clear();  

    nPho_ = 0;

    edm::Handle<edm::View<pat::Photon> > photonHandle;
    e.getByToken(photonCollection_, photonHandle);

    if (!photonHandle.isValid()) {
        edm::LogWarning("ggNtuplizer") << "no pat::Photons in event";
        return;
    }

    edm::Handle<reco::PhotonCollection> recoPhotonHandle;
    e.getByToken(recophotonCollection_, recoPhotonHandle);
    /*
    EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
    noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
    */
    for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho) {

        phoE_               .push_back(iPho->energy());
        phoCalibE_          .push_back(iPho->userFloat("ecalEnergyPostCorr"));
        phoEt_              .push_back(iPho->et());
        phoCalibEt_         .push_back(iPho->et()*iPho->userFloat("ecalEnergyPostCorr")/iPho->energy());
        phoSigmaE_          .push_back(iPho->userFloat("ecalEnergyErrPostCorr"));
        phoEta_             .push_back(iPho->eta());
        phoPhi_             .push_back(iPho->phi());
        phoSCE_             .push_back((*iPho).superCluster()->energy());
        phoSCRawE_          .push_back((*iPho).superCluster()->rawEnergy());
        phoESEnP1_          .push_back((*iPho).superCluster()->preshowerEnergyPlane1());
        phoESEnP2_          .push_back((*iPho).superCluster()->preshowerEnergyPlane2());
        phoSCEta_           .push_back((*iPho).superCluster()->eta());
        phoSCPhi_           .push_back((*iPho).superCluster()->phi());
        phoSCEtaWidth_      .push_back((*iPho).superCluster()->etaWidth());
        phoSCPhiWidth_      .push_back((*iPho).superCluster()->phiWidth());
        phoSCBrem_          .push_back((*iPho).superCluster()->phiWidth()/(*iPho).superCluster()->etaWidth());
        phohasPixelSeed_    .push_back((Int_t)iPho->hasPixelSeed());
        phoEleVeto_         .push_back((Int_t)iPho->passElectronVeto());
        phoR9_              .push_back(iPho->r9());
        phoHoverE_          .push_back(iPho->hadTowOverEm());
        phoESEffSigmaRR_    .push_back(lazyTool.eseffsirir(*((*iPho).superCluster())));
        phoPFChIso_         .push_back(iPho->chargedHadronIso()); 
        phoPFChPVIso_       .push_back(iPho->chargedHadronPFPVIso()); 
        phoPFPhoIso_        .push_back(iPho->photonIso());
        phoPFNeuIso_        .push_back(iPho->neutralHadronIso());
        phoPFChWorstIso_    .push_back(iPho->chargedHadronWorstVtxIso()); 
        phoPFChWorstVetoIso_.push_back(iPho->chargedHadronWorstVtxGeomVetoIso());
        phoEcalPFClusterIso_.push_back(iPho->ecalPFClusterIso());
        phoHcalPFClusterIso_.push_back(iPho->hcalPFClusterIso());
        phoIDMVA_           .push_back(iPho->userFloat("PhotonMVAEstimatorRunIIFall17v2Values"));  

        // VID decisions     
        UShort_t tmpphoIDbit = 0;        
        bool isPassLoose  = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
        if (isPassLoose)  setbit(tmpphoIDbit, 0);   
        bool isPassMedium = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-medium");
        if (isPassMedium) setbit(tmpphoIDbit, 1);    
        bool isPassTight  = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-tight");
        if (isPassTight)  setbit(tmpphoIDbit, 2);
        
        phoIDbit_.push_back(tmpphoIDbit);      

        // systematics for energy scale and resolution
        /*
        phoScale_stat_up_.push_back(iPho->userFloat("energyScaleStatUp"));
        phoScale_stat_dn_.push_back(iPho->userFloat("energyScaleStatDown"));
        phoScale_syst_up_.push_back(iPho->userFloat("energyScaleSystUp"));
        phoScale_syst_dn_.push_back(iPho->userFloat("energyScaleSystDown"));
        phoScale_gain_up_.push_back(iPho->userFloat("energyScaleGainUp"));
        phoScale_gain_dn_.push_back(iPho->userFloat("energyScaleGainDown"));
        phoResol_rho_up_ .push_back(iPho->userFloat("energySigmaRhoUp"));
        phoResol_rho_dn_ .push_back(iPho->userFloat("energySigmaRhoDown"));
        phoResol_phi_up_ .push_back(iPho->userFloat("energySigmaPhiUp"));
        phoResol_phi_dn_ .push_back(iPho->userFloat("energySigmaPhiDown"));
        */
        ///////////////////////////////SATURATED/UNSATURATED ///from ggFlash////
        DetId seed = (iPho->superCluster()->seed()->hitsAndFractions())[0].first;
        bool isBarrel = seed.subdetId() == EcalBarrel;
        const EcalRecHitCollection * rechits = (isBarrel?lazyTool.getEcalEBRecHitCollection():lazyTool.getEcalEERecHitCollection());
                
        EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
        if (theSeedHit != rechits->end()) {
        //std::cout<<"(*theSeedHit).time()"<<(*theSeedHit).time()<<"seed energy: "<<(*theSeedHit).energy()<<std::endl;  

        phoSeedTime_  .push_back((*theSeedHit).time());
        phoSeedEnergy_.push_back((*theSeedHit).energy());

        EBDetId dit   = theSeedHit->detid();
        phoSeedIEta_  .push_back(dit.ieta());
        phoSeedIPhi_  .push_back(dit.iphi());

        } else{
        phoSeedTime_  .push_back(-99.);
        phoSeedEnergy_.push_back(-99.);

        phoSeedIEta_  .push_back(-99.);
        phoSeedIPhi_  .push_back(-99.);
        }
        /*
        unsigned short nSaturated = 0, nLeRecovered = 0, nNeighRecovered = 0, nGain1 = 0, nGain6 = 0, nWeired = 0;
        int isSaturated       = 0;
        int isSaturated_gain6 = 0;
        
        UShort_t tmpxtalbit = 0;

        auto matrix5x5 = lazyTool.matrixDetId(seed,-2,+2,-2,+2);
        for (auto & deId : matrix5x5 ) {
        /// cout << "matrix " << deId.rawId() << endl;
        auto rh = rechits->find(deId);
        if( rh != rechits->end() ) {
        nSaturated += rh->checkFlag( EcalRecHit::kSaturated );
        nLeRecovered += rh->checkFlag( EcalRecHit::kLeadingEdgeRecovered );
        nNeighRecovered += rh->checkFlag( EcalRecHit::kNeighboursRecovered );
        nGain1 += rh->checkFlag( EcalRecHit::kHasSwitchToGain1 );
        nGain6 += rh->checkFlag( EcalRecHit::kHasSwitchToGain6 );
        nWeired += rh->checkFlag( EcalRecHit::kWeird ) || rh->checkFlag( EcalRecHit::kDiWeird );
        
        if( rh->checkFlag( EcalRecHit::kHasSwitchToGain1 ) && rh->checkFlag( EcalRecHit::kSaturated ) && !isSaturated){ //this is to fill only once, i.e. only if xtal has this, no need to check for other xtals

        setbit(tmpxtalbit, 0);
        isSaturated = 1;
        //break;
        }
        
        if( rh->checkFlag( EcalRecHit::kHasSwitchToGain6 ) && rh->checkFlag( EcalRecHit::kSaturated ) && !isSaturated_gain6){ //this is to fill only once, i.e. only if xtal has this, no need to check for other xtals

        setbit(tmpxtalbit, 1);
        isSaturated_gain6 = 1;
        //break;
        }
        
        }//if( rh != rechits->end() ) 
        
        if (nWeired>0) setbit(tmpxtalbit,2);      
        if (nGain6>0) setbit(tmpxtalbit,3); 

        }//for(auto & deId : matrix5x5 )
    
        phoxtalBits_.push_back(tmpxtalbit);
        */
        phoFiredSingleTrgs_     .push_back(matchSinglePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
        phoFiredDoubleTrgs_     .push_back(matchDoublePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
        phoFiredTripleTrgs_     .push_back(matchTriplePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
        phoFiredL1Trgs_         .push_back(matchL1TriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));

        phoSigmaIEtaIEtaFull5x5_ .push_back(iPho->full5x5_sigmaIetaIeta());
        phoSigmaIEtaIPhiFull5x5_ .push_back(iPho->full5x5_showerShapeVariables().sigmaIetaIphi);
        phoSigmaIPhiIPhiFull5x5_ .push_back(iPho->full5x5_showerShapeVariables().sigmaIphiIphi);
        phoE2x2Full5x5_          .push_back(lazyToolnoZS.e2x2(*((*iPho).superCluster()->seed())));
        phoE5x5Full5x5_          .push_back(iPho->full5x5_e5x5());
        phoR9Full5x5_            .push_back(iPho->full5x5_r9());
        phoMIPTotEnergy_         .push_back(iPho->mipTotEnergy());

        //phoSeedBCE_        .push_back((*iPho).superCluster()->seed()->energy());
        //phoSeedBCEta_      .push_back((*iPho).superCluster()->seed()->eta());
        
        //phoSeedTimeFull5x5_.push_back(lazyToolnoZS.SuperClusterSeedTime(*((*iPho).superCluster())));
        phoMIPChi2_        .push_back(iPho->mipChi2());
        phoMIPSlope_       .push_back(iPho->mipSlope());
        phoMIPIntercept_   .push_back(iPho->mipIntercept());
        phoMIPNhitCone_    .push_back(iPho->mipNhitCone());
        phoMIPIsHalo_      .push_back(iPho->mipIsHalo());    
        
        nPho_++;
    }

}

void ggNtuplizer::cleanupPhotons() {

}
