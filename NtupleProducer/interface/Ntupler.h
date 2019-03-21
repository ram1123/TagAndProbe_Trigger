#ifndef Ntupler_h
#define Ntupler_h

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//Electron related header files
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

//Vertex related header files
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// some common header files
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// Gen particles related header files
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// Muon related header files
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// c++ and root related header files
#include <regex>
#include "TString.h" 
#include "TTree.h"
#include "Math/VectorUtil.h"
#include <string> 
#include <algorithm>

// HLT related header files
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

// L1 related header files
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

//
// class declaration
//

class Ntupler : public edm::EDAnalyzer {
   public:
      explicit Ntupler(const edm::ParameterSet&);
      ~Ntupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      enum ElectronMatchType {UNMATCHED = 0, 
			  TRUE_PROMPT_ELECTRON, 
			  TRUE_ELECTRON_FROM_TAU,
			  TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      int matchToTruth(const edm::Ptr<reco::GsfElectron> el,  
		       const edm::Handle<edm::View<reco::GenParticle>>  &prunedGenParticles);
      void findFirstNonElectronMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);
      bool hasWZasMother(const reco::GenParticle  p)  ;
      static bool cmd(const reco::GenParticle & s1,const reco::GenParticle & s2);

  	// ----------member data ---------------------------
      // Data members that are the same for AOD and miniAOD
      edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
      edm::EDGetTokenT<GenEventInfoProduct> genEventInfoProduct_;

      // AOD case data members
      edm::EDGetToken electronsToken_;
      edm::EDGetToken muonsToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;

      // Trigger Tokens
      edm::EDGetToken triggerResultsToken_;
      edm::EDGetToken triggerSummaryToken_;
      edm::EDGetToken l1EGtoken_;
      edm::EDGetToken l1Muontoken_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales>triggerPrescale_; 
      HLTConfigProvider hltConfig; 
      std::vector<std::string> pathsToSave_;
      std::vector<std::string> filterToMatch_;
      std::string HLTprocess_;

      // MiniAOD case data members
      edm::EDGetToken electronsMiniAODToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxMiniAODToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsMiniAODToken_;
      edm::EDGetToken muonsMiniAODToken_;
  
     // VID decisions objects
     edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapLooseToken_;
     edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapMediumToken_;
     edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapTightToken_;
     edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapMVAnoIsoWP90Token_;
     edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapMVAnoIsoWP80Token_;
     edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapMVAnoIsoWPLooseToken_;
     edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapMVAIsoWP90Token_;
     edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapMVAIsoWP80Token_;
     edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapMVAIsoWPLooseToken_;
     edm::EDGetTokenT<edm::ValueMap<float> > eleMVAValuesMapTokenIso_;
     edm::EDGetTokenT<edm::ValueMap<float> > eleMVAValuesMapTokenNoIso_;

      // L1 Token
      edm::EDGetToken muToken;
  
     // Verbose output for ID
     bool isMC_;
     bool doEle_;
     bool doMuon_;

 
     uint32_t run_;
   
     // Vars for PVs
     Int_t pvNTracks_;
   
     // Vars for weight (can be negative)
     double genWeight_;
   
     // Vars for pile-up
     Int_t nPUTrue_;    // true pile-up
     Int_t nPU_;        // generated pile-up
     Int_t good_vertices_;        // generated pile-up
     Int_t nPV_;        // number of reconsrtucted primary vertices
     Float_t rho_;      // the rho variable
   
  // Trigger names and decision 
     std::vector<bool> triggerDecision;
     std::vector<std::string> triggerPath;
   
     // All Electron filters and variables
   
     std::vector<bool> passFilterEle32;
     std::vector<bool> passFilterEle35;
     std::vector<bool> passFilterEle23_12_leg1;
     std::vector<bool> passFilterEle23_12_leg2;
     std::vector<bool> passFilterMu12_Ele23_legEle;
     std::vector<bool> passFilterMu23_Ele12_legEle;
     std::vector<bool> L1EG_35 ;
     std::vector<bool> L1EG_23_12 ;
     std::vector<bool> passEleIdLoose_;
     std::vector<bool> passEleIdMedium_;
     std::vector<bool> passEleIdTight_;
     std::vector<bool> passMVAnoIsoWP90_;
     std::vector<bool> passMVAnoIsoWP80_;
     std::vector<bool> passMVAnoIsoWPLoose_;
     std::vector<bool> passMVAIsoWP90_;
     std::vector<bool> passMVAIsoWP80_;
     std::vector<bool> passMVAIsoWPLoose_;
     std::vector<float> valueMVAnoIso_;
     std::vector<float> valueMVAIso_;
        
     // all electron variables

     Int_t nElectrons_;
     std::vector<double> ele_pt_;
     std::vector<double> ele_eta_;
     std::vector<double> ele_etaSC_;
     std::vector<double> ele_phi_;
     std::vector<double> ele_phiSC_;
     std::vector<double> ele_energy_;
     std::vector<double> ele_energySC_;
     std::vector<double> ele_charge_;
     std::vector<double> ele_dEtaIn_;
     std::vector<double> ele_dEtaSeed_;
     std::vector<double> ele_dPhiIn_;
     std::vector<double> ele_hOverE_;
     std::vector<double> ele_full5x5_sigmaIetaIeta_;
     std::vector<double> ele_isoChargedHadrons_;
     std::vector<double> ele_isoNeutralHadrons_;
     std::vector<double> ele_isoPhotons_;
     std::vector<double> ele_isoChargedFromPU_;
     std::vector<double> ele_relCombIsoWithEA_;
     std::vector<double> ele_ooEmooP_;
     std::vector<double> ele_d0_;
     std::vector<double> ele_dz_;
     std::vector<double> ele_SIP_;
     std::vector<double> ele_dr03TkSumPt_;
     std::vector<double> ele_dr03EcalRecHitSumEt_;
     std::vector<double> ele_dr03HcalDepth1TowerSumEt_;
     std::vector<Int_t>  ele_expectedMissingInnerHits_;
     std::vector<Int_t>  ele_passConversionVeto_;     
     EffectiveAreas   effectiveAreas_;
   
     std::vector<bool> hasMatchedToZ;
 

    // All Muon variables 
     Int_t nMuons_;
     std::vector<double> mu_pt_;
     std::vector<double> mu_eta_;
     std::vector<double> mu_phi_;
     std::vector<double> mu_energy_;
     std::vector<int> mu_charge_;
     std::vector<int> mu_type_;
     std::vector<double> mu_d0_;
     std::vector<double> mu_dz_;
     std::vector<double> mu_SIP_;
     std::vector<double> mu_Chi2NDF_;
     std::vector<double> mu_InnerD0_;
     std::vector<double> mu_InnerDz_;
     std::vector<int>      mu_TrkLayers_;
     std::vector<int>      mu_PixelLayers_;
     std::vector<int>      mu_PixelHits_;
     std::vector<int>      mu_MuonHits_;
     std::vector<int>      mu_Stations_;
     std::vector<int>      mu_Matches_;
     std::vector<int>      mu_TrkQuality_;
     std::vector<double>    mu_IsoTrk_;
     std::vector<double>    mu_PFChIso_;
     std::vector<double>    mu_PFPhoIso_;
     std::vector<double>    mu_PFNeuIso_;
     std::vector<double>    mu_PFPUIso_;
     std::vector<double>    mu_PFChIso03_;
     std::vector<double>    mu_PFPhoIso03_;
     std::vector<double>    mu_PFNeuIso03_;
     std::vector<double>    mu_PFPUIso03_;
     std::vector<double>    mu_InnervalidFraction_;
     std::vector<double>    mu_segmentCompatibility_;
     std::vector<double>    mu_chi2LocalPosition_;
     std::vector<double>    mu_trkKink_;
     std::vector<double>    mu_BestTrkPtError_;
     std::vector<double>    mu_BestTrkPt_;
     std::vector<int>      mu_BestTrkType_;

     std::vector<int>      mu_CutBasedIdLoose_;
     std::vector<int>      mu_CutBasedIdMedium_;
     std::vector<int>      mu_CutBasedIdTight_;
     std::vector<int>      mu_CutBasedIdMediumPrompt_;
     std::vector<int>      mu_CutBasedIdGlobalHighPt_;
     std::vector<int>      mu_CutBasedIdTrkHighPt_;
     std::vector<int>      mu_PFIsoVeryLoose_;
     std::vector<int>      mu_PFIsoLoose_;
     std::vector<int>      mu_PFIsoMedium_;
     std::vector<int>      mu_PFIsoTight_;
     std::vector<int>      mu_PFIsoVeryTight_;
     std::vector<int>      mu_PFIsoVeryVeryTight_;
     std::vector<int>      mu_TrkIsoLoose_;
     std::vector<int>      mu_TrkIsoTight_;
     std::vector<int>      mu_SoftCutBasedId_;
     std::vector<int>      mu_MvaLoose_;
     std::vector<int>      mu_MvaMedium_;
     std::vector<int>      mu_MvaTight_;
     std::vector<int>      mu_MiniIsoLoose_;
     std::vector<int>      mu_MiniIsoMedium_;
     std::vector<int>      mu_MiniIsoTight_;
     std::vector<int>      mu_MiniIsoVeryTight_;
     std::vector<int>      mu_TriggerIdLoose_;
     std::vector<int>      mu_InTimeMuon_;
     std::vector<int>      mu_MultiIsoLoose_;
     std::vector<int>      mu_MultiIsoMedium_;
     std::vector<bool> passFilterIsoMu24;
     std::vector<bool> passFilterIsoMu27;
     std::vector<bool> passFilterMu50;
     std::vector<bool> passFilterMu17_Mu8_leg1;
     std::vector<bool> passFilterMu17_Mu8_leg2;
     std::vector<bool> passFilterMu17_Mu8_leg2_wL1;
     std::vector<bool> passFilterMu17_Mu8_IsoLeg;
     std::vector<bool> passFilterMu12_Ele23_legMu;
     std::vector<bool> passFilterMu12_Ele23_legMu_L10p5;
     std::vector<bool> passFilterMu12_Ele23_legMu_L10p3;
     std::vector<bool> passFilterMu12_HLTOnly;
     std::vector<bool> passFilterMu12_L1TOnly;
     std::vector<bool> passFilterMu23_Ele12_legMu;

     //Gen particles
     unsigned int   genParticles_n;
     std::vector<double> genElectron_pt;
     std::vector<double> genElectron_eta;
     std::vector<double> genElectron_phi;
     std::vector<double> genElectron_energy;
     std::vector<bool>   genElectron_fromZ;

     std::vector<double> genMuon_pt;
     std::vector<double> genMuon_eta;
     std::vector<double> genMuon_phi;
     std::vector<double> genMuon_energy;
     std::vector<bool>   genMuon_fromZ;

    // Tree decleration 
     TTree *tree_;
   
   };

#endif

