#include "TagAndProbe_Trigger/NtupleProducer/interface/Ntupler.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Ntupler::Ntupler(const edm::ParameterSet& iConfig):
    pathsToSave_(iConfig.getParameter<std::vector<std::string>>("pathsToSave" )),
    filterToMatch_(iConfig.getParameter<std::vector<std::string>>("filterToMatch" )),
    HLTprocess_(iConfig.getParameter<std::string>("HLTprocess" )),
    eleIdMapLooseToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapLoose"))),
    eleIdMapMediumToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapMedium"))),
    eleIdMapTightToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapTight"))),
    eleIdMapMVAnoIsoWP90Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVA90noIso"))),
    eleIdMapMVAnoIsoWP80Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVA80noIso"))),
    eleIdMapMVAIsoWP90Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVA90Iso"))),
    eleIdMapMVAIsoWP80Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVA80Iso"))),
    isMC_(iConfig.getParameter<bool>("isMC")),
    doEle_(iConfig.getParameter<bool>("doEle")),
    doMuon_(iConfig.getParameter<bool>("doMuon")),
    effectiveAreas_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath() )
{
    
    //
    // Prepare tokens for all input collections and objects
    //
    
    // Universal tokens for AOD and miniAOD
    genEventInfoProduct_ = consumes<GenEventInfoProduct>
    (iConfig.getParameter <edm::InputTag>
     ("genEventInfoProduct"));
    
    pileupToken_ = consumes<edm::View<PileupSummaryInfo> >
    (iConfig.getParameter <edm::InputTag>
     ("pileup"));
    
    rhoToken_    = consumes<double>
    (iConfig.getParameter <edm::InputTag>
     ("rho"));
    
    beamSpotToken_    = consumes<reco::BeamSpot>
    (iConfig.getParameter <edm::InputTag>
     ("beamSpot"));
    
    // AOD tokens
    electronsToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electrons"));
    
    vtxToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("vertices"));
    
    genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));
    
    conversionsToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversions"));
    
    triggerResultsToken_    = mayConsume< edm::TriggerResults >
    (iConfig.getParameter<edm::InputTag>
     ("triggerResultTag"));
    
    triggerSummaryToken_    = mayConsume< trigger::TriggerEvent >
    (iConfig.getParameter<edm::InputTag>
     ("triggerSummaryTag"));
    
    // MiniAOD tokens
    // For electrons, use the fact that pat::Electron can be cast into
    // GsfElectron
    electronsMiniAODToken_    = mayConsume<edm::View<pat::Electron> >
    (iConfig.getParameter<edm::InputTag>
     ("electronsMiniAOD"));
    
    vtxMiniAODToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("verticesMiniAOD"));
    
    genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));
    
    conversionsMiniAODToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversionsMiniAOD"));
    
    // Trigger
    triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>
    (iConfig.getParameter<edm::InputTag>
     ("objects"));
    
    triggerPrescale_ = consumes<pat::PackedTriggerPrescales>
    (iConfig.getParameter<edm::InputTag>
     ("prescale"));
    
    muonsMiniAODToken_    = mayConsume<edm::View<pat::Muon> >
    (iConfig.getParameter<edm::InputTag>
     ("muonsMiniAOD"));
    
    muToken    = consumes<BXVector<l1t::Muon>>(iConfig.getParameter<edm::InputTag>("muInputTag"));
    egToken    = consumes<BXVector<l1t::EGamma>>(iConfig.getParameter<edm::InputTag>("egInputTag"));
    
    //
    // Set up the ntuple structure
    //
    
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree> ("EventTree", "Event data");
    
    tree_->Branch("run",  &run_,  "run/I");
    
    tree_->Branch("pvNTracks", &pvNTracks_ , "pvNTracks/I");
    tree_->Branch("good_vertices",&good_vertices_, "good_vertices/I");
    tree_->Branch("nPV", &nPV_     , "nPV/I");
    tree_->Branch("nPU", &nPU_     , "nPU/I");
    tree_->Branch("nPUTrue", &nPUTrue_ , "nPUTrue/I");
    tree_->Branch("rho", &rho_ , "rho/F");
    
    tree_->Branch("genWeight"    ,  &genWeight_ , "genWeight/F");
    
    tree_->Branch("genParticles_n" ,  &genParticles_n);
    tree_->Branch("genElectron_pt" ,  &genElectron_pt);
    tree_->Branch("genElectron_eta" ,  &genElectron_eta);
    tree_->Branch("genElectron_phi" ,  &genElectron_phi);
    tree_->Branch("genElectron_energy" ,  &genElectron_energy);
    tree_->Branch("genElectron_fromZ" ,  &genElectron_fromZ);
    
    tree_->Branch("nEle"    ,  &nElectrons_ , "nEle/I");
    tree_->Branch("ele_pt"    ,  &ele_pt_    );
    tree_->Branch("ele_eta" ,  &ele_eta_ );
    tree_->Branch("ele_etaSC" ,  &ele_etaSC_ );
    tree_->Branch("ele_phi" ,  &ele_phi_ );
    tree_->Branch("ele_tricharge" ,  &ele_tricharge_ );
    tree_->Branch("ele_phiSC" ,  &ele_phiSC_ );
    tree_->Branch("ele_energy" ,  &ele_energy_ );
    tree_->Branch("ele_energySC" ,  &ele_energySC_ );
    tree_->Branch("ele_charge" ,  &ele_charge_ );
    tree_->Branch("ele_dEtaIn",  &ele_dEtaIn_);
    tree_->Branch("ele_dEtaSeed",  &ele_dEtaSeed_);
    tree_->Branch("ele_dPhiIn",  &ele_dPhiIn_);
    tree_->Branch("ele_hOverE",  &ele_hOverE_);
    tree_->Branch("ele_full5x5_sigmaIetaIeta", &ele_full5x5_sigmaIetaIeta_);
    tree_->Branch("ele_isoChargedHadrons"      , &ele_isoChargedHadrons_);
    tree_->Branch("ele_isoNeutralHadrons"      , &ele_isoNeutralHadrons_);
    tree_->Branch("ele_isoPhotons"             , &ele_isoPhotons_);
    tree_->Branch("ele_relCombIsoWithEA"       , &ele_relCombIsoWithEA_);
    tree_->Branch("ele_isoChargedFromPU"       , &ele_isoChargedFromPU_);
    tree_->Branch("ele_ooEmooP", &ele_ooEmooP_);
    tree_->Branch("ele_dr03TkSumPt", &ele_dr03TkSumPt_);
    tree_->Branch("ele_dr03EcalRecHitSumEt", &ele_dr03EcalRecHitSumEt_);
    tree_->Branch("ele_dr03HcalDepth1TowerSumEt", &ele_dr03HcalDepth1TowerSumEt_);
    tree_->Branch("ele_d0"     , &ele_d0_);
    tree_->Branch("ele_dz"     , &ele_dz_);
    tree_->Branch("ele_SIP"     , &ele_SIP_);
    tree_->Branch("ele_expectedMissingInnerHits", &ele_expectedMissingInnerHits_);
    tree_->Branch("ele_passConversionVeto", &ele_passConversionVeto_);
    
    //  tree_->Branch("isTrue"    , &isTrue_);
    
    tree_->Branch("passEleIdLoose"  ,  &passEleIdLoose_ );
    tree_->Branch("passEleIdMedium"  ,  &passEleIdMedium_ );
    tree_->Branch("passEleIdTight"  ,  &passEleIdTight_ );
    tree_->Branch("passMVAnoIsoWP90"  ,  &passMVAnoIsoWP90_ );
    tree_->Branch("passMVAnoIsoWP80"  ,  &passMVAnoIsoWP80_ );
    tree_->Branch("passMVAIsoWP90"  ,  &passMVAIsoWP90_ );
    tree_->Branch("passMVAIsoWP80"  ,  &passMVAIsoWP80_ );
    
    tree_->Branch("hasMatchedToZ" , &hasMatchedToZ);
    // Electron Trigger branch
    tree_->Branch("passL1EG10", &passL1EG10);
    tree_->Branch("passL1EG17", &passL1EG17);
    tree_->Branch("passL1EG23", &passL1EG23);
    tree_->Branch("passL1EG23Iso", &passL1EG23Iso);
    tree_->Branch("passL1EG20Iso", &passL1EG20Iso);
    tree_->Branch("triggerPath" ,  &triggerPath);
    tree_->Branch("triggerDecision" ,  &triggerDecision);
    tree_->Branch("passFilterEle32"           ,  &passFilterEle32);
    tree_->Branch("passFilterEle23_12_leg1"   ,  &passFilterEle23_12_leg1);
    tree_->Branch("passFilterEle23_12_leg2"   ,  &passFilterEle23_12_leg2);
    tree_->Branch("passFilterEle115"           ,  &passFilterEle115);
    tree_->Branch("passFilterEle50"           ,  &passFilterEle50);
    tree_->Branch("passFilterEle25"           ,  &passFilterEle25);
    tree_->Branch("passFilterEle27"           ,  &passFilterEle27);
    tree_->Branch("passFilterMu12_Ele23_legEle"   ,  &passFilterMu12_Ele23_legEle);
    tree_->Branch("passFilterMu23_Ele12_legEle"   ,  &passFilterMu23_Ele12_legEle);
    
    
    
    // Trigger objects
    
    // Muon variables
    if(doMuon_) {
        tree_->Branch("nMu", &nMuons_ , "nMu/I" );
        tree_->Branch("mu_pt",&mu_pt_);
        tree_->Branch("mu_eta",&mu_eta_);
        tree_->Branch("mu_phi",&mu_phi_);
        tree_->Branch("mu_energy",&mu_energy_);
        tree_->Branch("mu_charge",&mu_charge_);
        tree_->Branch("mu_type",&mu_type_);
        tree_->Branch("mu_d0",&mu_d0_);
        tree_->Branch("mu_dz",&mu_dz_);
        tree_->Branch("mu_SIP",&mu_SIP_);
        tree_->Branch("mu_Chi2NDF",&mu_Chi2NDF_);
        tree_->Branch("mu_InnerD0",&mu_InnerD0_);
        tree_->Branch("mu_InnerDz",&mu_InnerDz_);
        tree_->Branch("mu_TrkLayers",&mu_TrkLayers_);
        tree_->Branch("mu_PixelLayers",&mu_PixelLayers_);
        tree_->Branch("mu_PixelHits",&mu_PixelHits_);
        tree_->Branch("mu_MuonHits",&mu_MuonHits_);
        tree_->Branch("mu_Stations",&mu_Stations_);
        tree_->Branch("mu_Matches",&mu_Matches_);
        tree_->Branch("mu_TrkQuality",&mu_TrkQuality_);
        tree_->Branch("mu_IsoTrk",&mu_IsoTrk_);
        tree_->Branch("mu_PFChIso",&mu_PFChIso_);
        tree_->Branch("mu_PFPhoIso",&mu_PFPhoIso_);
        tree_->Branch("mu_PFNeuIso",&mu_PFNeuIso_);
        tree_->Branch("mu_PFPUIso",&mu_PFPUIso_);
        tree_->Branch("mu_PFCHIso03",&mu_PFChIso03_);
        tree_->Branch("mu_PFPhoIso03",&mu_PFPhoIso03_);
        tree_->Branch("mu_PFNeuIso03",&mu_PFNeuIso03_);
        tree_->Branch("mu_PFPUIso03",&mu_PFPUIso03_);
        tree_->Branch("mu_InnervalidFraction",&mu_InnervalidFraction_);
        tree_->Branch("mu_segmentCompatibility",&mu_segmentCompatibility_);
        tree_->Branch("mu_chi2LocalPosition",&mu_chi2LocalPosition_);
        tree_->Branch("mu_trkKink",&mu_trkKink_);
        tree_->Branch("mu_BestTrkPtError",&mu_BestTrkPtError_);
        tree_->Branch("mu_BestTrkPt",&mu_BestTrkPt_);
        tree_->Branch("mu_BestTrkType",&mu_BestTrkType_);
        tree_->Branch("mu_CutBasedIdLoose",&mu_CutBasedIdLoose_);
        tree_->Branch("mu_CutBasedIdMedium",&mu_CutBasedIdMedium_);
        tree_->Branch("mu_CutBasedIdTight",&mu_CutBasedIdTight_);
        tree_->Branch("mu_CutBasedIdMediumPrompt",&mu_CutBasedIdMediumPrompt_);
        tree_->Branch("mu_CutBasedIdGlobalHighPt",&mu_CutBasedIdGlobalHighPt_);
        tree_->Branch("mu_CutBasedIdTrkHighPt",&mu_CutBasedIdTrkHighPt_);
        tree_->Branch("mu_PFIsoVeryLoose",&mu_PFIsoVeryLoose_);
        tree_->Branch("mu_PFIsoLoose",&mu_PFIsoLoose_);
        tree_->Branch("mu_PFIsoMedium",&mu_PFIsoMedium_);
        tree_->Branch("mu_PFIsoTight",&mu_PFIsoTight_);
        tree_->Branch("mu_PFIsoVeryTight",&mu_PFIsoVeryTight_);
        tree_->Branch("mu_PFIsoVeryVeryTight",&mu_PFIsoVeryVeryTight_);
        tree_->Branch("mu_TrkIsoLoose",&mu_TrkIsoLoose_);
        tree_->Branch("mu_TrkIsoTight",&mu_TrkIsoTight_);
        tree_->Branch("mu_SoftCutBasedId",&mu_SoftCutBasedId_);
        tree_->Branch("mu_MvaLoose",&mu_MvaLoose_);
        tree_->Branch("mu_MvaMedium",&mu_MvaMedium_);
        tree_->Branch("mu_MvaTight",&mu_MvaTight_);
        tree_->Branch("mu_MiniIsoLoose",&mu_MiniIsoLoose_);
        tree_->Branch("mu_MiniIsoMedium",&mu_MiniIsoMedium_);
        tree_->Branch("mu_MiniIsoTight",&mu_MiniIsoTight_);
        tree_->Branch("mu_MiniIsoVeryTight",&mu_MiniIsoVeryTight_);
        tree_->Branch("mu_TriggerIdLoose",&mu_TriggerIdLoose_);
        tree_->Branch("mu_InTimeMuon",&mu_InTimeMuon_);
        tree_->Branch("mu_MultiIsoLoose",&mu_MultiIsoLoose_);
        tree_->Branch("mu_MultiIsoMedium",&mu_MultiIsoMedium_);
        
        tree_->Branch("passFilterIsoMu24" ,  &passFilterIsoMu24);
        tree_->Branch("passFilterMu17_Mu8_leg1" ,  &passFilterMu17_Mu8_leg1);
        tree_->Branch("passFilterMu17_Mu8_leg2" ,  &passFilterMu17_Mu8_leg2);
        tree_->Branch("passFilterMu17_Mu8_IsoLeg" ,  &passFilterMu17_Mu8_IsoLeg);
        tree_->Branch("passFilterMu12_Ele23_legMu" ,  &passFilterMu12_Ele23_legMu);
        tree_->Branch("passFilterMu23_Ele12_legMu" ,  &passFilterMu23_Ele12_legMu);
        
        tree_->Branch("passFilterMu12_Ele23_legMu_L10p5" ,  &passFilterMu12_Ele23_legMu_L10p5);
        tree_->Branch("passFilterMu12_Ele23_legMu_L10p3" ,  &passFilterMu12_Ele23_legMu_L10p3);
        tree_->Branch("passFilterMu23_Ele12_legMu_L10p5" ,  &passFilterMu23_Ele12_legMu_L10p5);
        tree_->Branch("passFilterMu23_Ele12_legMu_L10p3" ,  &passFilterMu23_Ele12_legMu_L10p3);
        tree_->Branch("passFilterMu12_L10p5" ,  &passFilterMu12_L10p5);
        tree_->Branch("passFilterMu12_L10p3" ,  &passFilterMu12_L10p3);
        tree_->Branch("passFilterMu23_L10p5" ,  &passFilterMu23_L10p5);
        tree_->Branch("passFilterMu23_L10p3" ,  &passFilterMu23_L10p3);
        
        
    }
}


Ntupler::~Ntupler()
{ 
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)   
}


//
// member functions
//

// ------------ method called for each event  ------------
void
Ntupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace std;
    using namespace edm;
    using namespace reco;
    
    
    TString ele_filters[9] ={ "hltEle32WPTightGsfTrackIsoFilter",
        "hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",
        "hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",
        "hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter",
        "hltEle50CaloIdVTGsfTrkIdTGsfDphiFilter",
        "hltDiEle25CaloIdLMWPMS2UnseededFilter",
        "hltDiEle27L1DoubleEGWPTightHcalIsoFilter",
        "hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",
        "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"
    };
    
    TString mu_filters[6] = { "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07",
        "hltL3fL1DoubleMu155fPreFiltered8",
        "hltL3fL1DoubleMu155fFiltered17",
        "hltDiMuon178RelTrkIsoFiltered0p4",
        "hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered12",
        "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"
    };
    
    if(isMC_)
    {  // Get gen weight info
        edm::Handle< GenEventInfoProduct > genWeightH;
        iEvent.getByToken(genEventInfoProduct_,genWeightH);
        genWeight_ = genWeightH->GenEventInfoProduct::weight();
        
        // Get Pileup info
        Handle<edm::View<PileupSummaryInfo> > pileupHandle;
        iEvent.getByToken(pileupToken_, pileupHandle);
        for( auto & puInfoElement : *pileupHandle)
        {
            if( puInfoElement.getBunchCrossing() == 0 )
            {
                nPU_    = puInfoElement.getPU_NumInteractions();
                nPUTrue_= puInfoElement.getTrueNumInteractions();
            }
        }
    }
    
    run_ = iEvent.id().run();
    
    // Get rho value
    edm::Handle< double > rhoH;
    iEvent.getByToken(rhoToken_,rhoH);
    rho_ = *rhoH;
    
    // Get the beam spot
    edm::Handle<reco::BeamSpot> theBeamSpot;
    iEvent.getByToken(beamSpotToken_,theBeamSpot);
    
    // Retrieve the collection of electrons from the event.
    // If we fail to retrieve the collection with the standard AOD
    // name, we next look for the one with the stndard miniAOD name.
    // We use exactly the same handle for AOD and miniAOD formats
    // since pat::Electron objects can be recast as reco::GsfElectron objects.
    edm::Handle<edm::View<pat::Electron> > electrons;
    bool isAOD = false;
    iEvent.getByToken(electronsMiniAODToken_, electrons);
    //  if( !electrons.isValid() ){
    //    isAOD = false;
    //    iEvent.getByToken(electronsMiniAODToken_,electrons);
    //  }
    
    // Get the MC collection
    reco::GenParticleCollection genElectrons;
    Handle<edm::View<reco::GenParticle> > genParticles;
    if(isMC_)
    {
        if( isAOD )
            iEvent.getByToken(genParticlesToken_,genParticles);
        else
            iEvent.getByToken(genParticlesMiniAODToken_,genParticles);
    }

    // Get PV
    edm::Handle<reco::VertexCollection> vertices;
    if( isAOD )
        iEvent.getByToken(vtxToken_, vertices);
    else
        iEvent.getByToken(vtxMiniAODToken_, vertices);
    
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &pv = vertices->front();
    nPV_    = vertices->size();
    good_vertices_ = 0;
    if (vertices.isValid())
        if (vertices->size() > 0)
            for (auto v : *vertices)
            {
                if (v.ndof() >= 4 && !v.isFake())
                    ++good_vertices_;
            }
    
    // NOTE FOR RUN 2 THE OLD SELECTION OF GOOD VERTICES BELOW IS DISCOURAGED
    // // Find the first vertex in the collection that passes
    // // good quality criteria
    // VertexCollection::const_iterator firstGoodVertex = vertices->end();
    // int firstGoodVertexIdx = 0;
    // for (VertexCollection::const_iterator vtx = vertices->begin();
    //      vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    //   // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    //   // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    //   bool isFake = vtx->isFake();
    //   if( !isAOD )
    //     isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    //   // Check the goodness
    //   if ( !isFake
    // 	 &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
    // 	 && fabs(vtx->position().Z())<=24.0) {
    //     firstGoodVertex = vtx;
    //     break;
    //   }
    // }
    
    // if ( firstGoodVertex==vertices->end() )
    //   return; // skip event if there are no good PVs
    
    // // Seems always zero. Not stored in miniAOD...?
    // pvNTracks_ = firstGoodVertex->nTracks();
    pvNTracks_ = pv.nTracks();
    
    
    // Get the conversions collection
    edm::Handle<reco::ConversionCollection> conversions;
    if(isAOD)
        iEvent.getByToken(conversionsToken_, conversions);
    else
        iEvent.getByToken(conversionsMiniAODToken_, conversions);
    
    // Get the electron ID data from the event stream.
    // Note: this implies that the VID ID modules have been run upstream.
    // If you need more info, check with the EGM group.
    edm::Handle<edm::ValueMap<bool> > loose_ele_id_decisions;
    edm::Handle<edm::ValueMap<bool> > medium_ele_id_decisions;
    edm::Handle<edm::ValueMap<bool> > tight_ele_id_decisions;
    edm::Handle<edm::ValueMap<bool> > eleMVAnoIsoWP90;
    edm::Handle<edm::ValueMap<bool> > eleMVAnoIsoWP80;
    edm::Handle<edm::ValueMap<bool> > eleMVAIsoWP90;
    edm::Handle<edm::ValueMap<bool> > eleMVAIsoWP80;
    
    iEvent.getByToken(eleIdMapLooseToken_ ,loose_ele_id_decisions);
    iEvent.getByToken(eleIdMapMediumToken_ ,medium_ele_id_decisions);
    iEvent.getByToken(eleIdMapTightToken_ ,tight_ele_id_decisions);
    iEvent.getByToken(eleIdMapMVAnoIsoWP90Token_ ,eleMVAnoIsoWP90);
    iEvent.getByToken(eleIdMapMVAnoIsoWP80Token_ ,eleMVAnoIsoWP80);
    iEvent.getByToken(eleIdMapMVAIsoWP90Token_ ,eleMVAIsoWP90);
    iEvent.getByToken(eleIdMapMVAIsoWP80Token_ ,eleMVAIsoWP80);
    
    // Get Triggers
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByToken(triggerResultsToken_, triggerResults);
    
    triggerPath.clear();
    triggerDecision.clear();
    
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);
    for(unsigned int iPath=0 ; iPath < pathsToSave_.size(); iPath++)
    {
        TString path = pathsToSave_.at(iPath);
        bool trigDec(false);
        size_t j;
        for (j=0; j < triggerResults->size(); j++)
        {
            if (TString(names.triggerName(j)).Contains(path))
            {
                if (triggerResults->accept(j))
                {
                    trigDec = true;
                }
            }
        }
        j=0;
        triggerPath.push_back( path.Data() );
        triggerDecision.push_back( trigDec );
        //	cout<<"path : "<<path<<"    , Decision : "<<triggerDecision[iPath]<<endl;
    }
    
    edm::Handle<trigger::TriggerEvent> triggerSummary;
    if(isAOD) iEvent.getByToken(triggerSummaryToken_, triggerSummary);
    trigger::TriggerObjectCollection allTriggerObjects;
    if(isAOD) allTriggerObjects = triggerSummary->getObjects();
    
    int numberOfFilters = filterToMatch_.size();
    trigger::TriggerObjectCollection *legObjects = new trigger::TriggerObjectCollection[numberOfFilters];
    // find the ref of the legs
    //
    if (isAOD)
    { 
        for (size_t iteFilter=0; iteFilter<filterToMatch_.size(); iteFilter++)
        {
            edm::InputTag filterTag = edm::InputTag(filterToMatch_.at(iteFilter), "", HLTprocess_);
            size_t filterIndex = (*triggerSummary).filterIndex(filterTag);
            if (filterIndex < (*triggerSummary).sizeFilters())
            { //check if the trigger object is present
                //save the trigger objects corresponding to muon leg
                //           cout<<"filterIndex : "<<filterIndex<<"   , filterName :  "<<(*triggerSummary).filterLabel(filterIndex)<<"  , filterTag : "<<filterTag<<endl;
                const trigger::Keys &keys = (*triggerSummary).filterKeys(filterIndex);
                for (size_t j = 0; j < keys.size(); j++)
                {
                    trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
                    legObjects[iteFilter].push_back(foundObject);
                }
            }
        }
    }
    
    Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    if(!isAOD) iEvent.getByToken(triggerObjects_, triggerObjects);
    
    if(!isAOD)
    {
        for (size_t iteFilter=0; iteFilter<filterToMatch_.size(); iteFilter++)
        {
            std::string filterTag = filterToMatch_.at(iteFilter);
            for ( pat::TriggerObjectStandAlone obj: *triggerObjects )
            {
                //        obj.unpackFilterLabels(filterToMatch_);
                //        obj.unpackPathNames(names);
                obj.unpackNamesAndLabels(iEvent,*triggerResults);
                if (obj.hasFilterLabel(filterTag))
                {
                    legObjects[iteFilter].push_back(obj);
                }
            }
        }
    }
    
    // Handle over L1-EG
    Handle<BXVector<l1t::EGamma>> L1EG;
    iEvent.getByToken(egToken,L1EG);
    
    if(doEle_)
    {
        // Loop over electrons
        nElectrons_ = 0;
        ele_pt_.clear();
        ele_etaSC_.clear();
        ele_phiSC_.clear();
        ele_eta_.clear();
        ele_phi_.clear();
        ele_tricharge_.clear();
        ele_energy_.clear();
        ele_energySC_.clear();
        ele_charge_.clear();
        ele_dEtaIn_.clear();
        ele_dEtaSeed_.clear();
        ele_dPhiIn_.clear();
        ele_hOverE_.clear();
        ele_full5x5_sigmaIetaIeta_.clear();
        ele_isoChargedHadrons_.clear();
        ele_isoNeutralHadrons_.clear();
        ele_isoPhotons_.clear();
        ele_relCombIsoWithEA_.clear();
        ele_isoChargedFromPU_.clear();
        ele_ooEmooP_.clear();
        ele_d0_.clear();
        ele_dz_.clear();
        ele_SIP_.clear();
        ele_dr03TkSumPt_.clear();
        ele_dr03EcalRecHitSumEt_.clear();
        ele_dr03HcalDepth1TowerSumEt_.clear();
        ele_expectedMissingInnerHits_.clear();
        ele_passConversionVeto_.clear();
        passEleIdLoose_.clear();
        passEleIdMedium_.clear();
        passEleIdTight_.clear();
        passMVAnoIsoWP90_.clear();
        passMVAnoIsoWP80_.clear();
        passMVAIsoWP90_.clear();
        passMVAIsoWP80_.clear();
        
        passL1EG10 .clear();
        passL1EG17 .clear();
        passL1EG23 .clear();
        passL1EG20Iso .clear();
        passL1EG23Iso .clear();
        passFilterEle32          .clear();
        passFilterEle115          .clear();
        passFilterEle50          .clear();
        passFilterEle27          .clear();
        passFilterEle25          .clear();
        passFilterEle23_12_leg1  .clear();
        passFilterEle23_12_leg2  .clear();
        passFilterMu12_Ele23_legEle.clear();
        passFilterMu23_Ele12_legEle.clear();
        
        for (size_t i = 0; i < electrons->size(); ++i)
        {
            const auto el = electrons->ptrAt(i);
            // for (const pat::Electron &el : *electrons)
            // Kinematics
            
            nElectrons_++;
            ele_pt_.push_back( el->pt() );
            ele_etaSC_.push_back( el->superCluster()->eta() );
            ele_phiSC_.push_back( el->superCluster()->phi() );
            ele_eta_.push_back( el->eta() );
            ele_phi_.push_back( el->phi() );
            ele_tricharge_.push_back( el->chargeInfo().isGsfCtfScPixConsistent );
            ele_energy_.push_back( el->energy() );
            ele_energySC_.push_back( el->superCluster()->energy() );
            ele_charge_.push_back( el->charge() );
            
            // L1 EGamma triggers
            
            float maxL1MatchedNorm = -1;
            float maxL1MatchedIso = -1;
            bool L1EG10(false), L1EG17(false), L1EG23(false), L1EG20Iso(false), L1EG23Iso(false);
            if (L1EG.isValid())
            {
                for(int ibx=L1EG->getFirstBX(); ibx<=L1EG->getLastBX();ibx++)
                {
                    for(std::vector<l1t::EGamma>::const_iterator L1eg = L1EG->begin(ibx); L1eg != L1EG->end(ibx); ++L1eg)
                    {
                        float L1EGPt = L1eg->pt();
                        float L1EGEta = L1eg->eta();
                        float L1EGPhi = L1eg->phi();
                        float L1EGiso = L1eg->hwIso();
                        
                        float delRL1_EG = deltaR(L1EGEta,L1EGPhi ,el->eta(),el->phi());
                        if (delRL1_EG < 0.5)
                        {
                            if(L1eg->pt() > maxL1MatchedNorm) maxL1MatchedNorm = L1eg->pt();
                            if(L1eg->hwIso() == 1 && L1eg->pt()>maxL1MatchedIso) maxL1MatchedIso = L1eg->pt();
                        }
                    }
                }
                if(maxL1MatchedNorm >= 10) L1EG10 = true;
                if(maxL1MatchedNorm >= 17) L1EG17 = true;
                if(maxL1MatchedNorm >= 23) L1EG23 = true;
                if(maxL1MatchedIso >= 20) L1EG20Iso = true;
                if(maxL1MatchedIso >= 23) L1EG23Iso = true;
            }
            
            passL1EG10.push_back(L1EG10);
            passL1EG17.push_back(L1EG17);
            passL1EG23.push_back(L1EG23);
            passL1EG20Iso.push_back(L1EG20Iso);
            passL1EG23Iso.push_back(L1EG23Iso);
            
            // Trigger matching
            bool filterEle32 = false;
            bool filterEle23_12_leg1 = false;
            bool filterEle23_12_leg2 = false;
            bool filterEle115 = false;
            bool filterEle50 = false;
            bool filterEle25 = false;
            bool filterEle27 = false;
            bool filterMu12_Ele23_legEle = false;
            bool filterMu23_Ele12_legEle = false;
            for (unsigned int iteTrigObj = 0 ; iteTrigObj < filterToMatch_.size() ; iteTrigObj++)
            {
                bool foundTheLeg = false;
                TString filter = filterToMatch_.at(iteTrigObj);
                for (unsigned int i = 0 ; i < legObjects[iteTrigObj].size() ; i++)
                {
                    float delR = deltaR(legObjects[iteTrigObj].at(i).eta(), legObjects[iteTrigObj].at(i).phi(),el->superCluster()->eta(),el->superCluster()->phi());
                    
                    if (delR<0.1)
                    {
                        foundTheLeg = true;break;
                    }
                }
                //       cout<<"filter : "<<ele_filters[4].Contain(filter)<<"    foundTheLeg : "<<foundTheLeg<<endl;
                //
                //
                if(ele_filters[0].Contains(filter) && foundTheLeg)  filterEle32 = true;
                if(ele_filters[1].Contains(filter) && foundTheLeg)  filterEle23_12_leg1 = true;
                if(ele_filters[2].Contains(filter) && foundTheLeg)  filterEle23_12_leg2 = true;
                if(ele_filters[3].Contains(filter) && foundTheLeg)  filterEle115 = true;
                if(ele_filters[4].Contains(filter) && foundTheLeg)  filterEle50 = true;
                if(ele_filters[5].Contains(filter) && foundTheLeg)  filterEle25 = true;
                if(ele_filters[6].Contains(filter) && foundTheLeg)  filterEle27 = true;
                if(ele_filters[7].Contains(filter) && foundTheLeg)  filterMu12_Ele23_legEle = true;
                if(ele_filters[8].Contains(filter) && foundTheLeg)  filterMu23_Ele12_legEle = true;
            }
            
            passFilterEle32          .push_back(filterEle32);
            passFilterEle23_12_leg1  .push_back(filterEle23_12_leg1);
            passFilterEle23_12_leg2  .push_back(filterEle23_12_leg2);
            passFilterEle115          .push_back(filterEle115);
            passFilterEle50          .push_back(filterEle50);
            passFilterEle25          .push_back(filterEle25);
            passFilterEle27          .push_back(filterEle27);
            passFilterMu12_Ele23_legEle  .push_back(filterMu12_Ele23_legEle);
            passFilterMu23_Ele12_legEle  .push_back(filterMu23_Ele12_legEle);
            
            // ID and matching
            ele_dEtaIn_.push_back( el->deltaEtaSuperClusterTrackAtVtx() );
            // Calculation of dEtaSeed is taken from VID (by HEEP folks)
            //   https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleDEtaInSeedCut.cc#L31-L32
            float dEtaSeedValue =
                        el->superCluster().isNonnull() && el->superCluster()->seed().isNonnull()
                        ?
                        el->deltaEtaSuperClusterTrackAtVtx()
                        - el->superCluster()->eta()
                        + el->superCluster()->seed()->eta()
                        : std::numeric_limits<float>::max();
            ele_dEtaSeed_.push_back( dEtaSeedValue );
            ele_dPhiIn_.push_back( el->deltaPhiSuperClusterTrackAtVtx() );
            ele_hOverE_.push_back( el->hadronicOverEm() );
            ele_full5x5_sigmaIetaIeta_.push_back( el->full5x5_sigmaIetaIeta() );
            // |1/E-1/p| = |1/E - EoverPinner/E| is computed below
            // The if protects against ecalEnergy == inf or zero
            // (always the case for miniAOD for electrons <5 GeV)
            if( el->ecalEnergy() == 0 )
            {
                //        printf("Electron energy is zero!\n");
                ele_ooEmooP_.push_back( 1e30 );
            } else if ( !std::isfinite(el->ecalEnergy()))
            {
                printf("Electron energy is not finite!\n");
                ele_ooEmooP_.push_back( 1e30 );
            } else
            {
                ele_ooEmooP_.push_back( fabs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() ) );
            }
            
            // Isolation
            GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
            // Compute individual PF isolations
            ele_isoChargedHadrons_.push_back( pfIso.sumChargedHadronPt );
            ele_isoNeutralHadrons_.push_back( pfIso.sumNeutralHadronEt );
            ele_isoPhotons_.push_back( pfIso.sumPhotonEt );
            ele_isoChargedFromPU_.push_back( pfIso.sumPUPt );
            
            // Compute combined relative PF isolation with the effective area correction for pile-up
            float abseta =  abs(el->superCluster()->eta());
            float eA = effectiveAreas_.getEffectiveArea(abseta);
            ele_relCombIsoWithEA_.push_back( ( pfIso.sumChargedHadronPt
                                              + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho_) )
                                            / el->pt() );
            
            // Impact parameter
            reco::GsfTrackRef theTrack = el->gsfTrack();
            ele_d0_.push_back( (-1) * theTrack->dxy(pv.position() ) );
            ele_dz_.push_back( theTrack->dz( pv.position() ) );
            // Conversion rejection
            ele_expectedMissingInnerHits_.push_back(el->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) );
            
            bool passConvVeto = !ConversionTools::hasMatchedConversion(*el,*conversions,theBeamSpot->position());
            
            ele_passConversionVeto_.push_back( (int) passConvVeto );
            //   ele_SIP_.push_back(fabs(el->dB(pat::Electron::PV3D))/el->edB(pat::Electron::PV3D) );
            //      ele_dr03TkSumPt_.push_back(el->dr03TkSumPt() );
            //     ele_dr03EcalRecHitSumEt_.push_back(el-> dr03EcalRecHitSumEt());
            //     ele_dr03HcalDepth1TowerSumEt_.push_back( el-> dr03HcalDepth1TowerSumEt());
            
            //
            // Look up and save the ID decisions
            //
            bool isPassEleIdLoose  = (*loose_ele_id_decisions)[el];
            bool isPassEleIdMedium  = (*medium_ele_id_decisions)[el];
            bool isPassEleIdTight  = (*tight_ele_id_decisions)[el];
            bool isPassMVAnoIsoWP90_  = (*eleMVAnoIsoWP90)[el];
            bool isPassMVAnoIsoWP80_  = (*eleMVAnoIsoWP80)[el];
            bool isPassMVAIsoWP90_  = (*eleMVAIsoWP90)[el];
            bool isPassMVAIsoWP80_  = (*eleMVAIsoWP80)[el];
            passEleIdLoose_.push_back  ( (int)isPassEleIdLoose  );
            passEleIdMedium_.push_back  ( (int)isPassEleIdMedium  );
            passEleIdTight_.push_back  ( (int)isPassEleIdTight  );
            passMVAnoIsoWP90_.push_back  ( (int)isPassMVAnoIsoWP90_ );
            passMVAnoIsoWP80_.push_back  ( (int)isPassMVAnoIsoWP80_ );
            passMVAIsoWP90_.push_back  ( (int)isPassMVAIsoWP90_ );
            passMVAIsoWP80_.push_back  ( (int)isPassMVAIsoWP80_ );
        }
    }

    if(isMC_)
    {
        genElectron_pt.clear();
        genElectron_eta.clear();
        genElectron_phi.clear();
        genElectron_energy.clear();
        genElectron_fromZ.clear();
        //  genMuon_pt.clear();
        //  genMuon_eta.clear();
        //  genMuon_phi.clear();
        //  genMuon_energy.clear();
        //  genMuon_fromZ.clear();
        genParticles_n = genParticles->size();
        for (unsigned int iteGen = 0 ; iteGen < genParticles_n ; iteGen++)
        {
            reco::GenParticle genPart = (*genParticles)[iteGen];
            reco::GenParticle genElectron;
            bool fromZ_ele = false;
            if (genPart.isPromptFinalState() && abs(genPart.pdgId())==11 && genPart.fromHardProcessFinalState())
            {
                // std::cout << "[DEBUG]: mother " << genPart.mother()->pdgId() << std::endl;
                genElectron = genPart;
                if(genElectron.pt()>5)
                {
                    genElectron_pt.push_back(genElectron.pt());
                    genElectron_eta.push_back(genElectron.eta());
                    genElectron_phi.push_back(genElectron.phi());
                    genElectron_energy.push_back(genElectron.energy());
                    // if (hasWZasMother(genElectron))fromZ_ele = true;
                    // We should not explicitly check the mother requirement. After requiring above 
                    // conditons, we are sure that the leptons are from Z-boson. 
                    // if we check mother then for many events it shows the mother of leptons as lepton itself.
                    // https://hypernews.cern.ch/HyperNews/CMS/get/generators/2802/1.html
                    fromZ_ele = true;
                    genElectron_fromZ.push_back(fromZ_ele);
                }
            }
            /*
             reco::GenParticle  genMuon;
             bool fromZ_mu = false;
             if(abs(genPart.pdgId())==13){
             genMuon = genPart;
             if(genMuon.pt()>5) {
             genMuon_pt.push_back(genMuon.pt());
             genMuon_eta.push_back(genMuon.eta());
             genMuon_phi.push_back(genMuon.phi());
             genMuon_energy.push_back(genMuon.energy());
             if (hasWZasMother(genMuon)) fromZ_mu=true;
             genMuon_fromZ.push_back(fromZ_mu);
             }
             }
             */
        }
        
    }
    
    // Muons collection starts
    
    edm::Handle<edm::View<pat::Muon> > muons;
    //  if(isAOD) iEvent.getByToken(muonsToken_,muons);
    iEvent.getByToken(muonsMiniAODToken_,muons);
    
    // Handle over L1-muon
    Handle<BXVector<l1t::Muon>> L1muons;
    iEvent.getByToken(muToken,L1muons);
    
    if(doMuon_)
    {
        nMuons_ = 0;
        mu_pt_.clear();
        mu_eta_.clear();
        mu_phi_.clear();
        mu_energy_.clear();
        mu_charge_.clear();
        mu_type_.clear();
        mu_d0_.clear();
        mu_dz_.clear();
        mu_SIP_.clear();
        mu_Chi2NDF_.clear();
        mu_InnerD0_.clear();
        mu_InnerDz_.clear();
        mu_TrkLayers_.clear();
        mu_PixelLayers_.clear();
        mu_PixelHits_.clear();
        mu_MuonHits_.clear();
        mu_Stations_.clear();
        mu_Matches_.clear();
        mu_TrkQuality_.clear();
        mu_IsoTrk_.clear();
        mu_PFChIso_.clear();
        mu_PFPhoIso_.clear();
        mu_PFNeuIso_.clear();
        mu_PFPUIso_.clear();
        mu_PFChIso03_.clear();
        mu_PFPhoIso03_.clear();
        mu_PFNeuIso03_.clear();
        mu_PFPUIso03_.clear();
        mu_InnervalidFraction_.clear();
        mu_segmentCompatibility_.clear();
        mu_chi2LocalPosition_.clear();
        mu_trkKink_.clear();
        mu_BestTrkPtError_.clear();
        mu_BestTrkPt_.clear();
        mu_BestTrkType_.clear();
        
        mu_CutBasedIdLoose_.clear();
        mu_CutBasedIdMedium_.clear();
        mu_CutBasedIdTight_.clear();
        mu_CutBasedIdMediumPrompt_.clear();
        mu_CutBasedIdGlobalHighPt_.clear();
        mu_CutBasedIdTrkHighPt_.clear();
        mu_PFIsoVeryLoose_.clear();
        mu_PFIsoLoose_.clear();
        mu_PFIsoMedium_.clear();
        mu_PFIsoTight_.clear();
        mu_PFIsoVeryTight_.clear();
        mu_PFIsoVeryVeryTight_.clear();
        mu_TrkIsoLoose_.clear();
        mu_TrkIsoTight_.clear();
        mu_SoftCutBasedId_.clear();
        mu_MvaLoose_.clear();
        mu_MvaMedium_.clear();
        mu_MvaTight_.clear();
        mu_MiniIsoLoose_.clear();
        mu_MiniIsoMedium_.clear();
        mu_MiniIsoTight_.clear();
        mu_MiniIsoVeryTight_.clear();
        mu_TriggerIdLoose_.clear();
        mu_InTimeMuon_.clear();
        mu_MultiIsoLoose_.clear();
        mu_MultiIsoMedium_.clear();
        
        passFilterIsoMu24       .clear();
        passFilterMu17_Mu8_leg1 .clear();
        passFilterMu17_Mu8_leg2 .clear();
        passFilterMu17_Mu8_IsoLeg .clear();
        passFilterMu12_Ele23_legMu.clear();
        passFilterMu23_Ele12_legMu.clear();
        
        passFilterMu12_Ele23_legMu_L10p5.clear();
        passFilterMu12_Ele23_legMu_L10p3.clear();
        passFilterMu23_Ele12_legMu_L10p5.clear();
        passFilterMu23_Ele12_legMu_L10p3.clear();
        passFilterMu12_L10p5.clear();
        passFilterMu12_L10p3.clear();
        passFilterMu23_L10p5.clear();
        passFilterMu23_L10p3.clear();
        
        for(unsigned i=0; i < muons->size();++i )
        {
            const auto mu = muons->ptrAt(i);
            
            bool L1Mu7 = false;
            bool L1Mu23 = false;
            
            float temp1 = 9999.0;
            float temp2 = 9999.0;
            float temp3 = 9999.0;
            float delRL1_7=9999.;
            float delRL1_23=9999.;
            float delRL1_23_L10p5=9999.;
            float delRL1_23_L10p3=9999.;
            float delRL1_7_L10p5=9999.;
            float delRL1_7_L10p3=9999.;
            
            // L1 Matching with check of quality cut
            if (L1muons.isValid())
            {
                for(int i=L1muons->getFirstBX(); i<=L1muons->getLastBX();i++)
                {
                    for(std::vector<l1t::Muon>::const_iterator L1mu = L1muons->begin(i); L1mu != L1muons->end(i); ++L1mu)
                    {
                        float L1MuPt = L1mu->pt();
                        float L1MuEta = L1mu->eta();
                        float L1MuPhi = L1mu->phi();
                        float L1MuQual = L1mu->hwQual();
                        if((L1MuPt >= 5 || L1MuPt >= 7) && (L1MuQual>=12 && L1MuQual<=15))
                        {
                            delRL1_7 = deltaR(L1MuEta,L1MuPhi ,mu->eta(),mu->phi());
                            if (delRL1_7<temp1) temp1 = delRL1_7;
                        }
                        if((L1MuPt >= 20 || L1MuPt >= 23) && (L1MuQual>=12 && L1MuQual<=15))
                        {
                            delRL1_23 = deltaR(L1MuEta,L1MuPhi ,mu->eta(),mu->phi());
                            if (delRL1_23<temp2) temp2 = delRL1_23;
                        }
                        if(temp1<0.3) delRL1_7_L10p3=temp1;
                        if(temp1<0.5) delRL1_7_L10p5=temp1;
                        if(temp2<0.3) delRL1_23_L10p3=temp2;
                        if(temp2<0.5) delRL1_23_L10p5=temp2;
                    }
                }
            }
            
            bool filterIsoMu24 = false;
            bool filterMu17_Mu8_Leg1 = false;
            bool filterMu17_Mu8_Leg2 = false;
            bool filterMu17_Mu8_IsoLeg = false;
            bool filterMu12_Ele23_legMu = false;
            bool filterMu23_Ele12_legMu = false;
            
            bool filterMu23_Ele12_legMu_L10p3 = false;
            bool filterMu23_Ele12_legMu_L10p5 = false;
            bool filterMu12_Ele23_legMu_L10p5 = false;
            bool filterMu12_Ele23_legMu_L10p3 = false;
            bool filterMu12_L1T0p5 = false;
            bool filterMu12_L1T0p3 = false;
            bool filterMu23_L1T0p5 = false;
            bool filterMu23_L1T0p3 = false;
            
            // Trigger matching
            for (unsigned int iteTrigObj = 0 ; iteTrigObj < filterToMatch_.size() ; iteTrigObj++)
            {
                bool foundTheLeg = false;
                TString filter = filterToMatch_.at(iteTrigObj);
                int iPass = 0;
                for (unsigned int i = 0 ; i < legObjects[iteTrigObj].size() ; i++)
                {
                    float delR = deltaR(legObjects[iteTrigObj].at(i).eta(), legObjects[iteTrigObj].at(i).phi(),mu->eta(),mu->phi());
                    if (delR<0.1)
                    {
                        foundTheLeg = true;iPass = i;break;
                    }
                }
                if(mu_filters[0].Contains(filter) && foundTheLeg)  filterIsoMu24 = true;
                // if(mu_filters[1].Contains(filter) && foundTheLeg)  {filterMu17_Mu8_Leg2 = true; }
                if(mu_filters[2].Contains(filter) && foundTheLeg)  {filterMu17_Mu8_Leg1 = true; }
                if(mu_filters[3].Contains(filter) && foundTheLeg)  {filterMu17_Mu8_Leg2 = true; }
                if(mu_filters[3].Contains(filter) && foundTheLeg)  {filterMu17_Mu8_IsoLeg = true; }
                if(mu_filters[3].Contains(filter) && foundTheLeg && legObjects[iteTrigObj].at(iPass).pt()>=12 && delRL1_7_L10p3 < 0.3 )  {filterMu12_Ele23_legMu_L10p3 = true; }
                if(mu_filters[3].Contains(filter) && foundTheLeg && legObjects[iteTrigObj].at(iPass).pt()>=12 && delRL1_7_L10p5 < 0.5 )  {filterMu12_Ele23_legMu_L10p5 = true; }
                if(mu_filters[3].Contains(filter) && foundTheLeg && legObjects[iteTrigObj].at(iPass).pt()>=23 &&  delRL1_23_L10p3 < 0.3 ) {filterMu23_Ele12_legMu_L10p3 = true; }
                if(mu_filters[3].Contains(filter) && foundTheLeg && legObjects[iteTrigObj].at(iPass).pt()>=23 &&  delRL1_23_L10p5 < 0.5 ) {filterMu23_Ele12_legMu_L10p5 = true; }
                
                if(delRL1_7_L10p5 < 0.5 )  {filterMu12_L1T0p5 = true;}
                if(delRL1_7_L10p3 < 0.3 )  {filterMu12_L1T0p3 = true;}
                if(delRL1_23_L10p5 < 0.5 )  {filterMu23_L1T0p5 = true;}
                if(delRL1_23_L10p3 < 0.3 )  {filterMu23_L1T0p3 = true;}
                
            }
            passFilterIsoMu24       .push_back(filterIsoMu24);
            passFilterMu17_Mu8_leg1 .push_back(filterMu17_Mu8_Leg1);
            passFilterMu17_Mu8_leg2 .push_back(filterMu17_Mu8_Leg2);
            passFilterMu17_Mu8_IsoLeg .push_back(filterMu17_Mu8_IsoLeg);
            
            passFilterMu12_Ele23_legMu_L10p5 .push_back(filterMu12_Ele23_legMu_L10p5);
            passFilterMu12_Ele23_legMu_L10p3 .push_back(filterMu12_Ele23_legMu_L10p3);
            passFilterMu23_Ele12_legMu_L10p5 .push_back(filterMu23_Ele12_legMu_L10p5);
            passFilterMu23_Ele12_legMu_L10p3 .push_back(filterMu23_Ele12_legMu_L10p3);
            
            passFilterMu12_L10p5 .push_back(filterMu12_L1T0p5);
            passFilterMu12_L10p3 .push_back(filterMu12_L1T0p3);
            passFilterMu23_L10p5 .push_back(filterMu23_L1T0p5);
            passFilterMu23_L10p3 .push_back(filterMu23_L1T0p3);
            
            nMuons_++;
            
            float muonIso = (mu->pfIsolationR04().sumChargedHadronPt + max(0., mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5*mu->pfIsolationR04().sumPUPt))/mu->pt();
            // Kinematics
            mu_pt_    .push_back(mu->pt());
            mu_energy_.push_back(mu->energy());
            mu_eta_   .push_back(mu->eta());
            mu_phi_   .push_back(mu->phi());
            mu_charge_.push_back(mu->charge());
            mu_d0_    .push_back(mu->muonBestTrack()->dxy(pv.position()));
            mu_dz_    .push_back(mu->muonBestTrack()->dz(pv.position()));
            mu_SIP_   .push_back(fabs(mu->dB(pat::Muon::PV3D))/mu->edB(pat::Muon::PV3D));
            const reco::TrackRef glbmu = mu->globalTrack();
            const reco::TrackRef innmu = mu->innerTrack();
            
            if (glbmu.isNull())
            {
                mu_Chi2NDF_ .push_back(-99.);
                mu_MuonHits_.push_back(-99);
            } else 
            {
                mu_Chi2NDF_.push_back(glbmu->normalizedChi2());
                mu_MuonHits_.push_back(glbmu->hitPattern().numberOfValidMuonHits());
            }
            
            if (innmu.isNull())
            {
                mu_InnerD0_     .push_back(-99.);
                mu_InnerDz_     .push_back(-99.);
            } else {
                mu_InnerD0_     .push_back(innmu->dxy(pv.position()));
                mu_InnerDz_     .push_back(innmu->dz(pv.position()));
            }
            
            mu_CutBasedIdLoose_. 	push_back(mu->passed(reco::Muon::CutBasedIdLoose));
            mu_CutBasedIdMedium_.	push_back(mu->passed(reco::Muon::CutBasedIdMedium));
            mu_CutBasedIdTight_.	push_back(mu->passed(reco::Muon::CutBasedIdTight));
            mu_CutBasedIdMediumPrompt_.push_back(mu->passed(reco::Muon::CutBasedIdMediumPrompt));
            mu_CutBasedIdGlobalHighPt_.push_back(mu->passed(reco::Muon::CutBasedIdGlobalHighPt));
            mu_CutBasedIdTrkHighPt_.	push_back(mu->passed(reco::Muon::CutBasedIdTrkHighPt));
            mu_PFIsoVeryLoose_.	push_back(mu->passed(reco::Muon::PFIsoVeryLoose));
            mu_PFIsoLoose_.		push_back(mu->passed(reco::Muon::PFIsoLoose));
            mu_PFIsoMedium_.		push_back(mu->passed(reco::Muon::PFIsoMedium));
            mu_PFIsoTight_.		push_back(mu->passed(reco::Muon::PFIsoTight));
            mu_PFIsoVeryTight_.	push_back(mu->passed(reco::Muon::PFIsoVeryTight));
            mu_TrkIsoLoose_.		push_back(mu->passed(reco::Muon::TkIsoLoose));
            mu_TrkIsoTight_.		push_back(mu->passed(reco::Muon::TkIsoTight));
            mu_SoftCutBasedId_.	push_back(mu->passed(reco::Muon::SoftCutBasedId));
            mu_MvaLoose_.		push_back(mu->passed(reco::Muon::MvaLoose));
            mu_MvaMedium_.		push_back(mu->passed(reco::Muon::MvaMedium));
            mu_MvaTight_.		push_back(mu->passed(reco::Muon::MvaTight));
            mu_MiniIsoLoose_.		push_back(mu->passed(reco::Muon::MiniIsoLoose));
            mu_MiniIsoMedium_.		push_back(mu->passed(reco::Muon::MiniIsoMedium));
            mu_MiniIsoTight_.		push_back(mu->passed(reco::Muon::MiniIsoTight));
            mu_MiniIsoVeryTight_.	push_back(mu->passed(reco::Muon::MiniIsoVeryTight));
            
        }
        
    }
    // write all information into the tree
    tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
Ntupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Ntupler::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
 void
 Ntupler::beginRun(edm::Run const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method called when ending the processing of a run  ------------
/*
 void
 Ntupler::endRun(edm::Run const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method called when starting to processes a luminosity block  ------------
/*
 void
 Ntupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
 void
 Ntupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ntupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

int Ntupler::matchToTruth(const edm::Ptr<reco::GsfElectron> el, 
                          const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles)
{
    //
    // Explicit loop and geometric matching method (advised by Josh Bendavid)
    //
    
    // Find the closest status 1 gen electron to the reco electron
    double dR = 999;
    const reco::Candidate *closestElectron = 0;
    for(size_t i=0; i<prunedGenParticles->size();i++)
    {
        const reco::Candidate *particle = &(*prunedGenParticles)[i];
        // Drop everything that is not electron or not status 1
        if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
            continue;
        //
        double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
        if( dRtmp < dR )
        {
            dR = dRtmp;
            closestElectron = particle;
        }
    }
    // See if the closest electron (if it exists) is close enough.
    // If not, no match found.
    if( !(closestElectron != 0 && dR < 0.1) )
    {
        return UNMATCHED;
    }
    
    //
    int ancestorPID = -999;
    int ancestorStatus = -999;
    findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);
    
    if( ancestorPID == -999 && ancestorStatus == -999 )
    {
        // No non-electron parent??? This should never happen.
        // Complain.
        printf("Ntupler: ERROR! Electron does not apper to have a non-electron parent\n");
        return UNMATCHED;
    }
    
    if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
        return TRUE_NON_PROMPT_ELECTRON;
    
    if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
        return TRUE_ELECTRON_FROM_TAU;
    
    // What remains is true prompt electrons
    return TRUE_PROMPT_ELECTRON;
}

bool Ntupler::cmd(const reco::GenParticle & s1, const reco::GenParticle & s2)
{  
    return s1.pt() > s2.pt();
}

void Ntupler::findFirstNonElectronMother(const reco::Candidate *particle,
                                         int &ancestorPID, int &ancestorStatus)
{
    
    if( particle == 0 )
    {
        printf("Ntupler: ERROR! null candidate pointer, this should never happen\n");
        return;
    }
    
    // Is this the first non-electron parent? If yes, return, otherwise
    // go deeper into recursion
    if( abs(particle->pdgId()) == 11 )
    {
        findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
    }else{
        ancestorPID = particle->pdgId();
        ancestorStatus = particle->status();
    }
    
    return;
}



bool
Ntupler::hasWZasMother(const reco::GenParticle  p)  
{
    bool foundZ = false;
    if (p.numberOfMothers()==0) return foundZ;
    const reco::Candidate  *part = (p.mother());
    // loop on the mother particles to check if is has a Z has mother
    //    while ((part->numberOfMothers()>0)) {
    //        const reco::Candidate  *MomPart =part->mother();
    //if ((fabs(MomPart->pdgId())>=22)&&(fabs(MomPart->pdgId())<=24)){
    if ((fabs(part->pdgId())==23))
    {
        foundZ = true;
        //    break;
    }
    //        part = MomPart;
    //    }
    return foundZ;
}



/*void Ntupler::printCutFlowResult(vid::CutFlowResult &cutflow){
 
 printf("    CutFlow name= %s    decision is %d\n",
 cutflow.cutFlowName().c_str(),
 (int) cutflow.cutFlowPassed());
 int ncuts = cutflow.cutFlowSize();
 printf(" Index                               cut name              isMasked    value-cut-upon     pass?\n");
 for(int icut = 0; icut<ncuts; icut++){
 printf("  %2d      %50s    %d        %f          %d\n", icut,
 cutflow.getNameAtIndex(icut).c_str(),
 (int)cutflow.isCutMasked(icut),
 cutflow.getValueCutUpon(icut),
 (int)cutflow.getCutResultByIndex(icut));
 }
 
 }
 */

//define this as a plug-in
DEFINE_FWK_MODULE(Ntupler);
