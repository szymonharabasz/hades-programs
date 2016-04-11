#include "hades.h"
#include "htool.h"
#include "htime.h"
#include "hphysicsconstants.h"
#include "hrootsource.h"
#include "hiterator.h"
#include "hloop.h"
#include "hhistmap.h"
#include "hdst.h"

#include "haddef.h"
#include "heventheader.h"
//#include "hgeantdef.h"

#include "hparticleanglecor.h"
#include "hparticledef.h"
#include "hparticlestructs.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hcategorymanager.h"
#include "hparticletracksorter.h"
#include "hparticlevertexfind.h"
#include "hparticleevtinfo.h"
#include "hparticlet0reco.h"
#include "htaskset.h"
//#include "hgeantkine.h"
#include "TMVA/Reader.h"

#include "hstart2hit.h"
#include "hstart2cal.h"
#include "hmdcdef.h"
#include "hmdctrackddef.h"
#include "hmdctrackgdef.h"
#include "horadef.h"
#include "horasimdef.h"
#include "hstartdef.h"
#include "richdef.h"
#include "rpcdef.h"
#include "showerdef.h"
#include "simulationdef.h"
#include "tofdef.h"
#include "walldef.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TCutG.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TNtuple.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <stdio.h>
#include "selectFunctions.h"
static const Float_t oAngleCut      = 9.;
static map<Int_t, Bool_t> sorterFlag;
static map<Int_t, Bool_t> richQaFlag;
static map<Int_t, Bool_t> betaFlag;
static map<Int_t, Bool_t> betaFlagTof;
static map<Int_t, Bool_t> betaFlagRpc;
static map<Int_t, Bool_t> showerFlag;
static map<Int_t, Bool_t> leptonFlagHC;
static map<Int_t, Bool_t> leptonFlagMLP[NWEIGHTS];
static map<Int_t, Bool_t> removedNearest;
static map<Int_t, Bool_t> removedIncomplete;
//------------------------------------------------------------------
#include __HELPERFC_FILE__
#include "yield_pars.h"
#include "fline.h"
#include "eventclassifier.h"

using namespace std;


//------------------------------------------------------------

static Bool_t selectLeptonsBeta(HParticleCand* pcand){
    //  selection function for lepton candidates.
    Bool_t selectEpEm = kFALSE;

    if(pcand->isFakeRejected()) return kFALSE;

    if(pcand->isFlagAND(5,
			Particle::kIsAcceptedHitRICH,
			Particle::kIsAcceptedHitInnerMDC,
			Particle::kIsAcceptedHitOuterMDC,
			Particle::kIsAcceptedHitMETA,
			Particle::kIsAcceptedRK)
      ) selectEpEm = kTRUE;

    return
	selectEpEm
	&& pcand->getBeta() > 0.9
	&& pcand->getMomentum() > 100
	&& pcand->getMomentum() < 2000
	&& pcand->getMetaMatchQuality() < 3.
	&& pcand->getChi2() < 100
        && HParticleTool::isGoodMetaCell(pcand,4,kTRUE)
        && gLoop->goodSector(pcand->getSector())
        && pcand->getSector() != 2 // don't use sector 2 anyway
        && !pcand->isAtAnyMdcEdge()
	;
}

//------------------------------------------------------------

Int_t dieleAna(TString inputlist, TString outfile, Int_t nev=-1)
{
    TH1::SetDefaultSumw2();
    setupMLPCuts("mlpmom_cutg_gen8.root");

    TH3F *p3DEffEle[6][NWEIGHTS+1]; // mult bins, MLP weights + HC
    TH3F *p3DEffPos[6][NWEIGHTS+1];
    TH3F *p3DAccEle[6][NWEIGHTS+1];
    TH3F *p3DAccPos[6][NWEIGHTS+1];
    readAccEffMatrices(p3DAccEle, p3DAccPos, p3DEffEle, p3DEffPos);

    HLoop* loop = new HLoop(kTRUE);  // kTRUE : create Hades  (needed to work with standard eventstructure)
    TString readCategories = "-*,+HParticleCand,+HParticleMdc,+HParticleEvtInfo,+HStart2Hit,+HStart2Cal,+HGeantKine,+HGeantMdc,+HGeantTof,+HGeantRpc,+HGeantShower";
    loop->addMultFiles(inputlist);

    if(!loop->setInput(readCategories)) { exit(1); }
    loop->printCategories();
    loop->readSectorFileList("FileListLepton.list");
    int sectors[6];

    HParticleCand* cand1;
    HCategory* candCat = (HCategory*)HCategoryManager::getCategory(catParticleCand);

    HHistMap hM(outfile.Data());
    hM.setSilentFail(kTRUE);

    HParticleAngleCor trackcor;
    trackcor.setDefaults("apr12");

    Float_t   _ringNP;
    Float_t   _ringAC;
    Float_t   _metaQa;
    Float_t   _richQa;
    Float_t   _ringHT;
    Float_t   _ringPM;
    Float_t   _beta;
    Float_t   _mdcdEdx;
    Float_t   _theta;
    Float_t   _showerDq;
    Float_t   _mom;
    Float_t   _tofdEdx;

    TString sys0_weights_file[NWEIGHTS] = { "sys0_weights/TMVAClassification_MLP_NPACmetaQaPMbetamdcdEdxthetashowerDqmom_lessstat.weights.xml",  //  0
                                            "sys0_weights/TMVAClassification_MLP_NPACmetaQaHTbetamdcdEdxthetashowerDqmom.weights.xml",           //  1
                                            "sys0_weights/TMVAClassification_MLP_NPACHTbetamdcdEdxthetashowerDqmom.weights.xml",                 //  2
                                            "sys0_weights/TMVAClassification_MLP_NPACbetamdcdEdxthetashowerDqmom.weights.xml",                   //  3
                                            "sys0_weights/TMVAClassification_MLP_NPACbetamdcdEdxshowerDqmom.weights.xml",                        //  4
                                            "sys0_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxthetashowerDqmom.weights.xml",                 //  5
                                            "sys0_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxshowerDqmom.weights.xml",                      //  6 <- minimal/optimal
                                            "sys0_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxshowerDqmom.weights.xml",                      //  7 = 6
                                            "sys0_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxshowerDqmom.weights.xml",                      //  8 = 6
                                            "sys0_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxshowerDqmom.weights.xml",                      //  9 = 6
                                            "sys0_weights/TMVAClassification_MLP.NPACrichQaPMbetamdcdEdxshowerDqmom_sim.weights.xml",            // 11 <- train on sim, have richQa
                                            "sys0_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxshowerDqmom_sim.weights.xml",                  // 10 <- train on sim, don't have richQa
                                            "sys0_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxshowerDqmom_simtrainonrich.weights.xml",       // 12 <- train on sim like on exp
    };

    TString sys1_weights_file[NWEIGHTS] = { "sys1_weights/TMVAClassification_MLP.NPACmetaQaPMbetamdcdEdxthetatofdEdx_beta085.weights.xml", //  0
                                            "sys1_weights/TMVAClassification_MLP_NPACmetaQaHTbetamdcdEdxthetatofdEdx_beta085.weights.xml", //  1
                                            "sys1_weights/TMVAClassification_MLP_NPACHTbetamdcdEdxthetatofdEdx.weights.xml",               //  2
                                            "sys1_weights/TMVAClassification_MLP_NPACbetamdcdEdxthetatofdEdx.weights.xml",                 //  3
                                            "sys1_weights/TMVAClassification_MLP_NPACbetamdcdEdxtofdEdx.weights.xml",                      //  4
                                            "sys1_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxthetatofdEdx_beta085.weights.xml",       //  5
                                            "sys1_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxtofdEdx_beta085.weights.xml",            //  6 <-
                                            "sys1_weights/TMVAClassification_MLP.NPACPMbetamdcdEdx.weights.xml",                           //  7 <- remove tofdEdx
                                            "sys1_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxtofdEdxmom_beta085.weights.xml",         //  8 <- add mom
                                            "sys1_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxtofdEdx_beta085_nobt.weights.xml",       //  9
                                            "sys1_weights/TMVAClassification_MLP.NPACrichQaPMbetamdcdEdxtofdEdx_sim.weights.xml",          // 10 <- train on sim, have richQa
                                            "sys1_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxtofdEdx_sim.weights.xml",                // 11 <- train on sim, don't have richQa
                                            "sys1_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxtofdEdx_simtrainonrich.weights.xml",     // 12 <- train on sim like on exp
    };
    bool useForPairs[NWEIGHTS][NWEIGHTS] = {  // x, cols, second index: sys1; y, rows, first index: sys0
        { 1,0,0,0,0,0,0,0,0,0,0,0,0, }, // 0
        { 0,0,0,0,0,0,0,0,0,0,0,0,0, }, // 1
        { 0,0,0,0,0,0,0,0,0,0,0,0,0, }, // 2
        { 0,0,0,0,0,0,0,0,0,0,0,0,0, }, // 3
        { 0,0,0,0,0,0,0,0,0,0,0,0,0, }, // 4
        { 0,0,0,0,0,0,0,0,0,0,0,0,0, }, // 5
        { 0,0,0,0,0,0,1,0,0,0,0,0,0, }, // 6
        { 0,0,0,0,0,0,0,0,0,0,0,0,0, }, // 7
        { 0,0,0,0,0,0,0,0,0,0,0,0,0, }, // 8
        { 0,0,0,0,0,0,0,0,0,0,0,0,0, }, // 9
        { 0,0,0,0,0,0,0,0,0,0,1,0,0, }, // 10
        { 0,0,0,0,0,0,0,0,0,0,0,0,0, }, // 11
        { 0,0,0,0,0,0,0,0,0,0,0,0,1, }, // 12

//        0 1 2 3 4 5 6 7 8 9101112          

    };
    for (int w = 0; w < NWEIGHTS; ++w) {
        reader[w][0] = new TMVA::Reader( "!Color:!Silent" );
        reader[w][1] = new TMVA::Reader( "!Color:!Silent" );

        if (sys0_weights_file[w].Contains("NP"))       reader[w][0]->AddVariable( "ringNP",   &_ringNP );
        if (sys1_weights_file[w].Contains("NP"))       reader[w][1]->AddVariable( "ringNP",   &_ringNP );
        if (sys0_weights_file[w].Contains("AC"))       reader[w][0]->AddVariable( "ringAC",   &_ringAC );
        if (sys1_weights_file[w].Contains("AC"))       reader[w][1]->AddVariable( "ringAC",   &_ringAC );
        if (sys0_weights_file[w].Contains("metaQa"))   reader[w][0]->AddVariable( "metaQa",   &_metaQa );
        if (sys1_weights_file[w].Contains("metaQa"))   reader[w][1]->AddVariable( "metaQa",   &_metaQa );
        if (sys0_weights_file[w].Contains("richQa"))   reader[w][0]->AddVariable( "richQa",   &_richQa );
        if (sys1_weights_file[w].Contains("richQa"))   reader[w][1]->AddVariable( "richQa",   &_richQa );
        if (sys0_weights_file[w].Contains("HT"))       reader[w][0]->AddVariable( "ringHT",   &_ringHT );
        if (sys1_weights_file[w].Contains("HT"))       reader[w][1]->AddVariable( "ringHT",   &_ringHT );
        if (sys0_weights_file[w].Contains("PM"))       reader[w][0]->AddVariable( "ringPM",   &_ringPM );
        if (sys1_weights_file[w].Contains("PM"))       reader[w][1]->AddVariable( "ringPM",   &_ringPM );
        if (sys0_weights_file[w].Contains("beta"))     reader[w][0]->AddVariable( "beta",     &_beta );
        if (sys1_weights_file[w].Contains("beta"))     reader[w][1]->AddVariable( "beta",     &_beta );
        if (sys0_weights_file[w].Contains("mdcdEdx"))  reader[w][0]->AddVariable( "mdcdEdx",  &_mdcdEdx );
        if (sys1_weights_file[w].Contains("mdcdEdx"))  reader[w][1]->AddVariable( "mdcdEdx",  &_mdcdEdx );
        if (sys0_weights_file[w].Contains("theta"))    reader[w][0]->AddVariable( "theta",    &_theta );
        if (sys1_weights_file[w].Contains("theta"))    reader[w][1]->AddVariable( "theta",    &_theta );
        if (sys0_weights_file[w].Contains("showerDq")) reader[w][0]->AddVariable( "showerDq", &_showerDq );
        if (sys1_weights_file[w].Contains("tofdEdx"))  reader[w][1]->AddVariable( "tofdEdx",  &_tofdEdx );
        if (sys0_weights_file[w].Contains("mom"))      reader[w][0]->AddVariable( "mom",      &_mom );
        if (sys1_weights_file[w].Contains("mom"))      reader[w][1]->AddVariable( "mom",      &_mom );

        reader[w][0]->BookMVA( "MLP", sys0_weights_file[w] );
        reader[w][1]->BookMVA( "MLP", sys1_weights_file[w] );
    }
    hM.addHistArray("TH1F","eventClass","eventClass_w%i","Event Class, weights %i",2000,0,2000,0,0,0,0,0,0,"","","","",NWEIGHTS);
    //--------------- end histo booking -----------------------------------------------------

    //--------------------------CONFIGURATION---------------------------------------------------
    //At begin of the program (outside the event loop)
    HParticleTrackSorter sorter;
    //sorter.setDebug();                                            // for debug
    //sorter.setPrintLevel(3);                                      // max prints
    //sorter.setRICHMatching(HParticleTrackSorter::kUseRKRICHWindow,4.); // select matching RICH-MDC for selectLeptons() function
    sorter.setIgnoreInnerMDC();                                   // do not reject Double_t inner MDC hits
    //sorter.setIgnoreOuterMDC();                                   // do not reject Double_t outer MDC hits
    //sorter.setIgnoreMETA();                                       // do not reject Double_t META hits
    //sorter.setIgnorePreviousIndex();                              // do not reject indices from previous selctions
    sorter.init();                                                  // get catgegory pointers etc...
    //    --------------------------------------------------------------------------------------------

    TStopwatch timer;
    timer.Reset();
    timer.Start();

    Int_t evtsInFile = loop->getEntries();
    if(nev < 0 || nev > evtsInFile ) nev = evtsInFile;

    for(Int_t i = 1; i < nev; i++)
    {
        //----------break if last event is reached-------------
        //if(!gHades->eventLoop(1)) break;
        if(loop->nextEvent(i) <= 0) { cout<<" end recieved "<<endl; break; } // last event reached
        HTool::printProgress(i,nev,1,"Analyze :");
        loop->getSectors(sectors);
        // Remove sector 2 from analysis and take this into account to characterize the event
        sectors[2] = 0;
        sectorCorr = chooseSectorCorr(sectors);
//        if (sectorCorr == kNoCorr) sectorCorr = kCorr5to6;

        TString tempFileName;
        Int_t dayOfYear;
        Int_t hour;
        Int_t min;
        if (loop->isNewFile(tempFileName)) {
            cout << "sec corr: " << sectorCorr << endl;
            TString type;
            Int_t year;
            Int_t sec;
            Int_t eb;
            currentFileName = tempFileName;
            isSimulation = tempFileName.Contains("Au_Au") || tempFileName.Contains("dilepton") || tempFileName.Contains("urqmd");
            isEmbedding  = tempFileName.Contains("e+") || tempFileName.Contains("e-") || tempFileName.Contains("pi0") || tempFileName.Contains("eta");
            HTime::splitFileName(HTime::stripFileName(tempFileName),type,year,dayOfYear,hour,min,sec,eb);
        }

        HParticleEvtInfo* evtinfo;
        evtinfo = HCategoryManager::getObject(evtinfo,catParticleEvtInfo,0);
        HVertex PrimVertexReco = ((HEventHeader*)(gHades->getCurrentEvent()->getHeader()))->getVertexReco(); 
        Float_t vx=PrimVertexReco.getX(); 
        Float_t vy=PrimVertexReco.getY(); 
        Float_t vz=PrimVertexReco.getZ(); 

        if (!isSimulation) {
            if (!evtinfo->isGoodEvent(kGoodTRIGGER)) continue;
            if (!evtinfo->isGoodEvent(kGoodSTART)) continue;
            if (!evtinfo->isGoodEvent(kGoodVertexCand)) continue;
            if (!evtinfo->isGoodEvent(kNoPileUpSTART)) continue;
            if (!evtinfo->isGoodEvent(kNoVETO)) continue;
            if (!evtinfo->isGoodEvent(kGoodSTARTVETO)) continue;
            if (!evtinfo->isGoodEvent(kGoodSTARTMETA)) continue;
            if (sectorCorr == kSkip) continue;
        }
	// Select meta multiplicity / event
        Int_t mult_meta = evtinfo->getSumTofMultCut() + evtinfo->getSumRpcMultHitCut();
        Int_t mult_bin = 5;
        if (mult_meta >  60 && mult_meta <=  88) mult_bin = 4; // most peripheral
        if (mult_meta >  88 && mult_meta <= 121) mult_bin = 3;
        if (mult_meta > 121 && mult_meta <= 160) mult_bin = 2;
        if (mult_meta > 160 && mult_meta <= 250) mult_bin = 1; // most central
        if (mult_meta > 250) mult_bin = 0;
        if (mult_bin == 5) continue;


        Int_t size = candCat->getEntries();

        sorter.cleanUp();
        sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);            // reset all flags for flags (0-28) ,reject,used,lepton
        sorter.fill(selectLeptonsBeta);   // fill only good leptons
        sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsLepton);

        Int_t size = candCat->getEntries();
	sorterFlag.clear();
        richQaFlag.clear();
        betaFlag.clear();
        betaFlagTof.clear();
        betaFlagRpc.clear();
        showerFlag.clear();
	leptonFlagHC.clear();
	removedNearest.clear();
	removedIncomplete.clear();
        MLPs.clear();
	for(Int_t j = 0; j < size; j ++){

            // Track candidates (Leptons AND Hadrons)
	    cand1 = HCategoryManager::getObject(cand1,candCat,j);

	    if(cand1->getSystemUsed()==-1)  continue;
            if(!cand1->isFlagBit(Particle::kIsLepton)) continue;
//            if (startFlag != 2) continue;

	    // ------ variables needed for mlp
	    trackcor.recalcSetEmission(cand1);  // changes values inside candidate
	    trackcor.realignRichRing(cand1);    //

	    _ringNP        = cand1->getRingNumPads();
	    _ringAC        = cand1->getRingAmplitude()/_ringNP;
	    _metaQa        = cand1->getMetaMatchQuality();
	    _ringHT        = cand1->getRingHouTra();
	    _ringPM        = cand1->getRingPatternMatrix();
	    _beta          = cand1->getBeta();
	    _mdcdEdx       = cand1->getMdcdEdx();
	    _theta         = cand1->getTheta();
	    _showerDq      = cand1->getShowerDeltaSum();
	    _mom           = cand1->getMomentum();
	    _TofdEdx       = cand1->getTofdEdx();
            Float_t phi    = cand1->getPhi();
            Int_t sys      = cand1->getSystemUsed();
	    if(cand1->getSystemUsed()==-1)  continue;
            if (!cand1->isFlagBit(Particle::kIsLepton)) continue;
            sorterFlag[j] = kTRUE;
            for (int w = 0; w < NWEIGHTS; ++w) {
                MLPs[w][j] = reader[w][sys]->EvaluateMVA( "MLP" );
            }
            Float_t weight = getEfficiencyFactor(p3DEffEle[0], p3DEffPos[0], _mom, _theta, phi, cand1->getCharge());
            if (richQa < 1.0) richQaFlag[j]=kTRUE;
            for (int w = 0; w < NWEIGHTS; ++w) {
                if (richQaFlag[j] && isGoodMLP(cand1,MLPs[w][j],w)) {
                    leptonFlagMLP[w][j]=kTRUE;
                }
            }

	    //----------------------------------------------------------

	    cand1->calc4vectorProperties(); // lepton by default

	    //*************************************************************
	}  // end pre-loop

        bool isEmptyEvent[NWEIGHTS];
        for (int w = 0; w < NWEIGHTS; ++w) {
            isEmptyEvent[w] = kTRUE;
            for(Int_t j = 0; j < size; j ++){
                if (!sorterFlag[j]) continue;
                cand1 = HCategoryManager::getObject(cand1,candCat,j);
                if (leptonFlagMLP[w][j]) {
                    isEmptyEvent[w] = kFALSE;
                }
            }
            if (!isEmptyEvent[w]) {
                Int_t eventClass = eventClassifierMult();
                hM.get("eventClass",w)->Fill(eventClass);
            }
        }

    } // end event loop

    sorter.finalize();
    timer.Stop();

    hM.writeHists();

    //delete myHades;
    cout<<"####################################################"<<endl;

    return 0;
}
