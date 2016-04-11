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

#include "mygeantbooker.h"
#include "hparticlebooker.h"
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

#include <tr1/memory>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <set>

using namespace std;
#include "selectFunctions.h"
//------------------------------------------------------------------
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
#include "helperFunctions_multbin_09_19_intermediate.h"
#include "yield_pars.h"
#include "fline.h"
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
	&& pcand->getMomentum() < 2000
	&& pcand->getMetaMatchQuality() < 3.
	&& pcand->getChi2() < 100
        && HParticleTool::isGoodMetaCell(pcand,4,kTRUE)
        && gLoop->goodSector(pcand->getSector())
        && !pcand->isAtAnyMdcEdge()
	;

}

Int_t dieleAna(TString inputlist, TString outfile, Int_t nev=-1, TString seqnumlist = "dummy")
{
    setupMassCuts("m2_cut.root");
    setupRichCuts("iso_newer.root");
    setupMLPCuts("mlpmom_cutg_gen8_new.root");
    setupSectorCorrs(kPM,"sector_factor_pm.root");
    setupSectorCorrs(kPP,"sector_factor_pp.root");
    setupSectorCorrs(kMM,"sector_factor_mm.root");
    setupPiotrFactor("piotrfactor_gen8_consistent.root");

    HLoop* loop = new HLoop(kTRUE);  // kTRUE : create Hades  (needed to work with standard eventstructure)
    TString readCategories = "-*,+HParticleCand,+HParticleEvtInfo,+HParticleBtRing,+HStart2Hit,+HStart2Cal,+HGeantKine,+HGeantMdc,+HGeantTof,+HGeantRpc,+HGeantShower,+HTofHit,+HRpcHit";
    if (inputlist.EndsWith(".list")) {
        loop->addFilesList(inputlist);
    }
    else {
        loop->addMultFiles(inputlist);
    }
    if(!loop->setInput(readCategories)) { exit(1); }
    loop->printCategories();
    loop->readSectorFileList("FileListLepton.list");
    int sectors[6];
    HParticleCand *cand1;
    HParticleCand *cand2;
    HCategory* candCat = (HCategory*)HCategoryManager::getCategory(catParticleCand);

    TF1 *yield_pos_fit[6];
    TF1 *yield_neg_fit[6];
    for (int sec = 0; sec < 6; ++sec) {
        yield_pos_fit[sec] = new TF1(Form("yield_pos_fit_sec%i",sec),fline,96,126,2);
        yield_neg_fit[sec] = new TF1(Form("yield_neg_fit_sec%i",sec),fline,96,126,2);
        yield_pos_fit[sec]->SetParameters(yield_pos_pars[sec]);
        yield_neg_fit[sec]->SetParameters(yield_neg_pars[sec]);
    }

    HHistMap hMSingle(outfile);
    hMSingle.setSilentFail(kTRUE);

    hMSingle.addHistArray("TH1F","hEpHour_cand",  "hEpHour_sec%i_cand",  "hEpHour_sec%i_cand", 744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hEpHour_ided",  "hEpHour_sec%i_ided",  "hEpHour_sec%i_ided", 744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hEmHour_cand",  "hEmHour_sec%i_cand",  "hEmHour_sec%i_cand", 744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hEmHour_ided",  "hEmHour_sec%i_ided",  "hEmHour_sec%i_ided", 744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hPipHour",      "hPipHour_sec%i",      "hPipHour_sec%i",744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hPimHour",      "hPimHour_sec%i",      "hPimHour_sec%i",744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hPi0Hour_cand", "hPi0Hour_sec%i_cand", "hPi0Hour_sec%i_cand",744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hPi0Hour_open", "hPi0Hour_sec%i_open", "hPi0Hour_sec%i_cand",744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHist("TH1F","hEvtHour", "hEvtHour",744*60,96,127);
    hMSingle.addHist("TH1F","hCounter", "hCounter",20,0,20);
    string currentFileName;
    string currentFileRun;

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

    TString sys0_weights_file = "sys0_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxshowerDqmom.weights.xml";     //  6 <- minimal/optimal
    TString sys1_weights_file = "sys1_weights/TMVAClassification_MLP.NPACPMbetamdcdEdxtofdEdx_beta085.weights.xml"; //  6 <-
    TMVA::Reader *reader[2];
    reader[0] = new TMVA::Reader( "!Color:!Silent" );
    reader[1] = new TMVA::Reader( "!Color:!Silent" );

    if (sys0_weights_file.Contains("NP"))       reader[0]->AddVariable( "ringNP",   &_ringNP );
    if (sys1_weights_file.Contains("NP"))       reader[1]->AddVariable( "ringNP",   &_ringNP );
    if (sys0_weights_file.Contains("AC"))       reader[0]->AddVariable( "ringAC",   &_ringAC );
    if (sys1_weights_file.Contains("AC"))       reader[1]->AddVariable( "ringAC",   &_ringAC );
    if (sys0_weights_file.Contains("metaQa"))   reader[0]->AddVariable( "metaQa",   &_metaQa );
    if (sys1_weights_file.Contains("metaQa"))   reader[1]->AddVariable( "metaQa",   &_metaQa );
    if (sys0_weights_file.Contains("richQa"))   reader[0]->AddVariable( "richQa",   &_richQa );
    if (sys1_weights_file.Contains("richQa"))   reader[1]->AddVariable( "richQa",   &_richQa );
    if (sys0_weights_file.Contains("HT"))       reader[0]->AddVariable( "ringHT",   &_ringHT );
    if (sys1_weights_file.Contains("HT"))       reader[1]->AddVariable( "ringHT",   &_ringHT );
    if (sys0_weights_file.Contains("PM"))       reader[0]->AddVariable( "ringPM",   &_ringPM );
    if (sys1_weights_file.Contains("PM"))       reader[1]->AddVariable( "ringPM",   &_ringPM );
    if (sys0_weights_file.Contains("beta"))     reader[0]->AddVariable( "beta",     &_beta );
    if (sys1_weights_file.Contains("beta"))     reader[1]->AddVariable( "beta",     &_beta );
    if (sys0_weights_file.Contains("mdcdEdx"))  reader[0]->AddVariable( "mdcdEdx",  &_mdcdEdx );
    if (sys1_weights_file.Contains("mdcdEdx"))  reader[1]->AddVariable( "mdcdEdx",  &_mdcdEdx );
    if (sys0_weights_file.Contains("theta"))    reader[0]->AddVariable( "theta",    &_theta );
    if (sys1_weights_file.Contains("theta"))    reader[1]->AddVariable( "theta",    &_theta );
    if (sys0_weights_file.Contains("showerDq")) reader[0]->AddVariable( "showerDq", &_showerDq );
    if (sys1_weights_file.Contains("tofdEdx"))  reader[1]->AddVariable( "tofdEdx",  &_tofdEdx );
    if (sys0_weights_file.Contains("mom"))      reader[0]->AddVariable( "mom",      &_mom );
    if (sys1_weights_file.Contains("mom"))      reader[1]->AddVariable( "mom",      &_mom );

    reader[0]->BookMVA( "MLP", sys0_weights_file );
    reader[1]->BookMVA( "MLP", sys1_weights_file );

    //------------------------------------------------------------------------------------------
    //--------------- Lepton pairs spectra -----------------------------------------------------
    //------------------------------------------------------------------------------------------
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
        SectorCorr sectorCorr = chooseSectorCorr(sectors);

        HParticleEvtInfo* evtinfo;
        evtinfo = HCategoryManager::getObject(evtinfo,catParticleEvtInfo,0);

        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(1, "all");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(2, "kGoodTRIGGER");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(3, "kGoodSTART");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(4, "kGoodVertexCand");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(5, "kNoPileUpSTART");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(6, "kNoVETO");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(7, "kGoodSTARTVETO");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(8, "kGoodSTARTMETA");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(9, "#geq 4 secs");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(10,"#geq 5 secs");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(11,"6 secs");
        for (int mb = 0; mb <= 4; ++mb) {
            hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(12+mb,Form("mult bin %i",mb));
        }
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(17,"mult underflow");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(18,"4 sec d0");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(19,"4 sec d1");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(20,"4 sec d2");

        hMSingle.get("hCounter")->Fill(0);
        if (!evtinfo->isGoodEvent(kGoodTRIGGER)) continue;
        hMSingle.get("hCounter")->Fill(1);

        if (!evtinfo->isGoodEvent(kGoodSTART)) continue;
        hMSingle.get("hCounter")->Fill(2);

        if (!evtinfo->isGoodEvent(kGoodVertexCand)) continue;
        hMSingle.get("hCounter")->Fill(3);
        //if (!evtinfo->isGoodEvent(kNoPileUpSTART)) continue;
        hMSingle.get("hCounter")->Fill(4);
        //if (!evtinfo->isGoodEvent(kNoVETO)) continue;
        hMSingle.get("hCounter")->Fill(5);
        //if (!evtinfo->isGoodEvent(kGoodSTARTVETO)) continue;
        hMSingle.get("hCounter")->Fill(6);
        //if (!evtinfo->isGoodEvent(kGoodSTARTMETA)) continue;
        hMSingle.get("hCounter")->Fill(7);
        //if (sectorCorr == kSkip) continue;
        hMSingle.get("hCounter")->Fill(8);
        //if (sectorCorr != kCorr5to6) continue;
        //if (sectorCorr == kNoCorr || sectorCorr == kCorr5to6) {
        //    hMSingle.get("hCounter")->Fill(9);
        //}
        //if (sectorCorr == kNoCorr) {
        //    hMSingle.get("hCounter")->Fill(10);
        //}
	// Select meta multiplicity / event
        Int_t mult_meta = evtinfo->getSumTofMultCut() + evtinfo->getSumRpcMultHitCut();
        Int_t mult_bin = 5;
        if (mult_meta >  58 && mult_meta <=  88) mult_bin = 4;
        if (mult_meta >  88 && mult_meta <= 121) mult_bin = 3;
        if (mult_meta > 121 && mult_meta <= 160) mult_bin = 2;
        if (mult_meta > 160 && mult_meta <= 240) mult_bin = 1;
        if (mult_meta > 240) mult_bin = 0;
        if (mult_bin == 5) continue;

        hMSingle.get("hCounter")->Fill(11+mult_bin);

        if (sectorCorr == kCorr4to6d0) {
            hMSingle.get("hCounter")->Fill(17);
        }
        if (sectorCorr == kCorr4to6d1) {
            hMSingle.get("hCounter")->Fill(18);
        }
        if (sectorCorr == kCorr4to6d2) {
            hMSingle.get("hCounter")->Fill(19);
        }

        TString tempFileName;
        Int_t dayOfYear;
        Int_t hour;
        Int_t min;
        if (loop->isNewFile(tempFileName)) {
            TString type;
            Int_t year;
            Int_t sec;
            Int_t eb;
            currentFileName = tempFileName;
            size_t poslastslash = currentFileName.rfind("/");
            size_t poshld = currentFileName.find(".hld");
            if (poslastslash < poshld) {
                currentFileRun = currentFileName.substr(poslastslash+1,poshld-poslastslash-1);
            }
            else {
                currentFileRun = currentFileName.substr(0,poshld);
            }
            HTime::splitFileName(HTime::stripFileName(tempFileName),type,year,dayOfYear,hour,min,sec,eb);
        }
        Int_t binHour     = (dayOfYear-96)*24*60+hour*60+min+1;
        Float_t valueHour = hMSingle.get("hEvtHour")->GetXaxis()->GetBinCenter(binHour);
        hMSingle.get("hEvtHour")->Fill(valueHour,1);
        Float_t yield_pos_weight[6];
        Float_t yield_neg_weight[6];
        for (int _s = 0; _s < 6; ++_s) {
            if (_s != 2 && _s != 3) {
                yield_pos_weight[_s] = yield_pos_fit[_s]->Eval(126)/yield_pos_fit[_s]->Eval(valueHour);
                yield_neg_weight[_s] = yield_neg_fit[_s]->Eval(126)/yield_neg_fit[_s]->Eval(valueHour);
            }
            else {
                yield_pos_weight[_s] = 1.;
                yield_neg_weight[_s] = 1.;
            }
        }

        Int_t size = candCat->getEntries();

        sorter.cleanUp();
        sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);            // reset all flags for flags (0-28) ,reject,used,lepton
        sorter.fill(selectLeptonsBeta);   // fill only good leptons
//        Int_t n_lepparcnd = sorter.selectBest(HParticleTrackSorter::kIsBestRKRKMETA,HParticleTrackSorter::kIsLepton);
        Int_t n_lepparcnd = sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsLepton);
	//------------------------------------------------------------------------
        Float_t MLPs[size];
	Bool_t sorterFlag[size];
	Bool_t leptonFlag[size];
	Bool_t removedRecursive[size];
	Bool_t removedRecursiveNearestSeg[5][size];
        angleNearestSeg.clear();
        angleNearestBoth.clear();
        angleNearestFull.clear();
        nearestIndexSeg.clear();
	removedNearest.clear();
	removedIncomplete.clear();
        doTrackCorr(trackcor, candCat);
        // SETTNG FLAGS
	for(Int_t j = 0; j < size; j ++){
	    sorterFlag[j]                     = kFALSE;
	    leptonFlag[j]                     = kFALSE;
            for (int cp = 2; cp <= 6; ++cp) {
                removedRecursiveNearestSeg[cp-2][j]     = kFALSE;
            }

	    cand1 = HCategoryManager::getObject(cand1,candCat,j);
            if(!loop->goodSector(cand1->getSector())) continue;
            if(cand1->getSystemUsed()==-1)  continue;

	    _ringNP         = cand1->getRingNumPads();
	    _ringAC         = cand1->getRingAmplitude()/_ringNP;
	    _metaQa         = cand1->getMetaMatchQuality();
	    _richQa         = cand1->getRichMatchingQuality();
	    _ringHT         = cand1->getRingHouTra();
	    _ringPM         = cand1->getRingPatternMatrix();
	    _beta           = cand1->getBeta();
	    _mdcdEdx        = cand1->getMdcdEdx();
	    _theta          = cand1->getTheta();
	    _metaQa         = cand1->getMetaMatchQuality();
	    _showerDq       = cand1->getShowerDeltaSum();
	    _mom            = cand1->getMomentum();
	    _tofdEdx        = cand1->getTofdEdx();
            Int_t sys       = cand1->getSystemUsed();

            if (!cand1->isFlagBit(kIsLepton)) continue;
            sorterFlag[j] = kTRUE;
	    cand1->calc4vectorProperties(); // lepton by default
            nearestIndexSeg[j] = findNearestNeighbor(cand1,angleNearestBoth,angleNearestSeg,angleNearestFull,candCat);

            MLPs[j] = reader[sys]->EvaluateMVA( "MLP" );
            if (_richQa < 1.0) richQaFlag[j]=kTRUE;
            if (richQaFlag[j] && isGoodMLP(cand1,MLPs[j])) {
                leptonFlag[j]=kTRUE;
                if (nearestIndexSeg[j] > -1) {
                    HParticleCand *nearestCand;
                    nearestCand = HCategoryManager::getObject(nearestCand,candCat,nearestIndexSeg[j]);
                    if (cand1->getMetaHitInd() == nearestCand->getMetaHitInd() && (nearestCand->getChi2() > 40 || nearestCand->getChi2() < 0)) {
                        removedIncomplete[j] = kTRUE;
                    }
                    if (!(angleNearestSeg[j] < 0 || angleNearestSeg[j] > 6)) {
                        removedNearest[j] = kTRUE;
                    }
                }
            }
        }
        // RECURSIVE CUT
        for(Int_t j = 0; j < size; j ++){ // first candidate
            if (!sorterFlag[j]) continue;
            cand1 = HCategoryManager::getObject(cand1,candCat,j);

            for(Int_t k = j+1; k < size; k ++){ // second candidate, k > j condition applied below
                if (!sorterFlag[k]) continue;
                cand2 = HCategoryManager::getObject(cand2,candCat,k);

                if (k > j) {

                    TLorentzVector dilep = (*cand1) + (*cand2);

                    if((cand1->getCharge()==1  && cand2->getCharge()==-1) ||
                            (cand1->getCharge()==-1 && cand2->getCharge()==1)) {
                        Float_t oAngle   = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), cand2->getPhi(), cand2->getTheta());
                        if (oAngle < 9) {
                            removedRecursive[j] = kTRUE;
                            removedRecursive[k] = kTRUE;

                            if (!removedIncomplete[j] && !removedIncomplete[k]) {
                                for (int cp = 2; cp <= 6; ++cp) {
                                    if ((angleNearestSeg[j] < 0 || angleNearestSeg[j] > cp) && (angleNearestSeg[k] < 0 || angleNearestSeg[k] > cp)) {
                                        removedRecursiveNearestSeg[cp-2][j] = kTRUE;
                                        removedRecursiveNearestSeg[cp-2][k] = kTRUE;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } // END RECURSIVE CUT
        // SINLE LEPTON FILLING HISTOGRAMS
        for(Int_t j = 0; j < size; j ++){
            cand1 = HCategoryManager::getObject(cand1,candCat,j);

            Int_t sec      = cand1->getSector();
            if (cand1->getPID() == 9) {
                hMSingle.get("hPimHour",sec)->Fill(valueHour);
            }
            if (cand1->getPID() == 0) {
                hMSingle.get("hPipHour",sec)->Fill(valueHour);
            }
            if(cand1->getSystemUsed()==-1)  continue;
            if(!cand1->isFlagBit(Particle::kIsLepton)) continue;
            Float_t weight_yield = 1;
            if (cand1->getCharge() == -1) {
                weight_yield = yield_neg_weight[sec];
                hMSingle.get("hEmHour_cand",sec)->Fill(valueHour,weight_yield);
            }
            if (cand1->getCharge() == 1) {
                weight_yield = yield_pos_weight[sec];
                hMSingle.get("hEpHour_cand",sec)->Fill(valueHour,weight_yield);
            }
            if (leptonFlag[j]) {
                if (cand1->getCharge() == -1) {
                    hMSingle.get("hEmHour_ided",sec)->Fill(valueHour,weight_yield);
                }
                if (cand1->getCharge() == 1) {
                    hMSingle.get("hEpHour_ided",sec)->Fill(valueHour,weight_yield);
                }
            }
        }
	// PAIRS FILLING HISTOGRAMS
	for(Int_t j = 0; j < size-1; j ++){ // first candidate
	    for(Int_t k = j+1; k < size; k ++){ // second candidate
                cand1 = HCategoryManager::getObject(cand1,candCat,j);
		cand2 = HCategoryManager::getObject(cand2,candCat,k);

                UInt_t flags;
                HParticleTool::setPairFlags(flags,cand2,cand1);
                if(!HParticleTool::evalPairsFlags(kPairCase1,flags)) continue;

                if (leptonFlag[j] && leptonFlag[k]) {
                    TLorentzVector dilep = (*cand1) + (*cand2);
                    Float_t oAngle   = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), cand2->getPhi(), cand2->getTheta());

                    if (dilep.M() < 150 && cand1->getSector() == cand2->getSector()) {
                        Float_t weight_yield_pi0 = 1;
                        if (cand1->getCharge() == cand2->getCharge()) {
                            if (cand1->getCharge() == 1) {
                                weight_yield_pi0 = yield_pos_weight[cand1->getSector()]*yield_pos_weight[cand1->getSector()];
                            }
                            if (cand1->getCharge() == -1) {
                                weight_yield_pi0 = yield_neg_weight[cand1->getSector()]*yield_neg_weight[cand1->getSector()];
                            }
                            hMSingle.get("hPi0Hour_cand",cand1->getSector())->Fill(valueHour,-weight_yield_pi0);
                        }
                        else {
                            weight_yield_pi0 = yield_neg_weight[cand1->getSector()]*yield_pos_weight[cand1->getSector()];
                            hMSingle.get("hPi0Hour_cand",cand1->getSector())->Fill(valueHour,weight_yield_pi0);
                        }
                    }
                    if(!removedRecursiveNearestSeg[4][j] && !removedRecursiveNearestSeg[4][k] && oAngle > 9) {
                        //cout << "DUPA " << dilep.M() << " " << cand1->getSector() << " " << cand2->getSector() << endl;
                        if (dilep.M() < 150 && cand1->getSector() == cand2->getSector()) {
                            Float_t weight_yield_pi0 = 1;
                            if (cand1->getCharge() == cand2->getCharge()) {
                                //cout << "DUPA 1" << endl;
                                if (cand1->getCharge() == 1) {
                                    weight_yield_pi0 = yield_pos_weight[cand1->getSector()]*yield_pos_weight[cand1->getSector()];
                                }
                                if (cand1->getCharge() == -1) {
                                    weight_yield_pi0 = yield_neg_weight[cand1->getSector()]*yield_neg_weight[cand1->getSector()];
                                }
                                hMSingle.get("hPi0Hour_open",cand1->getSector())->Fill(valueHour,-weight_yield_pi0);
                            }
                            else {
                                //cout << "DUPA 2" << endl;
                                weight_yield_pi0 = yield_neg_weight[cand1->getSector()]*yield_pos_weight[cand1->getSector()];
                                hMSingle.get("hPi0Hour_open",cand1->getSector())->Fill(valueHour,weight_yield_pi0);
                            }
                        }
                    }
                } // leptonFlag[j] && leptonFlag[k]
                //----------------------------------------------------------------------
            } // inner loop
        } // END PAIRS FILLING HISTOGRAMS
    } // end event loop

    sorter.finalize();
    timer.Stop();
    hMSingle.writeHists("nomap");
    cout<<"####################################################"<<endl;
    return 0;
}

