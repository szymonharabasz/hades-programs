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

#include "TLorentzVector.h"
#include "TObjString.h"
#include "TStopwatch.h"
#include "TObjArray.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TCutG.h"
#include "TH1F.h"
#include "TF1.h"

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
static map<Int_t, Bool_t> removedClosePartner;
static map<Int_t, Bool_t> closePartnerCandFlag;
//------------------------------------------------------------------
#include __HELPERFC_FILE__
#include "analyzekine.h"
#include "yield_pars.h"
#include "fline.h"
#include "eventrecord.h"

class abstract_functor {
    public: virtual bool operator()(HParticleCand *) const = 0;
};
class accept_all : public abstract_functor {
    public: bool operator() (HParticleCand *cand) const { return true; }
};
class accept_leptons : public abstract_functor {
    public: bool operator() (HParticleCand *cand) const { return leptonFlagHC[cand->getIndex()]; }
};
class same_inner_seg : public abstract_functor {
    public:
        same_inner_seg (Int_t ind) : index(ind) {}
        same_inner_seg (HParticleCand *cand) : index(cand->getInnerSegInd() ) {}
        bool operator() (HParticleCand *cand) const { return cand->getInnerSegInd() == index; }
    private:
        Int_t index;
};
class same_inner_same_outer_seg : public abstract_functor {
    public:
        same_inner_same_outer_seg (Int_t ind_inn, Int_t ind_out) : index_inn(ind_inn), index_out(ind_out) {}
        same_inner_same_outer_seg (HParticleCand *cand) : index_inn(cand->getInnerSegInd()), index_out(cand->getOuterSegInd()) {}
        bool operator() (HParticleCand *cand) const { return cand->getInnerSegInd() == index_inn && cand->getOuterSegInd() == index_out; }
    private:
        Int_t index_inn;
        Int_t index_out;
};
class same_inner_diff_outer_seg : public abstract_functor {
    public:
        same_inner_diff_outer_seg (Int_t ind_inn, Int_t ind_out) : index_inn(ind_inn), index_out(ind_out) {}
        same_inner_diff_outer_seg (HParticleCand *cand) : index_inn(cand->getInnerSegInd()), index_out(cand->getOuterSegInd()) {}
        bool operator() (HParticleCand *cand) const { return cand->getInnerSegInd() == index_inn && cand->getOuterSegInd() != index_out; }
    private:
        Int_t index_inn;
        Int_t index_out;
};
class same_inner_seg_leptons : public abstract_functor {
    public:
        same_inner_seg_leptons (Int_t ind) : index(ind) {}
        same_inner_seg_leptons (HParticleCand *cand) : index(cand->getInnerSegInd() ) {}
        bool operator() (HParticleCand *cand) const { return leptonFlagHC[cand->getIndex()] && cand->getInnerSegInd() == index; }
    private:
        Int_t index;
};
class same_inner_seg_leptons_like : public abstract_functor {
    public:
        same_inner_seg_leptons_like (Int_t ind, Int_t chg) : index(ind), charge(chg) {}
        same_inner_seg_leptons_like (HParticleCand *cand) : index(cand->getInnerSegInd()), charge(cand->getCharge()) {}
        bool operator() (HParticleCand *cand) const { return leptonFlagHC[cand->getIndex()] && cand->getInnerSegInd() == index && cand->getCharge() == charge; }
    private:
        Int_t index;
        Int_t charge;
};
class same_inner_seg_leptons_unlike : public abstract_functor {
    public:
        same_inner_seg_leptons_unlike (Int_t ind, Int_t chg) : index(ind), charge(chg) {}
        same_inner_seg_leptons_unlike (HParticleCand *cand) : index(cand->getInnerSegInd()), charge(cand->getCharge()) {}
        bool operator() (HParticleCand *cand) const { return leptonFlagHC[cand->getIndex()] && cand->getInnerSegInd() == index && cand->getCharge() != charge; }
    private:
        Int_t index;
        Int_t charge;
};
//------------------------------------------------------------
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
	&& pcand->getMomentum() < 1000
	&& pcand->getMetaMatchQuality() < 3.
	&& pcand->getChi2() < 100
        && HParticleTool::isGoodMetaCell(pcand,4,kTRUE)
        && gLoop->goodSector(pcand->getSector())
        && pcand->getSector() != 2 // don't use sector 2 anyway
        && !pcand->isAtAnyMdcEdge()
	;

}

static Bool_t selectHadrons(HParticleCand* pcand){
    //  selection function for lepton candidates.
    Bool_t selectEpEm = kFALSE;

    if(pcand->isFakeRejected()) return kFALSE;

    if(pcand->isFlagAND(4,
			Particle::kIsAcceptedHitInnerMDC,
			Particle::kIsAcceptedHitOuterMDC,
			Particle::kIsAcceptedHitMETA,
			Particle::kIsAcceptedRK)
      ) selectEpEm = kTRUE;

    return
	selectEpEm
	&& pcand->getBeta() > 0.
	&& pcand->getMomentum() < 2000
	&& pcand->getMetaMatchQuality() < 3.
	&& pcand->getChi2() < 100
        && HParticleTool::isGoodMetaCell(pcand,4,kTRUE)
        && pcand->getSector() != 2
        && !pcand->isAtAnyMdcEdge()
	;

}
//------------------------------------------------------------


Int_t dieleAna(TString inputlist, TString outfile, Int_t nev=-1, TString seqnumlist = "dummy")
{
    TH1::SetDefaultSumw2();
    setupMassCuts("m2_cut.root");
    setupRichCuts("iso_newer.root");
    setupMLPCuts("mlpmom_cutg_gen8_new.root");
    setupSectorCorrs(kPM,"sector_factor_pm.root");
    setupSectorCorrs(kPP,"sector_factor_pp.root");
    setupSectorCorrs(kMM,"sector_factor_mm.root");
    setupPiotrFactor("piotrfactor_gen8_consistent_w6momdepcuts.root");
//    setupPiotrFactor("piotrfactor_gen8.root");

    TH3F *p3DEffEle[6][NWEIGHTS+1]; // mult bins, MLP weights + HC
    TH3F *p3DEffPos[6][NWEIGHTS+1];
    TH3F *p3DAccEle[6][NWEIGHTS+1];
    TH3F *p3DAccPos[6][NWEIGHTS+1];
    readAccEffMatrices(p3DAccEle, p3DAccPos, p3DEffEle, p3DEffPos);

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

    HParticleCand* cand1;
    HParticleCand* cand2;

    HCategory* candCat = (HCategory*)HCategoryManager::getCategory(catParticleCand);

    TString outfileSingle = outfile;
    outfileSingle.ReplaceAll(".root","_single.root");
    HHistMap hMSingle(outfileSingle);
    hMSingle.setSilentFail(kTRUE);
    TString outfilePair = outfile;
    outfilePair.ReplaceAll(".root","_pair.root");
    HHistMap hMPair(outfilePair);
    hMPair.setSilentFail(kTRUE);
    TString outfileKine = outfile;
    outfileKine.ReplaceAll(".root","_kine.root");
    HHistMap hMKine(outfileKine);
    hMKine.setSilentFail(kTRUE);
    TString outfileEvents = outfile;
    outfileEvents.ReplaceAll(".root","_events.txt");
    TString outfileEventsHC = outfile;
    outfileEventsHC.ReplaceAll(".root","_events_hc.txt");

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
    //------------------------------------------------------------------

    //--------------- begin histo booking -----------------------------------------------------
    //------------------------------------------------------------------------------------------
    Int_t binsxmass  = 100;
    Double_t maxmom  = 2000;
    Double_t maxmass = 1000;

    const Int_t nbins = 26;
    Double_t xAxis1[nbins+1] = {0, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.110, 0.130, 0.150, 0.170, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 0.900, 1.};
//    const Int_t nbins = 200;
//    Double_t xAxis1[nbins+1];
//    for (int bb = 0; bb <= nbins; ++bb) {
//        xAxis1[bb] = 0.005*bb;
//    }

    //------------------------------------------------------------------------------------------
    TF1 *yield_pos_fit[6];
    TF1 *yield_neg_fit[6];
    for (int sec = 0; sec < 6; ++sec) {
        yield_pos_fit[sec] = new TF1(Form("yield_pos_fit_sec%i",sec),fline,96,126,2);
        yield_neg_fit[sec] = new TF1(Form("yield_neg_fit_sec%i",sec),fline,96,126,2);
        yield_pos_fit[sec]->SetParameters(yield_pos_pars[sec]);
        yield_neg_fit[sec]->SetParameters(yield_neg_pars[sec]);
    }

    hMSingle.addHist("TH1F","hEvtHour", "hEvtHour",744*60,96,127);
    hMSingle.addHistArray("TH1F","hEpHour_cand",  "hEpHour_sec%i_cand",  "hEpHour_sec%i_cand", 744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hEpHour_ided",  "hEpHour_sec%i_ided",  "hEpHour_sec%i_ided", 744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hEpHour_open",  "hEpHour_sec%i_open",  "hEpHour_sec%i_cand", 744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hEmHour_cand",  "hEmHour_sec%i_cand",  "hEmHour_sec%i_cand", 744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hEmHour_ided",  "hEmHour_sec%i_ided",  "hEmHour_sec%i_ided", 744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hEmHour_open",  "hEmHour_sec%i_open",  "hEmHour_sec%i_cand", 744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hPi0Hour_cand", "hPi0Hour_sec%i_cand", "hPi0Hour_sec%i_cand",744*60,96,127,0,0,0,0,0,0,"","","","",6);
    hMSingle.addHistArray("TH1F","hPi0Hour_open", "hPi0Hour_sec%i_open", "hPi0Hour_sec%i_cand",744*60,96,127,0,0,0,0,0,0,"","","","",6);

    hMSingle.addHist("TH1F","hCounter", "hCounter",20,0,20);
    hMSingle.addHist("TH1F","hCorr", "hCorr",20,0,20);

    hMSingle.addHist("TH1F","ep_angleNearestGoodLept","ep_angleNearestGoodLept",360,0,180);
    hMSingle.addHist("TH1F","ep_angleNearestLeptCand","ep_angleNearestLeptCand",360,0,180);
    hMSingle.addHist("TH1F","ep_angleNearestLeptCandNotGoodLept","ep_angleNearestLeptCandNotGoodLept",360,0,180);
    hMSingle.addHist("TH1F","ep_angleNearestOther","ep_angleNearestOther",360,0,180);

    hMSingle.addHist("TH1F","ep_anglesAllGoodLept","ep_anglesAllGoodLept",360,0,180);
    hMSingle.addHist("TH1F","ep_anglesAllLeptCand","ep_anglesAllLeptCand",360,0,180);
    hMSingle.addHist("TH1F","ep_anglesAllLeptCandNotGoodLept","ep_anglesAllLeptCandNotGoodLept",360,0,180);
    hMSingle.addHist("TH1F","ep_anglesAllOther","ep_anglesAllOther",360,0,180);

    hMSingle.addHist("TH1F","em_angleNearestGoodLept","em_angleNearestGoodLept",360,0,180);
    hMSingle.addHist("TH1F","em_angleNearestLeptCand","em_angleNearestLeptCand",360,0,180);
    hMSingle.addHist("TH1F","em_angleNearestLeptCandNotGoodLept","em_angleNearestLeptCandNotGoodLept",360,0,180);
    hMSingle.addHist("TH1F","em_angleNearestOther","em_angleNearestOther",360,0,180);

    hMSingle.addHist("TH1F","em_anglesAllGoodLept","em_anglesAllGoodLept",360,0,180);
    hMSingle.addHist("TH1F","em_anglesAllLeptCand","em_anglesAllLeptCand",360,0,180);
    hMSingle.addHist("TH1F","em_anglesAllLeptCandNotGoodLept","em_anglesAllLeptCandNotGoodLept",360,0,180);
    hMSingle.addHist("TH1F","em_anglesAllOther","em_anglesAllOther",360,0,180);

    hMSingle.addHist("TH2F","ep_angleNearestGoodLeptMDCdEdx","ep_angleNearestGoodLeptMDCdEdx",360,0,180,100,0,10);
    hMSingle.addHist("TH2F","ep_angleNearestLeptCandMDCdEdx","ep_angleNearestLeptCandMDCdEdx",360,0,180,100,0,10);
    hMSingle.addHist("TH2F","ep_angleNearestLeptCandNotGoodLeptMDCdEdx","ep_angleNearestLeptCandNotGoodLeptMDCdEdx",360,0,180,100,0,10);
    hMSingle.addHist("TH2F","ep_angleNearestOtherMDCdEdx","ep_angleNearestOtherMDCdEdx",360,0,180,100,0,10);

    hMSingle.addHist("TH2F","ep_anglesAllGoodLeptMDCdEdx","ep_anglesAllGoodLeptMDCdEdx",360,0,180,100,0,10);
    hMSingle.addHist("TH2F","ep_anglesAllLeptCandMDCdEdx","ep_anglesAllLeptCandMDCdEdx",360,0,180,100,0,10);
    hMSingle.addHist("TH2F","ep_anglesAllLeptCandNotGoodLeptMDCdEdx","ep_anglesAllLeptCandNotGoodLeptMDCdEdx",360,0,180,100,0,10);
    hMSingle.addHist("TH2F","ep_anglesAllOtherMDCdEdx","ep_anglesAllOtherMDCdEdx",360,0,180,100,0,10);

    hMSingle.addHist("TH2F","em_angleNearestGoodLeptMDCdEdx","em_angleNearestGoodLeptMDCdEdx",360,0,180,100,0,10);
    hMSingle.addHist("TH2F","em_angleNearestLeptCandMDCdEdx","em_angleNearestLeptCandMDCdEdx",360,0,180,100,0,10);
    hMSingle.addHist("TH2F","em_angleNearestLeptCandNotGoodLeptMDCdEdx","em_angleNearestLeptCandNotGoodLeptMDCdEdx",360,0,180,100,0,10);
    hMSingle.addHist("TH2F","em_angleNearestOtherMDCdEdx","em_angleNearestOtherMDCdEdx",360,0,180,100,0,10);

    hMSingle.addHist("TH2F","em_anglesAllGoodLeptMDCdEdx","em_anglesAllGoodLeptMDCdEdx",360,0,180,100,0,10);
    hMSingle.addHist("TH2F","em_anglesAllLeptCandMDCdEdx","em_anglesAllLeptCandMDCdEdx",360,0,180,100,0,10);
    hMSingle.addHist("TH2F","em_anglesAllLeptCandNotGoodLeptMDCdEdx","em_anglesAllLeptCandNotGoodLeptMDCdEdx",360,0,180,100,0,10);
    hMSingle.addHist("TH2F","em_anglesAllOtherMDCdEdx","em_anglesAllOtherMDCdEdx",360,0,180,100,0,10);

    //------------------------------------------------------------------------------------------
    //--------------- Lepton pairs spectra -----------------------------------------------------
    //------------------------------------------------------------------------------------------
    hMPair.addHistArray("TH1F","fill_counter", "fill_counter_w%i_cut%i_sgn%i", "fill_counter_w%i_cut%i_sgn%i",100,0,100,0,0,0,0,0,0,"","","","",4,3,3);
/*
    addHists(hMPair,"hcrossXZ",100,-100,50,100,-50,50);
    addHists(hMPair,"hcrossYZ",100,-100,50,100,-50,50);
    addHists(hMPair,"hclosApprXZ",100,-100,50,100,-50,50);
    addHists(hMPair,"hclosApprYZ",100,-100,50,100,-50,50);
    addHists(hMPair,"hvertXZ",100,-100,50,100,-50,50);
    addHists(hMPair,"hvertYZ",100,-100,50,100,-50,50);
    addHists(hMPair,"hminDistVertZ",100,-50,0,50,0,25);
    addHists(hMPair,"hph1ph2",90,0,360,90,0,360);
    addHists(hMPair,"hDthetaDphi",180,-360,360,180,-90,90);
    addHists(hMPair,"hDthetaDphiSinTheta",180,-360,360,180,-90,90);
    addHists(hMPair,"hphiDilPlane",90,-180,180,90,-180,180);
*/
//    addHists(hMPair,"hthetaAvgDil", 90,0,90,90,0,90);
//    addHists(hMPair,"hphiAvgDil", 360,-360,360,360,-360,360);
    addHists(hMPair,"hsqrtp1p2oa", binsxmass,0,maxmass,90,0,180);
    addHists(hMPair,"htheta1y", 45,0,90,100,0,2);
//    addHists(hMPair,"htheta2y", 45,0,90,100,0,2);
    addHists(hMPair,"hztheta1", 100,-80,20,45,0,90);
//    addHists(hMPair,"hztheta2", 100,-80,20,45,0,90);
    addHists(hMPair,"hthetaoa" ,90,0,90,100,0,180);

    addHists(hMPair,"hy" ,100,0,2);
    addHists(hMPair,"hpt",100,0,1000);
    addHists(hMPair,"hpt_0_200_",100,0,1000);
    addHists(hMPair,"hpt_200_400_",100,0,1000);
    addHists(hMPair,"hpt_400_600_",100,0,1000);
//    addHists(hMPair,"hp1p2diff",2000,-1000,1000);
//    addHists(hMPair,"hth1oAngle" ,45,0,90,90,0,180);
//    addHists(hMPair,"hth1oAngle" ,45,0,90,90,0,180);
//    addHists(hMPair,"hth2oAngle" ,45,0,90,90,0,180);
//    addHists(hMPair,"hth1mass" ,45,0,90,nbins,xAxis1);
//    addHists(hMPair,"hth2mass" ,45,0,90,nbins,xAxis1);
    addHists(hMPair,"hoAnglemass" ,90,0,180,nbins,xAxis1);
    addHists(hMPair,"hoAnglemassNoFactors" ,90,0,180,nbins,xAxis1);
//    addHists(hMPair,"hth1pt" ,45,0,90,120,0,1200);
//    addHists(hMPair,"hth2pt" ,45,0,90,120,0,1200);
    addHists(hMPair,"hoAnglept" ,90,0,180,120,0,1200);
    addHists(hMPair,"hoAngleptNoFactors" ,90,0,180,120,0,1200);
    addHists(hMPair,"hmasspt" ,nbins,xAxis1,120,0,1200);
    addHists(hMPair,"hmassptNoFactors" ,nbins,xAxis1,120,0,1200);
//    addHists(hMPair,"hth1y" ,45,0,90,100,0,2);
//    addHists(hMPair,"hth2y" ,45,0,90,100,0,2);
    addHists(hMPair,"hoAngley" ,90,0,180,100,0,2);
    addHists(hMPair,"hmassy" ,nbins,xAxis1,100,0,2);
    addHists(hMPair,"hpty" ,120,0,1200,100,0,2);
//    addHists(hMPair,"hp1th1", 100,0,1000,45,0,90);
//    addHists(hMPair,"hp2th2", 100,0,1000,45,0,90);
    addHists(hMPair,"hth1th2" ,45,0,90,45,0,90);
//    addHists(hMPair,"hpt1pt2" ,100,0,1000,100,0,1000);
    addHists(hMPair,"hp1p2LowM" ,100,0,1000,100,0,1000);
    addHists(hMPair,"hp1p2MidM" ,100,0,1000,100,0,1000);
    addHists(hMPair,"hp1p2HigM" ,100,0,1000,100,0,1000);
//    addHists(hMPair,"hy1y2" ,100,0,2,100,0,2);
    addHists(hMPair,"hzy" ,100,-80,20,100,0,2);
    addHists(hMPair,"hmtmassFine" ,100,0,1000,100,0,1);
    addHists(hMPair,"hmtmass" ,20,0,1000,nbins,xAxis1);
    addHists(hMPair,"hmass" ,nbins,xAxis1);
    addHists(hMPair,"hmassM" ,nbins,xAxis1);
    addHists(hMPair,"hmassNoFactors" ,nbins,xAxis1);
    addHists(hMPair,"hmassPiotrFactor" ,nbins,xAxis1);
    addHists(hMPair,"hmassSectorFactor" ,nbins,xAxis1);
    addHists(hMPair,"hmassNull" ,nbins,xAxis1);
    addHists(hMPair,"hmassConstWidth",100,0,1000);
    addHists(hMPair,"hoAngle" ,2000,0,200);
    addHists(hMPair,"hmass6sec" ,nbins,xAxis1);
    addHists(hMPair,"hoAngle6sec" ,2000,0,200);
    addHists(hMPair,"hmass5sec" ,nbins,xAxis1);
    addHists(hMPair,"hoAngle5sec" ,2000,0,200);
    addHists(hMPair,"hmass4secd0" ,nbins,xAxis1);
    addHists(hMPair,"hoAngle4secd0" ,2000,0,200);
    addHists(hMPair,"hmass4secd1" ,nbins,xAxis1);
    addHists(hMPair,"hoAngle4secd1" ,2000,0,200);
    addHists(hMPair,"hmass4secd2" ,nbins,xAxis1);
    addHists(hMPair,"hoAngle4secd2" ,2000,0,200);
    addHists(hMPair,"hmtLowM" ,20,0,1000);
    addHists(hMPair,"hmtMidM" ,20,0,1000);
    addHists(hMPair,"hmtHigM" ,20,0,1000);
    addHists(hMPair,"hmtMinusMinvLowM" ,20,0,1000);
    addHists(hMPair,"hmtMinusMinvMidM" ,20,0,1000);
    addHists(hMPair,"hmtMinusMinvHigM" ,20,0,1000);

    addHistsHC(hMPair,"hmtmassFine" ,100,0,1000,100,0,1);
    addHistsHC(hMPair,"hmtmass" ,20,0,1000,nbins,xAxis1);
    addHistsHC(hMPair,"hmasspt" ,nbins,xAxis1,120,0,1200);
    addHistsHC(hMPair,"hmass" ,nbins,xAxis1);
    addHistsHC(hMPair,"hmassConstWidth",100,0,1000);
    addHistsHC(hMPair,"hoAngle" ,2000,0,200);
    addHistsHC(hMPair,"hmtLowM" ,20,0,1000);
    addHistsHC(hMPair,"hmtMidM" ,20,0,1000);
    addHistsHC(hMPair,"hmtHigM" ,20,0,1000);
    addHistsHC(hMPair,"hmtMinusMinvLowM" ,20,0,1000);
    addHistsHC(hMPair,"hmtMinusMinvMidM" ,20,0,1000);
    addHistsHC(hMPair,"hmtMinusMinvHigM" ,20,0,1000);

    hMSingle.addHist("TH2F","thetaMomWon_all","thetaMomWon_all",200,-1000,1000,90,0,90);
    hMSingle.addHist("TH2F","thetaMomWon_sel","thetaMomWon_sel",200,-1000,1000,90,0,90);
    hMSingle.addHist("TH2F","thetaMomWon_mlp","thetaMomWon_mlp",200,-1000,1000,90,0,90);
    hMSingle.addHist("TH2F","thetaMomWon_pid","thetaMomWon_pid",200,-1000,1000,90,0,90);
    hMSingle.addHist("TH2F","ShowerMom_HC_sys0_cutRich", "ShowerMom_HC_sys0_cutRich",400,-1000,1000,200,-500,1500);
    hMSingle.addHist("TH2F","ShowerMom_HC_sys0_cutMassRich", "ShowerMom_HC_sys0_cutMassRich",400,-1000,1000,200,-500,1500);

//    addSingleHists(hMSingle,"_fromPi0");
    addSingleHists(hMSingle,"_noSorter");
//    addSingleHists(hMSingle,"_isUsed");
//    addSingleHists(hMSingle,"_hasShower");
//    addSingleHists(hMSingle,"_leadStop");
    addSingleHists(hMSingle,"",kTRUE);
    addSingleHists(hMSingle,"_cutRich");
    addSingleHists(hMSingle,"_cutRichShower",kTRUE);
    addSingleHists(hMSingle,"_cutRichShowerMass",kTRUE);
    addSingleHists(hMSingle,"_cutRichShowerMass_notIncomp",kTRUE);
    addSingleHists(hMSingle,"_cutRichShowerMass_notIncomp_notNearst",kTRUE);
    addSingleHists(hMSingle,"_cutRichShowerMass_notIncomp_notNearst_notRecur",kTRUE);
    addSingleHists(hMSingle,"_cutRichShowerMassCorr",kTRUE);
    addSingleHists(hMSingle,"_cutMLP");
    addSingleHists(hMSingle,"_cutMLPRich");
    addSingleHists(hMSingle,"_cutMLPRich_notIncomp");
    addSingleHists(hMSingle,"_cutMLPRich_notIncomp_notNearst");
    addSingleHists(hMSingle,"_cutMLPRich_notIncomp_notNearst_notRecur");
    addSingleHists(hMSingle,"_cutMLPRichCorr");

    hMSingle.addHist("TH2F","nepem_hc","nemep_hc",10,-0.5,9.5,10,-0.5,9.5);
    hMSingle.addHist("TH2F","vzx","vzx",200,-200,100,200,-50,50);
    hMSingle.addHist("TH2F","vzr","vzr",200,-200,100,200,0,100);

    TProfile *pmassNP_eff = new TProfile("pmassNP_eff","pmassNP_eff",100,0,1);
    TProfile *pmassPP_eff = new TProfile("pmassPP_eff","pmassPP_eff",100,0,1);
    TProfile *pmassNN_eff = new TProfile("pmassNN_eff","pmassNN_eff",100,0,1);
    pmassNP_eff->Sumw2();
    pmassPP_eff->Sumw2();
    pmassNN_eff->Sumw2();
    //--------------- end histo booking -----------------------------------------------------

    const Float_t shower_cut_pars[6] = {
        -16454.2, 242.44, 177.133
    };

    TF1 *shower_fnc = new TF1("shower_fnc", "[0]/(x-[1])+[2]",0,maxmom);
    for (int par = 0; par < 3; ++par) {
        shower_fnc->SetParameter(par, shower_cut_pars[par]);
    }
    //
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

    EventRecord<NWEIGHTS> record;
    EventRecord<1> record_hc;
    if (!seqnumlist.Contains("dummy")) {
        record.readRecord(seqnumlist.Data());
        TString seqnumlist_hc = seqnumlist;
        seqnumlist_hc.ReplaceAll(".txt","_hc.txt");
        record_hc.readRecord(seqnumlist_hc.Data());
    }

    TStopwatch timer;
    timer.Reset();
    timer.Start();

    Int_t evtsInFile = loop->getEntries();
    if(nev < 0 || nev > evtsInFile ) nev = evtsInFile;

    string currentFileName;
    string currentFileRun;
    bool isSimulation;
    bool isEmbedding;

    std::vector<std::unique_ptr<abstract_functor> > functors;
    std::vector<bool> uniqueInner;
    std::vector<Int_t> counts;
/* !!!
    uniqueInner.push_back(false); functors.push_back(unique_ptr<abstract_functor>(new accept_all()));
    uniqueInner.push_back(true);  functors.push_back(unique_ptr<abstract_functor>(new accept_all()));
    uniqueInner.push_back(false); functors.push_back(unique_ptr<abstract_functor>(new same_inner_seg(cand1)));
    uniqueInner.push_back(false); functors.push_back(unique_ptr<abstract_functor>(new same_inner_same_outer_seg(cand1)));
    uniqueInner.push_back(false); functors.push_back(unique_ptr<abstract_functor>(new same_inner_diff_outer_seg(cand1)));
    uniqueInner.push_back(false); functors.push_back(unique_ptr<abstract_functor>(new accept_leptons()));
    uniqueInner.push_back(false); functors.push_back(unique_ptr<abstract_functor>(new accept_leptons()));
    uniqueInner.push_back(false); functors.push_back(unique_ptr<abstract_functor>(new same_inner_seg_leptons(cand1)));
    uniqueInner.push_back(false); functors.push_back(unique_ptr<abstract_functor>(new same_inner_seg_leptons_like(cand1)));
    uniqueInner.push_back(false); functors.push_back(unique_ptr<abstract_functor>(new same_inner_seg_leptons_unlike(cand1)));
*/
    bool kineHistsNotAdded = kTRUE;
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
            size_t poslastslash = currentFileName.rfind("/");
            size_t poshld = currentFileName.find(".hld");
            if (poslastslash < poshld) {
                currentFileRun = currentFileName.substr(poslastslash+1,poshld-poslastslash-1);
            }
            else {
                currentFileRun = currentFileName.substr(0,poshld);
            }
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

        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(1,"all");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(2,"kGoodTRIGGER");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(3,"kGoodSTART");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(4,"kGoodVertexCand");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(5,"kNoPileUpSTART");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(6,"kNoVETO");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(7,"kGoodSTARTVETO");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(8,"kGoodSTARTMETA");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(9,"#geq 4 secs");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(10,"#geq 5 secs");
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(11,"6 secs");
        for (int mb = 0; mb <= 4; ++mb) {
            hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(12+mb,Form("mult bin %i",mb));
        }
        hMSingle.get("hCounter")->GetXaxis()->SetBinLabel(17,"mult underflow");

        hMSingle.get("hCounter")->Fill(0);
        if (!isSimulation) {
            if (!evtinfo->isGoodEvent(kGoodTRIGGER)) continue;
            hMSingle.get("hCounter")->Fill(1);

            if (!evtinfo->isGoodEvent(kGoodSTART)) continue;
            hMSingle.get("hCounter")->Fill(2);

            hMSingle.get("vzx")->Fill(vz,vx);
            hMSingle.get("vzr")->Fill(vz,TMath::Sqrt(vx*vx+vy*vy));
            if (vz < -50) continue;

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
//            if (sectorCorr != kCorr5to6) continue;
//            if (sectorCorr == kNoCorr || sectorCorr == kCorr5to6) continue;
//            hMSingle.get("hCounter")->Fill(9);
//            if (sectorCorr != kNoCorr) continue;
//            hMSingle.get("hCounter")->Fill(10);
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

/*
	// Select track multiplicity / event
        Int_t mult_sel = 0;
        Int_t good_sec = 0;
        for (int s = 0; s < 6; ++s) {
            if (sectors[s]) {
                mult_sel += evtinfo->getSelectedParticleCandMult(s);
                good_sec += 1;
            }
        }
        mult_sel *= 6./good_sec;

        Int_t mult_bin = 5;
        if (mult_sel >= 17 && mult_sel <  28) mult_bin = 4;
        if (mult_sel >= 28 && mult_sel <  44) mult_bin = 3;
        if (mult_sel >= 44 && mult_sel <  68) mult_bin = 2;
        if (mult_sel >= 68 && mult_sel < 140) mult_bin = 1;
        if (mult_sel >= 140) mult_bin = 0;
*/
        hMSingle.get("hCounter")->Fill(11+mult_bin);
        hMSingle.get("hCorr")->Fill(sectorCorr);

        Int_t binHour     = (dayOfYear-96)*24*60+hour*60+min+1;
        Float_t valueHour = hMSingle.get("hEvtHour")->GetXaxis()->GetBinCenter(binHour);
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

        hMSingle.get("hEvtHour")->Fill(valueHour,1);
        Float_t oaEmbed = -1;
        if (isSimulation) {
            if (kineHistsNotAdded) {
                addKineHists(hMKine);
                kineHistsNotAdded = kFALSE;
            }
            analyzeKine(hMKine,loop,oaEmbed,p3DEffEle[mult_bin][6],p3DEffPos[mult_bin][6],mult_bin);
        }

        Int_t size = candCat->getEntries();

        sorter.cleanUp();
        sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);            // reset all flags for flags (0-28) ,reject,used,lepton
        sorter.fill(selectLeptonsBeta);   // fill only good leptons
//        Int_t n_lepparcnd = sorter.selectBest(HParticleTrackSorter::kIsBestRKRKMETA,HParticleTrackSorter::kIsLepton);
        Int_t n_lepparcnd = sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsLepton);
        sorter.fill(selectHadrons);   // fill only good leptons
//        Int_t n_hadparcnd = sorter.selectBest(HParticleTrackSorter::kIsBestRKRKMETA,HParticleTrackSorter::kIsHadron);
        Int_t n_hadparcnd = sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsHadron);
	//------------------------------------------------------------------------
	Bool_t removedEmbedding[size];
	Bool_t removedRecursive[size];
	// clean vectors and index arrays
        angleNearestBoth.clear();
        angleNearestFull.clear();
        angleNearestSeg.clear();
        nearestIndexSeg.clear();
        countWithin.clear();
        countUniqueWithin.clear();
        countSameSegWithin.clear();
        countSameSegSameOutWithin.clear();
        countSameSegDiffOutWithin.clear();
        countLeptonsWithin.clear();
        countUniqueLeptonsWithin.clear();
        countSameSegLeptonsWithin.clear();
        countSameSegLeptonsLikeWithin.clear();
        countSameSegLeptonsUnlikeWithin.clear();

	sorterFlag.clear();
        richQaFlag.clear();
        betaFlag.clear();
        betaFlagTof.clear();
        betaFlagRpc.clear();
        showerFlag.clear();
	leptonFlagHC.clear();
	removedNearest.clear();
	removedIncomplete.clear();
	removedClosePartner.clear();
        closePartnerCandFlag.clear();
        for (int w = 0; w < NWEIGHTS; ++w) {
            leptonFlagMLP[w].clear();
        }

        int n_ep_cand = 0;
        int n_em_cand = 0;
        int n_ep_hc = 0;
        int n_em_hc = 0;
        int n_lp_hc[1][2] = {{ 0,0 }};
        int n_ep_mlp[NWEIGHTS][2];
        int n_em_mlp[NWEIGHTS][2];
        int n_lp_mlp[NWEIGHTS][2];
        int n_ep_mlp_notIncomp[NWEIGHTS][2];
        int n_em_mlp_notIncomp[NWEIGHTS][2];
        int n_lp_mlp_notIncomp[NWEIGHTS][2];
        int n_ep_mlp_notIncomp_notNearst[NWEIGHTS][2];
        int n_em_mlp_notIncomp_notNearst[NWEIGHTS][2];
        int n_lp_mlp_notIncomp_notNearst[NWEIGHTS][2];
        int n_ep_mlp_notIncomp_notNearst_notRecurs[NWEIGHTS][2];
        int n_em_mlp_notIncomp_notNearst_notRecurs[NWEIGHTS][2];
        int n_lp_mlp_notIncomp_notNearst_notRecurs[NWEIGHTS][2];
        for (int w = 0; w < NWEIGHTS; ++w) {
            n_ep_mlp[w][0] = 0;
            n_em_mlp[w][0] = 0;
            n_lp_mlp[w][0] = 0;
            n_ep_mlp_notIncomp[w][0] = 0;
            n_em_mlp_notIncomp[w][0] = 0;
            n_lp_mlp_notIncomp[w][0] = 0;
            n_ep_mlp_notIncomp_notNearst[w][0] = 0;
            n_em_mlp_notIncomp_notNearst[w][0] = 0;
            n_lp_mlp_notIncomp_notNearst[w][0] = 0;
            n_ep_mlp_notIncomp_notNearst_notRecurs[w][0] = 0;
            n_em_mlp_notIncomp_notNearst_notRecurs[w][0] = 0;
            n_lp_mlp_notIncomp_notNearst_notRecurs[w][0] = 0;
            n_ep_mlp[w][1] = 0;
            n_em_mlp[w][1] = 0;
            n_lp_mlp[w][1] = 0;
            n_ep_mlp_notIncomp[w][1] = 0;
            n_em_mlp_notIncomp[w][1] = 0;
            n_lp_mlp_notIncomp[w][1] = 0;
            n_ep_mlp_notIncomp_notNearst[w][1] = 0;
            n_em_mlp_notIncomp_notNearst[w][1] = 0;
            n_lp_mlp_notIncomp_notNearst[w][1] = 0;
            n_ep_mlp_notIncomp_notNearst_notRecurs[w][1] = 0;
            n_em_mlp_notIncomp_notNearst_notRecurs[w][1] = 0;
            n_lp_mlp_notIncomp_notNearst_notRecurs[w][1] = 0;
        }

        doTrackCorr(trackcor, candCat);
        // SETTNG FLAGS
	for(Int_t j = 0; j < size; j ++){
	    sorterFlag[j]           = kFALSE;
            richQaFlag[j]           = kFALSE;
            betaFlag[j]             = kFALSE;
            betaFlagTof[j]          = kFALSE;
            betaFlagRpc[j]          = kFALSE;
            showerFlag[j]           = kFALSE;
	    leptonFlagHC[j]         = kFALSE;
	    removedNearest[j]       = kFALSE;
	    removedIncomplete[j]    = kFALSE;
	    removedClosePartner[j]  = kFALSE;
            removedRecursive[j]     = kFALSE;
            removedEmbedding[j]     = kFALSE;
            closePartnerCandFlag[j] = kFALSE;
            for (int w = 0; w < NWEIGHTS; ++w) {
                leptonFlagMLP[w][j] = kFALSE;
            }
            // Track candidates (Leptons AND Hadrons)
	    cand1 = HCategoryManager::getObject(cand1,candCat,j);
            closePartnerCandFlag[j] = selectClosePairCands(cand1);

            if (isEmbedding) {
                HParticleCandSim *candSim1 = dynamic_cast<HParticleCandSim *>(cand1); 
                if (candSim1 != NULL) {
                    if (!((candSim1->getGeantPID() == 2 || candSim1->getGeantPID() == 3) && candSim1->getGeantParentPID() == -1)) {
                        removedEmbedding[j] = kTRUE;
                    }
                    if (candSim1->getGeantPID() != 2 && candSim1->getGeantPID() != 3) {
                        removedEmbedding[j] = kTRUE;
                    }
                }
            }

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
            Float_t richQa  = cand1->getRichMatchingQuality();
//            Float_t phi     = cand1->getPhi();
            Float_t m2      = _mom*_mom*(1-_beta*_beta)/(_beta*_beta);
            Int_t charge    = cand1->getCharge();
            Int_t sys       = cand1->getSystemUsed();
            Int_t sec       = cand1->getSector();

            if (!selectLeptonsBeta(cand1)) continue;
	    cand1->calc4vectorProperties(); // lepton by default
            nearestIndexSeg[j] = findNearestNeighbor(cand1,angleNearestBoth,angleNearestSeg,angleNearestFull,candCat);

/* !!!
            countWithinAngle(cand1,counts,functors,uniqueInner,4);

            countWithin[cand1->getIndex()] = counts[0];
            countUniqueWithin[cand1->getIndex()] = counts[1];
            countSameSegWithin[cand1->getIndex()] = counts[2];
            countSameSegSameOutWithin[cand1->getIndex()] = counts[3];
            countSameSegDiffOutWithin[cand1->getIndex()] = counts[4];
            countLeptonsWithin[cand1->getIndex()] = counts[5];
            countUniqueLeptonsWithin[cand1->getIndex()] = counts[6];
            countSameSegLeptonsWithin[cand1->getIndex()] = counts[7];
            countSameSegLeptonsLikeWithin[cand1->getIndex()] = counts[8];
            countSameSegLeptonsUnlikeWithin[cand1->getIndex()] = counts[9];
*/
/*
            vector <HParticleCand *> found;
            countWithin[cand1->getIndex()] = countWithinAngle(cand1,found,4,false,accept_all);
            countUniqueWithin[cand1->getIndex()] = countWithinAngle(cand1,found,4,true,accept_all);
            countSameSegWithin[cand1->getIndex()] = countWithinAngle(cand1,found,4,false,same_inner_seg(cand1));
            countSameSegSameOutWithin[cand1->getIndex()] = countWithinAngle(cand1,found,4,false,same_inner_same_outer_seg(cand1));
            countSameSegDiffOutWithin[cand1->getIndex()] = countWithinAngle(cand1,found,4,false,same_inner_diff_outer_seg(cand1));
            countLeptonsWithin[cand1->getIndex()] = countWithinAngle(cand1,found,4,false,accept_leptons);
            countUniqueLeptonsWithin[cand1->getIndex()] = countWithinAngle(cand1,found,4,true,accept_leptons);
            countSameSegLeptonsWithin[cand1->getIndex()] = countWithinAngle(cand1,found,4,false,same_inner_seg_leptons(cand1));
            countSameSegLeptonsLikeWithin[cand1->getIndex()] = countWithinAngle(cand1,found,4,false,same_inner_seg_leptons_like(cand1));
            countSameSegLeptonsUnlikeWithin[cand1->getIndex()] = countWithinAngle(cand1,found,4,false,same_inner_seg_leptons_unlike(cand1));
*/
            hMSingle.get("thetaMomWon_all")->Fill(_mom*cand1->getCharge(),_theta);
            if (!cand1->isFlagBit(Particle::kIsLepton)) continue;
            //if (cand1->getMomentum() < 100) continue;
            sorterFlag[j] = kTRUE;
            hMSingle.get("thetaMomWon_sel")->Fill(_mom*cand1->getCharge(),_theta);

            for (int w = 0; w < NWEIGHTS; ++w) {
                MLPs[w][j] = reader[w][sys]->EvaluateMVA( "MLP" );
            }
            ///////////////////////////////////////////////////// //////////////////////////////
            //                                                    
            // Hard cut flags                                     
            //
            ///////////////////////////////////////////////////////////////////////////////////
            // Ring matching
            //if (richQa < 3.5e-3*_mom + 1) richQaFlag[j]=kTRUE;
            //if (richQa < 1) richQaFlag[j]=kTRUE;

            bool cutrich;
            if (sys == 0) {
                cutrich = richQa < hcutRichQaSys0->Interpolate(_mom*charge);
            }
            else {
                cutrich = richQa < hcutRichQaSys1->Interpolate(_mom*charge);
            }
            if (cutrich) richQaFlag[j]=kTRUE;
            // beta
            Int_t charge_index;
            if (charge < 0) {
                charge_index = 0;
            }
            else {
                charge_index = 1;
            }
            if (sys==1) {
                if (m2 < hcutMassSys1->Interpolate(_mom*charge)) betaFlagTof[j] = kTRUE;
            }
            if (sys==0) {
                if (m2 < hcutMassSys0->Interpolate(_mom*charge)) betaFlagRpc[j] = kTRUE;
            }
            if (betaFlagTof[j] || betaFlagRpc[j]) betaFlag[j] = kTRUE;

            // Shower
            if (isGoodShower(cand1,shower_fnc)) showerFlag[j] = kTRUE;
            if (richQaFlag[j] && betaFlag[j] && showerFlag[j]) {
                leptonFlagHC[j]=kTRUE;
            }
            for (int w = 0; w < NWEIGHTS; ++w) {
                if (richQaFlag[j] && isGoodMLP(cand1,MLPs[w][j],w)) {
                    leptonFlagMLP[w][j]=kTRUE;
/*
                    for (int jj = 0; jj < size; ++jj) {
                        HParticleCand *candOther;
                        candOther = HCategoryManager::getObject(candOther,candCat,j);
                        if (jj != j && cand1->getMetaHitInd() == candOther->getMetaHitInd() && (candOther->getChi2() > 40 || candOther->getChi2() < 0)) {
                            removedIncomplete[j] = kTRUE;
                        }
                    }
*/
                    for (int jj = 0; jj < size; ++jj) {
                        if (closePartnerCandFlag[jj]) {
                            HParticleCand *candOther;
                            candOther = HCategoryManager::getObject(candOther,candCat,j);
                            Float_t angleToOther = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), candOther->getPhi(), candOther->getTheta());
                            if (angleToOther < 8) {
                                removedClosePartner[j] = kTRUE;
                            }
                        }
                    }
                    if (nearestIndexSeg[j] > -1) {
                        HParticleCand *nearestCand;
                        nearestCand = HCategoryManager::getObject(nearestCand,candCat,nearestIndexSeg[j]);

                        if (cand1->getMetaHitInd() == nearestCand->getMetaHitInd() && (nearestCand->getChi2() > 40 || nearestCand->getChi2() < 0)) {
                            removedIncomplete[j] = kTRUE;
                        }

                        if (/*!removedIncomplete[j] && */(!(angleNearestSeg[j] < 0 || angleNearestSeg[j] > 6))) {
                            removedNearest[j] = kTRUE;
                        }
                    }
                    
                }
            }
	    //*************************************************************
	}  // END SETTING FLAGS
        // RECURSIVE CUT
        for(Int_t j = 0; j < size; j ++){ // first candidate
            if (!sorterFlag[j]) continue;
            cand1 = HCategoryManager::getObject(cand1,candCat,j);

            for(Int_t k = j+1; k < size; k ++){ // second candidate, k > j condition applied below
                if (!sorterFlag[k]) continue;
                cand2 = HCategoryManager::getObject(cand2,candCat,k);

                //UInt_t flags;
                //HParticleTool::setPairFlags(flags,cand2,cand1);
                //if(!HParticleTool::evalPairsFlags(kPairCase1,flags)) continue;

                if (k > j && leptonFlagMLP[6][j] && leptonFlagMLP[6][k]) {

                    TLorentzVector dilep = (*cand1) + (*cand2);

                    if((cand1->getCharge()==1  && cand2->getCharge()==-1) ||
                            (cand1->getCharge()==-1 && cand2->getCharge()==1)) {
                        Float_t oAngle   = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), cand2->getPhi(), cand2->getTheta());
                        //if (!removedIncomplete[j] && !removedIncomplete[k]) {
                        //if (!removedNearest[j] && !removedNearest[k]) {
                        if (true) {
                            if (/*!removedRecursive[j] && !removedRecursive[k] && */oAngle < oAngleCut) {
                                removedRecursive[j] = kTRUE;
                                removedRecursive[k] = kTRUE;
                            }
                        }
                    }
/*
                    if (cand1->getCharge() == cand2->getCharge()) {
                        TLorentzVector dilep = (*cand1) + (*cand2);
                        if (dilep.Rapidity() > 1.2 && dilep.Perp() < 200) {
                            removedRecursive[j] = kTRUE;
                            removedRecursive[k] = kTRUE;
                        }
                    }
*/
                }
            }
        } // END RECURSIVE CUT
        // LOOP TO COUNT LEPTONS
	for(Int_t j = 0; j < size; j ++){
            if (isEmbedding && removedEmbedding[j]) continue;
            if (!sorterFlag[j]) continue;
	    cand1 = HCategoryManager::getObject(cand1,candCat,j);
            Int_t sys = cand1->getSystemUsed();

            if (cand1->getCharge() == 1) n_ep_cand++;
            if (cand1->getCharge() ==-1) n_em_cand++;
            if (leptonFlagHC[j]) {
                if (cand1->getCharge() == 1) n_ep_hc++;
                if (cand1->getCharge() ==-1) n_em_hc++;
                n_lp_hc[0][sys]++;
            }
            for (int w = 0; w < NWEIGHTS; ++w) {
                if (leptonFlagMLP[w][j]) {
                    if (cand1->getCharge() == 1) n_ep_mlp[w][sys]++;
                    if (cand1->getCharge() ==-1) n_em_mlp[w][sys]++;
                    n_lp_mlp[w][sys]++;
                    if (!removedIncomplete[j]) {
                        if (cand1->getCharge() == 1) n_ep_mlp_notIncomp[w][sys]++;
                        if (cand1->getCharge() ==-1) n_em_mlp_notIncomp[w][sys]++;
                        n_lp_mlp_notIncomp[w][sys]++;
                    }
                    if (!removedIncomplete[j] && !removedNearest[j]) {
                        if (cand1->getCharge() == 1) n_ep_mlp_notIncomp_notNearst[w][sys]++;
                        if (cand1->getCharge() ==-1) n_em_mlp_notIncomp_notNearst[w][sys]++;
                        n_lp_mlp_notIncomp_notNearst[w][sys]++;
                    }
                    if (!removedIncomplete[j] && !removedNearest[j] && !removedRecursive[j]) {
                        if (cand1->getCharge() == 1) n_ep_mlp_notIncomp_notNearst_notRecurs[w][sys]++;
                        if (cand1->getCharge() ==-1) n_em_mlp_notIncomp_notNearst_notRecurs[w][sys]++;
                        n_lp_mlp_notIncomp_notNearst_notRecurs[w][sys]++;
                    }
                }
            }
	    //*************************************************************
	}  // END LOOP TO COUNT LEPTONS
        if (n_ep_hc + n_em_hc > 1) {
            hMSingle.get("nepem_hc")->Fill(n_ep_hc,n_em_hc);
        }
        record.addEvent(currentFileRun, loop->getEventHeader()->getEventSeqNumber(), n_lp_mlp);
        record_hc.addEvent(currentFileRun, loop->getEventHeader()->getEventSeqNumber(), n_lp_hc);
        // END OF PREPARATIONS

        // SINLE LEPTON FILLING HISTOGRAMS
	for(Int_t j = 0; j < size; j ++){
            if (isEmbedding && removedEmbedding[j]) continue;
	    cand1 = HCategoryManager::getObject(cand1,candCat,j);

	    _theta         = cand1->getTheta();
	    _mom           = cand1->getMomentum();
	    Float_t phi    = cand1->getPhi();
	    _showerDq      = cand1->getShowerDeltaSum();
            Int_t sys      = cand1->getSystemUsed();
            Int_t sec      = cand1->getSector();

            if(cand1->getSystemUsed()==-1)  continue;
/*
            if(cand1->isFlagBit(Particle::kIsUsed)) {
                fillSingle(hMSingle,cand1,"_isUsed",mult_bin,nearestIndexSeg[j]);
                Float_t showerSum0 = cand1->getShowerSum0();
                Float_t showerSum1 = cand1->getShowerSum1();
                Float_t showerSum2 = cand1->getShowerSum2();
                if (showerSum0 > 0) {
                    fillSingle(hMSingle,cand1,"_hasShower",mult_bin,nearestIndexSeg[j]);
                }
                if (showerSum0 > 0 && showerSum2 == 0) {
                    fillSingle(hMSingle,cand1,"_leadStop",mult_bin,nearestIndexSeg[j]);
                }
            }
*/
            if(!selectLeptonsBeta(cand1)) continue;
            fillSingle(hMSingle,cand1,"_noSorter",mult_bin,nearestIndexSeg[j]);


            if(!cand1->isFlagBit(Particle::kIsLepton)) continue;
            if (cand1->getMomentum() < 100) continue;

            Float_t weight_yield = 1;
            if (cand1->getCharge() == -1) {
//                weight_yield = yield_neg_weight[sec];
                hMSingle.get("hEmHour_cand",sec)->Fill(valueHour,weight_yield);
            }
            if (cand1->getCharge() == 1) {
//                weight_yield = yield_pos_weight[sec];
                hMSingle.get("hEpHour_cand",sec)->Fill(valueHour,weight_yield);
            }

            if (isSimulation) {
                hMKine.get2("h_reco_PT",0)->Fill(phi,_theta);
                hMKine.get2("h_reco_MT",0)->Fill(_mom,_theta);
                hMKine.get2("h_reco_PM",0)->Fill(phi,_mom);
            }
            if (_showerDq != -1) {
                if (richQaFlag[j]) {
                    if (sys==0) { 
                        hMSingle.get("ShowerMom_HC_sys0_cutRich")->Fill(_mom*cand1->getCharge(),_showerDq);
                    }
                }
                if (richQaFlag[j] && betaFlag[j]) {
                    if (sys==0) { 
                        hMSingle.get("ShowerMom_HC_sys0_cutMassRich")->Fill(_mom*cand1->getCharge(),_showerDq);
                    }
                }
            }
            if (richQaFlag[j] && showerFlag[j]) {
                if ((sys == 0 && record_hc.checkEvent(currentFileRun, loop->getEventHeader()->getEventSeqNumber(), 0, 0)) ||
                    (sys == 1 && record_hc.checkEvent(currentFileRun, loop->getEventHeader()->getEventSeqNumber(), 0, 0))) {
                    fillSingle(hMSingle,cand1,"_cutRichShower",mult_bin,nearestIndexSeg[j],NWEIGHTS);
                }
            }
            if (richQaFlag[j] && showerFlag[j] && betaFlag[j]) {
                if ((sys == 0 && record_hc.checkEvent(currentFileRun, loop->getEventHeader()->getEventSeqNumber(), 0, 0)) ||
                    (sys == 1 && record_hc.checkEvent(currentFileRun, loop->getEventHeader()->getEventSeqNumber(), 0, 0))) {
                    fillSingle(hMSingle,cand1,"_cutRichShowerMass",mult_bin,nearestIndexSeg[j],NWEIGHTS);
                    if (cand1->getCharge() == -1) {
                        fillSingle(hMSingle,cand1,"_cutRichShowerMassCorr",mult_bin,nearestIndexSeg[j],NWEIGHTS,p3DEffEle[mult_bin][NWEIGHTS]);
                    }
                    if (cand1->getCharge() == 1) {
                        fillSingle(hMSingle,cand1,"_cutRichShowerMassCorr",mult_bin,nearestIndexSeg[j],NWEIGHTS,p3DEffPos[mult_bin][NWEIGHTS]);
                    }
                    if (!removedIncomplete[j]) {
                        fillSingle(hMSingle,cand1,"_cutRichShowerMass_notIncomp",mult_bin,nearestIndexSeg[j],NWEIGHTS);
                    }
                    if (!removedIncomplete[j] && !removedNearest[j]) {
                        fillSingle(hMSingle,cand1,"_cutRichShowerMass_notIncomp_notNearst",mult_bin,nearestIndexSeg[j],NWEIGHTS);
                    }
                    if (!removedIncomplete[j] && !removedNearest[j] && !removedRecursive[j]) {
                        fillSingle(hMSingle,cand1,"_cutRichShowerMass_notIncomp_notNearst_notRecur",mult_bin,nearestIndexSeg[j],NWEIGHTS);
                    }
                }
            }
            if (richQaFlag[j]) {
                if (isSimulation) {
                    hMKine.get2("h_reco_PT_richQa",0)->Fill(phi,_theta); 
                    hMKine.get2("h_reco_MT_richQa",0)->Fill(_mom,_theta);
                    hMKine.get2("h_reco_PM_richQa",0)->Fill(phi,_mom);
                }
            }
            fillSingle(hMSingle,cand1,"",mult_bin,nearestIndexSeg[j],NWEIGHTS);
            for (int w = 0; w < NWEIGHTS; ++w) {
                if (richQaFlag[j]) {
                    if ((sys == 0 && record.checkEvent(currentFileRun, loop->getEventHeader()->getEventSeqNumber(), w, 2)) ||
                        (sys == 1 && record.checkEvent(currentFileRun, loop->getEventHeader()->getEventSeqNumber(), 6, w))) {
                        fillSingle(hMSingle,cand1,"_cutRich",mult_bin,nearestIndexSeg[j],w);
                    }
                }
                if (leptonFlagMLP[w][j]) {
                    if (w == 6) {
//                        cout << weight_yield << endl;
                        if (cand1->getCharge() == -1) {
                            hMSingle.get("hEmHour_ided",sec)->Fill(valueHour,weight_yield);
                        }
                        if (cand1->getCharge() == 1) {
                            hMSingle.get("hEpHour_ided",sec)->Fill(valueHour,weight_yield);
                        }
                        Float_t angleNearestGoodLept = 999;
                        Float_t angleNearestLeptCand = 999;
                        Float_t angleNearestLeptCandNotGoodLept = 999;
                        Float_t angleNearestOther = 999;
                        
                        Int_t indNearestGoodLept = -1;
                        Int_t indNearestLeptCand = -1;
                        Int_t indNearestLeptCandNotGoodLept = -1;
                        Int_t indNearestOther = -1;
                        
                        for (Int_t jj = 0; jj < size; ++jj) {
                            if (jj == j) continue;
                            HParticleCand *candOther;
                            candOther = HCategoryManager::getObject(candOther,candCat,jj);
                            Float_t angleToOther = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), candOther->getPhi(), candOther->getTheta());
                            
                            if (leptonFlagMLP[6][jj]) {
                                if (angleToOther < angleNearestGoodLept) {
                                    angleNearestGoodLept = angleToOther;
                                    indNearestGoodLept = jj;
                                }
                                if (cand1->getCharge() > 0) {
                                    hMSingle.get("ep_anglesAllGoodLept")->Fill(angleToOther);
                                    hMSingle.get("ep_anglesAllGoodLeptMDCdEdx")->Fill(angleToOther,candOther->getMdcdEdx());
                                } else{
                                    hMSingle.get("em_anglesAllGoodLept")->Fill(angleToOther);
                                    hMSingle.get("em_anglesAllGoodLeptMDCdEdx")->Fill(angleToOther,candOther->getMdcdEdx());
                                }
                            }
                            if (closePartnerCandFlag[jj]) {
                                if (angleToOther < angleNearestLeptCand) {
                                    angleNearestLeptCand = angleToOther;
                                    indNearestLeptCand = jj;
                                }
                                if (cand1->getCharge() > 0) {
                                    hMSingle.get("ep_anglesAllLeptCand")->Fill(angleToOther);
                                    hMSingle.get("ep_anglesAllLeptCandMDCdEdx")->Fill(angleToOther,candOther->getMdcdEdx());
                                } else{
                                    hMSingle.get("em_anglesAllLeptCand")->Fill(angleToOther);
                                    hMSingle.get("em_anglesAllLeptCandMDCdEdx")->Fill(angleToOther,candOther->getMdcdEdx());
                                }
                            }
                            if (closePartnerCandFlag[jj] && !leptonFlagMLP[6][jj]) {
                                if (angleToOther < angleNearestLeptCandNotGoodLept) {
                                    angleNearestLeptCandNotGoodLept = angleToOther;
                                    indNearestLeptCandNotGoodLept = jj;
                                }
                                if (cand1->getCharge() > 0) {
                                    hMSingle.get("ep_anglesAllLeptCandNotGoodLept")->Fill(angleToOther);
                                    hMSingle.get("ep_anglesAllLeptCandNotGoodLeptMDCdEdx")->Fill(angleToOther,candOther->getMdcdEdx());
                                } else{
                                    hMSingle.get("em_anglesAllLeptCandNotGoodLept")->Fill(angleToOther);
                                    hMSingle.get("em_anglesAllLeptCandNotGoodLeptMDCdEdx")->Fill(angleToOther,candOther->getMdcdEdx());
                                }
                            }
                            if (!closePartnerCandFlag[jj] && !leptonFlagMLP[6][jj]) {
                                if (angleToOther < angleNearestOther) {
                                    angleNearestOther = angleToOther;
                                    indNearestOther = jj;
                                }
                                if (cand1->getCharge() > 0) {
                                    hMSingle.get("ep_anglesAllOther")->Fill(angleToOther);
                                    hMSingle.get("ep_anglesAllOtherMDCdEdx")->Fill(angleToOther,candOther->getMdcdEdx());
                                } else{
                                    hMSingle.get("em_anglesAllOther")->Fill(angleToOther);
                                    hMSingle.get("em_anglesAllOtherMDCdEdx")->Fill(angleToOther,candOther->getMdcdEdx());
                                }
                            }
                        }
                        if (cand1->getCharge() > 0) {
                            hMSingle.get("ep_angleNearestGoodLept")->Fill(angleNearestGoodLept);
                            hMSingle.get("ep_angleNearestLeptCand")->Fill(angleNearestLeptCand);
                            hMSingle.get("ep_angleNearestLeptCandNotGoodLept")->Fill(angleNearestLeptCandNotGoodLept);
                            hMSingle.get("ep_angleNearestOther")->Fill(angleNearestGoodLept);

                            HParticleCand *candOther;
                            candOther = HCategoryManager::getObject(candOther,candCat,indNearestGoodLept);
                            if (candOther) {
                                hMSingle.get("ep_angleNearestGoodLeptMDCdEdx")->Fill(angleNearestGoodLept,candOther->getMdcdEdx());
                            }
                            candOther = HCategoryManager::getObject(candOther,candCat,indNearestLeptCand);
                            if (candOther) {
                                hMSingle.get("ep_angleNearestLeptCandMDCdEdx")->Fill(angleNearestLeptCand,candOther->getMdcdEdx());
                            }
                            candOther = HCategoryManager::getObject(candOther,candCat,indNearestLeptCandNotGoodLept);
                            if (candOther) {
                                hMSingle.get("ep_angleNearestLeptCandNotGoodLeptMDCdEdx")->Fill(angleNearestLeptCandNotGoodLept,candOther->getMdcdEdx());
                            }
                            candOther = HCategoryManager::getObject(candOther,candCat,indNearestOther);
                            if (candOther) {
                                hMSingle.get("ep_angleNearestOtherMDCdEdx")->Fill(angleNearestGoodLept,candOther->getMdcdEdx());
                            }
                        } else{                                                               
                            hMSingle.get("em_angleNearestGoodLept")->Fill(angleNearestGoodLept);
                            hMSingle.get("em_angleNearestLeptCand")->Fill(angleNearestLeptCand);
                            hMSingle.get("em_angleNearestLeptCandNotGoodLept")->Fill(angleNearestLeptCandNotGoodLept);
                            hMSingle.get("em_angleNearestOther")->Fill(angleNearestGoodLept);

                            HParticleCand *candOther;
                            candOther = HCategoryManager::getObject(candOther,candCat,indNearestGoodLept);
                            if (candOther) {
                                hMSingle.get("em_angleNearestGoodLeptMDCdEdx")->Fill(angleNearestGoodLept,candOther->getMdcdEdx());
                            }
                            candOther = HCategoryManager::getObject(candOther,candCat,indNearestLeptCand);
                            if (candOther) {
                                hMSingle.get("em_angleNearestLeptCandMDCdEdx")->Fill(angleNearestLeptCand,candOther->getMdcdEdx());
                            }
                            candOther = HCategoryManager::getObject(candOther,candCat,indNearestLeptCandNotGoodLept);
                            if (candOther) {
                                hMSingle.get("em_angleNearestLeptCandNotGoodLeptMDCdEdx")->Fill(angleNearestLeptCandNotGoodLept,candOther->getMdcdEdx());
                            }
                            candOther = HCategoryManager::getObject(candOther,candCat,indNearestOther);
                            if (candOther) {
                                hMSingle.get("em_angleNearestOtherMDCdEdx")->Fill(angleNearestGoodLept,candOther->getMdcdEdx());
                            }
                        }
                    }
                    
                    if ((sys == 0 && record.checkEvent(currentFileRun, loop->getEventHeader()->getEventSeqNumber(), w, 2)) ||
                        (sys == 1 && record.checkEvent(currentFileRun, loop->getEventHeader()->getEventSeqNumber(), 6, w))) {
                        fillSingle(hMSingle,cand1,"_cutMLPRich",mult_bin,nearestIndexSeg[j],w);
                        if (!removedIncomplete[j]) {
                            fillSingle(hMSingle,cand1,"_cutMLPRich_notIncomp",mult_bin,nearestIndexSeg[j],w);
                        }
                        if (!removedIncomplete[j] && !removedNearest[j]) {
                            fillSingle(hMSingle,cand1,"_cutMLPRich_notIncomp_notNearst",mult_bin,nearestIndexSeg[j],w);
                        }
                        if (!removedIncomplete[j] && !removedNearest[j] && !removedRecursive[j]) {
                            fillSingle(hMSingle,cand1,"_cutMLPRich_notIncomp_notNearst_notRecur",mult_bin,nearestIndexSeg[j],w);
                        }
                        if (cand1->getCharge() == -1) {
                            fillSingle(hMSingle,cand1,"_cutMLPRichCorr",mult_bin,nearestIndexSeg[j],w,p3DEffEle[mult_bin][w]);
                            if (isSimulation) {
                                hMKine.get3("h_reco_MPT_em",w)->Fill(phi,_theta,_mom);
                            }
                        }
                        if (cand1->getCharge() == 1) {
                            fillSingle(hMSingle,cand1,"_cutMLPRichCorr",mult_bin,nearestIndexSeg[j],w,p3DEffPos[mult_bin][w]);
                            if (isSimulation) {
                                hMKine.get3("h_reco_MPT_ep",w)->Fill(phi,_theta,_mom);
                            }
                        }
                        if (isSimulation) {
                            hMKine.get2("h_reco_PT_mlp_richQa",w)->Fill(phi,_theta);
                            hMKine.get2("h_reco_MT_mlp_richQa",w)->Fill(_mom,_theta);
                            hMKine.get2("h_reco_PM_mlp_richQa",w)->Fill(phi,_mom);

                            Float_t weight = getEfficiencyFactor(p3DEffEle[mult_bin][w], p3DEffPos[mult_bin][w], _mom, _theta, phi, cand1->getCharge());
                            hMKine.get2("h_corr_PT_mlp_richQa",w)->Fill(phi,_theta,weight);
                            hMKine.get2("h_corr_MT_mlp_richQa",w)->Fill(_mom,_theta,weight);
                            hMKine.get2("h_corr_PM_mlp_richQa",w)->Fill(phi,_mom,weight);
                        }
                    }
                }
                if (isGoodMLP(cand1,MLPs[w][j],w)) {
                    if ((sys == 0 && record.checkEvent(currentFileRun, loop->getEventHeader()->getEventSeqNumber(), w, 2)) ||
                        (sys == 1 && record.checkEvent(currentFileRun, loop->getEventHeader()->getEventSeqNumber(), 6, w))) {
                        fillSingle(hMSingle,cand1,"_cutMLP",mult_bin,nearestIndexSeg[j],w);
                        if (isSimulation) {
                            hMKine.get2("h_reco_PT_mlp",w)->Fill(phi,_theta);
                            hMKine.get2("h_reco_MT_mlp",w)->Fill(_mom,_theta);
                            hMKine.get2("h_reco_PM_mlp",w)->Fill(phi,_mom);
                        }
                    }
                }
            }
            //*************************************************************
        }  // EBD SINGLE LEPTONS FILLING HISTOGRAMS

	// PAIRS FILLING HISTOGRAMS
        set<Int_t> usedIndices;
	for(Int_t j = 0; j < size-1; j ++){ // first candidate
	    for(Int_t k = j+1; k < size; k ++){ // second candidate
                cand1 = HCategoryManager::getObject(cand1,candCat,j);
		cand2 = HCategoryManager::getObject(cand2,candCat,k);
                Int_t sysj = cand1->getSystemUsed();
                Int_t sysk = cand2->getSystemUsed();

                Int_t chg = kPM;
                if (cand1->getCharge() == cand2->getCharge()) {
                    if (cand1->getCharge() > 0) chg = kPP;
                    if (cand1->getCharge() < 0) chg = kMM;
                }
                TH2F *pSectorCorr;
                if (sectorCorr != kNoCorr) {
                    pSectorCorr = hSectorCorr[chg][sectorCorr-1];
                }
                else {
                    pSectorCorr = NULL;
                }
                if (!sorterFlag[j]) continue;
                if (!sorterFlag[k]) continue;
                if (isEmbedding && removedEmbedding[j]) continue;
                if (isEmbedding && removedEmbedding[k]) continue;

                UInt_t flags;
                HParticleTool::setPairFlags(flags,cand2,cand1);
                if(!HParticleTool::evalPairsFlags(kPairCase1,flags)) continue;

                TLorentzVector dilep = (*cand1) + (*cand2);
                Float_t invM     = dilep.M();

                Float_t weight_yield_pair = 1;
                if (cand1->getCharge() == cand2->getCharge()) {
                    if (cand1->getCharge() == 1) {
                        weight_yield_pair = yield_pos_weight[cand1->getSector()]*yield_pos_weight[cand1->getSector()];
                    }
                    if (cand1->getCharge() == -1) {
                        weight_yield_pair = yield_neg_weight[cand1->getSector()]*yield_neg_weight[cand1->getSector()];
                    }
                    hMSingle.get("hPi0Hour_cand",cand1->getSector())->Fill(valueHour,-weight_yield_pair);
                }
                else {
                    weight_yield_pair = yield_neg_weight[cand1->getSector()]*yield_pos_weight[cand1->getSector()];
                }
                if (leptonFlagHC[j] && leptonFlagHC[k]) {

                    Float_t oAngle   = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), cand2->getPhi(), cand2->getTheta());

                    fillPairSpectraHC(hMPair,cand1,cand2,"",p3DEffEle[0][13],p3DEffPos[0][13],pSectorCorr,ALL_MULTS,NO_CUT,weight_yield_pair,4);
//                    fillPairSpectraHC(hMPair,cand1,cand2,"",p3DEffEle[mult_bin][13],p3DEffPos[mult_bin][13],pSectorCorr,ALL_MULTS,NO_CUT,weight_yield_pair,4);
                    if (mult_bin > 0 && mult_bin < 5) {
                        fillPairSpectraHC(hMPair,cand1,cand2,"",p3DEffEle[mult_bin][13],p3DEffPos[mult_bin][13],pSectorCorr,mult_bin,NO_CUT,weight_yield_pair,4);
                    }
                    if(oAngle > oAngleCut) {
                        fillPairSpectraHC(hMPair,cand1,cand2,"_oa9deg",p3DEffEle[0][13],p3DEffPos[0][13],pSectorCorr,ALL_MULTS,NO_CUT,weight_yield_pair,4);
//                        fillPairSpectraHC(hMPair,cand1,cand2,"_oa9deg",p3DEffEle[mult_bin][13],p3DEffPos[mult_bin][13],pSectorCorr,ALL_MULTS,NO_CUT,weight_yield_pair,4);
                        if (mult_bin > 0 && mult_bin < 5) {
                            fillPairSpectraHC(hMPair,cand1,cand2,"_oa9deg",p3DEffEle[mult_bin][13],p3DEffPos[mult_bin][13],pSectorCorr,mult_bin,NO_CUT,weight_yield_pair,4);
                        }
                    }
                    if(!removedRecursive[j] && !removedRecursive[k] && oAngle > oAngleCut) {
                        fillPairSpectraHC(hMPair,cand1,cand2,"_oa9degRecur",p3DEffEle[0][13],p3DEffPos[0][13],pSectorCorr,ALL_MULTS,NO_CUT,weight_yield_pair,4);
//                        fillPairSpectraHC(hMPair,cand1,cand2,"_oa9degRecur",p3DEffEle[mult_bin][13],p3DEffPos[mult_bin][13],pSectorCorr,ALL_MULTS,NO_CUT,weight_yield_pair,4);
                        if (mult_bin > 0 && mult_bin < 5) {
                            fillPairSpectraHC(hMPair,cand1,cand2,"_oa9degRecur",p3DEffEle[mult_bin][13],p3DEffPos[mult_bin][13],pSectorCorr,mult_bin,NO_CUT,weight_yield_pair,4);
                        }
                    }
                    /*
                    if (oAngle > 9 && invM < 150) {
                        if (usedIndices.find(cand1->getIndex()) == usedIndices.end()) {
                            fillSingle(hMSingle,cand1,"_fromPi0",mult_bin,nearestIndexSeg[j]);
                            usedIndices.insert(cand1->getIndex());
                        }
                        if (usedIndices.find(cand2->getIndex()) == usedIndices.end()) {
                            fillSingle(hMSingle,cand2,"_fromPi0",mult_bin,nearestIndexSeg[j]);
                            usedIndices.insert(cand2->getIndex());
                        }
                    }
*/

                } // leptonFlagHC[j] && leptonFlagHC[k]
                int w = 0;
                for (int w_sys0 = 0; w_sys0 < NWEIGHTS; ++w_sys0) {
                    for (int w_sys1 = 0; w_sys1 < NWEIGHTS; ++w_sys1) {
                        if (useForPairs[w_sys0][w_sys1]) {
                            bool leptonFlag_j, leptonFlag_k;
                            if (sysj == 0) {
                                leptonFlag_j = leptonFlagMLP[w_sys0][j];
                            }
                            else {
                                leptonFlag_j = leptonFlagMLP[w_sys1][j];
                            }
                            if (sysk == 0) {
                                leptonFlag_k = leptonFlagMLP[w_sys0][k];
                            }
                            else {
                                leptonFlag_k = leptonFlagMLP[w_sys1][k];
                            }
                            TString strW = TString("_w") + TString::Itoa(w_sys0,10) + "w" + TString::Itoa(w_sys1,10);
                            Int_t pfind = -1;
                            if (w_sys0 == 0  && w_sys1 == 0) {
                                pfind = 0;
                            }
                            if (w_sys0 == 6  && w_sys1 == 6) {
                                pfind = 1;
                            }
                            if (w_sys0 == 10 && w_sys1 == 10) {
                                pfind = 2;
                            }
                            if (w_sys0 == 12 && w_sys1 == 12) {
                                pfind = 3;
                            }
/*
                            pSectorCorr = NULL;
                            pfind = -1;
                            weight_yield_pair = 1;
*/
                            if (leptonFlag_j && leptonFlag_k && !removedIncomplete[j] && !removedIncomplete[k] && !removedNearest[j] && !removedNearest[k]/* && !removedClosePartner[j] && !removedClosePartner[k]*/) {
                                Float_t oAngle   = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), cand2->getPhi(), cand2->getTheta());

//                                fillPairSpectra(hMPair,cand1,cand2,strW,p3DEffEle[mult_bin][w],p3DEffPos[mult_bin][w],pSectorCorr,ALL_MULTS,NO_CUT,weight_yield_pair,pfind);

                                fillPairSpectra(hMPair,cand1,cand2,strW,p3DEffEle[0][w],p3DEffPos[0][w],pSectorCorr,ALL_MULTS,NO_CUT,weight_yield_pair,pfind);
                                if (mult_bin > 0 && mult_bin < 5) {
                                    fillPairSpectra(hMPair,cand1,cand2,strW,p3DEffEle[mult_bin][w],p3DEffPos[mult_bin][w],pSectorCorr,mult_bin,NO_CUT,weight_yield_pair,pfind);
                                }
                                if(oAngle > oAngleCut) {
//                                    fillPairSpectra(hMPair,cand1,cand2,strW+"_oa9deg",p3DEffEle[mult_bin][w],p3DEffPos[mult_bin][w],pSectorCorr,ALL_MULTS,NO_CUT,weight_yield_pair,pfind);
                                    fillPairSpectra(hMPair,cand1,cand2,strW+"_oa9deg",p3DEffEle[0][w],p3DEffPos[0][w],pSectorCorr,ALL_MULTS,NO_CUT,weight_yield_pair,pfind);
                                    if (mult_bin > 0 && mult_bin < 5) {
                                        fillPairSpectra(hMPair,cand1,cand2,strW+"_oa9deg",p3DEffEle[mult_bin][w],p3DEffPos[mult_bin][w],pSectorCorr,mult_bin,NO_CUT,weight_yield_pair,pfind);
                                    }
                                }
                                if(!removedRecursive[j] && !removedRecursive[k] && oAngle > oAngleCut) {
                                    //fillPairSpectra(hMPair,cand1,cand2,strW+"_oa9degRecur",p3DEffEle[mult_bin][w],p3DEffPos[mult_bin][w],pSectorCorr,ALL_MULTS,NO_CUT,weight_yield_pair,pfind);
                                    fillPairSpectra(hMPair,cand1,cand2,strW+"_oa9degRecur",p3DEffEle[0][w],p3DEffPos[0][w],pSectorCorr,ALL_MULTS,NO_CUT,weight_yield_pair,pfind);
                                    if (mult_bin > 0 && mult_bin < 5) {
                                        fillPairSpectra(hMPair,cand1,cand2,strW+"_oa9degRecur",p3DEffEle[mult_bin][w],p3DEffPos[mult_bin][w],pSectorCorr,mult_bin,NO_CUT,weight_yield_pair,pfind);
                                    }
                                    Float_t eff1 = getEfficiencyFactor(p3DEffEle[mult_bin][w], p3DEffPos[mult_bin][w], cand1->getMomentum(), cand1->getTheta(), cand1->getPhi(), cand1->getCharge());
                                    Float_t eff2 = getEfficiencyFactor(p3DEffEle[mult_bin][w], p3DEffPos[mult_bin][w], cand2->getMomentum(), cand2->getTheta(), cand2->getPhi(), cand2->getCharge());
                                    if (eff1 > 0.02 && eff2 > 0.02 && w_sys0 == 6 && w_sys1 == 6) {
                                        Float_t effPair = eff1*eff2;
                                        TLorentzVector dilep = *cand1 + *cand2;
                                        Float_t invM = dilep.M()/1000;
                                        if (cand1->getCharge() != cand2->getCharge()) {
                                            pmassNP_eff->Fill(invM, effPair);
                                        } else if (cand1->getCharge() == 1) {
                                            pmassPP_eff->Fill(invM, effPair);
                                        } else {
                                            pmassNN_eff->Fill(invM, effPair);
                                        }
                                    }
                                }
                            }
                            ++w;
                        }
                    }
                }
                //----------------------------------------------------------------------
            } // inner loop
        } // END PAIRS FILLING HISTOGRAMS

    } // end event loop

    sorter.finalize();
    timer.Stop();

    hMSingle.getFile()->cd();
    TMacro m1(__DIELEANA_FILE__);
    TMacro m2(__HELPERFC_FILE__);
    m1.Write();
    m2.Write();
    hMPair.getFile()->cd();
    pmassNP_eff->Write();
    pmassPP_eff->Write();
    pmassNN_eff->Write();

    hMSingle.writeHists("nomap");
    hMPair.writeHists("nomap");
    if (isSimulation) hMKine.writeHists("nomap");
    record.writeRecord(outfileEvents.Data());
    record_hc.writeRecord(outfileEventsHC.Data());

    
    cout<<"####################################################"<<endl;
    return 0;
}

