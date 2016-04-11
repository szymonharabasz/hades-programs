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
#include "/u/harabasz/myhydralibs-newmixer-generic/eventmixer/heventmixer.h"
#include "hparticlebooker.h"
#include "hparticleanglecor.h"
#include "hparticledef.h"
#include "hparticlestructs.h"
#include "hparticlepair.h"
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
static Int_t WEIGHT = 13;
static const Float_t oAngleCut      = 9.;
static map<Int_t, Bool_t> sorterFlag;
static map<Int_t, Bool_t> richQaFlag;
static map<Int_t, Bool_t> betaFlag;
static map<Int_t, Bool_t> betaFlagTof;
static map<Int_t, Bool_t> betaFlagRpc;
static map<Int_t, Bool_t> showerFlag;
static map<Int_t, Bool_t> leptonFlag;
static map<Int_t, Bool_t> leptonFlagHC;
static map<Int_t, Bool_t> leptonFlagMLP[NWEIGHTS];
static map<Int_t, Bool_t> removedNearest;
static map<Int_t, Bool_t> removedIncomplete;
//------------------------------------------------------------------
#include __HELPERFC_FILE__
#include "yield_pars.h"
#include "fline.h"
#include "eventclassifier.h"

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

//------------------------------------------------------------


Int_t dieleAna(TString inputlist, TString outfile, Int_t nev=-1, Int_t whichweight = 0)
{
    WEIGHT = whichweight;

    TH1::SetDefaultSumw2();
    setupRichCuts("iso_newer.root");
    setupMLPCuts("mlpmom_cutg_gen8_new.root");
    setupMassCuts("m2_cut.root");
    setupSectorCorrs(kPM,"sector_factor_pm.root");
    setupSectorCorrs(kPP,"sector_factor_pp.root");
    setupSectorCorrs(kMM,"sector_factor_mm.root");
    setupPiotrFactor("piotrfactor_gen8_consistent.root");
    setupPiotrFactor("piotrfactor_gen8.root");

    TH3F *p3DEffEle[6][NWEIGHTS+1]; // mult bins, MLP weights + HC
    TH3F *p3DEffPos[6][NWEIGHTS+1];
    TH3F *p3DAccEle[6][NWEIGHTS+1];
    TH3F *p3DAccPos[6][NWEIGHTS+1];
    readAccEffMatrices(p3DAccEle, p3DAccPos, p3DEffEle, p3DEffPos);

    TH1F *pEventClass;
    TH1F *pEventClass_recur;
    TFile *pEventClassFile;
//    pEventClassFile = new TFile("eventClass_mult_nonempty_4secmult_200fpj_wait.root");
//    pEventClassFile = new TFile("eventClass_target_mult_rplane_nonempty_4secmult.root");
    pEventClassFile = new TFile(TString("eventClass_target_mult_rplane_minmom_nmix_w")+TString::Itoa(WEIGHT,10)+".root");
    if (pEventClassFile) {
        pEventClass = (TH1F*)pEventClassFile->Get("eventClass");
        pEventClass_recur = (TH1F*)pEventClassFile->Get("eventClass_recur");
        if (pEventClass == NULL || pEventClass_recur == NULL) {
            Error("DrawFromNtuple constructor","Histogram not found in the event class file");
            exit (-1);
        }
    }    
    else {
        Error("DrawFromNtuple constructor","Event class file not found");
        exit (-1);
    }
    HLoop* loop = new HLoop(kTRUE);  // kTRUE : create Hades  (needed to work with standard eventstructure)
//    TString readCategories = "-*,+HParticleCand,+HParticleEvtInfo,+HStart2Hit,+HStart2Cal,+HGeantKine,+HGeantMdc,+HGeantTof,+HGeantRpc,+HGeantShower";
    TString readCategories = "-*,+HParticleCand,+HParticleEvtInfo,+HStart2Hit,+HStart2Cal";
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

    hM.addHist("TH1F","hEvtHour", "hEvtHour",744*60,96,127);
    hM.addHist("TH1F","hCounter", "hCounter",20,0,20);
    hM.addHist("TH1D","hClassFirstFilled","hClassFirstFilled",2000,0,2000);
    hM.addHist("TH1F","hClass",           "hClass",           2000,0,2000);

    //------------------------------------------------------------------------------------------
    //--------------- Lepton pairs spectra -----------------------------------------------------
    //------------------------------------------------------------------------------------------
/*
    addHists(hM,"hcrossXZ",100,-100,50,100,-50,50);
    addHists(hM,"hcrossYZ",100,-100,50,100,-50,50);
    addHists(hM,"hclosApprXZ",100,-100,50,100,-50,50);
    addHists(hM,"hclosApprYZ",100,-100,50,100,-50,50);
    addHists(hM,"hvertXZ",100,-100,50,100,-50,50);
    addHists(hM,"hvertYZ",100,-100,50,100,-50,50);
    addHists(hM,"hminDistVertZ",100,-50,0,100,-50,50);
    addHists(hM,"hph1ph2",90,0,360,90,0,360);
    addHists(hM,"hDthetaDphi",180,-360,360,180,-90,90);
    addHists(hM,"hDthetaDphiSinTheta",180,-360,360,180,-90,90);
    addHists(hM,"hphiDilPlane",90,-180,180,90,-180,180);
*/
//    addHists(hM,"hthetaAvgDil", 90,0,90,90,0,90);
//    addHists(hM,"hphiAvgDil", 360,-360,360,360,-360,360);
    addHists(hM,"hsqrtp1p2oa", binsxmass,0,maxmass,90,0,180);
    addHists(hM,"htheta1y", 45,0,90,100,0,2);
//    addHists(hM,"htheta2y", 45,0,90,100,0,2);
    addHists(hM,"hztheta1", 100,-80,20,45,0,90);
//    addHists(hM,"hztheta2", 100,-80,20,45,0,90);

    addHists(hM,"hy" ,100,0,2);
    addHists(hM,"hpt",100,0,1000);
//    addHists(hM,"hp1p2diff",2000,-1000,1000);
//    addHists(hM,"hth1oAngle" ,45,0,90,90,0,180);
//    addHists(hM,"hth1oAngle" ,45,0,90,90,0,180);
//    addHists(hM,"hth2oAngle" ,45,0,90,90,0,180);
//    addHists(hM,"hth1mass" ,45,0,90,nbins,xAxis1);
//    addHists(hM,"hth2mass" ,45,0,90,nbins,xAxis1);
    addHists(hM,"hoAnglemass" ,90,0,180,nbins,xAxis1);
    addHists(hM,"hoAnglemassNoFactors" ,90,0,180,nbins,xAxis1);
//    addHists(hM,"hth1pt" ,45,0,90,120,0,1200);
//    addHists(hM,"hth2pt" ,45,0,90,120,0,1200);
    addHists(hM,"hoAnglept" ,90,0,180,120,0,1200);
    addHists(hM,"hoAngleptNoFactors" ,90,0,180,120,0,1200);
    addHists(hM,"hmasspt" ,nbins,xAxis1,120,0,1200);
    addHists(hM,"hmassptNoFactors" ,nbins,xAxis1,120,0,1200);
//    addHists(hM,"hth1y" ,45,0,90,100,0,2);
//    addHists(hM,"hth2y" ,45,0,90,100,0,2);
    addHists(hM,"hoAngley" ,90,0,180,100,0,2);
    addHists(hM,"hmassy" ,nbins,xAxis1,100,0,2);
    addHists(hM,"hpty" ,120,0,1200,100,0,2);
//    addHists(hM,"hp1th1", 100,0,1000,45,0,90);
//    addHists(hM,"hp2th2", 100,0,1000,45,0,90);
    addHists(hM,"hth1th2" ,45,0,90,45,0,90);
//    addHists(hM,"hpt1pt2" ,100,0,1000,100,0,1000);
    addHists(hM,"hp1p2LowM" ,100,0,1000,100,0,1000);
    addHists(hM,"hp1p2MidM" ,100,0,1000,100,0,1000);
    addHists(hM,"hp1p2HigM" ,100,0,1000,100,0,1000);
//    addHists(hM,"hy1y2" ,100,0,2,100,0,2);
    addHists(hM,"hzy" ,100,-80,20,100,0,2);
    addHists(hM,"hmass" ,nbins,xAxis1);
    addHists(hM,"hmassNoFactors" ,nbins,xAxis1);
    addHists(hM,"hmassPiotrFactor" ,nbins,xAxis1);
    addHists(hM,"hmassSectorFactor" ,nbins,xAxis1);
    addHists(hM,"hmtmass" ,20,0,1000,nbins,xAxis1);
    addHists(hM,"hoAngle" ,2000,0,200);
    addHists(hM,"hoAngleNoFactors" ,200,0,200);
    addHists(hM,"hoAnglePiotrFactor" ,200,0,200);
    addHists(hM,"hoAngleSectorFactor" ,200,0,200);
    addHists(hM,"hmtLowM" ,20,0,1000);
    addHists(hM,"hmtMidM" ,20,0,1000);
    addHists(hM,"hmtHigM" ,20,0,1000);
    addHists(hM,"hmtMinusMinvLowM" ,20,0,1000);
    addHists(hM,"hmtMinusMinvMidM" ,20,0,1000);
    addHists(hM,"hmtMinusMinvHigM" ,20,0,1000);
    addHists(hM,"hmass4secd1" ,nbins,xAxis1);
    addHists(hM,"hoAngle4secd1" ,2000,0,200);
    addHists(hM,"hmass4secd2" ,nbins,xAxis1);
    addHists(hM,"hoAngle4secd2" ,2000,0,200);
    
    hM.addHist(new TH2F("hbeta1Sqrtp1p2NP","",100,0,1000,400,0,2));
    hM.addHist(new TH2F("hbeta2Sqrtp1p2NP","",100,0,1000,400,0,2));
    hM.addHist(new TH2F("hdBetaSqrtp1p2NP","",100,0,1000,400,-1,1));
    hM.addHist(new TH2F("hbeta1Sqrtp1p2PP","",100,0,1000,400,0,2));
    hM.addHist(new TH2F("hbeta2Sqrtp1p2PP","",100,0,1000,400,0,2));
    hM.addHist(new TH2F("hdBetaSqrtp1p2PP","",100,0,1000,400,-1,1));
    hM.addHist(new TH2F("hbeta1Sqrtp1p2NN","",100,0,1000,400,0,2));
    hM.addHist(new TH2F("hbeta2Sqrtp1p2NN","",100,0,1000,400,0,2));
    hM.addHist(new TH2F("hdBetaSqrtp1p2NN","",100,0,1000,400,-1,1));

    hM.addHist(new TH2F("hbeta1MassNP","",100,0,1000,400,0,2));
    hM.addHist(new TH2F("hbeta2MassNP","",100,0,1000,400,0,2));
    hM.addHist(new TH2F("hdBetaMassNP","",100,0,1000,400,-1,1));
    hM.addHist(new TH2F("hbeta1MassPP","",100,0,1000,400,0,2));
    hM.addHist(new TH2F("hbeta2MassPP","",100,0,1000,400,0,2));
    hM.addHist(new TH2F("hdBetaMassPP","",100,0,1000,400,-1,1));
    hM.addHist(new TH2F("hbeta1MassNN","",100,0,1000,400,0,2));
    hM.addHist(new TH2F("hbeta2MassNN","",100,0,1000,400,0,2));
    hM.addHist(new TH2F("hdBetaMassNN","",100,0,1000,400,-1,1));

    addHistsHC(hM,"hmass" ,nbins,xAxis1);
    addHistsHC(hM,"hmtmass" ,20,0,1000,nbins,xAxis1);
    addHistsHC(hM,"hmasspt" ,nbins,xAxis1,120,0,1200);
    addHistsHC(hM,"hoAngle" ,2000,0,200);
    addHistsHC(hM,"hmtLowM" ,20,0,1000);
    addHistsHC(hM,"hmtMidM" ,20,0,1000);
    addHistsHC(hM,"hmtHigM" ,20,0,1000);
    addHistsHC(hM,"hmtMinusMinvLowM" ,20,0,1000);
    addHistsHC(hM,"hmtMinusMinvMidM" ,20,0,1000);
    addHistsHC(hM,"hmtMinusMinvHigM" ,20,0,1000);

    addHists(hM,"hclassFirstFilled",2000,0,2000);

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

    HParticleEventMixer eventmixer;
    eventmixer.setPIDs(2,3,1);
    eventmixer.setBuffSize(40);
    eventmixer.setEventClassifier(eventClassifierMultTarget);

    HParticleEventMixer eventmixer_recur;
    eventmixer_recur.setPIDs(2,3,1);
    eventmixer_recur.setBuffSize(40);
    eventmixer_recur.setEventClassifier(eventClassifierMultTarget);

    TStopwatch timer;
    timer.Reset();
    timer.Start();

    Int_t evtsInFile = loop->getEntries();
    if(nev < 0 || nev > evtsInFile ) nev = evtsInFile;

    string currentFileName;
    bool isSimulation;
    bool isEmbedding;
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

        hM.get("hCounter")->GetXaxis()->SetBinLabel(1,"all");
        hM.get("hCounter")->GetXaxis()->SetBinLabel(2,"kGoodTRIGGER");
        hM.get("hCounter")->GetXaxis()->SetBinLabel(3,"kGoodSTART");
        hM.get("hCounter")->GetXaxis()->SetBinLabel(4,"kGoodVertexCand");
        hM.get("hCounter")->GetXaxis()->SetBinLabel(5,"kNoPileUpSTART");
        hM.get("hCounter")->GetXaxis()->SetBinLabel(6,"kNoVETO");
        hM.get("hCounter")->GetXaxis()->SetBinLabel(7,"kGoodSTARTVETO");
        hM.get("hCounter")->GetXaxis()->SetBinLabel(8,"kGoodSTARTMETA");
        hM.get("hCounter")->GetXaxis()->SetBinLabel(9,"#geq 4 secs");
        hM.get("hCounter")->GetXaxis()->SetBinLabel(10,"#geq 5 secs");
        hM.get("hCounter")->GetXaxis()->SetBinLabel(11,"6 secs");
        for (int mb = 0; mb <= 4; ++mb) {
            hM.get("hCounter")->GetXaxis()->SetBinLabel(12+mb,Form("mult bin %i",mb));
        }
        hM.get("hCounter")->GetXaxis()->SetBinLabel(17,"mult underflow");

        hM.get("hCounter")->Fill(0);
        if (!isSimulation) {
            if (!evtinfo->isGoodEvent(kGoodTRIGGER)) continue;
            hM.get("hCounter")->Fill(1);

            if (!evtinfo->isGoodEvent(kGoodSTART)) continue;
            hM.get("hCounter")->Fill(2);

            if (!evtinfo->isGoodEvent(kGoodVertexCand)) continue;
            hM.get("hCounter")->Fill(3);
            if (!evtinfo->isGoodEvent(kNoPileUpSTART)) continue;
            hM.get("hCounter")->Fill(4);
            if (!evtinfo->isGoodEvent(kNoVETO)) continue;
            hM.get("hCounter")->Fill(5);
            if (!evtinfo->isGoodEvent(kGoodSTARTVETO)) continue;
            hM.get("hCounter")->Fill(6);
            if (!evtinfo->isGoodEvent(kGoodSTARTMETA)) continue;
            hM.get("hCounter")->Fill(7);
            if (sectorCorr == kSkip) continue;
            hM.get("hCounter")->Fill(8);
//            if (sectorCorr != kCorr5to6) continue;
//            if (sectorCorr == kNoCorr || sectorCorr == kCorr5to6) continue;
//            hM.get("hCounter")->Fill(9);
//            if (sectorCorr != kNoCorr) continue;
//            hM.get("hCounter")->Fill(10);
        }
	// Select meta multiplicity / event
        Int_t mult_meta = evtinfo->getSumTofMultCut() + evtinfo->getSumRpcMultHitCut();
        Int_t mult_bin = 5;
        if (mult_meta >  60 && mult_meta <=  88) mult_bin = 4; // most peripheral
        if (mult_meta >  88 && mult_meta <= 121) mult_bin = 3;
        if (mult_meta > 121 && mult_meta <= 160) mult_bin = 2;
        if (mult_meta > 160 && mult_meta <= 250) mult_bin = 1; // most central
        if (mult_meta > 250) mult_bin = 0;
        //if (mult_bin == 5) continue;

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
        hM.get("hCounter")->Fill(11+mult_bin);

        Int_t size = candCat->getEntries();

        sorter.cleanUp();
        sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);            // reset all flags for flags (0-28) ,reject,used,lepton
        sorter.fill(selectLeptonsBeta);   // fill only good leptons
        Int_t n_lepparcnd = sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsLepton);
	//------------------------------------------------------------------------
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
	leptonFlag.clear();
	leptonFlagHC.clear();
	removedNearest.clear();
	removedIncomplete.clear();
	Bool_t removedRecursive[size];
	Bool_t removedRecursiveNearestSeg[5][size];
	Bool_t removedEmbedding[size];
        for (int w = 0; w < NWEIGHTS; ++w) {
            leptonFlagMLP[w].clear();
        }

        int n_ep_cand = 0;
        int n_em_cand = 0;
        int n_ep_hc = 0;
        int n_em_hc = 0;

        doTrackCorr(trackcor, candCat);
        // SETTNG FLAGS
	for(Int_t j = 0; j < size; j ++){
	    sorterFlag[j]                     = kFALSE;
            richQaFlag[j]                     = kFALSE;
            betaFlag[j]                       = kFALSE;
            betaFlagTof[j]                    = kFALSE;
            betaFlagRpc[j]                    = kFALSE;
            showerFlag[j]                     = kFALSE;
            leptonFlag[j]                     = kFALSE;
            leptonFlagHC[j]                   = kFALSE;
	    removedNearest[j]                 = kFALSE;
	    removedIncomplete[j]              = kFALSE;
            removedRecursive[j]               = kFALSE;
            removedEmbedding[j]               = kFALSE;
            for (int w = 0; w < NWEIGHTS; ++w) {
                leptonFlagMLP[w][j]           = kFALSE;
            }
            for (int cp = 2; cp <= 6; ++cp) {
                removedRecursiveNearestSeg[cp-2][j] = kFALSE;
            }
            // Track candidates (Leptons AND Hadrons)
	    cand1 = HCategoryManager::getObject(cand1,candCat,j);

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

            if (!cand1->isFlagBit(Particle::kIsLepton)) continue;
            sorterFlag[j] = kTRUE;

            for (int w = 0; w < NWEIGHTS; ++w) {
                MLPs[w][j] = reader[w][sys]->EvaluateMVA( "MLP" );
            }
            ///////////////////////////////////////////////////// //////////////////////////////
            //                                                    
            // Hard cut flags                                     
            //
            ///////////////////////////////////////////////////////////////////////////////////
            // Ring matching
            //if (richQa < 3.5e-3*_mom+1) richQaFlag[j]=kTRUE;
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
            if (WEIGHT == 13) {
                leptonFlag[j] = leptonFlagHC[j];
            }
            else {
                leptonFlag[j] = leptonFlagMLP[WEIGHT][j];
            }
	    //*************************************************************
	}  // END SETTING FLAGS
        // LOOP TO COUNT LEPTONS
	for(Int_t j = 0; j < size; j ++){
            if (!sorterFlag[j]) continue;
	    cand1 = HCategoryManager::getObject(cand1,candCat,j);

            if (cand1->getCharge() == 1) n_ep_cand++;
            if (cand1->getCharge() ==-1) n_em_cand++;
            if (leptonFlagHC[j]) {
                if (cand1->getCharge() == 1) n_ep_hc++;
                if (cand1->getCharge() ==-1) n_em_hc++;
            }
	    //*************************************************************
	}  // END LOOP TO COUNT LEPTONS
        // RECURSIVE CUT
        for(Int_t j = 0; j < size; j ++){ // first candidate
            if (!sorterFlag[j]) continue;
            cand1 = HCategoryManager::getObject(cand1,candCat,j);

            for(Int_t k = j+1; k < size; k ++){ // second candidate, k > j condition applied below
                if (!sorterFlag[k]) continue;
                cand2 = HCategoryManager::getObject(cand2,candCat,k);

//                UInt_t flags;
//                HParticleTool::setPairFlags(flags,cand2,cand1);
//                if(!HParticleTool::evalPairsFlags(kPairCase1,flags)) continue;

                if (k > j && leptonFlag[j] && leptonFlag[k]) {

                    TLorentzVector dilep = (*cand1) + (*cand2);

                    if((cand1->getCharge()==1  && cand2->getCharge()==-1) ||
                            (cand1->getCharge()==-1 && cand2->getCharge()==1)) {
                        Float_t oAngle   = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), cand2->getPhi(), cand2->getTheta());
                        if (/*!removedRecursive[j] && !removedRecursive[k] && */oAngle < oAngleCut) {
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
        } // end recursive removal of close pairs

        // Additional loop to fill vector
        vector<HParticleCand *> vep;
        vector<HParticleCand *> vem;
        vector<HParticleCand *> vep_recur;
        vector<HParticleCand *> vem_recur;
	for(Int_t j = 0; j < size; j ++){
            if (isEmbedding && removedEmbedding[j]) continue;
            if (!sorterFlag[j]) continue;
	    cand1 = HCategoryManager::getObject(cand1,candCat,j);

            if (leptonFlag[j]) {
                if (cand1->getCharge() == -1) {
                    vem.push_back(new HParticleCand(*cand1));
                }
                if (cand1->getCharge() == 1) {
                    vep.push_back(new HParticleCand(*cand1));
                }
//                    }
//                }
                if(!removedRecursiveNearestSeg[4][j]) {
                    if (cand1->getCharge() == -1) {
                        vem_recur.push_back(new HParticleCand(*cand1));
                    }
                    if (cand1->getCharge() == 1) {
                        vep_recur.push_back(new HParticleCand(*cand1));
                    }
                }
            }
            //*************************************************************
        }
        bool dontmix = false;
        bool dontmix_recur = false;

        //dontmix       = (vep.size()       < 1 || vem.size()       < 1);
        //dontmix_recur = (vep_recur.size() < 1 || vem_recur.size() < 1);
        //dontmix       = ((vep.size()       + vem.size()      ) < 2);
        //dontmix_recur = ((vep_recur.size() + vem_recur.size()) < 2);
        if (dontmix) {
            vep.clear();
            vem.clear();
        }
        if (dontmix_recur) {
            vep_recur.clear();
            vem_recur.clear();
        }

        eventmixer.nextEvent();
        eventmixer.addVector(vep,2);
        eventmixer.addVector(vem,3);
        vector<HParticlePair>& pairsVec = eventmixer.getMixedVector();

        eventmixer_recur.nextEvent();
        eventmixer_recur.addVector(vep_recur,2);
        eventmixer_recur.addVector(vem_recur,3);
        vector<HParticlePair>& pairsVecRecur =  eventmixer_recur.getMixedVector();

        Int_t eventClass = eventmixer.currentEventClass();
        Double_t eventClassMult = 0;
        eventClassMult = pEventClass->GetBinContent(pEventClass->FindBin(eventClass));
        Double_t eventClassMult_recur = 0;
        eventClassMult_recur = pEventClass_recur->GetBinContent(pEventClass->FindBin(eventClass));


//        for (int rphi = 1; rphi <= 13; ++rphi) {            
//            eventClassMult += pEventClass->GetBinContent(pEventClass->FindBin(16*rphi+eventClass));
//        }        
//        for (int target = 0; target <= 14; ++target) {            
//            eventClassMult += pEventClass->GetBinContent(pEventClass->FindBin(target+eventClass));
//        }        
//        for (int target = 0; target <= 14; ++target) {            
//            for (int rphi = 1; rphi <= 13; ++rphi) {            
//                eventClassMult += pEventClass->GetBinContent(pEventClass->FindBin(16*rphi+target+eventClass));
//            }        
//        }        
        Double_t eventClassWeight = 1;//pEventClass->Integral()/eventClassMult;
        Double_t eventClassWeight_recur = 1;//pEventClass_recur->Integral()/eventClassMult_recur;
//        cout << "eventClassWeight : " << eventClassWeight << endl;
        size = pairsVec.size();
//        if (vep.size() || vem.size() ) cout << vep.size() << " " << vem.size() << " " << eventClass << " " << eventmixer.dataStructure[eventClass].size() << " " << size << endl;
        TString strW = TString("_w") + TString::Itoa(WEIGHT,10) + "w" + TString::Itoa(WEIGHT,10);
        Int_t pfind = -1;
        if (WEIGHT == 0) {
            pfind = 0;
        }
        if (WEIGHT == 6) {
            pfind = 1;
        }
        if (WEIGHT == 10) {
            pfind = 2;
        }
        if (WEIGHT == 12) {
            pfind = 3;
        }
        for (Int_t j = 0; j < size; j ++) {
            HParticlePair& pair = pairsVec[j];

            cand1 = pair.getCand(0);
            cand2 = pair.getCand(1);
            //                if (cand1->getSector() != cand2->getSector()) continue;
            if (cand1->getCharge() == -1 && cand2->getCharge() == 1) {
                HParticleCand *temp = cand1;
                cand1 = cand2;
                cand2 = temp;
            }


            TLorentzVector dilep = (*cand1) + (*cand2);
            Float_t oAngle   = cand1->Angle(cand2->Vect())*TMath::RadToDeg();
//            if (oAngle < oAngleCut) continue;
//            if (oAngle > 130) continue;
            Int_t chg = kPM;
            if (cand1->getCharge() == cand2->getCharge()) {
                if (cand1->getCharge() > 0) chg = kPP;
                if (cand1->getCharge() < 0) chg = kMM;
            }
            TH2F *pSectorCorr = hSectorCorr[chg][sectorCorr-1];

            if (WEIGHT == 13) {
                fillPairSpectraHC(hM,cand1,cand2,"",p3DEffEle[0][13],p3DEffPos[0][13],pSectorCorr,ALL_MULTS,NO_CUT,eventClassWeight);
                if (mult_bin > 0 && mult_bin < 5) {
                    fillPairSpectraHC(hM,cand1,cand2,"",p3DEffEle[mult_bin][13],p3DEffPos[mult_bin][13],pSectorCorr,mult_bin,NO_CUT,eventClassWeight);
                }
                if(oAngle > oAngleCut) {
                    fillPairSpectraHC(hM,cand1,cand2,"_oa9deg",p3DEffEle[0][13],p3DEffPos[0][13],pSectorCorr,ALL_MULTS,NO_CUT,eventClassWeight);
                    if (mult_bin > 0 && mult_bin < 5) {
                        fillPairSpectraHC(hM,cand1,cand2,"_oa9deg",p3DEffEle[mult_bin][13],p3DEffPos[mult_bin][13],pSectorCorr,mult_bin,NO_CUT,eventClassWeight);
                    }
                }
            }
            else {
                fillPairSpectra(hM,cand1,cand2,strW,p3DEffEle[0][WEIGHT],p3DEffPos[0][WEIGHT],pSectorCorr,ALL_MULTS,0,eventClassWeight,pfind);
                if (mult_bin > 0 && mult_bin < 5) {
                    fillPairSpectra(hM,cand1,cand2,strW,p3DEffEle[mult_bin][WEIGHT],p3DEffPos[mult_bin][WEIGHT],pSectorCorr,mult_bin,0,eventClassWeight,pfind);
                }

                if (oAngle > oAngleCut) {
                    fillPairSpectra(hM,cand1,cand2,strW+"_oa9deg",p3DEffEle[0][WEIGHT],p3DEffPos[0][WEIGHT],pSectorCorr,ALL_MULTS,0,eventClassWeight,pfind);
                    if (mult_bin > 0 && mult_bin < 5) {
                        fillPairSpectra(hM,cand1,cand2,strW+"_oa9deg",p3DEffEle[mult_bin][WEIGHT],p3DEffPos[mult_bin][WEIGHT],pSectorCorr,mult_bin,0,eventClassWeight,pfind);
                    }
                }
            }

            Int_t bin = hM.get("hClassFirstFilled")->FindBin(eventClass);
            Double_t content = hM.get("hClassFirstFilled")->GetBinContent(bin);
            if (content == 0.) {
                hM.get("hClassFirstFilled")->Fill(eventClass,i);
            }
            hM.get("hClass")->Fill(eventClass);
        }
        size = pairsVecRecur.size();
        for (Int_t j = 0; j < size; j ++) {
            HParticlePair& pair = pairsVecRecur[j];

            cand1 = pair.getCand(0);
            cand2 = pair.getCand(1);
            //                if (cand1->getSector() != cand2->getSector()) continue;
            if (cand1->getCharge() == -1 && cand2->getCharge() == 1) {
                HParticleCand *temp = cand1;
                cand1 = cand2;
                cand2 = temp;
            }


            TLorentzVector dilep = (*cand1) + (*cand2);
            Float_t oAngle   = cand1->Angle(cand2->Vect())*TMath::RadToDeg();
            if (oAngle < oAngleCut) continue;
            Int_t chg = kPM;
            if (cand1->getCharge() == cand2->getCharge()) {
                if (cand1->getCharge() > 0) chg = kPP;
                if (cand1->getCharge() < 0) chg = kMM;
            }
            TH2F *pSectorCorr = hSectorCorr[chg][sectorCorr-1];

            if (WEIGHT == 13) {
                fillPairSpectraHC(hM,cand1,cand2,"_oa9degRecur",p3DEffEle[0][13],p3DEffPos[0][13],pSectorCorr,ALL_MULTS,NO_CUT,eventClassWeight_recur);
                if (mult_bin > 0 && mult_bin < 5) {
                    fillPairSpectraHC(hM,cand1,cand2,"_oa9degRecur",p3DEffEle[mult_bin][13],p3DEffPos[mult_bin][13],pSectorCorr,mult_bin,NO_CUT,eventClassWeight_recur);
                }
            }
            else {
                fillPairSpectra(hM,cand1,cand2,strW+"_oa9degRecur",p3DEffEle[0][WEIGHT],p3DEffPos[0][WEIGHT],pSectorCorr,ALL_MULTS,0,eventClassWeight_recur,pfind);
                if (mult_bin > 0 && mult_bin < 5) {
                    fillPairSpectra(hM,cand1,cand2,strW+"_oa9degRecur",p3DEffEle[mult_bin][WEIGHT],p3DEffPos[mult_bin][WEIGHT],pSectorCorr,mult_bin,0,eventClassWeight_recur,pfind);
                }

            }
        }
#define DELETE_MIX
#ifdef DELETE_MIX
        vector <HParticleCand *>* toDel = eventmixer.getObjectsToDelete();
        for (unsigned int ii = 0; ii < toDel->size(); ++ii) {
            delete toDel->at(ii);
        }
        toDel->clear();
        delete toDel;
        vector <HParticleCand *>* toDel_recur = eventmixer_recur.getObjectsToDelete();
        for (unsigned int ii = 0; ii < toDel_recur->size(); ++ii) {
            delete toDel_recur->at(ii);
        }
        toDel_recur->clear();
        delete toDel_recur;
#endif

    } // end event loop

    sorter.finalize();
    timer.Stop();

    hM.getFile()->cd();
    TMacro m1(__DIELEANA_FILE__);
    TMacro m2(__HELPERFC_FILE__);
    m1.Write();
    m2.Write();
    hM.writeHists("nomap");

    cout<<"####################################################"<<endl;
    return 0;
}
