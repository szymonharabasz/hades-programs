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
#include <set>
#include <stdio.h>

using namespace std;
//------------------------------------------------------------------
static const Float_t oAngleCut      = 9.;
//------------------------------------------------------------------
#include "selectFunctions.h"
#include "yield_pars.h"
#include "fline.h"
#include __HELPERFC_FILE__

//------------------------------------------------------------
Bool_t isTriggerPT3(){
    return gHades->getCurrentEvent()->getHeader()->isTBit(13);
    return kTRUE;
}

//------------------------------------------------------------
Bool_t selectStart() {
    HCategory *fCatStartHit = HCategoryManager::getCategory(catStart2Hit,1,"catStart2Hit");
    if(!((HStart2Hit*)fCatStartHit->getObject(0)))  return kFALSE;    //object not available (3%)
    if(((HStart2Hit*)fCatStartHit->getObject(0))->getCorrFlag() < 2) return kFALSE;  //No start time found (5%)   // besser 2
//    if(((HStart2Hit*)fCatStartHit->getObject(0))->getCorrFlag() == -1) return kFALSE;  //No start time found (5%)   // besser 2
    return kTRUE; // one ore two hits in Start (92%)
}

Int_t getStartFlag() {
    HCategory *fCatStartHit = HCategoryManager::getCategory(catStart2Hit,1,"catStart2Hit");
    if (fCatStartHit != NULL) {
        if(!((HStart2Hit*)fCatStartHit->getObject(0)))  return -2;    //object not available (3%)
        return ((HStart2Hit*)fCatStartHit->getObject(0))->getCorrFlag();  //No start time found (5%)   // besser 2
    }
    else {
        return -2;
    }
}

Bool_t isGoodStart() {
    HCategory* catStart  = (HCategory*)HCategoryManager::getCategory(catStart2Hit);
    HCategory* catStartcal  = (HCategory*)HCategoryManager::getCategory(catStart2Cal);
    HStart2Hit* start = (HStart2Hit*) catStart->getObject(0);
    if(start) {
        if(start->getCorrFlag()<0 ||start->getCorrFlag()==1) return kFALSE;
        Float_t startHitTime = start->getTime();
        if(catStartcal) {
            HStart2Cal* sHc=0;
            for(Int_t i=0;i< catStartcal->getEntries();i++) {
                if((sHc = (HStart2Cal *) catStartcal->getObject(i)) != NULL){
                    // find the corresponding start hits to the used
                    if( sHc->getModule()!=0 && sHc->getModule()!=1 ) continue;
                    for(Int_t j=1;j<=sHc->getMultiplicity()&&j<=sHc->getMaxMultiplicity();j++) {
                        if( fabs(sHc->getTime(j)-startHitTime)>0.25 && fabs(sHc->getTime(j)-startHitTime)<15.) {
                            return kFALSE;
                        }
                    }
                }
            }
        }
        return kTRUE;
    }else return kFALSE;
}

Bool_t goodVertex() { 
    HVertex PrimVertexReco = ((HEventHeader*)(gHades->getCurrentEvent()->getHeader()))->getVertexReco(); 
    //Float_t VertexX=PrimVertexReco.getX(); 
    //Float_t VertexY=PrimVertexReco.getY(); 
 
    if(PrimVertexReco.getZ()>0 || PrimVertexReco.getZ()<-59.0){return kFALSE;}        // get only events with Vertex Reconstruction (59 mm) 
    //if(sqrt((VertexX*VertexX)+(VertexY*VertexY)) > 4.  ){return kFALSE;}                 // in XYPlane better than 4mm 
 
    // get only events where Vertex Reconstruction not failed (Chi2, Iterations, SumOfWeights == -1 && ==0) 
    if(PrimVertexReco.getChi2()<0.01) {return kFALSE;}  
    if(PrimVertexReco.getIterations()<0.01) {return kFALSE;}  
    if(PrimVertexReco.getSumOfWeights()<0.01) {return kFALSE;} 
 
    // Vertex Reconstruction Quality with Good Chi2 (inefficient in low multipicity)  
    //if(PrimVertexReco.getChi2()>100) {return kFALSE;} 
   return kTRUE; 
}

Bool_t isOutlier() { 
     HCategory* fParticleEvtInfoCat =  (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo,kTRUE,"catParticleEvtInfo"); 
     if(!fParticleEvtInfoCat) { cout<<"No catParticleEvtInfo in input!"<<endl; return kTRUE;} 
     HParticleEvtInfo *event_info = (HParticleEvtInfo*)fParticleEvtInfoCat->getObject(0); 
     if(!event_info) {std::cout<<"No Event INFo"<<std::endl;return kTRUE;} 
 
     // pile-up events with bad tof&rpc timing 
     if((((event_info->getSumTofMult()+event_info->getSumRpcMult()))-(event_info->getSumTofMultCut()+event_info->getSumRpcMultCut()))>30){return kTRUE;} 
 
     //pile-up events 
     //if(event_info->getSumTofMultCut()+event_info->getSumRpcMultCut()<0.5*event_info->getSumSelectedParticleCandMult() ){return kTRUE;} 
 
     //"monster event" 
     if(event_info->getSumMdcSegFittedMult(1)/(Float_t)event_info->getSumSelectedParticleCandMult()>2.5 ){return kTRUE;} 
 
     return kFALSE; 
}

Bool_t selectVeto()
{
    HCategory *fCatStartCal = HCategoryManager::getCategory(catStart2Cal,1,"catStart2Cal");
    HCategory *fCatStartHit = HCategoryManager::getCategory(catStart2Hit,1,"catStart2Hit");

    Int_t nStartHits = fCatStartHit->getEntries();
    if(nStartHits != 1) {
	return kFALSE;
    }

    Int_t nStartCals = fCatStartCal->getEntries();

    Int_t mult_start = 0;
    Int_t mult_veto  = 0;

    for(Int_t n = 0; n < nStartCals; n++) {
	HStart2Cal *start_cal = (HStart2Cal*)fCatStartCal->getObject(n);
	if(start_cal->getModule() == 0
	   && TMath::Abs(start_cal->getTime(1)) < 10.
	  ) {
	    mult_start++;
	}
    }
    if(mult_start > 1) return kFALSE;


    for(Int_t n = 0; n < nStartCals; n++) {
	HStart2Cal *start_cal = (HStart2Cal*)fCatStartCal->getObject(n);
	if(start_cal->getModule() == 3
	   && TMath::Abs(start_cal->getTime(1)) < 10.
	  ) {
	    mult_veto++;
	}
    }
    if(mult_veto != 0) return kFALSE;
    return kTRUE;
}

//------------------------------------------------------------
Bool_t isGoodEventPT3(){
    return  (isTriggerPT3()&&selectStart());
    return kTRUE;
}

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
Float_t getAngleFactor(TH1F *p1DAngleFactor, Float_t oAngle) {
    Float_t factor = p1DAngleFactor->Interpolate(oAngle);
    if (factor > 0) return 1./factor;
    return 1.;
}
//------------------------------------------------------------
//------------------------------------------------------------


Int_t dieleAna(TString inputlist, TString outfile, Int_t nev=-1, TString seqnumlist = "dummy")
{
    TH1::SetDefaultSumw2();

    setupMassCuts("m2_cut.root");
    setupMassCuts("m2_cut.root");
    setupRichCuts("iso_newer.root");
    setupMLPCuts("mlpmom_cutg_gen8_new.root");
    setupSectorCorrs(kPM,"sector_factor_pm.root");
    setupSectorCorrs(kPP,"sector_factor_pp.root");
    setupSectorCorrs(kMM,"sector_factor_mm.root");
    setupPiotrFactor("piotrfactor_gen8_consistent_w6momdepcuts.root");

    TH3F *p3DEffEle[6][NWEIGHTS+1];
    TH3F *p3DEffPos[6][NWEIGHTS+1];
    TH3F *p3DAccEle[6][NWEIGHTS+1];
    TH3F *p3DAccPos[6][NWEIGHTS+1];
    readAccEffMatrices(p3DAccEle, p3DAccPos, p3DEffEle, p3DEffPos);

    HLoop* loop = new HLoop(kTRUE);  // kTRUE : create Hades  (needed to work with standard eventstructure)
    TString readCategories = "-*,+HParticleCand,+HParticleEvtInfo,+HStart2Hit,+HStart2Cal,+HGeantKine,+HGeantMdc,+HGeantTof,+HGeantRpc,+HGeantShower";
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
//    HCategory* geantTofCat = (HCategory*)HCategoryManager::getCategory(catTofGeantRaw);
//    HCategory* geantRpcCat = (HCategory*)HCategoryManager::getCategory(catRpcGeantRaw);

    TString outfileSingle = outfile;
    outfileSingle.ReplaceAll(".root","_single.root");
    HHistMap hMSingle(outfileSingle);
    hMSingle.setSilentFail(kTRUE);
    TString outfilePair = outfile;
    outfilePair.ReplaceAll(".root","_pair.root");
    HHistMap hMPair(outfilePair);
    hMPair.setSilentFail(kTRUE);
    TString outAsciiPM = outfile;
    outAsciiPM.ReplaceAll(".root","_pm.txt");
    ofstream out_pm(outAsciiPM.Data());
    TString outAsciiCB = outfile;
    outAsciiCB.ReplaceAll(".root","_cb.txt");
    ofstream out_cb(outAsciiCB.Data());

    HParticleAngleCor trackcor;
    trackcor.setDefaults("apr12");
//    trackcor.setUseEventXYVertex();

    Float_t   _ringNP;
    Float_t   _ringAC;
    Float_t   _metaQa;
    Float_t   _ringHT;
    Float_t   _ringPM;
    Float_t   _beta;
    Float_t   _mdcdEdx;
    Float_t   _theta;
    Float_t   _showerDq;
    Float_t   _mom;
    Float_t   _TofdEdx;

    Float_t mlp = -1000. ;

    TMVA::Reader *reader[2];
    reader[0] = new TMVA::Reader( "!Color:!Silent" );
    reader[1] = new TMVA::Reader( "!Color:!Silent" );

    reader[0]->AddVariable( "ringNP", &_ringNP );
    reader[1]->AddVariable( "ringNP", &_ringNP );
    reader[0]->AddVariable( "ringAC", &_ringAC );
    reader[1]->AddVariable( "ringAC", &_ringAC );
    reader[0]->AddVariable( "metaQa", &_metaQa );
    reader[1]->AddVariable( "metaQa", &_metaQa );
    reader[0]->AddVariable( "ringHT", &_ringHT );
    reader[1]->AddVariable( "ringHT", &_ringHT );
//    reader[0]->AddVariable( "ringPM", &_ringPM );
//    reader[1]->AddVariable( "ringPM", &_ringPM );
    reader[0]->AddVariable( "beta",  &_beta );
    reader[1]->AddVariable( "beta",  &_beta );
    reader[0]->AddVariable( "mdcdEdx",  &_mdcdEdx );
    reader[1]->AddVariable( "mdcdEdx",  &_mdcdEdx );
    reader[0]->AddVariable( "theta",  &_theta );
    reader[1]->AddVariable( "theta",  &_theta );
    reader[0]->AddVariable( "showerDq",  &_showerDq );
    reader[0]->AddVariable( "mom",  &_mom );
    reader[1]->AddVariable( "tofdEdx",  &_TofdEdx );

    TString sys0_weights_file = "sys0_weights/TMVAClassification_MLP_NPACmetaQaHTbetamdcdEdxthetashowerDqmom.weights.xml";
    TString sys1_weights_file = "sys1_weights/TMVAClassification_MLP_NPACmetaQaHTbetamdcdEdxthetatofdEdx_beta085.weights.xml";
    //TString sys0_weights_file = "TMVA_sys0_gen7_ACNPHTPMshowerbetadedxtoftheta_metaQa_EXP_lotstraindata.xml";
    //TString sys1_weights_file = "TMVA_sys1_gen7_ACNPHTPMshowerbetadedxtoftheta_metaQa_EXP_lotstraindata.xml";
    reader[0]->BookMVA( "MLP", sys0_weights_file );
    reader[1]->BookMVA( "MLP", sys1_weights_file );


    //------------------------------------------------------------------

    //--------------- begin histo booking -----------------------------------------------------
    //------------------------------------------------------------------------------------------
    Int_t binsxmass  = 100;
    Double_t maxmom  = 2000;
    Double_t maxmass = 1000;

//    const Int_t nbins = 26;
//    Double_t xAxis1[nbins+1] = {0, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.110, 0.130, 0.150, 0.170, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 0.900, 1.};
    const Int_t nbins = 1000;
    Double_t xAxis1[nbins+1];
    for (int bb = 0; bb <= nbins; ++bb) {
        xAxis1[bb] = 0.001*bb;
    }

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
    cout << "Buff size " << gROOT->GetBufferSize() << endl;

    //------------------------------------------------------------------------------------------
    //--------------- Lepton pairs spectra -----------------------------------------------------
    //------------------------------------------------------------------------------------------
/*
    addHists(hMPair,"hcrossXZ",100,-100,50,100,-50,50);
    addHists(hMPair,"hcrossYZ",100,-100,50,100,-50,50);
    addHists(hMPair,"hclosApprXZ",100,-100,50,100,-50,50);
    addHists(hMPair,"hclosApprYZ",100,-100,50,100,-50,50);
    addHists(hMPair,"hvertXZ",100,-100,50,100,-50,50);
    addHists(hMPair,"hvertYZ",100,-100,50,100,-50,50);
    addHists(hMPair,"hminDistVertZ",100,-50,0,100,-50,50);
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
//    addHists(hMPair,"hth1pt" ,45,0,90,120,0,1200);
//    addHists(hMPair,"hth2pt" ,45,0,90,120,0,1200);
    addHists(hMPair,"hoAnglept" ,90,0,180,120,0,1200);
    addHists(hMPair,"hmasspt" ,nbins,xAxis1,120,0,1200);
//    addHists(hMPair,"hth1y" ,45,0,90,100,0,2);
//    addHists(hMPair,"hth2y" ,45,0,90,100,0,2);
    addHists(hMPair,"hoAngley" ,90,0,180,100,0,2);
    addHists(hMPair,"hmassy" ,nbins,xAxis1,100,0,2);
    addHists(hMPair,"hpty" ,120,0,1200,100,0,2);
//    addHists(hMPair,"hp1th1", 100,0,1000,45,0,90);
//    addHists(hMPair,"hp2th2", 100,0,1000,45,0,90);
    addHists(hMPair,"hth1th2" ,45,0,90,45,0,90);
//    addHists(hMPair,"hpt1pt2" ,100,0,1000,100,0,1000);
    addHists(hMPair,"hp1p2" ,100,0,1000,100,0,1000);
    addHists(hMPair,"hp1p2LowM" ,100,0,1000,100,0,1000);
    addHists(hMPair,"hp1p2MidM" ,100,0,1000,100,0,1000);
    addHists(hMPair,"hp1p2HigM" ,100,0,1000,100,0,1000);
//    addHists(hMPair,"hy1y2" ,100,0,2,100,0,2);
    addHists(hMPair,"hzy" ,100,-80,20,100,0,2);
    addHists(hMPair,"hmass" ,nbins,xAxis1);
    addHists(hMPair,"hoAngle" ,2000,0,200);
    addHists(hMPair,"hmtLowM" ,20,0,1000);
    addHists(hMPair,"hmtMidM" ,20,0,1000);
    addHists(hMPair,"hmtHigM" ,20,0,1000);
    addHists(hMPair,"hmtMinusMinvLowM" ,20,0,1000);
    addHists(hMPair,"hmtMinusMinvMidM" ,20,0,1000);
    addHists(hMPair,"hmtMinusMinvHigM" ,20,0,1000);
    
    hMPair.addHistN("TH1F","hmassNP_HC", "hmassNP_HC" ,nbins,xAxis1);
    hMPair.addHistN("TH1F","hmassPP_HC", "hmassPP_HC" ,nbins,xAxis1);
    hMPair.addHistN("TH1F","hmassNN_HC", "hmassNN_HC" ,nbins,xAxis1);
    hMPair.addHist("TH1F","hoAngleNP_HC", "hoAngleNP_HC" ,200,0,200);
    hMPair.addHist("TH1F","hoAnglePP_HC", "hoAnglePP_HC" ,200,0,200);
    hMPair.addHist("TH1F","hoAngleNN_HC", "hoAngleNN_HC" ,200,0,200);
    hMPair.addHistN("TH1F","hmassNP_HC_oa9deg", "hmassNP_HC_oa9deg" ,nbins,xAxis1);
    hMPair.addHistN("TH1F","hmassPP_HC_oa9deg", "hmassPP_HC_oa9deg" ,nbins,xAxis1);
    hMPair.addHistN("TH1F","hmassNN_HC_oa9deg", "hmassNN_HC_oa9deg" ,nbins,xAxis1);
    hMPair.addHist("TH1F","hoAngleNP_HC_oa9deg", "hoAngleNP_HC_oa9deg" ,200,0,200);
    hMPair.addHist("TH1F","hoAnglePP_HC_oa9deg", "hoAnglePP_HC_oa9deg" ,200,0,200);
    hMPair.addHist("TH1F","hoAngleNN_HC_oa9deg", "hoAngleNN_HC_oa9deg" ,200,0,200);

    hMSingle.addHist("TH2F","thetaMomWon_all","thetaMomWon_all",200,-1000,1000,90,0,90);
    hMSingle.addHist("TH2F","thetaMomWon_sel","thetaMomWon_sel",200,-1000,1000,90,0,90);
    hMSingle.addHist("TH2F","thetaMomWon_mlp","thetaMomWon_mlp",200,-1000,1000,90,0,90);
    hMSingle.addHist("TH2F","thetaMomWon_pid","thetaMomWon_pid",200,-1000,1000,90,0,90);
    hMSingle.addHist("TH2F","ShowerMom_HC_sys0_cutRich", "ShowerMom_HC_sys0_cutRich",400,-1000,1000,200,-500,1500);
    hMSingle.addHist("TH2F","ShowerMom_HC_sys0_cutMassRich", "ShowerMom_HC_sys0_cutMassRich",400,-1000,1000,200,-500,1500);

    addSingleHists(hMSingle,"_noSorter");
    addSingleHists(hMSingle,"");
    addSingleHists(hMSingle,"_cutRich");
    addSingleHists(hMSingle,"_cutMLP");
    addSingleHists(hMSingle,"_cutMLPRich");
    addSingleHists(hMSingle,"_cutMLPRich_blob");
    addSingleHists(hMSingle,"_cutMLPRich_notIncomp");
//    addSingleHists(hMSingle,"_cutMLPRich_notNearst");
//    addSingleHists(hMSingle,"_cutMLPRich_notRecurs");
//    addSingleHists(hMSingle,"_cutMLPRich_notNearst_notRecurs");
//    addSingleHists(hMSingle,"_cutMLPRich_notIncomp_notRecurs");
    addSingleHists(hMSingle,"_cutMLPRich_notIncomp_notNearst");
    addSingleHists(hMSingle,"_cutMLPRich_notIncomp_notNearst_notRecurs");
    addSingleHists(hMSingle,"_cutRichShower");
    addSingleHists(hMSingle,"_cutRichShowerMass");

//    addSingleHists(hMSingle,"_2lepHC");
//    addSingleHists(hMSingle,"_2lepHC_pid");
    addSingleHists(hMSingle,"_2lepMLP");
    addSingleHists(hMSingle,"_2lepMLP_pid");
    addSingleHists(hMSingle,"_2lepMLP_pid_notIncomp");
    addSingleHists(hMSingle,"_2lepMLP_pid_notIncomp_notNearst");
    addSingleHists(hMSingle,"_2lepMLP_pid_notIncomp_notNearst_notRecurs");
//    addSingleHists(hMSingle,"_2lepMLP_mlp");
//    addSingleHists(hMSingle,"_2lepCand");
//    addSingleHists(hMSingle,"_2lepCand_pid");
//    addSingleHists(hMSingle,"_2lepCand_mlp");

    hMSingle.addHist("TH2F","nepem_pid","nemep_pid",10,-0.5,9.5,10,-0.5,9.5);
    hMSingle.addHist("TH2F","nepem_pid_notIncomp","nemep_pid_notIncomp",10,-0.5,9.5,10,-0.5,9.5);
    hMSingle.addHist("TH2F","nepem_pid_notIncomp_notNearst","nemep_pid_notIncomp_notNearst",10,-0.5,9.5,10,-0.5,9.5);
    hMSingle.addHist("TH2F","nepem_pid_notIncomp_notNearst_notRecurs","nemep_pid_notIncomp_notNearst_notRecurs",10,-0.5,9.5,10,-0.5,9.5);
    hMSingle.addHist("TH2F","nepem_hc","nemep_hc",10,-0.5,9.5,10,-0.5,9.5);
    hMSingle.addHist("TH2F","vzx","vzx",200,-200,100,200,-50,50);
    hMSingle.addHist("TH2F","vzr","vzr",200,-200,100,200,0,100);
    hMSingle.addHist("TH2F","outlier","outlier",100,-10,190,500,0,20);

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
        HTool::printProgress(i,nev,1,"Analyze pairs :");
        loop->getSectors(sectors);
        // Remove sector 2 from analysis and take this into account to characterize the event
        sectors[2] = 0;
        sectorCorr = chooseSectorCorr(sectors);

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
            isSimulation = tempFileName.Contains("Au_Au") || tempFileName.Contains("dilepton");
            isEmbedding  = tempFileName.Contains("embedding");
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
            if (!evtinfo->isGoodEvent(kNoPileUpSTART)) continue;
            hMSingle.get("hCounter")->Fill(4);
            if (!evtinfo->isGoodEvent(kNoVETO)) continue;
            hMSingle.get("hCounter")->Fill(5);
            if (!evtinfo->isGoodEvent(kGoodSTARTVETO)) continue;
            hMSingle.get("hCounter")->Fill(6);
            if (!evtinfo->isGoodEvent(kGoodSTARTMETA)) continue;
            hMSingle.get("hCounter")->Fill(7);
            if (sectorCorr == kSkip) continue;
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


        Int_t binHour     = (dayOfYear-96)*24*60+hour*60+min+1;
        Float_t valueHour = hMSingle.get("hEvtHour")->GetXaxis()->GetBinCenter(binHour);
        Float_t yield_pos_weight[6];
        Float_t yield_neg_weight[6];
        for (int _s = 0; _s < 6; ++_s) {
            yield_pos_weight[_s] = yield_pos_fit[_s]->Eval(109)/yield_pos_fit[_s]->Eval(valueHour);
            yield_neg_weight[_s] = yield_neg_fit[_s]->Eval(109)/yield_neg_fit[_s]->Eval(valueHour);
        }
        hMSingle.get("hEvtHour")->Fill(valueHour,1);
        //Float_t oaEmbed = -1;
	//------------------------------------------------------------------------
	// clean vectors and index arrays
	sorter.cleanUp();
	sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);            // reset all flags for flags (0-28) ,reject,used,lepton
	sorter.fill(selectLeptonsBeta);   // fill only good leptons
	sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsLepton);

	//----------------looping data-------------------------
        Int_t size = candCat->getEntries();

        angleNearestBoth.clear();
        angleNearestFull.clear();
        angleNearestSeg.clear();
        nearestIndexSeg.clear();

        Int_t nearestIndexSeg[size];
	Bool_t sorterFlag[size];
        Bool_t richQaFlag[size];
        Bool_t betaFlag[size];
        Bool_t betaFlagTof[size];
        Bool_t betaFlagRpc[size];
        Bool_t showerFlag[size];
	Bool_t leptonFlag[size];
	Bool_t leptonFlagHC[size];
	Bool_t removedNearest[size];
	Bool_t removedIncomplete[size];
	Bool_t removedRecursive[size];
	Bool_t removedEmbedding[size];
	Bool_t removedRecursiveNearestBoth[5][size];
	Bool_t removedRecursiveNearestSeg[5][size];
	Bool_t removedRecursiveNearestFull[5][size];
	Bool_t removedRecursiveNearestBothAbs[5][size];
	Bool_t removedRecursiveNearestSegAbs[5][size];
	Bool_t removedRecursiveNearestFullAbs[5][size];
        Float_t MLPs[size];

        int n_ep_cand = 0;
        int n_em_cand = 0;
        int n_ep_mlp = 0;
        int n_em_mlp = 0;
        int n_ep_mlp_notIncomp = 0;
        int n_em_mlp_notIncomp = 0;
        int n_ep_mlp_notIncomp_notNearst = 0;
        int n_em_mlp_notIncomp_notNearst = 0;
        int n_ep_mlp_notIncomp_notNearst_notRecurs = 0;
        int n_em_mlp_notIncomp_notNearst_notRecurs = 0;
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
            for (int cp = 2; cp <= 6; ++cp) {
                removedRecursiveNearestBoth[cp-2][j]    = kFALSE;
                removedRecursiveNearestSeg[cp-2][j]     = kFALSE;
                removedRecursiveNearestFull[cp-2][j]    = kFALSE;
                removedRecursiveNearestBothAbs[cp-2][j] = kFALSE;
                removedRecursiveNearestSegAbs[cp-2][j]  = kFALSE;
                removedRecursiveNearestFullAbs[cp-2][j] = kFALSE;
            }
            // Track candidates (Leptons AND Hadrons)
	    cand1 = HCategoryManager::getObject(cand1,candCat,j);
            if (isEmbedding) {
                HParticleCandSim *candSim1 = dynamic_cast<HParticleCandSim *>(cand1); 
                if (candSim1 != NULL) {
                    if (!((candSim1->getGeantPID() == 2 || candSim1->getGeantPID() == 3) && candSim1->getGeantParentPID() == -1)) {
                        removedEmbedding[j] = kTRUE;
                    }
                }
            }

            if(cand1->getSystemUsed()==-1)  continue;

	    _ringNP         = cand1->getRingNumPads();
	    _ringAC         = cand1->getRingAmplitude()/_ringNP;
	    _metaQa         = cand1->getMetaMatchQuality();
	    _ringHT         = cand1->getRingHouTra();
	    _ringPM         = cand1->getRingPatternMatrix();
	    _beta           = cand1->getBeta();
	    _mdcdEdx        = cand1->getMdcdEdx();
	    _theta          = cand1->getTheta();
	    _showerDq       = cand1->getShowerDeltaSum();
	    _mom            = cand1->getMomentum();
	    _TofdEdx        = cand1->getTofdEdx();
            Float_t richQa  = cand1->getRichMatchingQuality();
//            Float_t phi     = cand1->getPhi();
            Float_t m2      = _mom*_mom*(1-_beta*_beta)/(_beta*_beta);
            Int_t charge    = cand1->getCharge();
            Int_t sys       = cand1->getSystemUsed();
            Int_t sec       = cand1->getSector();

            if (!selectLeptonsBeta(cand1)) continue;

	    cand1->calc4vectorProperties(HPhysicsConstants::mass(3));
            nearestIndexSeg[j] = findNearestNeighbor(cand1,angleNearestBoth,angleNearestSeg,angleNearestFull,candCat);

            hMSingle.get("thetaMomWon_all")->Fill(_mom*cand1->getCharge(),_theta);
            if (!cand1->isFlagBit(Particle::kIsLepton)) continue;
            sorterFlag[j] = kTRUE;
            hMSingle.get("thetaMomWon_sel")->Fill(_mom*cand1->getCharge(),_theta);

	    mlp = reader[sys]->EvaluateMVA( "MLP" );
            MLPs[j] = mlp;
//            Float_t weight = getEfficiencyFactor(p3DEffEle[0], p3DEffPos[0], _mom, _theta, phi, cand1->getCharge());
//            if (weight < 0.05 && weight > 0) cout << "ERROR exp weight " << weight << endl;

	    //----------------------------------------------------------


//            if (angleNearestSeg[j] == 0) {
///                HParticleCand *candNearest;
//                candNearest = HCategoryManager::getObject(candNearest,catParticleCand,nearestIndexSeg[j]);
//                Float_t angNear = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), candNearest->getPhi(), candNearest->getTheta());
//                cout << "same inner: " << (cand1->getInnerSegInd() == candNearest->getInnerSegInd() )
//                    << ", same outer: " << (cand1->getOuterSegInd() == candNearest->getOuterSegInd() )
//                    << ", same meta: " << (cand1->getMetaHitInd() == candNearest->getMetaHitInd() )
//                    << ", rk's: " << cand1->getChi2() << " " << candNearest->getChi2() 
//                    << ", meta Qa's: " << cand1->getMetaMatchQuality() << " " << candNearest->getMetaMatchQuality() << endl;
//            }

            ///////////////////////////////////////////////////////////////////////////////////
            //
            // Hard cut flags
            //
            ///////////////////////////////////////////////////////////////////////////////////
            // Ring matching
//            bool cutrich;
//            if (sys == 0) {
//                cutrich = richQa < hcutRichQaSys0->Interpolate(_mom*charge);
//            }
//            else {
//                cutrich = richQa < hcutRichQaSys1->Interpolate(_mom*charge);
//            }
//            if (cutrich) richQaFlag[j]=kTRUE;
            if (richQa < 1.5) richQaFlag[j]=kTRUE;

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
            if (isGoodMLP(cand1,mlp)) {
                if (richQaFlag[j]/* && sec == 5*/)
                {
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
	    //*************************************************************
	}  // END SETTING FLAGS
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

                if (k > j/* && leptonFlag[j] && leptonFlag[k]*/) {

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
        } // END RECURSIVE CUT
        // LOOP TO COUNT LEPTONS
	for(Int_t j = 0; j < size; j ++){
            if (isEmbedding && removedEmbedding[j]) continue;
            if (!sorterFlag[j]) continue;
	    cand1 = HCategoryManager::getObject(cand1,candCat,j);

            if (cand1->getCharge() == 1) n_ep_cand++;
            if (cand1->getCharge() ==-1) n_em_cand++;
            if (leptonFlagHC[j]) {
                if (cand1->getCharge() == 1) n_ep_hc++;
                if (cand1->getCharge() ==-1) n_em_hc++;
            }
            if (leptonFlag[j]) {
                if (cand1->getCharge() == 1) n_ep_mlp++;
                if (cand1->getCharge() ==-1) n_em_mlp++;
                if (!removedIncomplete[j]) {
                    if (cand1->getCharge() == 1) n_ep_mlp_notIncomp++;
                    if (cand1->getCharge() ==-1) n_em_mlp_notIncomp++;
                }
                if (!removedIncomplete[j] && !removedNearest[j]) {
                    if (cand1->getCharge() == 1) n_ep_mlp_notIncomp_notNearst++;
                    if (cand1->getCharge() ==-1) n_em_mlp_notIncomp_notNearst++;
                }
                if (!removedIncomplete[j] && !removedNearest[j] && !removedRecursiveNearestSeg[4][j]) {
                    if (cand1->getCharge() == 1) n_ep_mlp_notIncomp_notNearst_notRecurs++;
                    if (cand1->getCharge() ==-1) n_em_mlp_notIncomp_notNearst_notRecurs++;
                }
            }
	    //*************************************************************
	}  // END LOOP TO COUNT LEPTONS
        if (n_ep_mlp + n_em_mlp > 1) {
            hMSingle.get("nepem_pid")->Fill(n_ep_mlp,n_em_mlp);
        }
        if (n_ep_mlp_notIncomp + n_em_mlp_notIncomp > 1) {
            hMSingle.get("nepem_pid_notIncomp")->Fill(n_ep_mlp_notIncomp,n_em_mlp_notIncomp);
        }
        if (n_ep_mlp_notIncomp_notNearst + n_em_mlp_notIncomp_notNearst > 1) {
            hMSingle.get("nepem_pid_notIncomp_notNearst")->Fill(n_ep_mlp_notIncomp_notNearst,n_em_mlp_notIncomp_notNearst);
        }
        if (n_ep_mlp_notIncomp_notNearst_notRecurs + n_em_mlp_notIncomp_notNearst_notRecurs > 1) {
            hMSingle.get("nepem_pid_notIncomp_notNearst_notRecurs")->Fill(n_ep_mlp_notIncomp_notNearst_notRecurs,n_em_mlp_notIncomp_notNearst_notRecurs);
        }
        if (n_ep_hc + n_em_hc > 1) {
            hMSingle.get("nepem_hc")->Fill(n_ep_hc,n_em_hc);
        }
        // END OF PREPARATIONS

        // SINLE LEPTON FILLING HISTOGRAMS
	for(Int_t j = 0; j < size; j ++){
            if (isEmbedding && removedEmbedding[j]) continue;
	    cand1 = HCategoryManager::getObject(cand1,candCat,j);

//	    _ringNP        = cand1->getRingNumPads();
//	    _ringAC        = cand1->getRingAmplitude()/_ringNP;
//	    _metaQa        = cand1->getMetaMatchQuality();
//	    _ringHT        = cand1->getRingHouTra();
//	    _ringPM        = cand1->getRingPatternMatrix();
//	    _beta          = cand1->getBeta();
//	    _mdcdEdx       = cand1->getMdcdEdx();
	    _theta         = cand1->getTheta();
	    _showerDq      = cand1->getShowerDeltaSum();
	    _mom           = cand1->getMomentum();
//	    _TofdEdx       = cand1->getTofdEdx();
//            Float_t phi    = cand1->getPhi();
            Int_t sys      = cand1->getSystemUsed();
            Int_t sec      = cand1->getSector();
	    if(cand1->getSystemUsed()==-1)  continue;
            if(!selectLeptonsBeta(cand1)) continue;
            fillSingle(hMSingle,cand1,"_noSorter",mult_bin,nearestIndexSeg[j]);
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
//            Float_t weight = getEfficiencyFactor(p3DEffEle[0], p3DEffPos[0], _mom, _theta, phi, cand1->getCharge());
//            if (weight < 0.05 && weight > 0) cout << "ERROR exp weight " << weight << endl;
            mlp = MLPs[j];

            fillSingle(hMSingle,cand1,"",mult_bin,nearestIndexSeg[j]);
            if (_showerDq != -1) {
                if (richQaFlag[j]) {
                    if     (sys==0) { 
                        hMSingle.get("ShowerMom_HC_sys0_cutRich")->Fill(_mom*cand1->getCharge(),_showerDq);
                    }
                }
                if (richQaFlag[j] && betaFlag[j]) {
                    if     (sys==0) { 
                        hMSingle.get("ShowerMom_HC_sys0_cutMassRich")->Fill(_mom*cand1->getCharge(),_showerDq);
                    }
                }
            }
            if (richQaFlag[j]) {
                fillSingle(hMSingle,cand1,"_cutRich",mult_bin,nearestIndexSeg[j]);
            }
            if (richQaFlag[j] && showerFlag[j]) {
                fillSingle(hMSingle,cand1,"_cutRichShower",mult_bin,nearestIndexSeg[j]);
            }
            if (richQaFlag[j] && showerFlag[j] && betaFlag[j]) {
                fillSingle(hMSingle,cand1,"_cutRichShowerMass",mult_bin,nearestIndexSeg[j]);
            }
            if (isGoodMLP(cand1,mlp)) {
                hMSingle.get("thetaMomWon_mlp")->Fill(_mom*cand1->getCharge(),_theta);
                fillSingle(hMSingle,cand1,"_cutMLP",mult_bin,nearestIndexSeg[j]);
            }
            if (leptonFlag[j])
            {
                if (cand1->getCharge() == -1) {
                    hMSingle.get("hEmHour_ided",sec)->Fill(valueHour,weight_yield);
                }
                if (cand1->getCharge() == 1) {
                    hMSingle.get("hEpHour_ided",sec)->Fill(valueHour,weight_yield);
                }
                fillSingle(hMSingle,cand1,"_cutMLPRich",mult_bin,nearestIndexSeg[j]);
                if (cand1->getTheta() > 75 && gHades->getCurrentEvent()->getHeader()->getVertexReco().getZ() < -50) {
                    fillSingle(hMSingle,cand1,"_cutMLPRich_blob",mult_bin,nearestIndexSeg[j]);
                }
                if (!removedIncomplete[j]) {
                    fillSingle(hMSingle,cand1,"_cutMLPRich_notIncomp",mult_bin,nearestIndexSeg[j]);
                }
/*
                if (!removedRecursiveNearestSeg[4][j]) {
                    fillSingle(hMSingle,cand1,"_cutMLPRich_notRecurs",mult_bin,nearestIndexSeg[j]);
                }
*/
/*
                if (!removedNearest[j]) {
                    fillSingle(hMSingle,cand1,"_cutMLPRich_notNearst",mult_bin,nearestIndexSeg[j]);
                }
*/
/*
                if (!removedNearest[j] && !removedRecursiveNearestSeg[4][j]) {
                    fillSingle(hMSingle,cand1,"_cutMLPRich_notNearst_notRecurs",mult_bin,nearestIndexSeg[j]);
                }
*/
/*
                if (!removedIncomplete[j] && !removedRecursiveNearestSeg[4][j]) {
                    fillSingle(hMSingle,cand1,"_cutMLPRich_notIncomp_notRecurs",mult_bin,nearestIndexSeg[j]);
                }
*/
                if (!removedIncomplete[j] && !removedNearest[j]) {
                    fillSingle(hMSingle,cand1,"_cutMLPRich_notIncomp_notNearst",mult_bin,nearestIndexSeg[j]);
                }
                if (!removedIncomplete[j] && !removedNearest[j] && !removedRecursiveNearestSeg[4][j]) {
                    fillSingle(hMSingle,cand1,"_cutMLPRich_notIncomp_notNearst_notRecurs",mult_bin,nearestIndexSeg[j]);
                    if (cand1->getCharge() == -1) {
                        hMSingle.get("hEmHour_open",sec)->Fill(valueHour,weight_yield);
                    }
                    if (cand1->getCharge() == 1) {
                        hMSingle.get("hEpHour_open",sec)->Fill(valueHour,weight_yield);
                    }
                }
                hMSingle.get("thetaMomWon_pid")->Fill(_mom*cand1->getCharge(),_theta);
            }
/*
            if (n_ep_hc + n_em_hc >= 2) {
                fillSingle(hMSingle,cand1,"_2lepHC",mult_bin,nearestIndexSeg[j]);
                if (leptonFlagHC[j]) {
                    fillSingle(hMSingle,cand1,"_2lepHC_pid",mult_bin,nearestIndexSeg[j]);
                }
            }
*/
/*
            if (n_ep_cand + n_em_cand >= 2) {
                fillSingle(hMSingle,cand1,"_2lepCand",mult_bin,nearestIndexSeg[j]);
                if (isGoodMLP(cand1,mlp)) {
                    fillSingle(hMSingle,cand1,"_2lepCand_mlp",mult_bin,nearestIndexSeg[j]);
                }
                if (leptonFlag[j]) {
                    fillSingle(hMSingle,cand1,"_2lepCand_pid",mult_bin,nearestIndexSeg[j]);
                }
            }
*/
            if (n_ep_mlp + n_em_mlp >= 2) {
                fillSingle(hMSingle,cand1,"_2lepMLP",mult_bin,nearestIndexSeg[j]);
/*
                if (isGoodMLP(cand1,mlp)) {
                    fillSingle(hMSingle,cand1,"_2lepMLP_mlp",mult_bin,nearestIndexSeg[j]);
                }
*/
                if (leptonFlag[j]) {
                    fillSingle(hMSingle,cand1,"_2lepMLP_pid",mult_bin,nearestIndexSeg[j]);
                }
            }
            if (n_ep_mlp_notIncomp + n_em_mlp_notIncomp >= 2) {
                if (leptonFlag[j]) {
                    if (!removedIncomplete[j]) {
                        fillSingle(hMSingle,cand1,"_2lepMLP_pid_notIncomp",mult_bin,nearestIndexSeg[j]);
                    }
                }
            }
            if (n_ep_mlp_notIncomp_notNearst + n_em_mlp_notIncomp_notNearst >= 2) {
                if (leptonFlag[j]) {
                    if (!removedIncomplete[j] && !removedRecursiveNearestSeg[4][j]) {
                        fillSingle(hMSingle,cand1,"_2lepMLP_pid_notIncomp_notNearst",mult_bin,nearestIndexSeg[j]);
                    }
                }
            }
            if (n_ep_mlp_notIncomp_notNearst_notRecurs + n_em_mlp_notIncomp_notNearst_notRecurs >= 2) {
                if (leptonFlag[j]) {
                    if (!removedIncomplete[j] && !removedNearest[j] && !removedRecursiveNearestSeg[4][j]) {
                        fillSingle(hMSingle,cand1,"_2lepMLP_pid_notIncomp_notNearst_notRecurs",mult_bin,nearestIndexSeg[j]);
                    }
                }
            }

            //*************************************************************
        }  // EBD SINGLE LEPTONS FILLING HISTOGRAMS

	// PAIRS FILLING HISTOGRAMS
	for(Int_t j = 0; j < size-1; j ++){ // first candidate
	    for(Int_t k = j+1; k < size; k ++){ // second candidate
                cand1 = HCategoryManager::getObject(cand1,candCat,j);
		cand2 = HCategoryManager::getObject(cand2,candCat,k);
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

                if (isEmbedding && removedEmbedding[j]) continue;
                if (!sorterFlag[j]) continue;
                if (isEmbedding && removedEmbedding[k]) continue;
                if (!sorterFlag[k]) continue;

                UInt_t flags;
                HParticleTool::setPairFlags(flags,cand2,cand1);
                if(!HParticleTool::evalPairsFlags(kPairCase1,flags)) continue;

                if (leptonFlagHC[j] && leptonFlagHC[k]) {

                    if (isSimulation) {
                        HParticleCandSim *candSim1 = dynamic_cast<HParticleCandSim *>(cand1);
                        HParticleCandSim *candSim2 = dynamic_cast<HParticleCandSim *>(cand2);
                        if (candSim1 != NULL && candSim2 != NULL) {
                            if (candSim1->isGhostTrack() || candSim2->isGhostTrack()) continue;
                            if (candSim1->getGeantParentPID() != -1 || candSim2->getGeantParentPID() != -1) continue;
                            if ( ! ( (candSim1->getGeantPID() == 2 && candSim1->getCharge() > 0) || (candSim1->getGeantPID() == 3 && candSim1->getCharge() < 0) ) )
                                continue;
                            if ( ! ( (candSim2->getGeantPID() == 2 && candSim2->getCharge() > 0) || (candSim2->getGeantPID() == 3 && candSim2->getCharge() < 0) ) )
                                continue;
                        }
                    }

                    TLorentzVector dilep = (*cand1) + (*cand2);
                    Float_t oAngle   = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), cand2->getPhi(), cand2->getTheta());
                    Float_t invM     = dilep.M();

                    fillPair(hMPair, "hmass%s_HC", cand1, cand2, invM/1000);
                    fillPair(hMPair, "hoAngle%s_HC", cand1, cand2, oAngle);

                } // leptonFlagHC[j] && leptonFlagHC[k]
                if (leptonFlag[j] && leptonFlag[k]) {
                    if (isSimulation) {
                        HParticleCandSim *candSim1 = dynamic_cast<HParticleCandSim *>(cand1);
                        HParticleCandSim *candSim2 = dynamic_cast<HParticleCandSim *>(cand2);
                        if (candSim1 != NULL && candSim2 != NULL) {
                            if (candSim1->isGhostTrack() || candSim2->isGhostTrack()) continue;
                            if (candSim1->getGeantParentPID() != -1 || candSim2->getGeantParentPID() != -1) continue;
                            if ( ! ( (candSim1->getGeantPID() == 2 && candSim1->getCharge() > 0) || (candSim1->getGeantPID() == 3 && candSim1->getCharge() < 0) ) )
                                continue;
                            if ( ! ( (candSim2->getGeantPID() == 2 && candSim2->getCharge() > 0) || (candSim2->getGeantPID() == 3 && candSim2->getCharge() < 0) ) )
                                continue;
                        }
                    }

                    TLorentzVector dilep = (*cand1) + (*cand2);
                    Float_t oAngle   = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), cand2->getPhi(), cand2->getTheta());
//                    if (oAngle > 130) continue;

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
                    fillPairSpectra(hMPair,cand1,cand2,"",p3DEffEle[0][0],p3DEffPos[0][0],pSectorCorr,ALL_MULTS,NO_CUT);
                    if (mult_bin > 0 && mult_bin < 5) {
                        fillPairSpectra(hMPair,cand1,cand2,"",p3DEffEle[mult_bin][0],p3DEffPos[mult_bin][0],pSectorCorr,mult_bin,NO_CUT);
                    }
                    if(oAngle > oAngleCut) {
                        fillPairSpectra(hMPair,cand1,cand2,"_oa9deg",p3DEffEle[0][0],p3DEffPos[0][0],pSectorCorr,ALL_MULTS,NO_CUT);
                        if (mult_bin > 0 && mult_bin < 5) {
                            fillPairSpectra(hMPair,cand1,cand2,"_oa9deg",p3DEffEle[mult_bin][0],p3DEffPos[mult_bin][0],pSectorCorr,mult_bin,NO_CUT);
                        }
                    }
                    if(!removedRecursiveNearestSeg[4][j] && !removedRecursiveNearestSeg[4][k] && oAngle > oAngleCut) {
                        fillPairSpectra(hMPair,cand1,cand2,"_oa9degRecur",p3DEffEle[0][0],p3DEffPos[0][0],pSectorCorr,ALL_MULTS,NO_CUT);
                        if (mult_bin > 0 && mult_bin < 5) {
                            fillPairSpectra(hMPair,cand1,cand2,"_oa9degRecur",p3DEffEle[mult_bin][0],p3DEffPos[mult_bin][0],pSectorCorr,mult_bin,NO_CUT);
                        }
                    }

                    if (!removedIncomplete[cand1->getIndex()] && !removedIncomplete[cand2->getIndex()]) {
                        fillPairSpectra(hMPair,cand1,cand2,"",p3DEffEle[0][0],p3DEffPos[0][0],pSectorCorr,ALL_MULTS,CUT_INCOMP);
                        if (mult_bin > 0 && mult_bin < 5) {
                            fillPairSpectra(hMPair,cand1,cand2,"",p3DEffEle[mult_bin][0],p3DEffPos[mult_bin][0],pSectorCorr,mult_bin,CUT_INCOMP);
                        }
                        if(oAngle > oAngleCut) {
                            fillPairSpectra(hMPair,cand1,cand2,"_oa9deg",p3DEffEle[0][0],p3DEffPos[0][0],pSectorCorr,ALL_MULTS,CUT_INCOMP);
                            if (mult_bin > 0 && mult_bin < 5) {
                                fillPairSpectra(hMPair,cand1,cand2,"_oa9deg",p3DEffEle[mult_bin][0],p3DEffPos[mult_bin][0],pSectorCorr,mult_bin,CUT_INCOMP);
                            }
                        }
                        if(!removedRecursiveNearestSeg[4][j] && !removedRecursiveNearestSeg[4][k] && oAngle > oAngleCut) {
                            if (dilep.M() < 150 && cand1->getSector() == cand2->getSector()) {
                                Float_t weight_yield_pi0 = 1;
                                if (cand1->getCharge() == cand2->getCharge()) {
                                    if (cand1->getCharge() == 1) {
                                        weight_yield_pi0 = yield_pos_weight[cand1->getSector()]*yield_pos_weight[cand1->getSector()];
                                    }
                                    if (cand1->getCharge() == -1) {
                                        weight_yield_pi0 = yield_neg_weight[cand1->getSector()]*yield_neg_weight[cand1->getSector()];
                                    }
                                    hMSingle.get("hPi0Hour_open",cand1->getSector())->Fill(valueHour,-weight_yield_pi0);
                                }
                                else {
                                    weight_yield_pi0 = yield_neg_weight[cand1->getSector()]*yield_pos_weight[cand1->getSector()];
                                    hMSingle.get("hPi0Hour_open",cand1->getSector())->Fill(valueHour,weight_yield_pi0);
                                }
                            }
                            fillPairSpectra(hMPair,cand1,cand2,"_oa9degRecur",p3DEffEle[0][0],p3DEffPos[0][0],pSectorCorr,ALL_MULTS,CUT_INCOMP);
                            if (mult_bin > 0 && mult_bin < 5) {
                                fillPairSpectra(hMPair,cand1,cand2,"_oa9degRecur",p3DEffEle[mult_bin][0],p3DEffPos[mult_bin][0],pSectorCorr,mult_bin,CUT_INCOMP);
                            }
                        }
                    }


                    if (!removedNearest[cand1->getIndex()] && !removedNearest[cand2->getIndex()] && !removedIncomplete[cand1->getIndex()] && !removedIncomplete[cand2->getIndex()]) {
                        fillPairSpectra(hMPair,cand1,cand2,"",p3DEffEle[0][0],p3DEffPos[0][0],pSectorCorr,ALL_MULTS,CUT_INCOMP_NEARST);
                        if (mult_bin > 0 && mult_bin < 5) {
                            fillPairSpectra(hMPair,cand1,cand2,"",p3DEffEle[mult_bin][0],p3DEffPos[mult_bin][0],pSectorCorr,mult_bin,CUT_INCOMP_NEARST);
                        }
                        if(oAngle > oAngleCut) {
                            fillPairSpectra(hMPair,cand1,cand2,"_oa9deg",p3DEffEle[0][0],p3DEffPos[0][0],pSectorCorr,ALL_MULTS,CUT_INCOMP_NEARST);
                            if (mult_bin > 0 && mult_bin < 5) {
                                fillPairSpectra(hMPair,cand1,cand2,"_oa9deg",p3DEffEle[mult_bin][0],p3DEffPos[mult_bin][0],pSectorCorr,mult_bin,CUT_INCOMP_NEARST);
                            }
                        }
                        if(!removedRecursiveNearestSeg[4][j] && !removedRecursiveNearestSeg[4][k] && oAngle > oAngleCut) {
                            if (cand1->getCharge() != cand2->getCharge()) {
//                                out_pm << cand1->getPhi() << " " << cand1->getTheta() << " " << cand1->getMomentum() << " " << cand2->getPhi() << " " << cand2->getTheta() << " " << cand2->getMomentum() << " "
//                                    << dilep.M() << " " << oAngle << endl;
                            }
                            else {
//                                out_cb << cand1->getPhi() << " " << cand1->getTheta() << " " << cand1->getMomentum() << " " << cand2->getPhi() << " " << cand2->getTheta() << " " << cand2->getMomentum() << " "
//                                    << dilep.M() << " " << oAngle << endl;
                            }
                            fillPairSpectra(hMPair,cand1,cand2,"_oa9degRecur",p3DEffEle[0][0],p3DEffPos[0][0],pSectorCorr,ALL_MULTS,CUT_INCOMP_NEARST);
                            if (mult_bin > 0 && mult_bin < 5) {
                                fillPairSpectra(hMPair,cand1,cand2,"_oa9degRecur",p3DEffEle[mult_bin][0],p3DEffPos[mult_bin][0],pSectorCorr,mult_bin,CUT_INCOMP_NEARST);
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
    hMPair.writeHists("nomap");
    out_pm.close();
    out_cb.close();
    cout<<"####################################################"<<endl;
    return 0;
}

