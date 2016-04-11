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
#include "TRandom.h"
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

#include "commonFunctions.h"

using namespace std;
//------------------------------------------------------------------
static const Float_t oAngleCut      = 9.;
//------------------------------------------------------------------
#include "eventclassifier.h"
#include "selectFunctions.h"
#include __HELPERFC_FILE__
//------------------------------------------------------------

Int_t dieleAna(TString inputlist, TString outfile, Int_t nev=-1, Int_t whichweight = 0)
{

    TH1::SetDefaultSumw2();
    TH3F *p3DEffEle[6][NWEIGHTS+1]; // mult bins, MLP weights + HC
    TH3F *p3DEffPos[6][NWEIGHTS+1];
    TH3F *p3DAccEle[6][NWEIGHTS+1];
    TH3F *p3DAccPos[6][NWEIGHTS+1];
    readAccEffMatrices(p3DAccEle, p3DAccPos, p3DEffEle, p3DEffPos);
    TH2F *smear_ele, *smear_pos;
    TFile *file_smear = new TFile("smearing_matrix.root","read");
    smear_ele = (TH2F*)file_smear->Get("smear_ele");
    smear_pos = (TH2F*)file_smear->Get("smear_pos");

    TRandom random;

/*
    TFile *pEffFile;

    pEffFile = new TFile("Input/EffMatrixMVA2RefAccNewCP_100Mio.root");
    if (pEffFile)
    {
	pEffFile->cd();
	for(Int_t i = 0 ; i < 5 ; i++){
            p3DEffEle[i][6] = (TH3F*) pEffFile->Get(Form("hHistEff3DMult%iNeg",i));
            p3DEffPos[i][6] = (TH3F*) pEffFile->Get(Form("hHistEff3DMult%iPos",i));
     //       p3DEffEle[i] = (TH3F*) pEffFile->Get("hHistEff3DNeg");
      //      p3DEffPos[i] = (TH3F*) pEffFile->Get("hHistEff3DPos");
        }
    }
    else
    {
	Error("DrawFromNtuple constructor","pointer to eff matrix file is NULL");
	for(Int_t i = 0 ; i < 5 ; i++){
	    p3DEffEle[i][6] = NULL;
	    p3DEffPos[i][6] = NULL;
        }
    }
*/

    TH1F *pEventClass;
    TH1F *pEventClass_recur;
    TFile *pEventClassFile;
//    pEventClassFile = new TFile("eventClass_mult_nonempty_4secmult_200fpj_wait.root");
//    pEventClassFile = new TFile("eventClass_target_mult_rplane_nonempty_4secmult.root");
    pEventClassFile = new TFile("eventClass_target_mult_rplane_minmom_nmix_w6.root");
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
    TString readCategories = "";
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

    HGeantKine *kine1;
    HGeantKine *kine2;

    HCategory* kineCat = (HCategory*)HCategoryManager::getCategory(catGeantKine);

    HHistMap hM(outfile.Data());
    hM.setSilentFail(kTRUE);

    //------------------------------------------------------------------
    //--------------- begin histo booking -----------------------------------------------------
    //------------------------------------------------------------------------------------------

    const Int_t nbins = 26;
    Double_t xAxis1[nbins+1] = {0, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.110, 0.130, 0.150, 0.170, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.700, 0.800, 0.900, 1.};

    hM.addHist(new TH1F(TString("hmassNP"),TString("hmassNP"),nbins,xAxis1));
    hM.addHist(new TH1F(TString("hmassPP"),TString("hmassPP"),nbins,xAxis1));
    hM.addHist(new TH1F(TString("hmassNN"),TString("hmassNN"),nbins,xAxis1));
    hM.addHist(new TH1F(TString("hoAngleNP"),TString("hoAngleNP"),2000,0,200));
    hM.addHist(new TH1F(TString("hoAnglePP"),TString("hoAnglePP"),2000,0,200));
    hM.addHist(new TH1F(TString("hoAngleNN"),TString("hoAngleNN"),2000,0,200));
    hM.addHist(new TH1F(TString("hyNP"),TString("hyNP"),100,0,2));
    hM.addHist(new TH1F(TString("hyPP"),TString("hyPP"),100,0,2));
    hM.addHist(new TH1F(TString("hyNN"),TString("hyNN"),100,0,2));
    hM.addHist(new TH1F(TString("hptNP"),TString("hptNP"),100,0,1000));
    hM.addHist(new TH1F(TString("hptPP"),TString("hptPP"),100,0,1000));
    hM.addHist(new TH1F(TString("hptNN"),TString("hptNN"),100,0,1000));
    hM.addHist(new TH2F(TString("hoAnglemassNP"),TString("hoAnglemassNP"),90,0,180,nbins,xAxis1));
    hM.addHist(new TH2F(TString("hoAnglemassPP"),TString("hoAnglemassPP"),90,0,180,nbins,xAxis1));
    hM.addHist(new TH2F(TString("hoAnglemassNN"),TString("hoAnglemassNN"),90,0,180,nbins,xAxis1));
    hM.addHist(new TH2F(TString("hoAngleptNP"),TString("hoAngleptNP"),90,0,180,120,0,1200));
    hM.addHist(new TH2F(TString("hoAngleptPP"),TString("hoAngleptPP"),90,0,180,120,0,1200));
    hM.addHist(new TH2F(TString("hoAngleptNN"),TString("hoAngleptNN"),90,0,180,120,0,1200));
    hM.addHist(new TH2F(TString("hmassptNP"),TString("hmassptNP"),nbins,xAxis1,120,0,1200));
    hM.addHist(new TH2F(TString("hmassptPP"),TString("hmassptPP"),nbins,xAxis1,120,0,1200));
    hM.addHist(new TH2F(TString("hmassptNN"),TString("hmassptNN"),nbins,xAxis1,120,0,1200));
    hM.addHist(new TH2F(TString("hoAngleyNP"),TString("hoAngleyNP"),90,0,180,100,0,2));
    hM.addHist(new TH2F(TString("hoAngleyPP"),TString("hoAngleyPP"),90,0,180,100,0,2));
    hM.addHist(new TH2F(TString("hoAngleyNN"),TString("hoAngleyNN"),90,0,180,100,0,2));
    hM.addHist(new TH2F(TString("hmassyNP"),TString("hmassyNP"),nbins,xAxis1,100,0,2));
    hM.addHist(new TH2F(TString("hmassyPP"),TString("hmassyPP"),nbins,xAxis1,100,0,2));
    hM.addHist(new TH2F(TString("hmassyNN"),TString("hmassyNN"),nbins,xAxis1,100,0,2));
    hM.addHist(new TH2F(TString("hptyNP"),TString("hptyNP"),120,0,1200,100,0,2));
    hM.addHist(new TH2F(TString("hptyPP"),TString("hptyPP"),120,0,1200,100,0,2));
    hM.addHist(new TH2F(TString("hptyNN"),TString("hptyNN"),120,0,1200,100,0,2));
    hM.addHist(new TH2F(TString("hth1th2NP"),TString("hth1th2NP"),90,0,90,90,0,90));
    hM.addHist(new TH2F(TString("hth1th2PP"),TString("hth1th2PP"),90,0,90,90,0,90));
    hM.addHist(new TH2F(TString("hth1th2NN"),TString("hth1th2NN"),90,0,90,90,0,90));
    hM.addHist(new TH2F(TString("hp1p2NP"),TString("hp1p2NP"),100,0,1100,100,0,1100));
    hM.addHist(new TH2F(TString("hp1p2PP"),TString("hp1p2PP"),100,0,1100,100,0,1100));
    hM.addHist(new TH2F(TString("hp1p2NN"),TString("hp1p2NN"),100,0,1100,100,0,1100));

    for (int i = 0; i < 5; ++i) {
        hM.addHist(new TH1F(TString("hmassNP_eff_multbin")+TString::Itoa(i,10),TString("hmassNP_eff_multbin")+TString::Itoa(i,10),nbins,xAxis1));
        hM.addHist(new TH1F(TString("hmassPP_eff_multbin")+TString::Itoa(i,10),TString("hmassPP_eff_multbin")+TString::Itoa(i,10),nbins,xAxis1));
        hM.addHist(new TH1F(TString("hmassNN_eff_multbin")+TString::Itoa(i,10),TString("hmassNN_eff_multbin")+TString::Itoa(i,10),nbins,xAxis1));
        hM.addHist(new TH1F(TString("hoAngleNP_eff_multbin")+TString::Itoa(i,10),TString("hoAngleNP_eff_multbin")+TString::Itoa(i,10),2000,0,200));
        hM.addHist(new TH1F(TString("hoAnglePP_eff_multbin")+TString::Itoa(i,10),TString("hoAnglePP_eff_multbin")+TString::Itoa(i,10),2000,0,200));
        hM.addHist(new TH1F(TString("hoAngleNN_eff_multbin")+TString::Itoa(i,10),TString("hoAngleNN_eff_multbin")+TString::Itoa(i,10),2000,0,200));
        hM.addHist(new TH1F(TString("hyNP_eff_multbin")+TString::Itoa(i,10),TString("hyNP_eff_multbin")+TString::Itoa(i,10),100,0,2));
        hM.addHist(new TH1F(TString("hyPP_eff_multbin")+TString::Itoa(i,10),TString("hyPP_eff_multbin")+TString::Itoa(i,10),100,0,2));
        hM.addHist(new TH1F(TString("hyNN_eff_multbin")+TString::Itoa(i,10),TString("hyNN_eff_multbin")+TString::Itoa(i,10),100,0,2));
        hM.addHist(new TH1F(TString("hptNP_eff_multbin")+TString::Itoa(i,10),TString("hptNP_eff_multbin")+TString::Itoa(i,10),100,0,1000));
        hM.addHist(new TH1F(TString("hptPP_eff_multbin")+TString::Itoa(i,10),TString("hptPP_eff_multbin")+TString::Itoa(i,10),100,0,1000));
        hM.addHist(new TH1F(TString("hptNN_eff_multbin")+TString::Itoa(i,10),TString("hptNN_eff_multbin")+TString::Itoa(i,10),100,0,1000));
        hM.addHist(new TH2F(TString("hoAnglemassNP_eff_multbin")+TString::Itoa(i,10),TString("hoAnglemassNP_eff_multbin")+TString::Itoa(i,10),90,0,180,nbins,xAxis1));
        hM.addHist(new TH2F(TString("hoAnglemassPP_eff_multbin")+TString::Itoa(i,10),TString("hoAnglemassPP_eff_multbin")+TString::Itoa(i,10),90,0,180,nbins,xAxis1));
        hM.addHist(new TH2F(TString("hoAnglemassNN_eff_multbin")+TString::Itoa(i,10),TString("hoAnglemassNN_eff_multbin")+TString::Itoa(i,10),90,0,180,nbins,xAxis1));
        hM.addHist(new TH2F(TString("hoAngleptNP_eff_multbin")+TString::Itoa(i,10),TString("hoAngleptNP_eff_multbin")+TString::Itoa(i,10),90,0,180,120,0,1200));
        hM.addHist(new TH2F(TString("hoAngleptPP_eff_multbin")+TString::Itoa(i,10),TString("hoAngleptPP_eff_multbin")+TString::Itoa(i,10),90,0,180,120,0,1200));
        hM.addHist(new TH2F(TString("hoAngleptNN_eff_multbin")+TString::Itoa(i,10),TString("hoAngleptNN_eff_multbin")+TString::Itoa(i,10),90,0,180,120,0,1200));
        hM.addHist(new TH2F(TString("hmassptNP_eff_multbin")+TString::Itoa(i,10),TString("hmassptNP_eff_multbin")+TString::Itoa(i,10),nbins,xAxis1,120,0,1200));
        hM.addHist(new TH2F(TString("hmassptPP_eff_multbin")+TString::Itoa(i,10),TString("hmassptPP_eff_multbin")+TString::Itoa(i,10),nbins,xAxis1,120,0,1200));
        hM.addHist(new TH2F(TString("hmassptNN_eff_multbin")+TString::Itoa(i,10),TString("hmassptNN_eff_multbin")+TString::Itoa(i,10),nbins,xAxis1,120,0,1200));
        hM.addHist(new TH2F(TString("hoAngleyNP_eff_multbin")+TString::Itoa(i,10),TString("hoAngleyNP_eff_multbin")+TString::Itoa(i,10),90,0,180,100,0,2));
        hM.addHist(new TH2F(TString("hoAngleyPP_eff_multbin")+TString::Itoa(i,10),TString("hoAngleyPP_eff_multbin")+TString::Itoa(i,10),90,0,180,100,0,2));
        hM.addHist(new TH2F(TString("hoAngleyNN_eff_multbin")+TString::Itoa(i,10),TString("hoAngleyNN_eff_multbin")+TString::Itoa(i,10),90,0,180,100,0,2));
        hM.addHist(new TH2F(TString("hmassyNP_eff_multbin")+TString::Itoa(i,10),TString("hmassyNP_eff_multbin")+TString::Itoa(i,10),nbins,xAxis1,100,0,2));
        hM.addHist(new TH2F(TString("hmassyPP_eff_multbin")+TString::Itoa(i,10),TString("hmassyPP_eff_multbin")+TString::Itoa(i,10),nbins,xAxis1,100,0,2));
        hM.addHist(new TH2F(TString("hmassyNN_eff_multbin")+TString::Itoa(i,10),TString("hmassyNN_eff_multbin")+TString::Itoa(i,10),nbins,xAxis1,100,0,2));
        hM.addHist(new TH2F(TString("hptyNP_eff_multbin")+TString::Itoa(i,10),TString("hptyNP_eff_multbin")+TString::Itoa(i,10),120,0,1200,100,0,2));
        hM.addHist(new TH2F(TString("hptyPP_eff_multbin")+TString::Itoa(i,10),TString("hptyPP_eff_multbin")+TString::Itoa(i,10),120,0,1200,100,0,2));
        hM.addHist(new TH2F(TString("hptyNN_eff_multbin")+TString::Itoa(i,10),TString("hptyNN_eff_multbin")+TString::Itoa(i,10),120,0,1200,100,0,2));
        hM.addHist(new TH2F(TString("hth1th2NP_eff_multbin")+TString::Itoa(i,10),TString("hth1th2NP_eff_multbin")+TString::Itoa(i,10),90,0,90,90,0,90));
        hM.addHist(new TH2F(TString("hth1th2PP_eff_multbin")+TString::Itoa(i,10),TString("hth1th2PP_eff_multbin")+TString::Itoa(i,10),90,0,90,90,0,90));
        hM.addHist(new TH2F(TString("hth1th2NN_eff_multbin")+TString::Itoa(i,10),TString("hth1th2NN_eff_multbin")+TString::Itoa(i,10),90,0,90,90,0,90));
        hM.addHist(new TH2F(TString("hp1p2NP_eff_multbin")+TString::Itoa(i,10),TString("hp1p2NP_eff_multbin")+TString::Itoa(i,10),100,0,1100,100,0,1100));
        hM.addHist(new TH2F(TString("hp1p2PP_eff_multbin")+TString::Itoa(i,10),TString("hp1p2PP_eff_multbin")+TString::Itoa(i,10),100,0,1100,100,0,1100));
        hM.addHist(new TH2F(TString("hp1p2NN_eff_multbin")+TString::Itoa(i,10),TString("hp1p2NN_eff_multbin")+TString::Itoa(i,10),100,0,1100,100,0,1100));
    }

    //--------------- end histo booking -----------------------------------------------------

    HGenericEventMixer<HGeantKine> eventmixer;
    eventmixer.setPIDs(2,3,1);
    eventmixer.setBuffSize(80);
    //eventmixer.setBuffSize(2);
    HGenericEventMixer<HGeantKine> eventmixer_eff[5];
    for (int mb = 0; mb < 5; ++mb) {
        eventmixer_eff[mb].setPIDs(2,3,1);
        eventmixer_eff[mb].setBuffSize(80);
        //eventmixer_eff[mb].setBuffSize(2);
    }

    TStopwatch timer;
    timer.Reset();
    timer.Start();

    Float_t impB = -1.;
    Float_t impB_bins[]  = {9.3, 8.1, 6.6, 4.7, 0.};

    Int_t evtsInFile = loop->getEntries();
    if(nev < 0 || nev > evtsInFile ) nev = evtsInFile;

    for(Int_t i = 1; i < nev; i++)
    {
        //----------break if last event is reached-------------
        //if(!gHades->eventLoop(1)) break;
        if(loop->nextEvent(i) <= 0) { cout<<" end recieved "<<endl; break; } // last event reached
        HTool::printProgress(i,nev,1,"Analyze :");
        loop->getSectors(sectors);

	HPartialEvent *fSimul        = ((HRecEvent*)gHades->getCurrentEvent())->getPartialEvent(catSimul);
	HGeantHeader *fSubHeader = (HGeantHeader*)(fSimul->getSubHeader());
	impB = fSubHeader->getImpactParameter();

	Int_t multbin = 0;
        if (impB >= impB_bins[4] && impB <= impB_bins[3]) {multbin=1;} // most central
	if (impB >  impB_bins[3] && impB <= impB_bins[2]) {multbin=2;}
	if (impB >  impB_bins[2] && impB <= impB_bins[1]) {multbin=3;}
	if (impB >  impB_bins[1] && impB <= impB_bins[0]) {multbin=4;} // most peripheral
	if (impB >  impB_bins[0]) {multbin=5;}
/*
	HParticleEvtInfo* evtinfo;
	evtinfo = HCategoryManager::getObject(evtinfo,catParticleEvtInfo,0);

	Int_t multbin = 0;
	Int_t mult_meta = evtinfo->getSumTofMultCut() + evtinfo->getSumRpcMultHitCut();
        if (mult_meta >  60 && mult_meta <=  88) multbin = 4; // most peripheral
        if (mult_meta >  88 && mult_meta <= 121) multbin = 3;
        if (mult_meta > 121 && mult_meta <= 160) multbin = 2;
        if (mult_meta > 160 && mult_meta <= 250) multbin = 1; // most central
        if (mult_meta > 250) multbin = 5;
*/

        if (multbin == 0 || multbin == 5) continue;

        Int_t size = kineCat->getEntries();

        // Additional loop to fill vector
        vector<HGeantKine *> vep;
        vector<HGeantKine *> vem;
        vector<HGeantKine *> vep_eff;
        vector<HGeantKine *> vem_eff;
        vector<HGeantKine *> vep_eff_multbin;
        vector<HGeantKine *> vem_eff_multbin;
	for(Int_t j = 0; j < size; j ++){
	    kine1 = HCategoryManager::getObject(kine1,kineCat,j);
            Float_t vx,vy,vz;
            kine1->getVertex(vx,vy,vz);
	    Float_t vr = TMath::Sqrt(vx*vx+vy*vy);
            if (vz < -60 || vz > 0) continue;
            if (vr > 2) continue;
            Int_t mamaNum = kine1->getParentTrack();

            if (kine1->isInAcceptance()) {
                if (kine1->getTotalMomentum() > 100 && kine1->getTotalMomentum() < 1000) {
                    Float_t px,py,pz;
                    kine1->getMomentum(px,py,pz);
                    TH2F *smear_matr;
                    if (kine1->getID() == 2) {
                        smear_matr = smear_pos;
                    } else {
                        smear_matr = smear_ele;
                    }
                    Float_t mom_ideal = kine1->getTotalMomentum();
                    Float_t mom_reco = smear(mom_ideal,smear_matr,random);
                    Float_t reco_over_ideal = mom_reco / mom_ideal;
                    kine1->setMomentum(px*reco_over_ideal,py*reco_over_ideal,pz*reco_over_ideal);
                    TLorentzVector vec;
                    HParticleTool::getTLorentzVector(kine1,vec,kine1->getID());

                    Float_t mom = vec.Vect().Mag();
                    Float_t the = vec.Theta()*TMath::RadToDeg();
                    Float_t phi = (vec.Phi()+TMath::Pi())*TMath::RadToDeg();
                    Float_t chg = (kine1->getID() == 2) ? 1 : -1;

                    if (kine1->getID() == 3) {
                        vem.push_back(new HGeantKine(*kine1));
                    }
                    if (kine1->getID() == 2) {
                        vep.push_back(new HGeantKine(*kine1));
                    }
                    Float_t eff = 1./getEfficiencyFactor(p3DEffEle[0][6],p3DEffPos[0][6],mom_ideal,the,phi,chg,false,false); // don't debug, don't check min value
                    if (isinf(eff) || isnan(eff)) eff = 0.;
                    if (random.Uniform(1) > eff) {
                        if (kine1->getID() == 3) {
                            vem_eff.push_back(new HGeantKine(*kine1));
                        }
                        if (kine1->getID() == 2) {
                            vep_eff.push_back(new HGeantKine(*kine1));
                        }
                    }
                    Float_t eff_multbin = 1./getEfficiencyFactor(p3DEffEle[multbin][6],p3DEffPos[multbin][6],mom,the,phi,chg,false,false);
                    if (isinf(eff_multbin) || isnan(eff_multbin)) eff_multbin = 0.;
                    if (random.Uniform(1) > eff_multbin) {
                        if (kine1->getID() == 3) {
                            vem_eff_multbin.push_back(new HGeantKine(*kine1));
                        }
                        if (kine1->getID() == 2) {
                            vep_eff_multbin.push_back(new HGeantKine(*kine1));
                        }
                    }
                }
            }
        }

        eventmixer.nextEvent();
        eventmixer.addVector(vep,2);
        eventmixer.addVector(vem,3);
        vector<pair<HGeantKine *, HGeantKine* > >& pairsVec_acc = eventmixer.getMixedVector();

        eventmixer_eff[0].nextEvent();
        eventmixer_eff[0].addVector(vep_eff,2);
        eventmixer_eff[0].addVector(vem_eff,3);
        vector<pair<HGeantKine *, HGeantKine* > >& pairsVec_eff = eventmixer_eff[0].getMixedVector();

        eventmixer_eff[multbin].nextEvent();
        eventmixer_eff[multbin].addVector(vep_eff,2);
        eventmixer_eff[multbin].addVector(vem_eff,3);
        vector<pair<HGeantKine *, HGeantKine* > >& pairsVec_eff_multbin = eventmixer_eff[multbin].getMixedVector();

        for (int imix = 0; imix < 3; ++imix) {
            vector<pair<HGeantKine *, HGeantKine* > > pairsVec;
            TString suffix;
            switch (imix) {
                case 0: 
                    pairsVec = pairsVec_acc;
                    suffix = TString(""); 
                    break;
                case 1: 
                    pairsVec = pairsVec_eff;
                    suffix = TString("_eff_multbin0");
                    break;
                case 2: 
                    pairsVec = pairsVec_eff_multbin;
                    suffix = TString("_eff_multbin")+TString::Itoa(multbin,10);
                    break;
            }
            size = pairsVec.size();

            for (Int_t j = 0; j < size; j ++) {
                pair<HGeantKine*,HGeantKine*>& pair = pairsVec[j];

                kine1 = pair.first;
                kine2 = pair.second;
                TLorentzVector vec1, vec2;
                HParticleTool::getTLorentzVector(kine1,vec1,kine1->getID());
                HParticleTool::getTLorentzVector(kine2,vec2,kine2->getID());

                Float_t mom1 = vec1.Vect().Mag();
                Float_t the1 = vec1.Theta()*TMath::RadToDeg();
                Float_t mom2 = vec2.Vect().Mag();
                Float_t the2 = vec2.Theta()*TMath::RadToDeg();

                TLorentzVector dilep = vec1 + vec2;
                Float_t oAngle = vec1.Angle(vec2.Vect())*TMath::RadToDeg();
                Float_t mass   = dilep.M()/1000;
                Float_t pt     = dilep.Perp();
                Float_t y      = dilep.Rapidity();
                Float_t mbinw  = hM.get("hmassNP")->GetBinWidth(hM.get("hmassNP")->FindBin(mass));

                if (oAngle < 9) continue;

                TString chg = "NP";
                if (kine1->getID() == kine2->getID()) {
                    if (kine1->getID() == 2) chg = "PP";
                    if (kine1->getID() == 3) chg = "NN";
                }

                hM.get( TString("hmass")      +chg+suffix)->Fill(mass,       1./mbinw);
                hM.get( TString("hoAngle")    +chg+suffix)->Fill(oAngle      );
                hM.get( TString("hy")         +chg+suffix)->Fill(y           );
                hM.get( TString("hpt")        +chg+suffix)->Fill(pt          );
                hM.get2(TString("hoAnglemass")+chg+suffix)->Fill(oAngle,mass,1./mbinw);
                hM.get2(TString("hoAnglept")  +chg+suffix)->Fill(oAngle,pt   );
                hM.get2(TString("hmasspt")    +chg+suffix)->Fill(mass,pt,    1./mbinw);
                hM.get2(TString("hoAngley")   +chg+suffix)->Fill(oAngle,y    );
                hM.get2(TString("hmassy")     +chg+suffix)->Fill(mass,y,     1./mbinw);
                hM.get2(TString("hpty")       +chg+suffix)->Fill(pt,y        );
                hM.get2(TString("hth1th2")    +chg+suffix)->Fill(the1,the2   );
                hM.get2(TString("hp1p2")      +chg+suffix)->Fill(mom1,mom2   );
            }
        }
//#define DELETE_MIX
#ifdef DELETE_MIX
        vector <HGeantKine *>* toDel = eventmixer.getObjectsToDelete();
        for (unsigned int ii = 0; ii < toDel->size(); ++ii) {
            delete toDel->at(ii);
        }
        toDel->clear();
        delete toDel;
        vector <HGeantKine *>* toDel_eff = eventmixer_eff[0].getObjectsToDelete();
        for (unsigned int ii = 0; ii < toDel_eff->size(); ++ii) {
            delete toDel_eff->at(ii);
        }
        toDel_eff->clear();
        delete toDel_eff;
        vector <HGeantKine *>* toDel_eff_multbin = eventmixer_eff[multbin].getObjectsToDelete();
        for (unsigned int ii = 0; ii < toDel_eff_multbin->size(); ++ii) {
            delete toDel_eff_multbin->at(ii);
        }
        toDel_eff_multbin->clear();
        delete toDel_eff_multbin;
#endif

    } // end event loop

    timer.Stop();

    hM.getFile()->cd();
    TMacro m1(__DIELEANA_FILE__);
    m1.Write();
    hM.writeHists("nomap");

    cout<<"####################################################"<<endl;
    return 0;
}
