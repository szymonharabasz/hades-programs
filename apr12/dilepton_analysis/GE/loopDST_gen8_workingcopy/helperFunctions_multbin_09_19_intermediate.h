#include "hparticlebtring.h"
#include "TSystem.h"
#include "TRandom.h"

const int NSIGNS = 3;
const int NMULTS = 5;
const int NCUTS = 4;
const TString signs[NSIGNS] = { "NP", "PP", "NN" };
const TString mults[NMULTS] = { "", "_multbin1", "_multbin2", "_multbin3", "_multbin4" };
const TString cuts[NCUTS]  = { 
    "_w0w0",  "_w6w6",  "_w10w10",  "_w12w12",
};

enum PairCharge { kPM, kPP, kMM, };
enum SectorCorr { kNoCorr, kCorr5to6, kCorr4to6d0, kCorr4to6d1, kCorr4to6d2, kSkip };

static const int ALL_MULTS = 0;
static const int NO_CUT = 0;
static const int CUT_INCOMP = 1;
static const int CUT_INCOMP_NEARST = 2;
static map<Int_t, Int_t> countWithin;
static map<Int_t, Int_t> countUniqueWithin;
static map<Int_t, Int_t> countSameSegWithin;
static map<Int_t, Int_t> countSameSegSameOutWithin;
static map<Int_t, Int_t> countSameSegDiffOutWithin;
static map<Int_t, Int_t> countLeptonsWithin;
static map<Int_t, Int_t> countUniqueLeptonsWithin;
static map<Int_t, Int_t> countSameSegLeptonsWithin;
static map<Int_t, Int_t> countSameSegLeptonsLikeWithin;
static map<Int_t, Int_t> countSameSegLeptonsUnlikeWithin;
static TH2F *hSectorCorr[3][kSkip-1];
static TH1F *hPiotrFactor[5]; // weights 0,6,10,12,HC
static SectorCorr sectorCorr;

static HParticleBooker booker;

Float_t myDotProduct(HParticleCand *, HParticleCand *);
Float_t myOpeningAngle(HParticleCand *, HParticleCand *);

TH3F *p3DEffEle_common[NMULTS];
TH3F *p3DEffPos_common[NMULTS];

TH3F *p3DAccEle_common[NMULTS];
TH3F *p3DAccPos_common[NMULTS];

void setupMassCuts(TString filename) {
    TFile *fileCutMass = new TFile(filename);
    if (fileCutMass) {
        hcutMassSys0 = (TH1F*)fileCutMass->Get("m2_cut_sys0");
        hcutMassSys1 = (TH1F*)fileCutMass->Get("m2_cut_sys1");
        if (!hcutMassSys0 || !hcutMassSys1) {
            cout << "No histogram with mass cuts found in file. Exitting..." << endl;
            exit(-1);
        }
    }
    else {
        cout << "No file with mass cuts found. Exitting..." << endl;
        exit(-1);
    }
}

void setupRichCuts(TString filename) {
    TFile *fileCutRichQa = new TFile(filename);
    if (fileCutRichQa) {
        //hcutRichQaSys0 = (TH1F*)fileCutRichQa->Get("hdeltaRichQaP_mlp_sys0_purity94");
        //hcutRichQaSys1 = (TH1F*)fileCutRichQa->Get("hdeltaRichQaP_mlp_sys1_purity94");
        hcutRichQaSys0 = (TH1F*)fileCutRichQa->Get("RichQaMom_cutMLP_sys0_w6_purity89");
        hcutRichQaSys1 = (TH1F*)fileCutRichQa->Get("RichQaMom_cutMLP_sys1_w6_purity89");
        if (!hcutRichQaSys0 || !hcutRichQaSys1) {
            cout << "No histogram with RICH QA cuts found in file. Exitting..." << endl;
            exit(-1);
        }
    }
    else {
        cout << "No file with RICH QA cuts found. Exitting..." << endl;
        exit(-1);
    }
}

void setupMLPCuts(TString filename) {
    TFile *fileCutMLP = new TFile(filename);
    if (fileCutMLP) {
        for (int w = 0; w < NWEIGHTS; ++w) {
            mlplocmin_sys0[w] = (TH1F*)fileCutMLP->Get(TString("mlplocmin_MLPMom_cutRich_sys0_w") + TString::Itoa(w,10));
            mlplocmin_sys1[w] = (TH1F*)fileCutMLP->Get(TString("mlplocmin_MLPMom_cutRich_sys1_w") + TString::Itoa(w,10));
            if (w >= 9) {
                mlplocmin_sys0[w] = (TH1F*)fileCutMLP->Get(TString("mlplocmin_MLPMom_cutRich_sys0_w") + TString::Itoa(8,10));
                mlplocmin_sys1[w] = (TH1F*)fileCutMLP->Get(TString("mlplocmin_MLPMom_cutRich_sys1_w") + TString::Itoa(8,10));
            }
            if (!mlplocmin_sys0[w] || !mlplocmin_sys1[w]) {
                cout << "No histogram with MLP cuts found in file. Exitting..." << endl;
                exit(-1);
            }
        }
    }
    else {
        cout << "No file with MLP cuts found. Exitting..." << endl;
        exit(-1);
    }
}

void setupSectorCorrs(PairCharge chg, TString filename) {
    TFile *f = new TFile(filename);
    // kCorr5to6 is 1 in enum
//    TString histname[5] = { "hoa_6sec", "hoa_5sec", "hoa_4sec_d0", "hoa_4sec_d1", "hoa_4sec_d2" };
    TString histname[5] = { "", "hthetaoa_6sec_0_0_f0", "hthetaoa_6sec_1_0_f0", "hthetaoa_6sec_2_0_f0", "hthetaoa_6sec_3_0_f0" };
    
    for (int corr = kCorr5to6; corr < kSkip; ++corr) {
        hSectorCorr[chg][corr-1] = (TH2F*)f->Get(histname[corr]);
    }
}

void setupPiotrFactor(TString filename) {
    TFile *f = new TFile(filename);
    hPiotrFactor[0] = (TH1F*)f->Get("hoAngleNP_eff1_w0");
    hPiotrFactor[1] = (TH1F*)f->Get("hoAngleNP_eff1_w6");
    hPiotrFactor[2] = (TH1F*)f->Get("hoAngleNP_eff1_w10");
    hPiotrFactor[3] = (TH1F*)f->Get("hoAngleNP_eff1_w12");
    hPiotrFactor[4] = (TH1F*)f->Get("hoAngleNP_eff1_w13");
    if (hPiotrFactor[4] == NULL) {
        hPiotrFactor[4] = (TH1F*)f->Get("hoAngleNP_eff1_w0");
    }
    if (hPiotrFactor[0] == NULL) {
        Error("setupPiotrFactor","pointer to histogram hoAngleNP_eff1_w0 is NULL");
    }
    if (hPiotrFactor[1] == NULL) {
        Error("setupPiotrFactor","pointer to histogram hoAngleNP_eff1_w6 is NULL");
    }
    if (hPiotrFactor[2] == NULL) {
        Error("setupPiotrFactor","pointer to histogram hoAngleNP_eff1_w12 is NULL");
    }
    if (hPiotrFactor[3] == NULL) {
        Error("setupPiotrFactor","pointer to histogram hoAngleNP_eff1_w12 is NULL");
    }
}

SectorCorr chooseSectorCorr(int sectors[6]) {
    vector<int> bad;
    for (int s = 0; s < 6; ++s) {
        if (!sectors[s]) {
            bad.push_back(s);
        }
    }
    if (bad.size() == 0) return kNoCorr;
    if (bad.size() == 1) return kCorr5to6;
    if (bad.size() > 2)  return kSkip;
    if (abs(bad[1]-bad[0]) == 1 || abs(bad[1] - bad[0]) == 5) return kCorr4to6d0;
    if (abs(bad[1]-bad[0]) == 2 || abs(bad[1] - bad[0]) == 4) return kCorr4to6d1;
    if (abs(bad[1]-bad[0]) == 3) return kCorr4to6d2;
    return kSkip; // shpuld not happen
}
/*
void readAccEffMatrices(TH3F *p3DAccEle[5][NWEIGHTS+1], TH3F *p3DAccPos[5][NWEIGHTS+1], TH3F *p3DEffEle[5][NWEIGHTS+1], TH3F *p3DEffPos[5][NWEIGHTS+1]) {
    TFile *pEleEffFile;
    TFile *pPosEffFile;
    TFile *pEleAccFile;
    TFile *pPosAccFile;

    TString filesForWeights[2][9] = {
        {
            "matricesEff_embUrQMD_Phi180Theta90Mom200_eMinus_apr12_gen8_wszystkieWagi_BestRK.root",
            "matricesEff_embUrQMD_Phi180Theta90Mom200_eMinus_apr12_gen8_HC.root",
        },
        {
            "matricesEff_embUrQMD_Phi180Theta90Mom200_ePlus_apr12_gen8_wszystkieWagi_BestRK.root",
            "matricesEff_embUrQMD_Phi180Theta90Mom200_ePlus_apr12_gen8_HC.root",
        }
    };

    for (int w = 0; w < NWEIGHTS; ++w) { // MLP + HC
        for (int mlt = 0; mlt < 5; ++mlt) {
            p3DEffEle[mlt][w] = NULL;
            p3DEffPos[mlt][w] = NULL;
            p3DAccEle[mlt][w] = NULL;
            p3DAccPos[mlt][w] = NULL;
        }
        for (int mlt = 0; mlt < 1; ++mlt) {
            pEleEffFile = new TFile(filesForWeights[0][0]);
            if (pEleEffFile)
            {
                p3DEffEle[mlt][w] = (TH3F*) pEleEffFile->Get(TString("hEffIdeal_eMinus_")+TString::Itoa(w,10));
            }
            else
            {
                Error("readAccEffMatrices","pointer to eff matrix file is NULL");
            }

            pPosEffFile = new TFile(filesForWeights[1][0]);
            if (pPosEffFile)
            {
                p3DEffPos[mlt][w] = (TH3F*) pPosEffFile->Get(TString("hEffIdeal_ePlus_")+TString::Itoa(w,10));
            }
            else
            {
                Error("readAccEffMatrices","pointer to eff matrix file is NULL");
            }

            // *********************** Acceptance corrections *********************************************************

            pEleAccFile = new TFile(filesForWeights[0][0]);
            if (pEleAccFile)
            {
                p3DAccEle[mlt][w] = (TH3F*) pEleAccFile->Get("hacc");
            }
            else
            {
                Error("readAccEffMatrices","pointer to eff matrix file is NULL");
            }

            pPosAccFile = new TFile(filesForWeights[1][0]);
            if (pPosAccFile)
            {
                p3DAccPos[mlt][w] = (TH3F*) pPosAccFile->Get("hacc");
            }
            else
            {
                Error("readAccEffMatrices","pointer to eff matrix file is NULL");
            }
        }
    }
    for (int mlt = 0; mlt < 5; ++mlt) {
        p3DEffEle[mlt][NWEIGHTS] = NULL;
        p3DEffPos[mlt][NWEIGHTS] = NULL;
        p3DAccEle[mlt][NWEIGHTS] = NULL;
        p3DAccPos[mlt][NWEIGHTS] = NULL;
    }
    for (int mlt = 0; mlt < 1; ++mlt) {
        pEleEffFile = new TFile(filesForWeights[0][1]);
        if (pEleEffFile)
        {
            p3DEffEle[mlt][NWEIGHTS] = (TH3F*) pEleEffFile->Get("eff3DeMinus");
        }
        else
        {
            Error("readAccEffMatrices","pointer to eff matrix file is NULL");
        }

        pPosEffFile = new TFile(filesForWeights[1][1]);
        if (pPosEffFile)
        {
            p3DEffPos[mlt][NWEIGHTS] = (TH3F*) pPosEffFile->Get("eff3DePlus");
        }
        else
        {
            Error("readAccEffMatrices","pointer to eff matrix file is NULL");
        }

        pEleAccFile = new TFile(filesForWeights[0][1]);
        if (pEleAccFile)
        {
            p3DAccEle[mlt][NWEIGHTS] = (TH3F*) pEleAccFile->Get("acc3DeMinus");
        }
        else
        {
            Error("readAccEffMatrices","pointer to eff matrix file is NULL");
        }

        pPosAccFile = new TFile(filesForWeights[1][1]);
        if (pPosAccFile)
        {
            p3DAccPos[mlt][NWEIGHTS] = (TH3F*) pPosAccFile->Get("acc3DePlus");
        }
        else
        {
            Error("readAccEffMatrices","pointer to eff matrix file is NULL");
        }
    }

    p3DEffEle_common = p3DEffEle[0][6];
    p3DEffPos_common = p3DEffPos[0][6];
}
*/
void readAccEffMatrices(TH3F *p3DAccEle[5][NWEIGHTS+1], TH3F *p3DAccPos[5][NWEIGHTS+1], TH3F *p3DEffEle[5][NWEIGHTS+1], TH3F *p3DEffPos[5][NWEIGHTS+1]) {
    TFile *pEleEffFile;
    TFile *pPosEffFile;
    TFile *pEleAccFile;
    TFile *pPosAccFile;
    TString filesForWeights[2][NWEIGHTS+1] = {
        {
            "matricesEff_embUrQMD_Electron_isBestRk_waga0_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga0_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga0_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga0_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga0_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga0_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            //"matricesEff_embUrQMDDelta_Electron_isBestRk_waga6_fillMode2_richqag_MLPqag_test.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga6_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga6_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga6_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga6_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga10_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga10_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga12_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga13_fillMode2_fullStat_allSec_richQa1_MLP06.root",
        },
        {
            "matricesEff_embUrQMD_Positron_isBestRk_waga0_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga0_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga0_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga0_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga0_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga0_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            //"matricesEff_embUrQMDDelta_Positron_isBestRk_waga6_fillMode2_richqag_MLPqag_test.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga6_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga6_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga6_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga6_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga10_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga10_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga12_fillMode2_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga13_fillMode2_fullStat_allSec_richQa1_MLP06.root",
        }
    };
/*
    TString filesForWeights[2][NWEIGHTS+1] = {
        {
            "matricesEff_embUrQMD_Electron_isBestRk_waga0_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga0_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga0_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga0_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga0_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga0_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_goodbins.root",
//            "matricesEff_embData_Electron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_MLP06.root",
//            "matricesEff_embData_Electron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_MLP06.root",
//            "matricesEff_embData_Electron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_MLP06.root",
//            "matricesEff_embData_Electron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga10_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga10_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga12_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Electron_isBestRk_waga13_fillMode3_fullStat_allSec_richQa1_goodbins.root",
        },
        {
            "matricesEff_embUrQMD_Positron_isBestRk_waga0_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga0_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga0_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga0_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga0_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga0_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_goodbins.root",
//            "matricesEff_embData_Positron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_MLP06.root",
//            "matricesEff_embData_Positron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_MLP06.root",
//            "matricesEff_embData_Positron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_MLP06.root",
//            "matricesEff_embData_Positron_isBestRk_waga6_fillMode3_fullStat_allSec_richQa1_MLP06.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga10_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga10_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga12_fillMode3_fullStat_allSec_richQa1_goodbins.root",
            "matricesEff_embUrQMD_Positron_isBestRk_waga13_fillMode3_fullStat_allSec_richQa1_goodbins.root",
        }
    };
*/
    for (int _mlt = 0; _mlt < 6; ++_mlt) {
        int mlt;
        if      (_mlt == 0) mlt = 0;
        else if (_mlt == 5) mlt = 5;
        else                mlt = 5 - _mlt;
        for (int w = 0; w < NWEIGHTS+1; ++w) { // MLP + HC
            // Multiplicity bins are numbered in opposite order in the macro for generating matrices than in this one
            p3DEffEle[mlt][w] = NULL;
            p3DEffPos[mlt][w] = NULL;
            p3DAccEle[mlt][w] = NULL;
            p3DAccPos[mlt][w] = NULL;

            pEleEffFile = new TFile(filesForWeights[0][w]);
            if (pEleEffFile)
            {
                p3DEffEle[mlt][w] = (TH3F*) pEleEffFile->Get(TString("heff3D_Electron_mult")+TString::Itoa(_mlt,10));
            }
            else
            {
                Error("readAccEffMatrices","pointer to eff matrix file is NULL");
            }

            pPosEffFile = new TFile(filesForWeights[1][w]);
            if (pPosEffFile)
            {
                p3DEffPos[mlt][w] = (TH3F*) pPosEffFile->Get(TString("heff3D_Positron_mult")+TString::Itoa(_mlt,10));
            }
            else
            {
                Error("readAccEffMatrices","pointer to eff matrix file is NULL");
            }

            cout << "pointers " << mlt << " " << w <<": " << p3DEffEle[mlt][w] << " " << p3DEffPos[mlt][w] << endl;
            // *********************** Acceptance corrections *********************************************************

            //pEleAccFile = new TFile("matricesEffSingle_electrons_richQa_mlp_fullStat_16072014.root");
            //pEleAccFile = new TFile("matricesEff_embData_Phi180Theta45Mom200_Ele_apr12_v0_newPur_isMLP_isRICH_chargeMinus1.root");
            pEleAccFile = new TFile(filesForWeights[0][w]);
            if (pEleAccFile)
            {
                p3DAccEle[mlt][w] = (TH3F*) pEleAccFile->Get(TString("hacc3D_Electron_mult")+TString::Itoa(_mlt,10));
            }
            else
            {
                Error("readAccEffMatrices","pointer to eff matrix file is NULL");
            }

            //pPosAccFile = new TFile("matricesEffSingle_electrons_richQa_mlp_fullStat_16072014.root");
            //pPosAccFile = new TFile("matricesEff_embData_Phi180Theta45Mom200_Pos_apr12_v0_newPur_isMLP_isRICH_chargePlus1.root");
            pPosAccFile = new TFile(filesForWeights[1][w]);
            if (pPosAccFile)
            {
                p3DAccPos[mlt][w] = (TH3F*) pPosAccFile->Get(TString("hacc3D_Positron_mult")+TString::Itoa(_mlt,10));
            }
            else
            {
                Error("readAccEffMatrices","pointer to eff matrix file is NULL");
            }
        }
        p3DEffEle_common[mlt] = p3DEffEle[mlt][6];
        p3DEffPos_common[mlt] = p3DEffPos[mlt][6];
        p3DAccEle_common[mlt] = p3DAccEle[mlt][6];
        p3DAccPos_common[mlt] = p3DAccPos[mlt][6];
        cout << "pointers common " << mlt << ": " << p3DEffEle_common[mlt] << " " << p3DEffPos_common[mlt] << endl;
    }
    p3DEffEle[5][NWEIGHTS] = p3DEffEle[4][NWEIGHTS];
    p3DEffPos[5][NWEIGHTS] = p3DEffPos[4][NWEIGHTS];
    p3DAccEle[5][NWEIGHTS] = p3DAccEle[4][NWEIGHTS];
    p3DAccPos[5][NWEIGHTS] = p3DAccPos[4][NWEIGHTS];
}

void doTrackCorr(HParticleAngleCor &trackcor, HCategory *candCat) {
    Int_t size = candCat->getEntries();
    for(Int_t j = 0; j < size; j ++){
        HParticleCand *cand;
        cand = HCategoryManager::getObject(cand,candCat,j);
        trackcor.recalcSetEmission(cand);  // changes values inside candidate
        trackcor.realignRichRing(cand);    //
    }
}

void fillVectorFromKine(TLorentzVector& v, HGeantKine* kine)
{
    // fills vector v with one leg of a pair from slot in indextable with kine information

    Float_t xmom,ymom,zmom;
    kine->getMomentum(xmom,ymom,zmom);
    Double_t mass  =HPhysicsConstants::mass(kine->getID());
    Double_t energy=TMath::Sqrt(mass*mass + xmom*xmom + ymom*ymom + zmom*zmom);
    v.SetPxPyPzE(xmom,ymom,zmom,energy);

}

Float_t transformPhi(Float_t Phi)
{

    Float_t dPhi;

    if( (dPhi = TMath::RadToDeg() * Phi) < 0.0 )
	dPhi += 360.0;

    return dPhi;
}

Bool_t checkPhi(Float_t phi, Float_t range)
{
    if(TMath::Abs(fmod(phi,60.f)) < range ||
       TMath::Abs(fmod(phi,60.f) > 60-range))
	return kFALSE;
    else
	return kTRUE;
}

Float_t getPiotrFactor(TH1F *pPiotrCorr, Float_t oa)
{
    if (pPiotrCorr == NULL) return 1.;
    //if (oa > 20.) pPiotrCorr->Interpolate(20);;
    if (oa > 20.) return 1.;
    Float_t fCorr = 1.;

    if (pPiotrCorr) {
        Float_t oa_min = pPiotrCorr->GetXaxis()->GetBinCenter(1);
        Float_t oa_maoa = 20.;

        if (oa > oa_min && oa < oa_maoa) {
            fCorr = pPiotrCorr->Interpolate(oa);
        }
        else {
            fCorr = pPiotrCorr->GetBinContent(pPiotrCorr->FindBin(oa));
        }
    }
    return fCorr;
}

Float_t getSectorFactor(TH1F *pSectorCorr, Float_t x)
{
    if (pSectorCorr == NULL) return 1.;
    Float_t fCorr = 1.;

    if (pSectorCorr) {
        Float_t x_min = pSectorCorr->GetXaxis()->GetBinCenter(1);
        Float_t x_max = pSectorCorr->GetXaxis()->GetBinCenter(pSectorCorr->GetNbinsX());

        if (x > x_min && x < x_max) {
            fCorr = pSectorCorr->Interpolate(x);
        }
        else {
            fCorr = pSectorCorr->GetBinContent(pSectorCorr->FindBin(x));
        }
    }
    return fCorr;
}

Float_t getSectorFactor(TH2F *pSectorCorr, Float_t x, Float_t y)
{
    if (pSectorCorr == NULL) return 1.;
    Float_t fCorr = 1.;

    if (pSectorCorr) {
        Float_t x_min = pSectorCorr->GetXaxis()->GetBinCenter(1);
        Float_t x_max = pSectorCorr->GetXaxis()->GetBinCenter(pSectorCorr->GetNbinsX());

        Float_t y_min = pSectorCorr->GetYaxis()->GetBinCenter(1);
        Float_t y_max = pSectorCorr->GetYaxis()->GetBinCenter(pSectorCorr->GetNbinsY());

        if (x > x_min && x < x_max &&
                y > y_min && y < y_max ) {
            fCorr = pSectorCorr->Interpolate(x,y);

        }
        else {
            fCorr = pSectorCorr->
                GetBinContent(pSectorCorr->FindBin(x,y));
        }
    }
    return fCorr;
}

Float_t getEfficiencyFactor(TH3F *p3DEffEle, TH3F *p3DEffPos, TH3F *p3DEffEle_common, TH3F *p3DEffPos_common, Float_t mom, Float_t theta, Float_t phi, Float_t charge, bool debug = false, bool checkmin = true)
{
    Float_t fEff1 = 1.;

    TH3F *ptr3DEffPos = p3DEffPos;
    TH3F *ptr3DEffEle = p3DEffEle;
/*
    Float_t p_min = ptr3DEffEle->GetZaxis()->GetBinCenter(1);
    Float_t p_max = ptr3DEffEle->GetZaxis()->GetBinCenter(ptr3DEffEle->GetNbinsZ());

    Float_t theta_min = ptr3DEffEle->GetYaxis()->GetBinCenter(1);
    Float_t theta_max = ptr3DEffEle->GetYaxis()->GetBinCenter(ptr3DEffEle->GetNbinsY());

    Float_t phi_min = ptr3DEffEle->GetXaxis()->GetBinCenter(1);
    Float_t phi_max = ptr3DEffEle->GetXaxis()->GetBinCenter(ptr3DEffEle->GetNbinsX());
*/
    if (charge==1) // positron
    {
	if (ptr3DEffPos) fEff1 = ptr3DEffPos->
	    GetBinContent(ptr3DEffPos->FindBin(phi,theta,mom));
	else fEff1 = 1.;
    }

    else if (charge==-1) // electron
    {
	if (ptr3DEffEle) fEff1 = ptr3DEffEle->
	    GetBinContent(ptr3DEffEle->FindBin(phi,theta,mom));
	else fEff1 = 1.;
    }

//    // TO CHECK WHAT IS THE INFLUENCE OF 20 PERCENT EFFICIENCY LOSS
//    fEff1 *= 1.2;

    if (debug) cout << " " << mom << " " << " " << theta << " " << phi << " " << fEff1 << endl;

    Float_t fCorrection = 1./fEff1;

    Float_t fEff2 = 1.;

    // ---
    ptr3DEffPos = p3DEffPos_common;
    ptr3DEffEle = p3DEffEle_common;
    if (charge==1) // positron
    {
	if (ptr3DEffPos) fEff2 = ptr3DEffPos->
	    GetBinContent(ptr3DEffPos->FindBin(phi,theta,mom));
	else fEff2 = 1.;
    }

    else if (charge==-1) // electron
    {
	if (ptr3DEffEle) fEff2 = ptr3DEffEle->
	    GetBinContent(ptr3DEffEle->FindBin(phi,theta,mom));
	else fEff2 = 1.;
    }

//    // TO CHECK WHAT IS THE INFLUENCE OF 20 PERCENT EFFICIENCY LOSS
//    fEff1 *= 1.2;
    if ((p3DEffEle == p3DEffEle_common || p3DEffPos == p3DEffPos_common) && fEff1 != fEff2) {
        cout << "something wrong fEff1: " << fEff1 << " fEff2: " << fEff2 << endl;
        gSystem->StackTrace();
        cout << "---" << endl;
        cout << "---" << endl;
        cout << "---" << endl;
    }


    if (fCorrection >= 1. && (fEff2 > 0.05 || !checkmin) && !isinf(fCorrection)) {
	return fCorrection;
    }
    else
    {
	return 0.;
    }

}

Float_t getEfficiencyFactor(TH3F *p3DEffEle, TH3F *p3DEffPos, Float_t mom, Float_t theta, Float_t phi, Float_t charge, bool debug = false, bool checkmin = true)
{
    Float_t fEff1 = 1.;

    TH3F *ptr3DEffPos = p3DEffPos;
    TH3F *ptr3DEffEle = p3DEffEle;
/*
    Float_t p_min = ptr3DEffEle->GetZaxis()->GetBinCenter(1);
    Float_t p_max = ptr3DEffEle->GetZaxis()->GetBinCenter(ptr3DEffEle->GetNbinsZ());

    Float_t theta_min = ptr3DEffEle->GetYaxis()->GetBinCenter(1);
    Float_t theta_max = ptr3DEffEle->GetYaxis()->GetBinCenter(ptr3DEffEle->GetNbinsY());

    Float_t phi_min = ptr3DEffEle->GetXaxis()->GetBinCenter(1);
    Float_t phi_max = ptr3DEffEle->GetXaxis()->GetBinCenter(ptr3DEffEle->GetNbinsX());
*/
    if (charge==1) // positron
    {
	if (ptr3DEffPos) fEff1 = ptr3DEffPos->
	    GetBinContent(ptr3DEffPos->FindBin(phi,theta,mom));
	else fEff1 = 1.;
    }

    else if (charge==-1) // electron
    {
	if (ptr3DEffEle) fEff1 = ptr3DEffEle->
	    GetBinContent(ptr3DEffEle->FindBin(phi,theta,mom));
	else fEff1 = 1.;
    }

//    // TO CHECK WHAT IS THE INFLUENCE OF 20 PERCENT EFFICIENCY LOSS
//    fEff1 *= 1.2;

    if (debug) cout << " " << mom << " " << " " << theta << " " << phi << " " << fEff1 << endl;

    Float_t fCorrection = 1./fEff1;


    if (fCorrection >= 1. && (fEff1 > 0.05 || !checkmin)) {

	return fCorrection;
    }
    else
    {
	return 0.;
    }

}

Float_t getAcceptanceFactor(TH3F *p3DAccEle, TH3F *p3DAccPos, Float_t mom, Float_t theta, Float_t phi, Float_t charge, bool debug = false)
{
    Float_t fAcc1 = 1.;

    TH3F *ptr3DAccPos = p3DAccPos;
    TH3F *ptr3DAccEle = p3DAccEle;
/*
    Float_t p_min = ptr3DAccEle->GetZaxis()->GetBinCenter(1);
    Float_t p_max = ptr3DAccEle->GetZaxis()->GetBinCenter(ptr3DAccEle->GetNbinsZ());

    Float_t theta_min = ptr3DAccEle->GetYaxis()->GetBinCenter(1);
    Float_t theta_max = ptr3DAccEle->GetYaxis()->GetBinCenter(ptr3DAccEle->GetNbinsY());

    Float_t phi_min = ptr3DAccEle->GetXaxis()->GetBinCenter(1);
    Float_t phi_max = ptr3DAccEle->GetXaxis()->GetBinCenter(ptr3DAccEle->GetNbinsX());
*/
    if (charge==1) // positron
    {
	if (ptr3DAccPos) fAcc1 = ptr3DAccPos->
	    GetBinContent(ptr3DAccPos->FindBin(phi,theta,mom));
	else fAcc1 = 1.;
    }

    else if (charge==-1) // electron
    {
	if (ptr3DAccEle) fAcc1 = ptr3DAccEle->
	    GetBinContent(ptr3DAccEle->FindBin(phi,theta,mom));
	else fAcc1 = 1.;
    }

//    // TO CHECK WHAT IS THE INFLUENCE OF 20 PERCENT EFFICIENCY LOSS
//    fAcc1 *= 1.2;

    if (debug) cout << " " << mom << " " << " " << theta << " " << phi << " " << fAcc1 << endl;

    Float_t fCorrection = 1./fAcc1;


    if (fCorrection >= 1. && fAcc1 > 0.05) {

	return fCorrection;
    }
    else
    {
	return 0.;
    }

}

// something that can be indexed and contains floats
template <typename T>
Int_t findNearestNeighbor(HParticleCand *cand1, T &angleNearestBoth, T &angleNearestSeg, T &angleNearestFull, HCategory *candCat) {
    Int_t nearestIndex = -1;
    Int_t j = cand1->getIndex();
    Int_t size = candCat->getEntries();
    angleNearestBoth[j] = -999;
    angleNearestSeg[j] = -999;
    angleNearestFull[j] = -999;
    HParticleCand *cand2;

    for(Int_t k = 0; k < size; k ++){ // second candidate, k > j condition applied below
        if (k == j) continue;
        cand2 = HCategoryManager::getObject(cand2,candCat,k);
//        cand2->calc4vectorProperties();
        Float_t mdcdEdx2 = cand2->getMdcdEdx();
        Bool_t notLeptonlike = mdcdEdx2 > 3;

        if (!cand2->isFakeRejected(-1) && cand2->getChi2() > 100 /* && !(cand1->getMetaHitInd() == cand2->getMetaHitInd() && cand1->getChi2() < cand2->getChi2())*/ && cand2->getChi2() != cand1->getChi2()/* && cand2->getOuterSegInd() == -1*/) {
            Float_t newAngleNearest = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), cand2->getPhi(), cand2->getTheta());

            if (angleNearestSeg[j] == -999 || fabs(newAngleNearest) < fabs(angleNearestSeg[j]) ) {
                angleNearestSeg[j] = newAngleNearest;

                if (cand2->getRichInd() == -1 || notLeptonlike) {
                    angleNearestSeg[j] *= -1;
                }
                nearestIndex = k;
            }
        }
    }
    return nearestIndex;
}


template <typename FuncType>
void countWithinAngle(HParticleCand *cand1, vector<Int_t>& counts, vector<FuncType>& conditions, vector<bool> uniqueInner, Float_t maxAngle) {
    if (conditions.size() != uniqueInner.size()) {
        cerr << "Error: conditions.size() != uniqueInner.size(), this should not happen, returning..." << endl;
    }
    counts.clear();
    HCategory *candCat = HCategoryManager::getCategory(catParticleCand);
    Int_t count[conditions.size()];
    Int_t j = cand1->getIndex();
    Int_t size = candCat->getEntries();
    HParticleCand *cand2;
    set<Int_t> usedInnerInd[conditions.size()];
    set<Float_t> usedChi2[conditions.size()];

    for (unsigned int m = 0; m < conditions.size(); ++m) {
        usedChi2[m].insert(cand1->getChi2());
        usedInnerInd[m].insert(cand1->getInnerSegInd());
    }

    for(Int_t k = 0; k < size; k ++){ // second candidate, k > j condition applied below
        if (k == j) continue;
        cand2 = HCategoryManager::getObject(cand2,candCat,k);

        for (unsigned int m = 0; m < conditions.size(); ++m) {
            if ((*conditions[m])(cand2) && (cand2->getChi2() == -1 || usedChi2[m].find(cand2->getChi2()) == usedChi2[m].end()) && 
                    (!uniqueInner[m] || cand2->getInnerSegInd() == -1|| usedInnerInd[m].find(cand2->getInnerSegInd()) == usedInnerInd[m].end())) {
                Float_t angle = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), cand2->getPhi(), cand2->getTheta());
                if (angle <= maxAngle) {
                    usedChi2[m].insert(cand2->getChi2());
                    usedInnerInd[m].insert(cand2->getInnerSegInd());
                    ++count[m];
                }
            }
        }
    }
    for (unsigned int m = 0; m < conditions.size(); ++m) {
        counts.push_back(count[m]);
    }
}

template <typename FuncType>
Int_t countWithinAngle(HParticleCand *cand1, vector<HParticleCand *>& found, Float_t maxAngle, bool uniqueInner, FuncType condition) {
    HCategory *candCat = HCategoryManager::getCategory(catParticleCand);
    Int_t count = 0;
    found.clear();
    Int_t j = cand1->getIndex();
    Int_t size = candCat->getEntries();
    HParticleCand *cand2;
    set<Int_t> usedInnerInd;
    set<Float_t> usedChi2;
    usedChi2.insert(cand1->getChi2());
    usedInnerInd.insert(cand1->getInnerSegInd());

    for(Int_t k = 0; k < size; k ++){ // second candidate, k > j condition applied below
        if (k == j) continue;
        cand2 = HCategoryManager::getObject(cand2,candCat,k);

        if (condition(cand2) && (cand2->getChi2() == -1 || usedChi2.find(cand2->getChi2()) == usedChi2.end()) && 
                (!uniqueInner || cand2->getInnerSegInd() == -1|| usedInnerInd.find(cand2->getInnerSegInd()) == usedInnerInd.end())) {
            Float_t angle = HParticleTool::getOpeningAngle(cand1->getPhi(), cand1->getTheta(), cand2->getPhi(), cand2->getTheta());
            if (angle <= maxAngle) {
                usedChi2.insert(cand2->getChi2());
                usedInnerInd.insert(cand2->getInnerSegInd());
                found.push_back(cand2);
                ++count;
            }
        }
    }
    return count;
}


void addHistsHC(HHistMap &hM, TString prefix, Int_t nbinx, Double_t x1, Double_t x2, Int_t nbiny=0, Double_t y1=0, Double_t y2=0) { 
    TString type = "TH2F";
    if (nbiny == 0) type = "TH1F";
    for (int i = 0; i < NSIGNS; ++i) {
        for (int b = 0; b < NMULTS; ++b) {
            hM.addHist(type,prefix       +signs[i]+"_HC"+mults[b],               prefix       +signs[i]+mults[b],               nbinx,x1,x2,nbiny,y1,y2);
            hM.addHist(type,prefix       +signs[i]+"_HC"+mults[b]+"_oa9deg",     prefix       +signs[i]+mults[b]+"_oa9deg",     nbinx,x1,x2,nbiny,y1,y2);
            hM.addHist(type,prefix       +signs[i]+"_HC"+mults[b]+"_oa9degRecur",prefix       +signs[i]+mults[b]+"_oa9degRecur",nbinx,x1,x2,nbiny,y1,y2);
            hM.addHist(type,prefix+"Corr"+signs[i]+"_HC"+mults[b],               prefix+"Corr"+signs[i]+mults[b],               nbinx,x1,x2,nbiny,y1,y2);
            hM.addHist(type,prefix+"Corr"+signs[i]+"_HC"+mults[b]+"_oa9deg",     prefix+"Corr"+signs[i]+mults[b]+"_oa9deg",     nbinx,x1,x2,nbiny,y1,y2);
            hM.addHist(type,prefix+"Corr"+signs[i]+"_HC"+mults[b]+"_oa9degRecur",prefix+"Corr"+signs[i]+mults[b]+"_oa9degRecur",nbinx,x1,x2,nbiny,y1,y2);
        }
    }
}

void addHistsHC(HHistMap &hM, TString prefix, Int_t nbinx, Double_t *x) { 
    for (int i = 0; i < NSIGNS; ++i) {
        for (int b = 0; b < NMULTS; ++b) {
            hM.addHist(new TH1F(prefix       +signs[i]+"_HC"+mults[b],               prefix       +signs[i]+mults[b],               nbinx,x));
            hM.addHist(new TH1F(prefix       +signs[i]+"_HC"+mults[b]+"_oa9deg",     prefix       +signs[i]+mults[b]+"_oa9deg",     nbinx,x));
            hM.addHist(new TH1F(prefix       +signs[i]+"_HC"+mults[b]+"_oa9degRecur",prefix       +signs[i]+mults[b]+"_oa9degRecur",nbinx,x));
            hM.addHist(new TH1F(prefix+"Corr"+signs[i]+"_HC"+mults[b],               prefix+"Corr"+signs[i]+mults[b],               nbinx,x));
            hM.addHist(new TH1F(prefix+"Corr"+signs[i]+"_HC"+mults[b]+"_oa9deg",     prefix+"Corr"+signs[i]+mults[b]+"_oa9deg",     nbinx,x));
            hM.addHist(new TH1F(prefix+"Corr"+signs[i]+"_HC"+mults[b]+"_oa9degRecur",prefix+"Corr"+signs[i]+mults[b]+"_oa9degRecur",nbinx,x));
        }
    }
}

void addHistsHC(HHistMap &hM, TString prefix, Int_t nbinx, Double_t *x, Int_t nbiny, Double_t y1, Double_t y2) { 
    for (int i = 0; i < NSIGNS; ++i) {
        for (int b = 0; b < NMULTS; ++b) {
            hM.addHist(new TH2F(prefix       +signs[i]+"_HC"+mults[b],               prefix       +signs[i]+mults[b],               nbinx,x,nbiny,y1,y2));
            hM.addHist(new TH2F(prefix       +signs[i]+"_HC"+mults[b]+"_oa9deg",     prefix       +signs[i]+mults[b]+"_oa9deg",     nbinx,x,nbiny,y1,y2));
            hM.addHist(new TH2F(prefix       +signs[i]+"_HC"+mults[b]+"_oa9degRecur",prefix       +signs[i]+mults[b]+"_oa9degRecur",nbinx,x,nbiny,y1,y2));
            hM.addHist(new TH2F(prefix+"Corr"+signs[i]+"_HC"+mults[b],               prefix+"Corr"+signs[i]+mults[b],               nbinx,x,nbiny,y1,y2));
            hM.addHist(new TH2F(prefix+"Corr"+signs[i]+"_HC"+mults[b]+"_oa9deg",     prefix+"Corr"+signs[i]+mults[b]+"_oa9deg",     nbinx,x,nbiny,y1,y2));
            hM.addHist(new TH2F(prefix+"Corr"+signs[i]+"_HC"+mults[b]+"_oa9degRecur",prefix+"Corr"+signs[i]+mults[b]+"_oa9degRecur",nbinx,x,nbiny,y1,y2));
        }
    }
}

void addHistsHC(HHistMap &hM, TString prefix, Int_t nbinx, Double_t x1, Double_t x2, Int_t nbiny, Double_t *y) { 
    for (int i = 0; i < NSIGNS; ++i) {
        for (int b = 0; b < NMULTS; ++b) {
            hM.addHist(new TH2F(prefix       +signs[i]+"_HC"+mults[b],               prefix       +signs[i]+mults[b],               nbinx,x1,x2,nbiny,y));
            hM.addHist(new TH2F(prefix       +signs[i]+"_HC"+mults[b]+"_oa9deg",     prefix       +signs[i]+mults[b]+"_oa9deg",     nbinx,x1,x2,nbiny,y));
            hM.addHist(new TH2F(prefix       +signs[i]+"_HC"+mults[b]+"_oa9degRecur",prefix       +signs[i]+mults[b]+"_oa9degRecur",nbinx,x1,x2,nbiny,y));
            hM.addHist(new TH2F(prefix+"Corr"+signs[i]+"_HC"+mults[b],               prefix+"Corr"+signs[i]+mults[b],               nbinx,x1,x2,nbiny,y));
            hM.addHist(new TH2F(prefix+"Corr"+signs[i]+"_HC"+mults[b]+"_oa9deg",     prefix+"Corr"+signs[i]+mults[b]+"_oa9deg",     nbinx,x1,x2,nbiny,y));
            hM.addHist(new TH2F(prefix+"Corr"+signs[i]+"_HC"+mults[b]+"_oa9degRecur",prefix+"Corr"+signs[i]+mults[b]+"_oa9degRecur",nbinx,x1,x2,nbiny,y));
        }
    }
}

void addHists(HHistMap &hM, TString prefix, Int_t nbinx, Double_t x1, Double_t x2, Int_t nbiny=0, Double_t y1=0, Double_t y2=0) { 
    TString type = "TH2F";
    if (nbiny == 0) type = "TH1F";
    for (int i = 0; i < NSIGNS; ++i) {
        for (int b = 0; b < NMULTS; ++b) {
            for (int c = 0; c < NCUTS; ++c) {
                hM.addHist(type,prefix       +signs[i]+mults[b]+cuts[c],               prefix       +signs[i]+mults[b]+cuts[c],               nbinx,x1,x2,nbiny,y1,y2);
                hM.addHist(type,prefix       +signs[i]+mults[b]+cuts[c]+"_oa9deg",     prefix       +signs[i]+mults[b]+cuts[c]+"_oa9deg",     nbinx,x1,x2,nbiny,y1,y2);
                hM.addHist(type,prefix       +signs[i]+mults[b]+cuts[c]+"_oa9degRecur",prefix       +signs[i]+mults[b]+cuts[c]+"_oa9degRecur",nbinx,x1,x2,nbiny,y1,y2);
                hM.addHist(type,prefix+"Corr"+signs[i]+mults[b]+cuts[c],               prefix+"Corr"+signs[i]+mults[b]+cuts[c],               nbinx,x1,x2,nbiny,y1,y2);
                hM.addHist(type,prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9deg",     prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9deg",     nbinx,x1,x2,nbiny,y1,y2);
                hM.addHist(type,prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9degRecur",prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9degRecur",nbinx,x1,x2,nbiny,y1,y2);
            }
        }
    }
}

void addHists(HHistMap &hM, TString prefix, Int_t nbinx, Double_t *x) { 
    for (int i = 0; i < NSIGNS; ++i) {
        for (int b = 0; b < NMULTS; ++b) {
            for (int c = 0; c < NCUTS; ++c) {
                hM.addHist(new TH1F(prefix       +signs[i]+mults[b]+cuts[c],               prefix       +signs[i]+mults[b]+cuts[c],               nbinx,x));
                hM.addHist(new TH1F(prefix       +signs[i]+mults[b]+cuts[c]+"_oa9deg",     prefix       +signs[i]+mults[b]+cuts[c]+"_oa9deg",     nbinx,x));
                hM.addHist(new TH1F(prefix       +signs[i]+mults[b]+cuts[c]+"_oa9degRecur",prefix       +signs[i]+mults[b]+cuts[c]+"_oa9degRecur",nbinx,x));
                hM.addHist(new TH1F(prefix+"Corr"+signs[i]+mults[b]+cuts[c],               prefix+"Corr"+signs[i]+mults[b]+cuts[c],               nbinx,x));
                hM.addHist(new TH1F(prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9deg",     prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9deg",     nbinx,x));
                hM.addHist(new TH1F(prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9degRecur",prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9degRecur",nbinx,x));
            }
        }
    }
}

void addHists(HHistMap &hM, TString prefix, Int_t nbinx, Double_t *x, Int_t nbiny, Double_t y1, Double_t y2) { 
    for (int i = 0; i < NSIGNS; ++i) {
        for (int b = 0; b < NMULTS; ++b) {
            for (int c = 0; c < NCUTS; ++c) {
                hM.addHist(new TH2F(prefix       +signs[i]+mults[b]+cuts[c],               prefix       +signs[i]+mults[b]+cuts[c],               nbinx,x,nbiny,y1,y2));
                hM.addHist(new TH2F(prefix       +signs[i]+mults[b]+cuts[c]+"_oa9deg",     prefix       +signs[i]+mults[b]+cuts[c]+"_oa9deg",     nbinx,x,nbiny,y1,y2));
                hM.addHist(new TH2F(prefix       +signs[i]+mults[b]+cuts[c]+"_oa9degRecur",prefix       +signs[i]+mults[b]+cuts[c]+"_oa9degRecur",nbinx,x,nbiny,y1,y2));
                hM.addHist(new TH2F(prefix+"Corr"+signs[i]+mults[b]+cuts[c],               prefix+"Corr"+signs[i]+mults[b]+cuts[c],               nbinx,x,nbiny,y1,y2));
                hM.addHist(new TH2F(prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9deg",     prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9deg",     nbinx,x,nbiny,y1,y2));
                hM.addHist(new TH2F(prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9degRecur",prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9degRecur",nbinx,x,nbiny,y1,y2));
            }
        }
    }
}

void addHists(HHistMap &hM, TString prefix, Int_t nbinx, Double_t x1, Double_t x2, Int_t nbiny, Double_t *y) { 
    for (int i = 0; i < NSIGNS; ++i) {
        for (int b = 0; b < NMULTS; ++b) {
            for (int c = 0; c < NCUTS; ++c) {
                hM.addHist(new TH2F(prefix       +signs[i]+mults[b]+cuts[c],               prefix       +signs[i]+mults[b]+cuts[c],               nbinx,x1,x2,nbiny,y));
                hM.addHist(new TH2F(prefix       +signs[i]+mults[b]+cuts[c]+"_oa9deg",     prefix       +signs[i]+mults[b]+cuts[c]+"_oa9deg",     nbinx,x1,x2,nbiny,y));
                hM.addHist(new TH2F(prefix       +signs[i]+mults[b]+cuts[c]+"_oa9degRecur",prefix       +signs[i]+mults[b]+cuts[c]+"_oa9degRecur",nbinx,x1,x2,nbiny,y));
                hM.addHist(new TH2F(prefix+"Corr"+signs[i]+mults[b]+cuts[c],               prefix+"Corr"+signs[i]+mults[b]+cuts[c],               nbinx,x1,x2,nbiny,y));
                hM.addHist(new TH2F(prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9deg",     prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9deg",     nbinx,x1,x2,nbiny,y));
                hM.addHist(new TH2F(prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9degRecur",prefix+"Corr"+signs[i]+mults[b]+cuts[c]+"_oa9degRecur",nbinx,x1,x2,nbiny,y));
            }
        }
    }
}

void fillPair(HHistMap &hM, const char *pattern, HParticleCand *cand1, HParticleCand *cand2, Float_t x, Float_t y, Float_t weight, bool likeSignOnly = kFALSE, bool dontCorrectNP = kFALSE) {
    char histname[128];
    if((1==cand1->getCharge() && -1==cand2->getCharge()) ||
            (1==cand2->getCharge() && -1==cand1->getCharge())) {
        if (likeSignOnly) {
            return;
        }
        else {
            sprintf(histname, pattern, "NP");
        }
        if (dontCorrectNP) weight = 1;
    }
    if(1==cand1->getCharge()  &&  1==cand2->getCharge()) {
        sprintf(histname, pattern, "PP");
    }
    if(-1==cand1->getCharge() && -1==cand2->getCharge()) {
        sprintf(histname, pattern, "NN");
    }
//    cout << histname << endl;
//    return;
    if (hM.get(histname) != NULL) {
        TString strHistname(histname);
        if (strHistname.Contains("mass")) {
            Int_t foundBinX = hM.get(histname)->GetXaxis()->FindBin(x);
            Float_t binWidthX = hM.get(histname)->GetXaxis()->GetBinWidth(foundBinX);
            Int_t foundBinY = hM.get(histname)->GetYaxis()->FindBin(y);
            Float_t binWidthY = hM.get(histname)->GetYaxis()->GetBinWidth(foundBinY);
            Float_t binArea = binWidthX*binWidthY;
            hM.get2(histname)->Fill(x,y,weight/binArea);
        }
        else {
            hM.get2(histname)->Fill(x,y,weight);
        }
    }
    else {
//        printf("fillPair(... , Float x, Float y, ...): histogram %s not found\n",histname);
    }
}

void fillPair(HHistMap &hM, const char *pattern, HParticleCand *cand1, HParticleCand *cand2, Float_t value, Float_t weight = 1., bool likeSignOnly = kFALSE) {
    TString strPattern(pattern);
    char histname[128];
//    if (TString(pattern).Contains("hmass%s")) {
//        cout << "test: " << pattern << ", " << value << endl;
//    }
    if((1==cand1->getCharge() && -1==cand2->getCharge()) ||
       (1==cand2->getCharge() && -1==cand1->getCharge())) {
        if (likeSignOnly) {
            return;
        }
        else {
            sprintf(histname, pattern, "NP");
        }
    }
    if(1==cand1->getCharge()  &&  1==cand2->getCharge()) {
        sprintf(histname, pattern, "PP");
    }
    if(-1==cand1->getCharge() && -1==cand2->getCharge()) {
        sprintf(histname, pattern, "NN");
    }
    if (hM.get(histname) != NULL) {
        TString strHistname(histname);
        if (strHistname.Contains("hmass") && !strHistname.Contains("hmassy") && !strHistname.Contains("hmasspt")) {
            Int_t foundBin = hM.get(histname)->FindBin(value);
            Float_t binWidth = hM.get(histname)->GetBinWidth(foundBin);
            hM.get(histname)->Fill(value,weight/binWidth);
        }
        else {
            hM.get(histname)->Fill(value,weight);
        }
    }
    else {
//        printf("fillPair(... , Float value, ...): histogram %s not found\n",histname);
    }
}

void fillPairSpectraHC(HHistMap &hM, HParticleCand *cand1, HParticleCand *cand2, TString suffix, 
        TH3F *p3DEffEle, TH3F *p3DEffPos, TH2F *pSectorCorr, Int_t mult_bin, Int_t cut, Float_t eventWeight = 1., Int_t pfind = -1) {

    HVertex PrimVertexReco = ((HEventHeader*)(gHades->getCurrentEvent()->getHeader()))->getVertexReco(); 
    Float_t vz=PrimVertexReco.getZ(); 
    TLorentzVector dilep = (*cand1) + (*cand2);
    Float_t oAngle   = cand1->Angle(cand2->Vect())*TMath::RadToDeg();
    Float_t invM     = dilep.M();
    Float_t pt       = dilep.Perp();
    Float_t rapidity = dilep.Rapidity();
    Float_t mt       = TMath::Sqrt(pt*pt+invM*invM);
    Float_t mt_w     = 1./TMath::Sqrt(mt*mt*mt);

    Float_t mom1    = cand1->getMomentum();
    Float_t theta1  = cand1->getTheta();
    Float_t phi1    = cand1->getPhi();
    Float_t pt1     = cand1->Perp();
    Float_t y1      = cand1->Rapidity();

    Float_t mom2    = cand2->getMomentum();
    Float_t theta2  = cand2->getTheta();
    Float_t phi2    = cand2->getPhi();
    Float_t pt2     = cand2->Perp();
    Float_t y2      = cand2->Rapidity();

    HCategory* fParticleEvtInfoCat =  (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo,kTRUE,"catParticleEvtInfo"); 
    HParticleEvtInfo *event_info = (HParticleEvtInfo*)fParticleEvtInfoCat->getObject(0); 
    Float_t phiRPl = event_info->getRPlanePhi();
    Float_t phiAvg = 0.5*(phi1+phi2);
    Float_t phiDil = dilep.Phi()*TMath::RadToDeg();
    if (phi1 > 270 && phi2 > 270) {
        phiAvg -=360;
    }
    if ((phi1 < 90 && phi2 > 270) || (phi2 < 90 && phi1 > 270)) {
        phiAvg -=180;
    }

    Float_t weight1 = getEfficiencyFactor(p3DEffEle, p3DEffPos, p3DEffEle_common[mult_bin], p3DEffPos_common[mult_bin], mom1, theta1, phi1, cand1->getCharge());
    Float_t weight2 = getEfficiencyFactor(p3DEffEle, p3DEffPos, p3DEffEle_common[mult_bin], p3DEffPos_common[mult_bin], mom2, theta2, phi2, cand2->getCharge());
    Float_t acc1 = getAcceptanceFactor(p3DAccEle_common[mult_bin], p3DAccPos_common[mult_bin], mom1, theta1, phi1, cand1->getCharge());
    Float_t acc2 = getAcceptanceFactor(p3DAccEle_common[mult_bin], p3DAccPos_common[mult_bin], mom2, theta2, phi2, cand2->getCharge());
//    weight1*=acc1;
//    weight2*=acc2;
    TString cuts[3]  = { "", "_notIncomp", "_notIncomp_notNearst" };
    if (!suffix.Contains("Sorter")) {
        if (mult_bin > 0 && mult_bin < 5) {
            suffix = TString("_multbin")+TString::Itoa(mult_bin,10)+cuts[cut]+suffix;
        }
        else {
            suffix = cuts[cut]+suffix;
        }
    }
    Float_t weight = weight1 * weight2;
    if (weight1 < 0.05 || weight2 < 0.05) {
        weight = 0.;
    }
    Float_t factorSector;
    Float_t factorPiotr;

    if (pSectorCorr != NULL) {
        factorSector = getSectorFactor(pSectorCorr, 0.5*(theta1+theta2), oAngle);
    }
    else {
        factorSector = 1.;
    }
    if (pfind > -1 && hPiotrFactor[pfind] != NULL && cand1->getSector() == cand2->getSector()) {
        factorPiotr = getPiotrFactor(hPiotrFactor[pfind], oAngle);
    }
    else {
        factorPiotr = 1.;
    }

//    factorPiotr = 1.;

    Float_t factorSectorPiotr = factorSector/factorPiotr;
    fillPair(hM, TString("hmtmassFine%s_HC")+suffix, cand1, cand2, mt-invM,invM/1000,mt_w*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmtmass%s_HC")+suffix, cand1, cand2, mt-invM,invM/1000,mt_w*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmasspt%s_HC")+suffix, cand1, cand2, invM/1000,pt,eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmass%s_HC")+suffix, cand1, cand2, invM/1000,eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hmassConstWidth%s_HC")+suffix, cand1, cand2, invM/1000,eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hoAngle%s_HC")+suffix, cand1, cand2, oAngle,eventWeight*factorSectorPiotr);
    if (invM < 150) {
        fillPair(hM, TString("hmtMinusMinvLowM%s_HC")+suffix, cand1, cand2, mt-invM, mt_w*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtLowM%s_HC")+suffix, cand1, cand2, mt, mt_w*eventWeight*factorSectorPiotr);
    }
    else if (invM < 550) {
        fillPair(hM, TString("hmtMinusMinvMidM%s_HC")+suffix, cand1, cand2, mt-invM, mt_w*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtMidM%s_HC")+suffix, cand1, cand2, mt, mt_w*eventWeight*factorSectorPiotr);
    }
    else if (invM < 1000) {
        fillPair(hM, TString("hmtMinusMinvHigM%s_HC")+suffix, cand1, cand2, mt-invM, mt_w*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtHigM%s_HC")+suffix, cand1, cand2, mt, mt_w*eventWeight*factorSectorPiotr);
    }
    fillPair(hM, TString("hmtmassFineCorr%s_HC")+suffix, cand1, cand2, mt-invM,invM/1000,mt_w*weight*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmtmassCorr%s_HC")+suffix, cand1, cand2, mt-invM,invM/1000,mt_w*weight*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmassptCorr%s_HC")+suffix, cand1, cand2, invM/1000,pt,weight*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmassCorr%s_HC")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hmassConstWidthCorr%s_HC")+suffix, cand1, cand2, invM, weight*eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hoAngleCorr%s_HC")+suffix, cand1, cand2, oAngle, weight*eventWeight*factorSectorPiotr);
    if (invM < 150) {
        fillPair(hM, TString("hmtMinusMinvLowMCorr%s_HC")+suffix, cand1, cand2, mt-invM, weight*mt_w*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtLowMCorr%s_HC")+suffix, cand1, cand2, mt, weight*mt_w*eventWeight*factorSectorPiotr);
    }
    else if (invM < 550) {
        fillPair(hM, TString("hmtMinusMinvMidMCorr%s_HC")+suffix, cand1, cand2, mt-invM, weight*mt_w*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtMidMCorr%s_HC")+suffix, cand1, cand2, mt, weight*mt_w*eventWeight*factorSectorPiotr);
    }
    else if (invM < 1000) {
        fillPair(hM, TString("hmtMinusMinvHigMCorr%s_HC")+suffix, cand1, cand2, mt-invM, weight*mt_w*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtHigMCorr%s_HC")+suffix, cand1, cand2, mt, weight*mt_w*eventWeight*factorSectorPiotr);
    }
}

void fillPairSpectra(HHistMap &hM, HParticleCand *cand1, HParticleCand *cand2, TString suffix, 
        TH3F *p3DEffEle, TH3F *p3DEffPos, TH2F *pSectorCorr, Int_t mult_bin, Int_t cut, Float_t eventWeight = 1., Int_t pfind = -1) {

    HVertex PrimVertexReco = ((HEventHeader*)(gHades->getCurrentEvent()->getHeader()))->getVertexReco(); 
    Float_t vz=PrimVertexReco.getZ(); 
    TLorentzVector dilep = (*cand1) + (*cand2);
    Float_t oAngle   = cand1->Angle(cand2->Vect())*TMath::RadToDeg();
    Float_t invM     = dilep.M();
    Float_t pt       = dilep.Perp();
    Float_t rapidity = dilep.Rapidity();
    Float_t mt       = TMath::Sqrt(pt*pt+invM*invM);
    Float_t mt_w     = 1./TMath::Sqrt(mt*mt*mt);

    Float_t mom1    = cand1->getMomentum();
    Float_t theta1  = cand1->getTheta();
    Float_t phi1    = cand1->getPhi();
    Float_t pt1     = cand1->Perp();
    Float_t y1      = cand1->Rapidity();

    Float_t mom2    = cand2->getMomentum();
    Float_t theta2  = cand2->getTheta();
    Float_t phi2    = cand2->getPhi();
    Float_t pt2     = cand2->Perp();
    Float_t y2      = cand2->Rapidity();
    
    HCategory* fParticleEvtInfoCat =  (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo,kTRUE,"catParticleEvtInfo"); 
    HParticleEvtInfo *event_info = (HParticleEvtInfo*)fParticleEvtInfoCat->getObject(0); 
    Float_t phiRPl = event_info->getRPlanePhi();
    Float_t phiAvg = 0.5*(phi1+phi2);
    Float_t phiDil = dilep.Phi()*TMath::RadToDeg();
    if (phi1 > 270 && phi2 > 270) {
        phiAvg -=360;
    }
    if ((phi1 < 90 && phi2 > 270) || (phi2 < 90 && phi1 > 270)) {
        phiAvg -=180;
    }

    Float_t weight1 = getEfficiencyFactor(p3DEffEle, p3DEffPos, p3DEffEle_common[mult_bin], p3DEffPos_common[mult_bin], mom1, theta1, phi1, cand1->getCharge());
    Float_t weight2 = getEfficiencyFactor(p3DEffEle, p3DEffPos, p3DEffEle_common[mult_bin], p3DEffPos_common[mult_bin], mom2, theta2, phi2, cand2->getCharge());
    Float_t acc1 = getAcceptanceFactor(p3DAccEle_common[mult_bin], p3DAccPos_common[mult_bin], mom1, theta1, phi1, cand1->getCharge());
    Float_t acc2 = getAcceptanceFactor(p3DAccEle_common[mult_bin], p3DAccPos_common[mult_bin], mom2, theta2, phi2, cand2->getCharge());
//    weight1*=acc1;
//    weight2*=acc2;
    TString cuts[3]  = { "", "_notIncomp", "_notIncomp_notNearst" };
    if (!suffix.Contains("Sorter")) {
        if (mult_bin > 0 && mult_bin < 5) {
            suffix = TString("_multbin")+TString::Itoa(mult_bin,10)+cuts[cut]+suffix;
        }
        else {
            suffix = cuts[cut]+suffix;
        }
    }
    Float_t weight = weight1 * weight2;
    if (weight1 < 0.05 || weight2 < 0.05) {
        weight = 0.;
//        return; // TEST: what happens if you reject low-eff even in RAW?
    }
/*
    HGeomVector base1, base2, dir1, dir2;
    HParticleTool::calcSegVector(cand1->getZ(),cand1->getR(),cand1->getPhi(),cand1->getTheta(),base1,dir1);
    HParticleTool::calcSegVector(cand2->getZ(),cand2->getR(),cand2->getPhi(),cand2->getTheta(),base2,dir2);
    HGeomVector cross   = HParticleTool::calculateCrossPoint(base1,dir1,base2,dir2);
    HGeomVector pClAp   = HParticleTool::calculatePointOfClosestApproach(base1,dir1,base2,dir2);
    HGeomVector verAn   = HParticleTool::calcVertexAnalytical(base1,dir1,base2,dir2);
    Double_t minDist    = HParticleTool::calculateMinimumDistance(base1,dir1,base2,dir2);
    if (cross.X() == -20000) {
        cout << dir1.X() << " " << dir1.Y() << " " << dir1.Z() << endl;
        cout << dir2.X() << " " << dir2.Y() << " " << dir2.Z() << endl;
        cout << "---" << endl;
    }
*/
    Float_t factorSector;
    Float_t factorPiotr;

    if (pSectorCorr != NULL) {
        factorSector = getSectorFactor(pSectorCorr, 0.5*(theta1+theta2), oAngle);
    }
    else {
        factorSector = 1.;
    }
    if (pfind > -1 && hPiotrFactor[pfind] != NULL && cand1->getSector() == cand2->getSector()) {
        factorPiotr = getPiotrFactor(hPiotrFactor[pfind], oAngle);
    }
    else {
        factorPiotr = 1.;
    }


//    factorPiotr = 1.;

    Float_t factorSectorPiotr = factorSector/factorPiotr;
//    cout << "pfind " << pfind << " oAngle " << oAngle << " pSectorCorr " << pSectorCorr << " factorSector " << factorSector << " factorPiotr " << factorPiotr << endl;
//    Float_t factorSectorPiotr = getSectorFactor(pSectorCorr, oAngle);
//    Float_t factorSectorPiotr = getSectorFactor(pSectorCorr, oAngle, 0.5*(theta1+theta2));
//    if (pt > 1000) return;
/*
    fillPair(hM, TString("hcrossXZ%s")+suffix, cand1, cand2, cross.getZ(), cross.getX());
    fillPair(hM, TString("hcrossYZ%s")+suffix, cand1, cand2, cross.getZ(), cross.getY());
    fillPair(hM, TString("hclosApprXZ%s")+suffix, cand1, cand2, pClAp.getZ(), pClAp.getX());
    fillPair(hM, TString("hclosApprYZ%s")+suffix, cand1, cand2, pClAp.getZ(), pClAp.getY());
    fillPair(hM, TString("hvertXZ%s")+suffix, cand1, cand2, verAn.getZ(), verAn.getX());
    fillPair(hM, TString("hvertYZ%s")+suffix, cand1, cand2, verAn.getZ(), verAn.getY());
    fillPair(hM, TString("hminDistVertZ%s")+suffix, cand1, cand2, vz, minDist);
    fillPair(hM,TString("hph1ph2%s")+suffix, cand1, cand2, phi1, phi2,eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hDthetaDphi%s")        +suffix, cand1, cand2,  cand1->getPhi()-cand2->getPhi(),                           cand1->getTheta()-cand2->getTheta(),eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hDthetaDphiSinTheta%s")+suffix, cand1, cand2, (cand1->getPhi()-cand2->getPhi())*TMath::Sin(dilep.Theta()),cand1->getTheta()-cand2->getTheta(),eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hphiDilPlane%s")       +suffix, cand1, cand2, phiDil, phiRPl, eventWeight*factorSectorPiotr);
//    fillPair(hM,TString("hphiAvgDil%s")+suffix, cand1, cand2, phiDil, phiAvg,eventWeight*factorSectorPiotr);
//    fillPair(hM,TString("hthetaAvgDil%s")+suffix, cand1, cand2, dilep.Theta()*TMath::RadToDeg(), 0.5*(cand1->getTheta()+cand2->getTheta()),eventWeight*factorSectorPiotr);
*/
    fillPair(hM, TString("hthetaoa%s")+suffix, cand1, cand2, .5*(theta1+theta2), oAngle,eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hsqrtp1p2oa%s")+suffix, cand1, cand2, TMath::Sqrt(mom1*mom2), oAngle,eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hoAnglemass%s")+suffix, cand1, cand2, oAngle,invM/1000,eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hoAnglemassNoFactors%s")+suffix, cand1, cand2, oAngle,invM/1000,eventWeight,false,false);
    fillPair(hM, TString("hoAnglept%s")+suffix, cand1, cand2, oAngle,pt,eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hoAngleptNoFactors%s")+suffix, cand1, cand2, oAngle,pt,eventWeight,false,false);
    fillPair(hM, TString("hmtmassFine%s")+suffix, cand1, cand2, mt-invM,invM/1000,mt_w*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmtmass%s")+suffix, cand1, cand2, mt-invM,invM/1000,mt_w*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmasspt%s")+suffix, cand1, cand2, invM/1000,pt,eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmassptNoFactors%s")+suffix, cand1, cand2, invM/1000,pt,eventWeight,false,false);
    fillPair(hM, TString("hoAngley%s")+suffix, cand1, cand2, oAngle,rapidity,eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmassy%s")+suffix, cand1, cand2, invM/1000,rapidity,eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hpty%s")+suffix, cand1, cand2, pt,rapidity,eventWeight*factorSectorPiotr,false,false);
    fillPair(hM,TString("hzy%s")+suffix, cand1, cand2, vz, rapidity,eventWeight*factorSectorPiotr);
    if (cand1->getCharge()==-1 && cand2->getCharge()==1) {
        fillPair(hM,TString("htheta1y%s")+suffix, cand2, cand1, theta2, rapidity,eventWeight*factorSectorPiotr);
//        fillPair(hM,TString("htheta2y%s")+suffix, cand2, cand1, theta1, rapidity,eventWeight*factorSectorPiotr);
        fillPair(hM,TString("hztheta1%s")+suffix, cand2, cand1, vz, theta2,eventWeight*factorSectorPiotr);
//        fillPair(hM,TString("hztheta2%s")+suffix, cand2, cand1, vz, theta1,eventWeight*factorSectorPiotr);
//        fillPair(hM, TString("hth1oAngle%s")+suffix, cand2, cand1, theta2,oAngle,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2oAngle%s")+suffix, cand2, cand1, theta1,oAngle,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth1mass%s")+suffix, cand2, cand1, theta2,invM/1000,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2mass%s")+suffix, cand2, cand1, theta1,invM/1000,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth1pt%s")+suffix, cand2, cand1, theta2,pt,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2pt%s")+suffix, cand2, cand1, theta1,pt,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth1y%s")+suffix, cand2, cand1, theta1,rapidity,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2y%s")+suffix, cand2, cand1, theta2,rapidity,eventWeight*factorSectorPiotr,false,false);
        fillPair(hM, TString("hth1th2%s")+suffix, cand2, cand1, theta2, theta1,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hp1th1%s")+suffix, cand2, cand1, mom2, theta2,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hp2th2%s")+suffix, cand2, cand1, mom1, theta1,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM,TString("hpt1pt2%s")+suffix, cand2, cand1, pt2, pt1,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM,TString("hy1y2%s")+suffix, cand2, cand1, y2, y1,eventWeight*factorSectorPiotr);
//        fillPair(hM, TString("hp1p2diff%s")+suffix, cand2, cand1, mom2-mom1,eventWeight*factorSectorPiotr);
    }
    else {
        fillPair(hM,TString("htheta1y%s")+suffix, cand1, cand2, theta1, rapidity,eventWeight*factorSectorPiotr);
//        fillPair(hM,TString("htheta2y%s")+suffix, cand1, cand2, theta2, rapidity,eventWeight*factorSectorPiotr);
        fillPair(hM,TString("hztheta1%s")+suffix, cand1, cand2, vz, theta1,eventWeight*factorSectorPiotr);
//        fillPair(hM,TString("hztheta2%s")+suffix, cand1, cand2, vz, theta2,eventWeight*factorSectorPiotr);
//        fillPair(hM, TString("hth1oAngle%s")+suffix, cand1, cand2, theta1,oAngle,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2oAngle%s")+suffix, cand1, cand2, theta2,oAngle,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth1mass%s")+suffix, cand1, cand2, theta1,invM/1000,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2mass%s")+suffix, cand1, cand2, theta2,invM/1000,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth1pt%s")+suffix, cand1, cand2, theta1,pt,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2pt%s")+suffix, cand1, cand2, theta2,pt,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth1y%s")+suffix, cand1, cand2, theta1,rapidity,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2y%s")+suffix, cand1, cand2, theta2,rapidity,eventWeight*factorSectorPiotr,false,false);
        fillPair(hM, TString("hth1th2%s")+suffix, cand1, cand2, theta1, theta2,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hp1th1%s")+suffix, cand1, cand2, mom1, theta1,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hp2th2%s")+suffix, cand1, cand2, mom2, theta2,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM,TString("hpt1pt2%s")+suffix, cand1, cand2, pt1, pt2,eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM,TString("hy1y2%s")+suffix, cand1, cand2, y1, y2,eventWeight*factorSectorPiotr);
//        fillPair(hM, TString("hp1p2diff%s")+suffix, cand1, cand2, mom1-mom2,eventWeight*factorSectorPiotr);
    }
    if (weight == 0) {
        fillPair(hM, TString("hmassNull%s")+suffix, cand1, cand2, invM/1000,eventWeight*factorSectorPiotr);
    }
    fillPair(hM, TString("hmass%s")+suffix, cand1, cand2, invM/1000,eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hmassM%s")+suffix, cand1, cand2, invM/1000,invM*eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hmassConstWidth%s")+suffix, cand1, cand2, invM,eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hmassSectorFactor%s")+suffix,    cand1, cand2, invM/1000, eventWeight*factorSector);
    fillPair(hM, TString("hmassPiotrFactor%s")+suffix,    cand1, cand2, invM/1000, eventWeight/factorPiotr);
    fillPair(hM, TString("hmassNoFactors%s")+suffix,    cand1, cand2, invM/1000, eventWeight);
    fillPair(hM, TString("hoAngle%s")+suffix, cand1, cand2, oAngle,eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hoAngleSectorFactor")+suffix,    cand1, cand2, oAngle, eventWeight*factorSector);
    fillPair(hM, TString("hoAnglePiotrFactor")+suffix,    cand1, cand2, oAngle, eventWeight/factorPiotr);
    fillPair(hM, TString("hoAngleNoFactors%s")+suffix,    cand1, cand2, oAngle, eventWeight);
    if (sectorCorr == kNoCorr) {
        fillPair(hM, TString("hmass6sec%s")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hoAngle6sec%s")+suffix, cand1, cand2, oAngle, weight*eventWeight*factorSectorPiotr);
    }
    if (sectorCorr == kCorr5to6) {
        fillPair(hM, TString("hmass5sec%s")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hoAngle5sec%s")+suffix, cand1, cand2, oAngle, weight*eventWeight*factorSectorPiotr);
    }
    if (sectorCorr == kCorr4to6d0) {
        fillPair(hM, TString("hmass4secd0%s")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hoAngle4secd0%s")+suffix, cand1, cand2, oAngle, weight*eventWeight*factorSectorPiotr);
    }
    if (sectorCorr == kCorr4to6d1) {
        fillPair(hM, TString("hmass4secd1%s")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hoAngle4secd1%s")+suffix, cand1, cand2, oAngle, weight*eventWeight*factorSectorPiotr);
    }
    if (sectorCorr == kCorr4to6d2) {
        fillPair(hM, TString("hmass4secd2%s")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hoAngle4secd2%s")+suffix, cand1, cand2, oAngle, weight*eventWeight*factorSectorPiotr);
    }
    fillPair(hM, TString("hy%s")+suffix, cand1, cand2, rapidity,eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hpt%s")+suffix, cand1, cand2, pt,eventWeight*factorSectorPiotr);
    if (invM < 150) {
        fillPair(hM,TString("hp1p2LowM%s")+suffix, cand2, cand1, mom2, mom1,eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtMinusMinvLowM%s")+suffix, cand1, cand2, mt-invM, mt_w*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtLowM%s")+suffix, cand1, cand2, mt, mt_w*eventWeight*factorSectorPiotr);
    }
    else if (invM < 550) {
        fillPair(hM,TString("hp1p2MidM%s")+suffix, cand2, cand1, mom2, mom1,eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtMinusMinvMidM%s")+suffix, cand1, cand2, mt-invM, mt_w*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtMidM%s")+suffix, cand1, cand2, mt, mt_w*eventWeight*factorSectorPiotr);
    }
    else if (invM < 1000) {
        fillPair(hM,TString("hp1p2HigM%s")+suffix, cand2, cand1, mom2, mom1,eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtMinusMinvHigM%s")+suffix, cand1, cand2, mt-invM, mt_w*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtHigM%s")+suffix, cand1, cand2, mt, mt_w*eventWeight*factorSectorPiotr);
    }
    if (invM < 200) {
        fillPair(hM, TString("hpt_0_200_%s")+suffix, cand1, cand2, pt,eventWeight*factorSectorPiotr);
    } else if (invM < 400) {
        fillPair(hM, TString("hpt_200_400_%s")+suffix, cand1, cand2, pt,eventWeight*factorSectorPiotr);
    } else if (invM < 600) {
        fillPair(hM, TString("hpt_400_600_%s")+suffix, cand1, cand2, pt,eventWeight*factorSectorPiotr);
    }
    //---
    fillPair(hM,TString("hsqrtp1p2oaCorr%s")+suffix, cand1, cand2, TMath::Sqrt(mom1*mom2), oAngle,weight*eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hoAnglemassCorr%s")+suffix, cand1, cand2, oAngle,invM/1000,weight*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hoAnglemassNoFactorsCorr%s")+suffix, cand1, cand2, oAngle,invM/1000,weight*eventWeight,false,false);
    fillPair(hM, TString("hoAngleptCorr%s")+suffix, cand1, cand2, oAngle,pt,weight*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hoAngleptNoFactorsCorr%s")+suffix, cand1, cand2, oAngle,pt,weight*eventWeight,false,false);
    fillPair(hM, TString("hmtmassFineCorr%s")+suffix, cand1, cand2, mt-invM,invM/1000,mt_w*weight*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmtmassCorr%s")+suffix, cand1, cand2, mt-invM,invM/1000,mt_w*weight*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmassptCorr%s")+suffix, cand1, cand2, invM/1000,pt,weight*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmassptNoFactorsCorr%s")+suffix, cand1, cand2, invM/1000,pt,weight*eventWeight,false,false);
    fillPair(hM, TString("hoAngleyCorr%s")+suffix, cand1, cand2, oAngle,rapidity,weight*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hmassyCorr%s")+suffix, cand1, cand2, invM/1000,rapidity,weight*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM, TString("hptyCorr%s")+suffix, cand1, cand2, pt,rapidity,weight*eventWeight*factorSectorPiotr,false,false);
    fillPair(hM,TString("hzyCorr%s")+suffix, cand1, cand2, vz, rapidity,weight*eventWeight*factorSectorPiotr);
    if (cand1->getCharge()==-1 && cand2->getCharge()==1) {
        fillPair(hM,TString("htheta1yCorr%s")+suffix, cand2, cand1, theta2, rapidity,weight*eventWeight*factorSectorPiotr);
//        fillPair(hM,TString("htheta2yCorr%s")+suffix, cand2, cand1, theta1, rapidity,weight*eventWeight*factorSectorPiotr);
        fillPair(hM,TString("hztheta1Corr%s")+suffix, cand2, cand1, vz, theta2,weight*eventWeight*factorSectorPiotr);
//        fillPair(hM,TString("hztheta2Corr%s")+suffix, cand2, cand1, vz, theta1,weight*eventWeight*factorSectorPiotr);
//        fillPair(hM, TString("hth1oAngleCorr%s")+suffix, cand2, cand1, theta2,oAngle,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2oAngleCorr%s")+suffix, cand2, cand1, theta1,oAngle,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth1massCorr%s")+suffix, cand2, cand1, theta2,invM/1000,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2massCorr%s")+suffix, cand2, cand1, theta1,invM/1000,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth1ptCorr%s")+suffix, cand2, cand1, theta2,pt,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2ptCorr%s")+suffix, cand2, cand1, theta1,pt,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth1yCorr%s")+suffix, cand2, cand1, theta1,rapidity,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2yCorr%s")+suffix, cand2, cand1, theta2,rapidity,weight*eventWeight*factorSectorPiotr,false,false);
        fillPair(hM, TString("hth1th2Corr%s")+suffix, cand2, cand1, theta2, theta1,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hp1th1Corr%s")+suffix, cand2, cand1, mom2, theta2,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hp2th2Corr%s")+suffix, cand2, cand1, mom1, theta1,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM,TString("hpt1pt2Corr%s")+suffix, cand2, cand1, pt2, pt1,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM,TString("hy1y2Corr%s")+suffix, cand2, cand1, y2, y1,weight*eventWeight*factorSectorPiotr);
    }
    else {
        fillPair(hM,TString("htheta1yCorr%s")+suffix, cand1, cand2, theta1, rapidity,weight*eventWeight*factorSectorPiotr);
//        fillPair(hM,TString("htheta2yCorr%s")+suffix, cand1, cand2, theta2, rapidity,weight*eventWeight*factorSectorPiotr);
        fillPair(hM,TString("hztheta1Corr%s")+suffix, cand1, cand2, vz, theta1,weight*eventWeight*factorSectorPiotr);
//        fillPair(hM,TString("hztheta2Corr%s")+suffix, cand1, cand2, vz, theta2,weight*eventWeight*factorSectorPiotr);
//        fillPair(hM, TString("hth1oAngleCorr%s")+suffix, cand1, cand2, theta1,oAngle,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2oAngleCorr%s")+suffix, cand1, cand2, theta2,oAngle,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth1massCorr%s")+suffix, cand1, cand2, theta1,invM/1000,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2massCorr%s")+suffix, cand1, cand2, theta2,invM/1000,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth1ptCorr%s")+suffix, cand1, cand2, theta1,pt,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2ptCorr%s")+suffix, cand1, cand2, theta2,pt,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth1yCorr%s")+suffix, cand1, cand2, theta1,rapidity,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hth2yCorr%s")+suffix, cand1, cand2, theta2,rapidity,weight*eventWeight*factorSectorPiotr,false,false);
        fillPair(hM, TString("hth1th2Corr%s")+suffix, cand1, cand2, theta1, theta2,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hp1th1Corr%s")+suffix, cand1, cand2, mom1, theta1,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM, TString("hp2th2Corr%s")+suffix, cand1, cand2, mom2, theta2,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM,TString("hpt1pt2Corr%s")+suffix, cand1, cand2, pt1, pt2,weight*eventWeight*factorSectorPiotr,false,false);
//        fillPair(hM,TString("hy1y2Corr%s")+suffix, cand1, cand2, y1, y2,weight*eventWeight*factorSectorPiotr);
    }

    fillPair(hM, TString("hmassCorr%s")+suffix,             cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hmassMCorr%s")+suffix,             cand1, cand2, invM/1000, invM*weight*eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hmassNoFactorsCorr%s")+suffix,    cand1, cand2, invM/1000, weight*eventWeight);
    fillPair(hM, TString("hmassPiotrFactorCorr%s")+suffix,  cand1, cand2, invM/1000, weight*eventWeight/factorPiotr);
    fillPair(hM, TString("hmassSectorFactorCorr%s")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSector);
    fillPair(hM, TString("hmassConstWidthCorr%s")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);

    fillPair(hM, TString("hoAngleCorr%s")+suffix,             cand1, cand2, oAngle, weight*eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hoAngleNoFactorsCorr%s")+suffix,    cand1, cand2, oAngle, weight*eventWeight);
    fillPair(hM, TString("hoAnglePiotrFactorCorr%s")+suffix,  cand1, cand2, oAngle, weight*eventWeight/factorPiotr);
    fillPair(hM, TString("hoAngleSectorFactorCorr%s")+suffix, cand1, cand2, oAngle, weight*eventWeight*factorSector);
    if (sectorCorr == kNoCorr) {
        fillPair(hM, TString("hmass6secCorr%s")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hoAngle6secCorr%s")+suffix, cand1, cand2, oAngle, weight*eventWeight*factorSectorPiotr);
    }
    if (sectorCorr == kCorr5to6) {
        fillPair(hM, TString("hmass5secCorr%s")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hoAngle5secCorr%s")+suffix, cand1, cand2, oAngle, weight*eventWeight*factorSectorPiotr);
    }
    if (sectorCorr == kCorr4to6d0) {
        fillPair(hM, TString("hmass4secd0Corr%s")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hoAngle4secd0Corr%s")+suffix, cand1, cand2, oAngle, weight*eventWeight*factorSectorPiotr);
    }
    if (sectorCorr == kCorr4to6d1) {
        fillPair(hM, TString("hmass4secd1Corr%s")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hoAngle4secd1Corr%s")+suffix, cand1, cand2, oAngle, weight*eventWeight*factorSectorPiotr);
    }
    if (sectorCorr == kCorr4to6d2) {
        fillPair(hM, TString("hmass4secd2Corr%s")+suffix, cand1, cand2, invM/1000, weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hoAngle4secd2Corr%s")+suffix, cand1, cand2, oAngle, weight*eventWeight*factorSectorPiotr);
    }
    fillPair(hM, TString("hyCorr%s")+suffix, cand1, cand2, rapidity,weight*eventWeight*factorSectorPiotr);
    fillPair(hM, TString("hptCorr%s")+suffix, cand1, cand2, pt,weight*eventWeight*factorSectorPiotr);
    if (invM < 150) {
        fillPair(hM,TString("hp1p2LowMCorr%s")+suffix, cand2, cand1, mom2, mom1,weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtMinusMinvLowMCorr%s")+suffix, cand1, cand2, mt-invM, weight*mt_w*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtLowMCorr%s")+suffix, cand1, cand2, mt, weight*mt_w*eventWeight*factorSectorPiotr);
    }
    else if (invM < 550) {
        fillPair(hM,TString("hp1p2MidMCorr%s")+suffix, cand2, cand1, mom2, mom1,weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtMinusMinvMidMCorr%s")+suffix, cand1, cand2, mt-invM, weight*mt_w*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtMidMCorr%s")+suffix, cand1, cand2, mt, weight*mt_w*eventWeight*factorSectorPiotr);
    }
    else if (invM < 1000) {
        fillPair(hM,TString("hp1p2HigMCorr%s")+suffix, cand2, cand1, mom2, mom1,weight*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtMinusMinvHigMCorr%s")+suffix, cand1, cand2, mt-invM, weight*mt_w*eventWeight*factorSectorPiotr);
        fillPair(hM, TString("hmtHigMCorr%s")+suffix, cand1, cand2, mt, weight*mt_w*eventWeight*factorSectorPiotr);
    }
    if (invM < 200) {
        fillPair(hM, TString("hpt_0_200_Corr%s")+suffix, cand1, cand2, pt,weight*eventWeight*factorSectorPiotr);
    } else if (invM < 400) {
        fillPair(hM, TString("hpt_200_400_Corr%s")+suffix, cand1, cand2, pt,weight*eventWeight*factorSectorPiotr);
    } else if (invM < 600) {
        fillPair(hM, TString("hpt_400_600_Corr%s")+suffix, cand1, cand2, pt,weight*eventWeight*factorSectorPiotr);
    }
}

inline void fillSingle(HHistMap &hM, HParticleCand *cand, const char *whatkind,int mlt = 0,Int_t nearestIndex = -1,Int_t w = 0, TH3F *p3DEff = NULL) {
    Float_t weight = 1.;
    if (p3DEff != NULL) {
        weight = getEfficiencyFactor(p3DEff, p3DEff, p3DEffEle_common[0], p3DEffPos_common[0], cand->getMomentum(), cand->getTheta(), cand->getPhi(), cand->getCharge());
    }
    Int_t sys = cand->getSystemUsed();
    Int_t sec = cand->getSector();
    Int_t ind = cand->getIndex();
    if (sys > -1) {
        if (w != NWEIGHTS) {
            if (cand->getCharge() == 1) {
                hM.get2(Form("ThetaMomEp%s_sys%i_sec",   whatkind,sys),cand->getSector(),w)->Fill(cand->getMomentum(),cand->getTheta(),weight);
            }
            else if (cand->getCharge() == -1) {
                hM.get2(Form("ThetaMomEm%s_sys%i_sec",   whatkind,sys),cand->getSector(),w)->Fill(cand->getMomentum(),cand->getTheta(),weight);
            }
            hM.get2(Form("MLPMom%s_sys%i",           whatkind,sys),w)->Fill(cand->getMomentum()*cand->getCharge(),MLPs[w][cand->getIndex()],weight);
            hM.get2(Form("BetaMom%s_sys%i",          whatkind,sys),w)->Fill(cand->getMomentum()*cand->getCharge(),cand->getBeta(),weight);
            hM.get2(Form("BetaMom%s_sys%i_sec",      whatkind,sys),cand->getSector(),w)->Fill(cand->getMomentum()*cand->getCharge(),cand->getBeta(),weight);
            hM.get2(Form("Mom%s_sys%i",              whatkind,sys),w)->Fill(cand->getMomentum()*cand->getCharge(),weight);
            hM.get2(Form("RichQaMom%s_sys%i",        whatkind,sys),w)->Fill(cand->getMomentum()*cand->getCharge(),cand->getRichMatchingQuality(),weight);
        }
        else {
            if (cand->getCharge() == 1) {
                hM.get2(Form("ThetaMomEp%s_sys%i_sec",   whatkind,sys),cand->getSector())->Fill(cand->getMomentum(),cand->getTheta(),weight);
            }
            else if (cand->getCharge() == -1) {
                hM.get2(Form("ThetaMomEm%s_sys%i_sec",   whatkind,sys),cand->getSector())->Fill(cand->getMomentum(),cand->getTheta(),weight);
            }
            hM.get2(Form("MLPMom%s_sys%i",           whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),MLPs[0][cand->getIndex()],weight);
            hM.get2(Form("BetaMom%s_sys%i",          whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getBeta(),weight);
            hM.get2(Form("BetaMom%s_sys%i_sec",      whatkind,sys),cand->getSector())->Fill(cand->getMomentum()*cand->getCharge(),cand->getBeta(),weight);
            hM.get2(Form("Mom%s_sys%i",              whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),weight);
            hM.get2(Form("RichQaMom%s_sys%i",        whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getRichMatchingQuality(),weight);
        }
        if (w == 0) {
            if (cand->getCharge() == -1) {
                hM.get2(Form("MetaQaMomSecEm%s_sys%i",whatkind,sys))->Fill(cand->getMomentum() + 1000*(sec%3), cand->getMetaMatchQuality() + 10*(sec<3),weight);
                hM.get2(Form("ThetaPhiEm%s_sys%i",   whatkind,sys))->Fill(cand->getPhi(), cand->getTheta(),weight);
                if (mlt > 0) {
                    hM.get2(Form("MetaQaMomSecEm%s_sys%i_mult%i",whatkind,sys,mlt))->Fill(cand->getMomentum() + 1000*(sec%3), cand->getMetaMatchQuality() + 10*(sec<3),weight);
                    hM.get2(Form("ThetaPhiEm%s_sys%i_mult%i",   whatkind,sys,mlt))->Fill(cand->getPhi(), cand->getTheta(),weight);
                }
            }
            if (cand->getCharge() == 1) {
                hM.get2(Form("MetaQaMomSecEp%s_sys%i",whatkind,sys))->Fill(cand->getMomentum() + 1000*(sec%3), cand->getMetaMatchQuality() + 10*(sec<3),weight);
                hM.get2(Form("ThetaPhiEp%s_sys%i",   whatkind,sys))->Fill(cand->getPhi(), cand->getTheta());
                if (mlt > 0) {
                    hM.get2(Form("MetaQaMomSecEp%s_sys%i_mult%i",whatkind,sys,mlt))->Fill(cand->getMomentum() + 1000*(sec%3), cand->getMetaMatchQuality() + 10*(sec<3),weight);
                    hM.get2(Form("ThetaPhiEp%s_sys%i_mult%i",   whatkind,sys,mlt))->Fill(cand->getPhi(), cand->getTheta(),weight);
                }
            }
            bool lowTh = cand->getTheta() < 45;
            bool secCenter = true;
            Float_t phi = cand->getPhi();
            if (phi < 15  && phi > 345) secCenter = false;
            if (phi > 45  && phi < 75)  secCenter = false;
            if (phi > 105 && phi < 135) secCenter = false;
            if (phi > 165 && phi < 195) secCenter = false;
            if (phi > 225 && phi < 255) secCenter = false;
            if (phi > 285 && phi < 315) secCenter = false;
            
            if (lowTh && secCenter)   hM.get2(Form("MDCdEdxMomLowThSide%s_sys%i", whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getMdcdEdx(),weight);
            if (!lowTh && secCenter)  hM.get2(Form("MDCdEdxMomHigThSide%s_sys%i", whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getMdcdEdx(),weight);
            if (lowTh && !secCenter)  hM.get2(Form("MDCdEdxMomLowThCent%s_sys%i", whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getMdcdEdx(),weight);
            if (!lowTh && !secCenter) hM.get2(Form("MDCdEdxMomHigThCent%s_sys%i", whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getMdcdEdx(),weight);


            hM.get2(Form("MDCdEdxMom%s_sys%i",       whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getMdcdEdx(),weight);
            hM.get2(Form("MassMom%s_sys%i",          whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getMass(),weight);
            hM.get2(Form("TOFdEdxMom%s_sys%i",       whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getTofdEdx(),weight);
            hM.get2(Form("hdeltaRichQaNormP%s_sys%i",whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getRichMatchingQualityNorm(),weight);
            hM.get2(Form("hdeltaRichQaP%s_sys%i",    whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getRichMatchingQuality(),weight);
            hM.get2(Form("hdeltaThetaP%s_sys%i",     whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getDeltaTheta(),weight);
            hM.get2(Form("hdeltaPhiP%s_sys%i",       whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getDeltaPhi(),weight);

            HParticleCand *nearestCand = NULL;

            nearestCand = HCategoryManager::getObject(nearestCand,catParticleCand,nearestIndex);
            Float_t angNear = angleNearestSeg[ind];
            if (nearestCand != NULL) {
                if (nearestCand->getMetaHitInd() == cand->getMetaHitInd()) {
                    hM.get2(Form("deltaMomMomNearestSeg%s_sys%i_sameMeta",whatkind,sys))->Fill(cand->getMomentum(),cand->getMomentum() - nearestCand->getMomentum(),weight);
                    hM.get2(Form("deltaMomAngleNearestSeg%s_sys%i_sameMeta",whatkind,sys))->Fill(cand->getMomentum() - nearestCand->getMomentum(),angNear,weight);
                    hM.get2(Form("angleNearestSeg%s_sys%i_sameMeta",       whatkind,sys))->Fill(angNear,weight);
                    hM.get2(Form("deltaMomNearestSeg%s_sys%i_sameMeta",    whatkind,sys))->Fill(cand->getMomentum() - nearestCand->getMomentum(),weight);
                }
                else {
                    hM.get2(Form("deltaMomMomNearestSeg%s_sys%i_notSameMeta",whatkind,sys))->Fill(cand->getMomentum(),angNear,weight);
                    hM.get2(Form("deltaMomAngleNearestSeg%s_sys%i_notSameMeta",whatkind,sys))->Fill(cand->getMomentum() - nearestCand->getMomentum(),angNear,weight);
                    hM.get2(Form("angleNearestSeg%s_sys%i_notSameMeta",    whatkind,sys))->Fill(angNear,weight);
                    hM.get2(Form("deltaMomNearestSeg%s_sys%i_notSameMeta", whatkind,sys))->Fill(cand->getMomentum() - nearestCand->getMomentum(),weight);
                }
            }

            vector<HParticleCand *> candsSameMeta;
            Int_t nSameMeta = 0;
            Int_t nSameMeta_angleLt4 = 0;
            Int_t nSameMeta_angleLt1 = 0;
            booker.getSameMeta(cand,candsSameMeta);
            for (size_t k = 0; k < candsSameMeta.size(); ++k) {
                if (cand->getChi2() != candsSameMeta[k]->getChi2()) {
                    Float_t angle = HParticleTool::getOpeningAngle(cand->getPhi(), cand->getTheta(), candsSameMeta[k]->getPhi(), candsSameMeta[k]->getTheta());
                    if (angle < 4) {
                        ++nSameMeta_angleLt4;
                    }
                    if (angle < 1) {
                        ++nSameMeta_angleLt1;
                    }
                    ++nSameMeta;
                    hM.get(Form("Angle2SameMeta%s",   whatkind))->Fill(angle,weight);
                }
            }
            hM.get(Form("nSameMeta%s_sys%i",        whatkind,sys))->Fill(nSameMeta,weight);
            hM.get(Form("nSameMetaAngleLt4%s_sys%i",whatkind,sys))->Fill(nSameMeta_angleLt4,weight);
            hM.get(Form("nSameMetaAngleLt1%s_sys%i",whatkind,sys))->Fill(nSameMeta_angleLt1,weight);
        }
        if (w == NWEIGHTS && !TString(whatkind).IsNull()) {
            bool lowTh = cand->getTheta() < 45;
            bool secCenter = true;
            Float_t phi = cand->getPhi();
            if (phi < 15  && phi > 345) secCenter = false;
            if (phi > 45  && phi < 75)  secCenter = false;
            if (phi > 105 && phi < 135) secCenter = false;
            if (phi > 165 && phi < 195) secCenter = false;
            if (phi > 225 && phi < 255) secCenter = false;
            if (phi > 285 && phi < 315) secCenter = false;

            if (lowTh && secCenter)   hM.get2(Form("MDCdEdxMomLowThSide%s_sys%i", whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getMdcdEdx(),weight);
            if (!lowTh && secCenter)  hM.get2(Form("MDCdEdxMomHigThSide%s_sys%i", whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getMdcdEdx(),weight);
            if (lowTh && !secCenter)  hM.get2(Form("MDCdEdxMomLowThCent%s_sys%i", whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getMdcdEdx(),weight);
            if (!lowTh && !secCenter) hM.get2(Form("MDCdEdxMomHigThCent%s_sys%i", whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getMdcdEdx(),weight);


            hM.get2(Form("MDCdEdxMom%s_sys%i",       whatkind,sys))->Fill(cand->getMomentum()*cand->getCharge(),cand->getMdcdEdx(),weight);
        }
    }
    if (w == 0) {
        hM.get(Form("TracksWithin4deg%s",whatkind))->Fill(countWithin[cand->getIndex()],weight);
        hM.get(Form("TracksUniqueWithin4deg%s",whatkind))->Fill(countUniqueWithin[cand->getIndex()],weight);
        hM.get(Form("TracksSameSegWithin4deg%s",whatkind))->Fill(countSameSegWithin[cand->getIndex()],weight);
        hM.get(Form("TracksSameSegSameOutWithin4deg%s",whatkind))->Fill(countSameSegSameOutWithin[cand->getIndex()],weight);
        hM.get(Form("TracksSameSegDiffOutWithin4deg%s",whatkind))->Fill(countSameSegDiffOutWithin[cand->getIndex()],weight);
        hM.get(Form("TracksLeptonsWithin4deg%s",whatkind))->Fill(countLeptonsWithin[cand->getIndex()],weight);
        hM.get(Form("TracksUniqueLeptonsWithin4deg%s",whatkind))->Fill(countUniqueLeptonsWithin[cand->getIndex()],weight);
        hM.get(Form("TracksSameSegLeptonsWithin4deg%s",whatkind))->Fill(countSameSegLeptonsWithin[cand->getIndex()],weight);
        hM.get(Form("TracksSameSegLeptonsLikeWithin4deg%s",whatkind))->Fill(countSameSegLeptonsLikeWithin[cand->getIndex()],weight);
        hM.get(Form("TracksSameSegLeptonsUnlikeWithin4deg%s",whatkind))->Fill(countSameSegLeptonsUnlikeWithin[cand->getIndex()],weight);
        hM.get2(Form("ThetaPhi%s",whatkind))->Fill(cand->getPhi(),cand->getTheta(),weight);
        hM.get2(Form("ThetaVz%s",whatkind))->Fill(gHades->getCurrentEvent()->getHeader()->getVertexReco().getZ(),cand->getTheta(),weight);
//        hM.get2(Form("hringAmpTheta%s",whatkind))->Fill(cand->getTheta(),cand->getRingAmplitude(),weight);
        hM.get(Form("hringNP%s",whatkind),cand->getSector())->Fill(cand->getRingNumPads());
        hM.get(Form("hringAC%s",whatkind),cand->getSector())->Fill(cand->getRingAmplitude()/cand->getRingNumPads());
        hM.get(Form("hringRC%s",whatkind),cand->getSector())->Fill(cand->getRingCentroid());
        hM.get(Form("hringPM%s",whatkind),cand->getSector())->Fill(cand->getRingPatternMatrix());
        hM.get(Form("hringHT%s",whatkind),cand->getSector())->Fill(cand->getRingHouTra());
    }
}

void addSingleHists(HHistMap &hM, const char *whatkind, bool hc = kFALSE) {
    for (int sys = 0; sys <= 1; ++sys) {
        for (int mlt = 1; mlt <= 5; ++mlt) {
            hM.addHist("TH2F",Form("MetaQaMomSecEm%s_sys%i_mult%i",whatkind,sys,mlt),   Form("MetaQaMomSysEm%s_sys%i_mult%i",whatkind,sys,mlt),   600,    0,3000,200,0,20);
            hM.addHist("TH2F",Form("MetaQaMomSecEp%s_sys%i_mult%i",whatkind,sys,mlt),   Form("MetaQaMomSysEp%s_sys%i_mult%i",whatkind,sys,mlt),   600,    0,3000,200,0,20);
            hM.addHist("TH2F",Form("ThetaPhiEm%s_sys%i_mult%i",whatkind,sys,mlt),       Form("ThetaPhiEm%s_sys%i_mult%i",whatkind,sys,mlt),    360,    0, 360, 90,0,90);
            hM.addHist("TH2F",Form("ThetaPhiEp%s_sys%i_mult%i",whatkind,sys,mlt),       Form("ThetaPhiEp%s_sys%i_mult%i",whatkind,sys,mlt),    360,    0, 360, 90,0,90);
        }
        if (!hc) {
            hM.addHistArray("TH2F",Form("MLPMom%s_sys%i",whatkind,sys),      Form("MLPMom%s_sys%i_w%%i",whatkind,sys),         "",200,-1000,1000,200,0,1,0,0,0,"","","","",NWEIGHTS);
            hM.addHistArray("TH2F",Form("BetaMom%s_sys%i",whatkind,sys),     Form("BetaMom%s_sys%i_w%%i",whatkind,sys),        "",200,-1000,1000,600,0.,1.2,0,0,0,"","","","",NWEIGHTS);
            hM.addHistArray("TH2F",Form("RichQaMom%s_sys%i",whatkind,sys),   Form("RichQaMom%s_sys%i_w%%i",whatkind,sys),      "",200,-1000,1000,200,0.,10,0,0,0,"","","","",NWEIGHTS);
            hM.addHistArray("TH1F",Form("Mom%s_sys%i",whatkind,sys),         Form("Mom%s_sys%i_w%%i",whatkind,sys),            "",200,-1000,1000,0,0,0,0,0,0,"","","","",NWEIGHTS);
            hM.addHistArray("TH2F",Form("BetaMom%s_sys%i_sec",whatkind,sys), Form("BetaMom%s_sys%i_sec%%i_w%%i",whatkind,sys), "",200,-1000,1000,100,0.,1.2,0,0,0,"","","","",6,NWEIGHTS);
            hM.addHistArray("TH2F",Form("ThetaMomEp%s_sys%i_sec",whatkind,sys), Form("ThetaMomEp%s_sys%i_sec%%i_w%%i",whatkind,sys), "",200,0,1000,45,0.,90,0,0,0,"","","","",6,NWEIGHTS);
            hM.addHistArray("TH2F",Form("ThetaMomEm%s_sys%i_sec",whatkind,sys), Form("ThetaMomEm%s_sys%i_sec%%i_w%%i",whatkind,sys), "",200,0,1000,45,0.,90,0,0,0,"","","","",6,NWEIGHTS);
        }
        else {
            hM.addHist("TH2F",Form("MLPMom%s_sys%i",whatkind,sys),       "",200,-1000,1000,200,0,1);
            hM.addHist("TH2F",Form("BetaMom%s_sys%i",whatkind,sys),      "",200,-1000,1000,600,0.,1.2);
            hM.addHist("TH1F",Form("Mom%s_sys%i",whatkind,sys),          "",200,-1000,1000,0,0,0);
            hM.addHist("TH2F",Form("RichQaMom%s_sys%i",whatkind,sys),    "",200,-1000,1000,200,0,10);
            hM.addHistArray("TH2F",Form("BetaMom%s_sys%i_sec",whatkind,sys),      Form("BetaMom%s_sys%i_sec%%i",whatkind,sys),    "",200,-1000,1000,100,0.,1.2,0,0,0,"","","","",6);
            hM.addHistArray("TH2F",Form("ThetaMomEp%s_sys%i_sec",whatkind,sys),   Form("ThetaMomEp%s_sys%i_sec%%i",whatkind,sys), "",200,0,1000,45,0.,90,0,0,0,"","","","",6);
            hM.addHistArray("TH2F",Form("ThetaMomEm%s_sys%i_sec",whatkind,sys),   Form("ThetaMomEm%s_sys%i_sec%%i",whatkind,sys), "",200,0,1000,45,0.,90,0,0,0,"","","","",6);
        }

        hM.addHist("TH2F",Form("MetaQaMomSecEm%s_sys%i",whatkind,sys),      Form("MetaQaMomSysEm%s_sys%i",whatkind,sys),   600,    0,3000,200,0,20);
        hM.addHist("TH2F",Form("MetaQaMomSecEp%s_sys%i",whatkind,sys),      Form("MetaQaMomSysEp%s_sys%i",whatkind,sys),   600,    0,3000,200,0,20);
        hM.addHist("TH2F",Form("ThetaPhiEm%s_sys%i",whatkind,sys),          Form("ThetaPhiEm%s_sys%i",whatkind,sys),    360,    0, 360, 90,0, 90);
        hM.addHist("TH2F",Form("ThetaPhiEp%s_sys%i",whatkind,sys),          Form("ThetaPhiEp%s_sys%i",whatkind,sys),    360,    0, 360, 90,0, 90);
        hM.addHist("TH2F",Form("MassMom%s_sys%i",whatkind,sys),             Form("MassMom%s_sys%i",whatkind,sys),          200,-1000,1000,400,-200,1800);
        hM.addHist("TH2F",Form("MDCdEdxMom%s_sys%i",whatkind,sys),          Form("MDCdEdxMom%s_sys%i",whatkind,sys),       200,-1000,1000,200,0,10);
        hM.addHist("TH2F",Form("MDCdEdxMomLowThCent%s_sys%i",whatkind,sys), Form("MDCdEdxMomLowThCent%s_sys%i",whatkind,sys),       200,-1000,1000,200,0,10);
        hM.addHist("TH2F",Form("MDCdEdxMomLowThSide%s_sys%i",whatkind,sys), Form("MDCdEdxMomLowThSide%s_sys%i",whatkind,sys),       200,-1000,1000,200,0,10);
        hM.addHist("TH2F",Form("MDCdEdxMomHigThCent%s_sys%i",whatkind,sys), Form("MDCdEdxMomHigThCent%s_sys%i",whatkind,sys),       200,-1000,1000,200,0,10);
        hM.addHist("TH2F",Form("MDCdEdxMomHigThSide%s_sys%i",whatkind,sys), Form("MDCdEdxMomHigThSide%s_sys%i",whatkind,sys),       200,-1000,1000,200,0,10);
        hM.addHist("TH2F",Form("TOFdEdxMom%s_sys%i",whatkind,sys),          Form("TOFdEdxMom%s_sys%i",whatkind,sys),       200,-1000,1000,200,0,10);
        hM.addHist("TH2F",Form("hdeltaThetaP%s_sys%i",whatkind,sys),        Form("hdeltaThetaP%s_sys%i",whatkind,sys),      60,-1000,1000,100,-10,10);
        hM.addHist("TH2F",Form("hdeltaPhiP%s_sys%i",whatkind,sys),          Form("hdeltaPhiP%s_sys%i",whatkind,sys),        60,-1000,1000,100,-10,10);
        hM.addHist("TH2F",Form("hdeltaRichQaP%s_sys%i",whatkind,sys),       Form("hdeltaRichQaP%s_sys%i",whatkind,sys),    200,-1000,1000,100,0,10);
        hM.addHist("TH2F",Form("hdeltaRichQaNormP%s_sys%i",whatkind,sys),   Form("hdeltaRichQaNormP%s_sys%i",whatkind,sys),200,-1000,1000,100,0,10);

        hM.addHist("TH2F",Form("deltaMomMomNearestSeg%s_sys%i_sameMeta",     whatkind,sys), Form("deltaMomMomNearestSeg%s_sys%i_sameMeta",      whatkind,sys), 200,  0,1000, 200,-10,10);
        hM.addHist("TH2F",Form("deltaMomMomNearestSeg%s_sys%i_notSameMeta",  whatkind,sys), Form("deltaMomMomNearestSeg%s_sys%i_notSameMeta"  , whatkind,sys), 200,  0,1000, 200,-10,10);
        hM.addHist("TH2F",Form("deltaMomAngleNearestSeg%s_sys%i_sameMeta",   whatkind,sys), Form("deltaMomAngleNearestSeg%s_sys%i_sameMeta",    whatkind,sys), 200,-10,10, 200,-10,10);
        hM.addHist("TH2F",Form("deltaMomAngleNearestSeg%s_sys%i_notSameMeta",whatkind,sys), Form("deltaMomAngleNearestSeg%s_sys%i_notSameMeta", whatkind,sys), 200,-10,10, 200,-10,10);
        hM.addHist("TH1F",Form("angleNearestSeg%s_sys%i_sameMeta",           whatkind,sys), Form("angleNearestSeg%s_sys%i_sameMeta",            whatkind,sys), 200,-10,10);
        hM.addHist("TH1F",Form("angleNearestSeg%s_sys%i_notSameMeta",        whatkind,sys), Form("angleNearestSeg%s_sys%i_notSameMeta"   ,      whatkind,sys), 200,-10,10);
        hM.addHist("TH1F",Form("deltaMomNearestSeg%s_sys%i_sameMeta",        whatkind,sys), Form("deltaMomNearestSeg%s_sys%i_sameMeta",         whatkind,sys), 200,-10,10);
        hM.addHist("TH1F",Form("deltaMomNearestSeg%s_sys%i_notSameMeta",     whatkind,sys), Form("deltaMomNearestSeg%s_sys%i_notSameMeta",      whatkind,sys), 200,-10,10);

        hM.addHist("TH1F",Form("nSameMeta%s_sys%i",                          whatkind,sys), Form("nSameMeta%s_sys%i",                           whatkind,sys), 10,0,10);
        hM.addHist("TH1F",Form("nSameMetaAngleLt4%s_sys%i",                  whatkind,sys), Form("nSameMetaAngleLt4%s_sys%i",                   whatkind,sys), 10,0,10);
        hM.addHist("TH1F",Form("nSameMetaAngleLt1%s_sys%i",                  whatkind,sys), Form("nSameMetaAngleLt1%s_sys%i",                   whatkind,sys), 10,0,10);
    }
    hM.addHist("TH2F",Form("ThetaPhi%s",      whatkind),   Form("ThetaPhi%s",      whatkind),  180,    0, 360, 90,  0,90);
    hM.addHist("TH2F",Form("ThetaVz%s",       whatkind),   Form("ThetaVz%s",       whatkind),  200,  -80,  20, 90,  0,90);
    hM.addHist("TH1F",Form("TracksWithin4deg%s",whatkind),   Form("TracksWithin4deg%s",whatkind),  40,  0,  40);
    hM.addHist("TH1F",Form("TracksUniqueWithin4deg%s",whatkind),   Form("TracksUniqueWithin4deg%s",whatkind),  40,  0,  40);
    hM.addHist("TH1F",Form("TracksSameSegWithin4deg%s",whatkind),   Form("TracksSameSegWithin4deg%s",whatkind),  40,  0,  40);
    hM.addHist("TH1F",Form("TracksSameSegSameOutWithin4deg%s",whatkind),   Form("TracksSameSegSameOutWithin4deg%s",whatkind),  40,  0,  40);
    hM.addHist("TH1F",Form("TracksSameSegDiffOutWithin4deg%s",whatkind),   Form("TracksSameSegDiffOutWithin4deg%s",whatkind),  40,  0,  40);
    hM.addHist("TH1F",Form("TracksLeptonsWithin4deg%s",whatkind),   Form("TracksLeptonsWithin4deg%s",whatkind),  40,  0,  40);
    hM.addHist("TH1F",Form("TracksUniqueLeptonsWithin4deg%s",whatkind),   Form("TracksUniqueLeptonsWithin4deg%s",whatkind),  40,  0,  40);
    hM.addHist("TH1F",Form("TracksSameSegLeptonsWithin4deg%s",whatkind),   Form("TracksSameSegLeptonsWithin4deg%s",whatkind),  40,  0,  40);
    hM.addHist("TH1F",Form("TracksSameSegLeptonsLikeWithin4deg%s",whatkind),   Form("TracksSameSegLeptonsLikeWithin4deg%s",whatkind),  40,  0,  40);
    hM.addHist("TH1F",Form("TracksSameSegLeptonsUnlikeWithin4deg%s",whatkind),   Form("TracksSameSegLeptonsUnlikeWithin4deg%s",whatkind),  40,  0,  40);
    hM.addHist("TH1F",Form("Angle2SameMeta%s",whatkind),   Form("Angle2SameMeta%s",whatkind),  100,  0,  10);
//    hM.addHist("TH2F",Form("hringAmpTheta%s",whatkind),"Number of pads vs. #theta %i",90,0,90,1000,0,10000);
    hM.addHistArray("TH1F",Form("hringNP%s",whatkind),Form("ringNP%s_sec%%i",whatkind),"Number of pads, sector %i",50,0,50,0,0,0,0,0,0,"","","","",6);
    hM.addHistArray("TH1F",Form("hringAC%s",whatkind),Form("ringAC%s_sec%%i",whatkind),"Average charge, sector %i",100,0,400,0,0,0,0,0,0,"","","","",6);
    hM.addHistArray("TH1F",Form("hringRC%s",whatkind),Form("ringRC%s_sec%%i",whatkind),"Ring centroid, sector %i",100,0,5,0,0,0,0,0,0,"","","","",6);
    hM.addHistArray("TH1F",Form("hringPM%s",whatkind),Form("ringPM%s_sec%%i",whatkind),"Pattern matrix, sector %i",100,0,1800,0,0,0,0,0,0,"","","","",6);
    hM.addHistArray("TH1F",Form("hringHT%s",whatkind),Form("ringHT%s_sec%%i",whatkind),"Hough transform, sector %i",100,0,4000,0,0,0,0,0,0,"","","","",6);
}

Float_t myDotProduct(HParticleCand *cand1, HParticleCand *cand2) {
    static const Float_t degToRad = 1.74532925199432955e-02;
    Float_t phi1 = cand1->getPhi()*degToRad;
    Float_t phi2 = cand2->getPhi()*degToRad;
    Float_t theta1 = cand1->getTheta()*degToRad;
    Float_t theta2 = cand2->getTheta()*degToRad;
    Float_t x1 = sin(theta1)*cos(phi1);
    Float_t y1 = sin(theta1)*sin(phi1);
    Float_t z1 = cos(theta1);
    Float_t x2 = sin(theta2)*cos(phi2);
    Float_t y2 = sin(theta2)*sin(phi2);
    Float_t z2 = cos(theta2);
    Float_t dot2 = x1*x2+y1*y2+z1*z2;
    return dot2;

}

Float_t myOpeningAngle(HParticleCand *cand1, HParticleCand *cand2) {
    Float_t dot = myDotProduct(cand1,cand2);
    return acos(dot)*TMath::RadToDeg();    
}

static Bool_t selectClosePairCands(HParticleCand* cand){
    const static Float_t betaThr[2] ={0.94F, 0.93F};
    //  selection function for lepton candidates.
    Int_t system = cand->getSystemUsed();
    Float_t beta = cand->getBeta();
    if( beta < betaThr[system] || beta > 1.15F) return kFALSE;
    Float_t metaThr = 2.2F;
    Float_t momThr = 540.0F;
    Float_t betaFac = 1.06F;
    if(system == 0) {
	metaThr= 2.0F;
	momThr = 600.0F;
	betaFac = 1.05F;
    }
    if(cand->getMetaMatchQuality() > metaThr) return kFALSE;
    Float_t mom =  cand->getMomentum();
    // meta match Qa cut && time-of-flight cut
    if(mom > momThr) return kFALSE;
    Float_t momNorm = mom/139.6F;
    Float_t betaPion = sqrt(momNorm*momNorm/(1.0F + momNorm*momNorm));
    if(beta > betaFac*betaPion) return kTRUE;
    return kFALSE;
}

Float_t smear(Float_t value, TH2F *matrix, TRandom &random) {
    Int_t bin1 = matrix->GetYaxis()->FindBin(value-50);
    Int_t bin2 = matrix->GetYaxis()->FindBin(value+50);
    TH1D *px = matrix->ProjectionX("_px",bin1,bin2);
    value += px->GetRandom();
    delete px;
    return value;
}
