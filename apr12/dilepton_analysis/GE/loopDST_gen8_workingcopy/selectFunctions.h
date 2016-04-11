//------------------------------------------------------------
//-- NEEDED GLOBAL POINTERS ----------------------------------
//------------------------------------------------------------
static const int NWEIGHTS = 13;
static TH1F *mlplocmin_sys0[NWEIGHTS];
static TH1F *mlplocmin_sys1[NWEIGHTS];
static TH1F *hcutRichQaSys0 = NULL;
static TH1F *hcutRichQaSys1 = NULL;
static TH1F *hcutMassSys0 = NULL;
static TH1F *hcutMassSys1 = NULL;
static map<Int_t, Float_t> MLPs[NWEIGHTS];
static map<Int_t, Float_t> angleNearestBoth;
static map<Int_t, Float_t> angleNearestFull;
static map<Int_t, Float_t> angleNearestSeg;
static map<Int_t, Int_t>   nearestIndexSeg;
static TMVA::Reader *reader[NWEIGHTS][2];

template <typename T>
Int_t findNearestNeighbor(HParticleCand *cand1, T &angleNearestBoth, T &angleNearestSeg, T &angleNearestFull, HCategory *candCat);

Bool_t isGoodMLP(HParticleCand *cand, Float_t MLP, Int_t w = 0) {
    TH1F *hcutg;
    if (cand->getSystemUsed() == 0) {
        hcutg = mlplocmin_sys0[w];
    }
    else {
        hcutg = mlplocmin_sys1[w];
    }
    Float_t x = cand->getMomentum() * cand->getCharge();
    Float_t minMLP;
    if (x < hcutg->GetBinCenter(1)) {
        minMLP = hcutg->GetBinContent(1);
    }
    else if (x < hcutg->GetBinCenter(hcutg->GetNbinsX())) {
        minMLP = hcutg->Interpolate(x);
    }
    else {
        minMLP = hcutg->GetBinContent(hcutg->GetNbinsX());
    }
    //minMLP = 0.6;
    return MLP > minMLP;
}

Bool_t selectStrict(HParticleCand *cand1) {
    if(cand1->getSystemUsed()==-1) return kFALSE;
    if(!cand1->isFlagBit(Particle::kIsLepton)) return kFALSE;
    HCategory* candCat = (HCategory*)HCategoryManager::getCategory(catParticleCand);
    Int_t nearestIndex = findNearestNeighbor(cand1,angleNearestBoth,angleNearestSeg,angleNearestFull,candCat);


    if (nearestIndex > -1) {
        HParticleCand *nearestCand = NULL;
        nearestCand = HCategoryManager::getObject(nearestCand,candCat,nearestIndex);
        Float_t nearestAngle = HParticleTool::getOpeningAngle(cand1->getPhi(),cand1->getTheta(),nearestCand->getPhi(),nearestCand->getTheta());
        if (cand1->getMetaHitInd() == nearestCand->getMetaHitInd() && nearestCand->getChi2() > 40) {
            return kFALSE;
        }
        if (cand1->getMetaHitInd() == nearestCand->getMetaHitInd() && nearestCand->getChi2() < 0) {
            return kFALSE;
        }
        if (nearestAngle > 0 && nearestAngle < 6) {
            return kFALSE;
        }
    }


    Float_t mom  = cand1->getMomentum();
    Float_t mlp  = MLPs[0][cand1->getIndex()];
    Int_t charge = cand1->getCharge();
    Int_t sys    = cand1->getSystemUsed();
    Int_t richQa = cand1->getRichMatchingQuality();

    bool cutrich;
    if (sys == 0) {
        cutrich = richQa < hcutRichQaSys0->Interpolate(mom*charge);
    }
    else {
        cutrich = richQa < hcutRichQaSys1->Interpolate(mom*charge);
    }
//    if (cand1->getTheta() < 25 || cand1->getTheta() > 75) {
//        return kFALSE;
//    }
//    if (charge > 0 && isGoodMLP(cand1,mlp) && cutrich) {
    if (isGoodMLP(cand1,mlp) && richQa < 2) {
            return kTRUE;
    }
    return kFALSE;

}

Bool_t selectPosStrict(HParticleCand *cand1) {
    return (cand1->getCharge() > 0 && selectStrict(cand1));
}

Bool_t selectNegStrict(HParticleCand *cand1) {
    return (cand1->getCharge() < 0 && selectStrict(cand1));
}

Bool_t selectPosStrict5degBoth(HParticleCand *cand1) {
    bool angleNearestCond = angleNearestBoth[cand1->getIndex()] < 0 || angleNearestBoth[cand1->getIndex()] > 5;
    return selectPosStrict(cand1) && angleNearestCond;
}

Bool_t selectNegStrict5degBoth(HParticleCand *cand1) {
    bool angleNearestCond = angleNearestBoth[cand1->getIndex()] < 0 || angleNearestBoth[cand1->getIndex()] > 5;
    return selectNegStrict(cand1) && angleNearestCond;
}

Bool_t selectPosStrict5degAbsBoth(HParticleCand *cand1) {
    bool angleNearestCond = fabs(angleNearestBoth[cand1->getIndex()]) > 5;
    return selectPosStrict(cand1) && angleNearestCond;
}

Bool_t selectNegStrict5degAbsBoth(HParticleCand *cand1) {
    bool angleNearestCond = fabs(angleNearestBoth[cand1->getIndex()]) > 5;
    return selectNegStrict(cand1) && angleNearestCond;
}

Bool_t selectPosStrict7degBoth(HParticleCand *cand1) {
    bool angleNearestCond = angleNearestBoth[cand1->getIndex()] < 0 || angleNearestBoth[cand1->getIndex()] > 7;
    return selectPosStrict(cand1) && angleNearestCond;
}

Bool_t selectNegStrict7degBoth(HParticleCand *cand1) {
    bool angleNearestCond = angleNearestBoth[cand1->getIndex()] < 0 || angleNearestBoth[cand1->getIndex()] > 7;
    return selectNegStrict(cand1) && angleNearestCond;
}

Bool_t selectPosStrict7degAbsBoth(HParticleCand *cand1) {
    bool angleNearestCond = fabs(angleNearestBoth[cand1->getIndex()]) > 7;
    return selectPosStrict(cand1) && angleNearestCond;
}

Bool_t selectNegStrict7degAbsBoth(HParticleCand *cand1) {
    bool angleNearestCond = fabs(angleNearestBoth[cand1->getIndex()]) > 7;
    return selectNegStrict(cand1) && angleNearestCond;
}

Bool_t selectPosStrict5degFull(HParticleCand *cand1) {
    bool angleNearestCond = angleNearestFull[cand1->getIndex()] < 0 || angleNearestFull[cand1->getIndex()] > 5;
    return selectPosStrict(cand1) && angleNearestCond;
}

Bool_t selectNegStrict5degFull(HParticleCand *cand1) {
    bool angleNearestCond = angleNearestFull[cand1->getIndex()] < 0 || angleNearestFull[cand1->getIndex()] > 5;
    return selectNegStrict(cand1) && angleNearestCond;
}

Bool_t selectPosStrict5degAbsFull(HParticleCand *cand1) {
    bool angleNearestCond = fabs(angleNearestFull[cand1->getIndex()]) > 5;
    return selectPosStrict(cand1) && angleNearestCond;
}

Bool_t selectNegStrict5degAbsFull(HParticleCand *cand1) {
    bool angleNearestCond = fabs(angleNearestFull[cand1->getIndex()]) > 5;
    return selectNegStrict(cand1) && angleNearestCond;
}

Bool_t selectPosStrict7degFull(HParticleCand *cand1) {
    bool angleNearestCond = angleNearestFull[cand1->getIndex()] < 0 || angleNearestFull[cand1->getIndex()] > 7;
    return selectPosStrict(cand1) && angleNearestCond;
}

Bool_t selectNegStrict7degFull(HParticleCand *cand1) {
    bool angleNearestCond = angleNearestFull[cand1->getIndex()] < 0 || angleNearestFull[cand1->getIndex()] > 7;
    return selectNegStrict(cand1) && angleNearestCond;
}

Bool_t selectPosStrict7degAbsFull(HParticleCand *cand1) {
    bool angleNearestCond = fabs(angleNearestFull[cand1->getIndex()]) > 7;
    return selectPosStrict(cand1) && angleNearestCond;
}

Bool_t selectNegStrict7degAbsFull(HParticleCand *cand1) {
    bool angleNearestCond = fabs(angleNearestFull[cand1->getIndex()]) > 7;
    return selectNegStrict(cand1) && angleNearestCond;
}

Bool_t selectPosStrict5degSeg(HParticleCand *cand1) {
    bool angleNearestCond = angleNearestSeg[cand1->getIndex()] < 0 || angleNearestSeg[cand1->getIndex()] > 5;
    return selectPosStrict(cand1) && angleNearestCond;
}

Bool_t selectNegStrict5degSeg(HParticleCand *cand1) {
    bool angleNearestCond = angleNearestSeg[cand1->getIndex()] < 0 || angleNearestSeg[cand1->getIndex()] > 5;
    return selectNegStrict(cand1) && angleNearestCond;
}

Bool_t selectPosStrict5degAbsSeg(HParticleCand *cand1) {
    bool angleNearestCond = fabs(angleNearestSeg[cand1->getIndex()]) > 5;
    return selectPosStrict(cand1) && angleNearestCond;
}

Bool_t selectNegStrict5degAbsSeg(HParticleCand *cand1) {
    bool angleNearestCond = fabs(angleNearestSeg[cand1->getIndex()]) > 5;
    return selectNegStrict(cand1) && angleNearestCond;
}

Bool_t selectPosStrict7degSeg(HParticleCand *cand1) {
    bool angleNearestCond = angleNearestSeg[cand1->getIndex()] < 0 || angleNearestSeg[cand1->getIndex()] > 7;
    return selectPosStrict(cand1) && angleNearestCond;
}

Bool_t selectNegStrict7degSeg(HParticleCand *cand1) {
    bool angleNearestCond = angleNearestSeg[cand1->getIndex()] < 0 || angleNearestSeg[cand1->getIndex()] > 7;
    return selectNegStrict(cand1) && angleNearestCond;
}

Bool_t selectPosStrict7degAbsSeg(HParticleCand *cand1) {
    bool angleNearestCond = fabs(angleNearestSeg[cand1->getIndex()]) > 7;
    return selectPosStrict(cand1) && angleNearestCond;
}

Bool_t selectNegStrict7degAbsSeg(HParticleCand *cand1) {
    bool angleNearestCond = fabs(angleNearestSeg[cand1->getIndex()]) > 7;
    return selectNegStrict(cand1) && angleNearestCond;
}

Bool_t isGoodShower(HParticleCand* cand,TF1* showerF) {

    // return kTRUE if good shower or no shower at all
    // or momentum < 400
    if (cand->getSystemUsed() == 1) return kTRUE;

    Double_t a = showerF->GetParameter(0);
    Double_t b = showerF->GetParameter(1);
    Double_t c = showerF->GetParameter(2);

    Double_t fSumDiff = cand->getShowerSum1() + cand->getShowerSum2() - cand->getShowerSum0();
    Double_t fMom     = cand->getMomentum();
    Double_t fMom0    = -a/(1000+c) + b;
    Double_t value    = showerF->Eval(fMom);
    Double_t thr;
    if (showerF->GetExpFormula().Contains("([0]/(x-[1]))+[2]"))
        thr = fMom > fMom0 ? value : -1000;
    else
        thr = value;
    return  fSumDiff > thr;
}
