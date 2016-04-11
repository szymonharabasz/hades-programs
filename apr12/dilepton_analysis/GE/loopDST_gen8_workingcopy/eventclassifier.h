Int_t getMultBin() {
    HParticleEvtInfo* evtinfo = NULL;
    evtinfo = HCategoryManager::getObject(evtinfo,catParticleEvtInfo,0);
    Int_t mult_meta = evtinfo->getSumTofMultCut() + evtinfo->getSumRpcMultHitCut();
    Int_t mult_bin = 5;
    if (mult_meta >  60 && mult_meta <=  88) mult_bin = 4; // most peripheral
    if (mult_meta >  88 && mult_meta <= 121) mult_bin = 3;
    if (mult_meta > 121 && mult_meta <= 160) mult_bin = 2;
    if (mult_meta > 160 && mult_meta <= 250) mult_bin = 1; // most central
    if (mult_meta > 250) mult_bin = 0;

    return mult_bin;
}

Int_t getMultBinCand() {
    HParticleEvtInfo* evtinfo = NULL;
    evtinfo = HCategoryManager::getObject(evtinfo,catParticleEvtInfo,0);
    Int_t mult_sel = 0;
    if (evtinfo) {
        mult_sel += evtinfo->getSelectedParticleCandMult(0);
        mult_sel += evtinfo->getSelectedParticleCandMult(1);
        mult_sel += evtinfo->getSelectedParticleCandMult(4);
        mult_sel += evtinfo->getSelectedParticleCandMult(5);
    }
    mult_sel *= 1.5;
    Int_t mult_bin = 5;
    if (mult_sel >= 17 && mult_sel <  28) mult_bin = 4;
    if (mult_sel >= 28 && mult_sel <  44) mult_bin = 3;
    if (mult_sel >= 44 && mult_sel <  68) mult_bin = 2;
    if (mult_sel >= 68 && mult_sel < 140) mult_bin = 1;
    if (mult_sel >= 140) mult_bin = 0;
    return mult_bin;
}

Int_t getTargetBin() {
    Float_t vz = gHades->getCurrentEvent()->getHeader()->getVertex().getZ();
    Int_t target_bin;
    if (vz < -52.0) target_bin = 0;
    else if (vz < -49.0) target_bin = 1;
    else if (vz < -45.0) target_bin = 2;
    else if (vz < -41.0) target_bin = 3;
    else if (vz < -37.0) target_bin = 4;
    else if (vz < -34.0) target_bin = 5;
    else if (vz < -31.0) target_bin = 6;
    else if (vz < -28.0) target_bin = 7;
    else if (vz < -23.0) target_bin = 8;
    else if (vz < -20.0) target_bin = 9;
    else if (vz < -16.0) target_bin = 10;
    else if (vz < -12.0) target_bin = 11;
    else if (vz <  -9.0) target_bin = 12;
    else if (vz <  -5.0) target_bin = 13;
    else target_bin = 14;
    return target_bin;
}

Int_t getRPhiBin() {
    HCategory* fParticleEvtInfoCat =  (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo,kTRUE,"catParticleEvtInfo"); 
    HParticleEvtInfo *event_info = (HParticleEvtInfo*)fParticleEvtInfoCat->getObject(0); 
    Float_t rphi = event_info->getRPlanePhi();
    Int_t rphi_bin;
    if (rphi < -200) rphi_bin = 1;
    else if (rphi < -150) rphi_bin = 2;
    else if (rphi < -120) rphi_bin = 3;
    else if (rphi <  -90) rphi_bin = 4;
    else if (rphi <  -60) rphi_bin = 5;
    else if (rphi <  -30) rphi_bin = 6;
    else if (rphi <    0) rphi_bin = 7;
    else if (rphi <   30) rphi_bin = 8;
    else if (rphi <   60) rphi_bin = 9;
    else if (rphi <   90) rphi_bin = 10;
    else if (rphi <  120) rphi_bin = 11;
    else if (rphi <  150) rphi_bin = 12;
    else rphi_bin = 13;
    return rphi_bin;
}

Int_t eventClassifierMultRPhiTargetFewer() {
    Int_t mult_bin   = getMultBin();
    Int_t target_bin = getTargetBin()/5;
    Int_t rphi_bin   = getRPhiBin()/4;

    return 300*mult_bin + 16*rphi_bin + target_bin;
}

Int_t eventClassifierMultRPhiTarget() {
    Int_t mult_bin   = getMultBin();
    Int_t target_bin = getTargetBin();
    Int_t rphi_bin   = getRPhiBin();

    return 300*mult_bin + 16*rphi_bin + target_bin;
}

Int_t eventClassifierMultTarget() {
    Int_t mult_bin   = getMultBin();
    Int_t target_bin = getTargetBin();

    return 300*mult_bin + target_bin;
}

Int_t eventClassifierMultRPhi() {
    Int_t mult_bin   = getMultBin();
    Int_t rphi_bin   = getRPhiBin();

    return 300*mult_bin + 16*rphi_bin;
}

Int_t eventClassifierMult() {
    Int_t mult_bin   = getMultBin();

    return 300*mult_bin;
}



