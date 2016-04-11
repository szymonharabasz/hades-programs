void addKineHists(HHistMap &hM) {
    Int_t binsxmom  = 100;
    Int_t binsxmass = 100;
    Float_t maxmass = 1000;

    hM.addHistArray("TH2F","nhitCellTof","nhitCellTof_%i","nhitCellTof_%i",10,0,10,64,0,64,0,0,0,"","","","",10 );
    hM.addHistArray("TH2F","nhitCellRpc","nhitCellRpc_%i","nhitCellRpc_%i",10,0,10,32,0,32,0,0,0,"","","","",10 );
    hM.addHist("TH1F","hnKine","hnKine",10,0,10);
    hM.addHist("TH2F","h_acc_PT","h_acc_PT",180,0,360,45,0,90);
    hM.addHist("TH2F","h_acc_MT","h_acc_MT",100,0,1000,45,0,90);
    hM.addHist("TH2F","h_acc_PM","h_acc_PM",180,0,360,100,0,1000);
    hM.addHistArray("TH3F","h_acc_MPT_ep","h_acc_MPT_ep_%i","h_acc_MPT_ep_%i",180,0,360,45,0,90,100,0,1000,"","","","",7);
    hM.addHistArray("TH3F","h_acc_MPT_em","h_acc_MPT_em_%i","h_acc_MPT_em_%i",180,0,360,45,0,90,100,0,1000,"","","","",7);

    hM.addHist("TH2F","h_acc_PT_oa","h_acc_PT_oa",180,0,360,45,0,90);
    hM.addHist("TH2F","h_acc_MT_oa","h_acc_MT_oa",100,0,1000,45,0,90);
    hM.addHist("TH2F","h_acc_PM_oa","h_acc_PM_oa",180,0,360,100,0,1000);

    hM.addHistArray("TH2F","h_reco_PT","h_reco_PT_w%i","h_reco_PT_w%i",180,0,360,45,0,90,0,0,0,"","","","",1);
    hM.addHistArray("TH2F","h_reco_MT","h_reco_MT_w%i","h_reco_MT_w%i",100,0,1000,45,0,90,0,0,0,"","","","",1);
    hM.addHistArray("TH2F","h_reco_PM","h_reco_PM_w%i","h_reco_PM_w%i",180,0,360,100,0,1000,0,0,0,"","","","",1);
    hM.addHistArray("TH3F","h_reco_MPT_ep","h_reco_MPT_ep_w%i","h_reco_MPT_ep_w%i",180,0,360,45,0,90,100,0,1000,"","","","",NWEIGHTS);
    hM.addHistArray("TH3F","h_reco_MPT_em","h_reco_MPT_em_w%i","h_reco_MPT_em_w%i",180,0,360,45,0,90,100,0,1000,"","","","",NWEIGHTS);

    hM.addHistArray("TH2F","h_reco_PT_richQa","h_reco_PT_richQa_w%i","h_reco_PT_richQa_w%i",180,0,360,45,0,90,0,0,0,"","","","",1);
    hM.addHistArray("TH2F","h_reco_MT_richQa","h_reco_MT_richQa_w%i","h_reco_MT_richQa_w%i",100,0,1000,45,0,90,0,0,0,"","","","",1);
    hM.addHistArray("TH2F","h_reco_PM_richQa","h_reco_PM_richQa_w%i","h_reco_PM_richQa_w%i",180,0,360,100,0,1000,0,0,0,"","","","",1);

    hM.addHistArray("TH2F","h_reco_PT_mlp","h_reco_PT_mlp_w%i","h_reco_PT_mlp_w%i",180,0,360,45,0,90,0,0,0,"","","","",NWEIGHTS);
    hM.addHistArray("TH2F","h_reco_MT_mlp","h_reco_MT_mlp_w%i","h_reco_MT_mlp_w%i",100,0,1000,45,0,90,0,0,0,"","","","",NWEIGHTS);
    hM.addHistArray("TH2F","h_reco_PM_mlp","h_reco_PM_mlp_w%i","h_reco_PM_mlp_w%i",180,0,360,100,0,1000,0,0,0,"","","","",NWEIGHTS);

    hM.addHistArray("TH2F","h_reco_PT_mlp_richQa","h_reco_PT_mlp_richQa_w%i","h_reco_PT_mlp_richQa_w%i",180,0,360,45,0,90,0,0,0,"","","","",NWEIGHTS);
    hM.addHistArray("TH2F","h_reco_MT_mlp_richQa","h_reco_MT_mlp_richQa_w%i","h_reco_MT_mlp_richQa_w%i",100,0,1000,45,0,90,0,0,0,"","","","",NWEIGHTS);
    hM.addHistArray("TH2F","h_reco_PM_mlp_richQa","h_reco_PM_mlp_richQa_w%i","h_reco_PM_mlp_richQa_w%i",180,0,360,100,0,1000,0,0,0,"","","","",NWEIGHTS);

    hM.addHistArray("TH2F","h_corr_PT_mlp_richQa","h_corr_PT_mlp_richQa_w%i","h_corr_PT_mlp_richQa_w%i",180,0,360,45,0,90,0,0,0,"","","","",NWEIGHTS);
    hM.addHistArray("TH2F","h_corr_MT_mlp_richQa","h_corr_MT_mlp_richQa_w%i","h_corr_MT_mlp_richQa_w%i",100,0,1000,45,0,90,0,0,0,"","","","",NWEIGHTS);
    hM.addHistArray("TH2F","h_corr_PM_mlp_richQa","h_corr_PM_mlp_richQa_w%i","h_corr_PM_mlp_richQa_w%i",180,0,360,100,0,1000,0,0,0,"","","","",NWEIGHTS);

    hM.addHistArray("TH2F","hmassptNP_kine", "hmassptNP_kine_parent%i", "hmassptNP_kine_parent%i" ,binsxmass,0,maxmass,binsxmom,0,maxmass,0,0,0,"","","","hmasspt_kine",8);
    hM.addHistArray("TH2F","hmassptPP_kine", "hmassptPP_kine_parent%i", "hmassptPP_kine_parent%i" ,binsxmass,0,maxmass,binsxmom,0,maxmass,0,0,0,"","","","hmasspt_kine",8);
    hM.addHistArray("TH2F","hmassptNN_kine", "hmassptNN_kine_parent%i", "hmassptNN_kine_parent%i" ,binsxmass,0,maxmass,binsxmom,0,maxmass,0,0,0,"","","","hmasspt_kine",8);

    hM.addHistArray("TH1F","hmassNP_kine", "hmassNP_kine_parent%i", "hmassNP_kine_parent%i" ,binsxmass,0,maxmass,0,0,0,0,0,0,"","","","hmass_kine",8);
    hM.addHistArray("TH1F","hmassPP_kine", "hmassPP_kine_parent%i", "hmassPP_kine_parent%i" ,binsxmass,0,maxmass,0,0,0,0,0,0,"","","","hmass_kine",8);
    hM.addHistArray("TH1F","hmassNN_kine", "hmassNN_kine_parent%i", "hmassNN_kine_parent%i" ,binsxmass,0,maxmass,0,0,0,0,0,0,"","","","hmass_kine",8);
    hM.addHistArray("TH1F","hoAngleNP_kine", "hoAngleNP_kine_parent%i", "hoAngleNP_kine_parent%i" ,200,0,200,0,0,0,0,0,0,"","","","hoAngle_kine",8);
    hM.addHistArray("TH1F","hoAnglePP_kine", "hoAnglePP_kine_parent%i", "hoAnglePP_kine_parent%i" ,200,0,200,0,0,0,0,0,0,"","","","hoAngle_kine",8);
    hM.addHistArray("TH1F","hoAngleNN_kine", "hoAngleNN_kine_parent%i", "hoAngleNN_kine_parent%i" ,200,0,200,0,0,0,0,0,0,"","","","hoAngle_kine",8);

    hM.addHistArray("TH1F","hmassNP_kine_oa9deg", "hmassNP_kine_oa9deg_parent%i", "hmassNP_kine_oa9deg_parent%i" ,binsxmass,0,maxmass,0,0,0,0,0,0,"","","","hmass_kine",8);
    hM.addHistArray("TH1F","hmassPP_kine_oa9deg", "hmassPP_kine_oa9deg_parent%i", "hmassPP_kine_oa9deg_parent%i" ,binsxmass,0,maxmass,0,0,0,0,0,0,"","","","hmass_kine",8);
    hM.addHistArray("TH1F","hmassNN_kine_oa9deg", "hmassNN_kine_oa9deg_parent%i", "hmassNN_kine_oa9deg_parent%i" ,binsxmass,0,maxmass,0,0,0,0,0,0,"","","","hmass_kine",8);
    hM.addHistArray("TH1F","hoAngleNP_kine_oa9deg", "hoAngleNP_kine_oa9deg_parent%i", "hoAngleNP_kine_oa9deg_parent%i" ,200,0,200,0,0,0,0,0,0,"","","","hoAngle_kine",8);
    hM.addHistArray("TH1F","hoAnglePP_kine_oa9deg", "hoAnglePP_kine_oa9deg_parent%i", "hoAnglePP_kine_oa9deg_parent%i" ,200,0,200,0,0,0,0,0,0,"","","","hoAngle_kine",8);
    hM.addHistArray("TH1F","hoAngleNN_kine_oa9deg", "hoAngleNN_kine_oa9deg_parent%i", "hoAngleNN_kine_oa9deg_parent%i" ,200,0,200,0,0,0,0,0,0,"","","","hoAngle_kine",8);
}

void analyzeKine(HHistMap &hM, HLoop *loop, Float_t &oaEmbed, TH3F *p3DEffEle, TH3F *p3DEffPos, Int_t mult_bin) {
    HGeantKine *kine1, *kine2;
    HCategory* kineCat = (HCategory*)HCategoryManager::getCategory(catGeantKine);
    if (kineCat == NULL) return;
    MyGeantBooker booker(kineCat);
    HIterator* iter1   = NULL;
    HIterator* iter2   = NULL;
    if(kineCat)
    {
        iter1   =(HIterator *)kineCat->MakeIterator("native");
        iter2   =(HIterator *)kineCat->MakeIterator("native");
    }
    Float_t imp = loop->getGeantHeader()->getImpactParameter();
    Int_t imp_bin = imp;
    if (imp_bin < 0) imp_bin = 0;
    if (imp_bin > 9) imp_bin = 9;

    map<HGeantTof *, int> tofHits = booker.countTofHits();
    map<HGeantTof *, int>::iterator itTof = tofHits.begin();
    for ( ; itTof != tofHits.end(); ++itTof) {
        HGeantTof *hit = itTof->first;
        hM.get("nhitCellTof",imp_bin)->Fill(itTof->second,(hit->getModule()-14)*8+hit->getCell());
    }
    map<HGeantRpc *, int> rpcHits = booker.countRpcHits();
    map<HGeantRpc *, int>::iterator itRpc = rpcHits.begin();
    for ( ; itRpc != rpcHits.end(); ++itRpc) {
        HGeantRpc *hit = itRpc->first;
        hM.get("nhitCellRpc",imp_bin)->Fill(itRpc->second,hit->getCell());
    }

    Int_t kineCatEntries = kineCat->getEntries();
    if (iter1 != NULL) {
        iter1->Reset();
        while ((kine1 = (HGeantKine*)iter1->Next())) {
            hM.get("hnKine")->Fill(0);
            Int_t parenttrack, medium, creation;
            kine1->getCreator(parenttrack,medium,creation);
            if (parenttrack != 0 && kine1->isPrimary()) cout << parenttrack << " " << kine1->isPrimary() << endl;

            TVector3 v1;
            Float_t px1, py1, pz1;
            kine1->getMomentum(px1,py1,pz1);
            v1.SetXYZ(px1,py1,pz1);
            Float_t p1     = kine1->getTotalMomentum();
            Float_t phi1   = transformPhi(v1.Phi());
            Float_t theta1 = v1.Theta() * TMath::RadToDeg();
            Bool_t acc1    = kine1->isInAcceptance(4,4,4,4,1,1);
            Int_t sec1     = 4;
            if      (phi1 < 60)  sec1 = 5;
            else if (phi1 < 120) sec1 = 0;
            else if (phi1 < 180) sec1 = 1;
            else if (phi1 < 240) sec1 = 2;
            else if (phi1 < 300) sec1 = 3;
            Float_t eff_res = getEfficiencyFactor(p3DEffEle, p3DEffPos, p1, theta1, phi1, -1);
            Float_t weight1 = (eff_res != 0.);
            hM.get("hnKine")->Fill(1);
            if (p1 < 100 || p1 > 1000) continue;
            hM.get("hnKine")->Fill(2);
            if (theta1 < 15 || theta1 > 85) continue;
            hM.get("hnKine")->Fill(3);
            if (phi1 > 180 && phi1 < 240) continue;
            hM.get("hnKine")->Fill(4);
            if (!checkPhi(phi1,4)) continue;
            hM.get("hnKine")->Fill(5);

            TLorentzVector vec1;
            fillVectorFromKine(vec1,kine1);

            if (acc1 && kine1->getGeneratorInfo() < 0 && eff_res > 0.05) {
                if (kine1->getID() == 3) {
                    hM.get("hnKine")->Fill(6);
                    hM.get2("h_acc_PT")->Fill(phi1,theta1);
                    hM.get2("h_acc_MT")->Fill(p1,theta1);
                    hM.get2("h_acc_PM")->Fill(p1,phi1);
                    hM.get3("h_acc_MPT_em",mult_bin)->Fill(phi1,theta1,p1);
                    hM.get3("h_acc_MPT_em",6)->Fill(phi1,theta1,p1);
                }
                else if (kine1->getID() == 2) {
                    hM.get2("h_acc_PT")->Fill(phi1,theta1);
                    hM.get2("h_acc_MT")->Fill(p1,theta1);
                    hM.get2("h_acc_PM")->Fill(p1,phi1);
                    hM.get3("h_acc_MPT_ep",mult_bin)->Fill(phi1,theta1,p1);
                    hM.get3("h_acc_MPT_ep",6)->Fill(phi1,theta1,p1);
                }
            }
/*
            for (int k = 0; k < kineCatEntries; ++k) {
                kine2 = HCategoryManager::getObject(kine2,kineCat,k);
                TVector3 v2;
                Float_t px2, py2, pz2;
                kine2->getMomentum(px2,py2,pz2);
                v2.SetXYZ(px2,py2,pz2);
                Float_t p2     = kine2->getTotalMomentum();
                Float_t phi2   = transformPhi(v2.Phi());
                Float_t theta2 = v2.Theta() * TMath::RadToDeg();
                Bool_t acc2    = kine2->isInAcceptance(4,4,4,4,1,1);
                Int_t sec2     = 4;
                TLorentzVector vec2;
                fillVectorFromKine(vec2,kine2);

                TLorentzVector dilep = vec1 + vec2;

                Float_t oa = v2.Angle(v1) * TMath::RadToDeg();
                Float_t minv = dilep.M();
                Float_t pt = dilep.Perp();
                bool sameparent = HGeantKine::getParent( kine1->getTrack() ) == HGeantKine::getParent( kine2->getTrack() );
                bool primary = kine1->isPrimary() && kine2->isPrimary();
                bool oacond = oa > oAngleCut;
                bool ids32 = kine1->getID() == 3 && kine2->getID() == 2;

                if (kine1->getGeneratorInfo2() < 0 || kine2->getGeneratorInfo2() < 0) continue;
                if (primary && ids32) {
                    oaEmbed = oa;
                    if (acc1 && oacond && eff_res > 0.05) {
                        hM.get("h_acc_PT_oa")->Fill(phi1,theta1);
                        hM.get("h_acc_MT_oa")->Fill(p1,theta1);
                        hM.get("h_acc_PM_oa")->Fill(p1,phi1);
                    }


                    if      (phi2 < 60)  sec2 = 5;
                    else if (phi2 < 120) sec2 = 0;
                    else if (phi2 < 180) sec2 = 1;
                    else if (phi2 < 240) sec2 = 2;
                    else if (phi2 < 300) sec2 = 3;

                    Int_t parent_id_ind = -1;
                    if (primary || sameparent) {
                        HGeantKine *parent = HGeantKine::getParent(kine1->getTrack());
                        int parent_id = -1;
                        if (parent != NULL) {
                            parent_id = parent->getID();
                        }
                        switch (parent_id) {
                            case -1:  parent_id_ind = 1; break;
                            case 1:   parent_id_ind = 2; break;
                            case 7:   parent_id_ind = 3; break;
                            case 17:  parent_id_ind = 4; break;
                            default : parent_id_ind = 5; break;
                        }
                    }

                    //                            if (p2 < 100 || p2 > 1000) continue;
                    if (theta2 < 15 || theta2 > 85) continue;
                    //                            if (phi2 > 180 && phi2 < 240) continue;
                    if (!checkPhi(phi2,4)) continue;
                    if (!(acc1 && acc2)) continue;;
                    if (oa < 0) continue;
                    //                        if (sec1 != sec2) continue;

                    Float_t weight2 = (getEfficiencyFactor(p3DEffEle, p3DEffPos, p2, theta2, phi2, 1) != 0.);
                    Float_t weight = weight1 * weight2;
                    //                            if (getEfficiencyFactor(p3DEffEle[0], p3DEffPos[0], p2, theta2, phi2, 1) < 0.05 && weight > 0) cout << "ERROR WITH EFF : " << getEfficiencyFactor(p4DEffEle[0], p3DEffPos[0], p2, theta2, phi2, 1) << endl;

                    if (parent_id_ind > -1) {
                        if ((kine1->getID() == 3 && kine2->getID() == 2)) {
                            if (oa > 180) cout << "Too big opening angle: " << oa << endl; 
                            hM.get2("hmassptNP_kine",parent_id_ind)->Fill(minv,pt,weight);
                            hM.get("hmassNP_kine",parent_id_ind)->Fill(minv,weight);
                            hM.get("hoAngleNP_kine",parent_id_ind)->Fill(oa,weight);
                            hM.get2("hmassptNP_kine",0)->Fill(minv,pt,weight);
                            hM.get("hmassNP_kine",0)->Fill(minv,weight);
                            hM.get("hoAngleNP_kine",0)->Fill(oa,weight);
                        }
                        if (kine1->getID() == 2 && kine2->getID() == 2) {
                            hM.get2("hmassptPP_kine",parent_id_ind)->Fill(minv,pt,weight);
                            hM.get("hmassPP_kine",parent_id_ind)->Fill(minv,weight);
                            hM.get("hoAnglePP_kine",parent_id_ind)->Fill(oa,weight);
                            hM.get2("hmassptPP_kine",0)->Fill(minv,pt,weight);
                            hM.get("hmassPP_kine",0)->Fill(minv,weight);
                            hM.get("hoAnglePP_kine",0)->Fill(oa,weight);
                        }
                        if (kine1->getID() == 3 && kine2->getID() == 3) {
                            hM.get2("hmassptNN_kine",parent_id_ind)->Fill(minv,pt,weight);
                            hM.get("hmassNN_kine",parent_id_ind)->Fill(minv,weight);
                            hM.get("hoAngleNN_kine",parent_id_ind)->Fill(oa,weight);
                            hM.get2("hmassptNN_kine",0)->Fill(minv,pt,weight);
                            hM.get("hmassNN_kine",0)->Fill(minv,weight);
                            hM.get("hoAngleNN_kine",0)->Fill(oa,weight);
                        }
                        if (oa > oAngleCut) {
                            if ((kine1->getID() == 3 && kine2->getID() == 2)) {
                                hM.get("hmassNP_kine_oa9deg",parent_id_ind)->Fill(minv,weight);
                                hM.get("hoAngleNP_kine_oa9deg",parent_id_ind)->Fill(oa,weight);
                                hM.get("hmassNP_kine_oa9deg",0)->Fill(minv,weight);
                                hM.get("hoAngleNP_kine_oa9deg",0)->Fill(oa,weight);
                            }
                            if (kine1->getID() == 2 && kine2->getID() == 2) {
                                hM.get("hmassPP_kine_oa9deg",parent_id_ind)->Fill(minv,weight);
                                hM.get("hoAnglePP_kine_oa9deg",parent_id_ind)->Fill(oa,weight);
                                hM.get("hmassPP_kine_oa9deg",0)->Fill(minv,weight);
                                hM.get("hoAnglePP_kine_oa9deg",0)->Fill(oa,weight);
                            }
                            if (kine1->getID() == 3 && kine2->getID() == 3) {
                                hM.get("hmassNN_kine_oa9deg",parent_id_ind)->Fill(minv,weight);
                                hM.get("hoAngleNN_kine_oa9deg",parent_id_ind)->Fill(oa,weight);
                                hM.get("hmassNN_kine_oa9deg",0)->Fill(minv,weight);
                                hM.get("hoAngleNN_kine_oa9deg",0)->Fill(oa,weight);
                            }
                        }
                    }
                    else {
                        if ((kine1->getID() == 3 && kine2->getID() == 2)) {
                            hM.get("hmassNP_kine",6)->Fill(minv,weight);
                            hM.get("hoAngleNP_kine",6)->Fill(oa,weight);
                        }
                        if (kine1->getID() == 2 && kine2->getID() == 2) {
                            hM.get("hmassPP_kine",6)->Fill(minv,weight);
                            hM.get("hoAnglePP_kine",6)->Fill(oa,weight);
                        }
                        if (kine1->getID() == 3 && kine2->getID() == 3) {
                            hM.get("hmassNN_kine",6)->Fill(minv,weight);
                            hM.get("hoAngleNN_kine",6)->Fill(oa,weight);
                        }
                        if (oa > oAngleCut) {
                            if ((kine1->getID() == 3 && kine2->getID() == 2)) {
                                hM.get("hmassNP_kine_oa9deg",6)->Fill(minv,weight);
                                hM.get("hoAngleNP_kine_oa9deg",6)->Fill(oa,weight);
                            }
                            if (kine1->getID() == 2 && kine2->getID() == 2) {
                                hM.get("hmassPP_kine_oa9deg",6)->Fill(minv,weight);
                                hM.get("hoAnglePP_kine_oa9deg",6)->Fill(oa,weight);
                            }
                            if (kine1->getID() == 3 && kine2->getID() == 3) {
                                hM.get("hmassNN_kine_oa9deg",6)->Fill(minv,weight);
                                hM.get("hoAngleNN_kine_oa9deg",6)->Fill(oa,weight);
                            }
                        }
                    }
                    if ((kine1->getID() == 3 && kine2->getID() == 2)) {
                        hM.get("hmassNP_kine",7)->Fill(minv,weight);
                        hM.get("hoAngleNP_kine",7)->Fill(oa,weight);
                    }
                    if (kine1->getID() == 2 && kine2->getID() == 2) {
                        hM.get("hmassPP_kine",7)->Fill(minv,weight);
                        hM.get("hoAnglePP_kine",7)->Fill(oa,weight);
                    }
                    if (kine1->getID() == 3 && kine2->getID() == 3) {
                        hM.get("hmassNN_kine",7)->Fill(minv,weight);
                        hM.get("hoAngleNN_kine",7)->Fill(oa,weight);
                    }
                    if (oa > oAngleCut) {
                        if ((kine1->getID() == 3 && kine2->getID() == 2)) {
                            hM.get("hmassNP_kine_oa9deg",7)->Fill(minv,weight);
                            hM.get("hoAngleNP_kine_oa9deg",7)->Fill(oa,weight);
                        }
                        if (kine1->getID() == 2 && kine2->getID() == 2) {
                            hM.get("hmassPP_kine_oa9deg",7)->Fill(minv,weight);
                            hM.get("hoAnglePP_kine_oa9deg",7)->Fill(oa,weight);
                        }
                        if (kine1->getID() == 3 && kine2->getID() == 3) {
                            hM.get("hmassNN_kine_oa9deg",7)->Fill(minv,weight);
                            hM.get("hoAngleNN_kine_oa9deg",7)->Fill(oa,weight);
                        }
                    }
                }
            }
*/
        }
    }
}
