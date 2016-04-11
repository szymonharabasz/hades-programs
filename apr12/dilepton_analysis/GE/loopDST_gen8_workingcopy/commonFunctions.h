#ifndef commonFunctions_H
#define commonFunctions_H

Float_t getEfficiencyFactorKine(TH3F *p3DEffEle, TH3F *p3DEffPos, Float_t mom, Float_t theta, Float_t phi, Float_t charge, Float_t corrMax = 1)
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

    fEff1 = fEff1*corrMax;
    //Float_t fCorrection = 1./fEff1;

    if (fEff1 > 0.05) {
	return fEff1;
    }
    else
    {
	return 0.;
    }

}

Float_t getEffWeightKine(Float_t phi1, Float_t theta1, Float_t p1, Int_t id1, Float_t phi2, Float_t theta2, Float_t p2, Int_t id2 ,TH3F *p3DEffEle, TH3F *p3DEffPos)
{
    Int_t charge[2];
    if(id1 == 2)
	charge[0] = 1;
    else
	charge[0] = -1;

    if(id2 == 2)
	charge[1] = 1;
    else
	charge[1] = -1;

   // p1   = p1*1000;
   // p2   = p2*1000;
   // theta1 = theta1*TMath::RadToDeg();             // polar angle in deg for my matrices
   // theta2 = theta2*TMath::RadToDeg();             // polar angle in deg for my matrices
    //phi1   = (phi1+TMath::Pi())*TMath::RadToDeg(); // azimuthal angle in deg for my matrices
    //phi2   = (phi2+TMath::Pi())*TMath::RadToDeg(); // azimuthal angle in deg for my matrices


    //Efficiency correction
    Float_t weight1 = -1;
    Float_t weight2 = -1;
	weight1 = getEfficiencyFactorKine(p3DEffEle, p3DEffPos, p1, theta1, phi1, charge[0],1);
	weight2 = getEfficiencyFactorKine(p3DEffEle, p3DEffPos, p2, theta2, phi2, charge[1],1);
	return weight1 * weight2;
}

#endif
