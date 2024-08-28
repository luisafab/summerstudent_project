#include <ROOT/RDataFrame.hxx>

// calculating the armenteros-podolanski plot for given data and saving new variables in root file 
// string file : name of file
// string tree : name of tree
// bool MC : true : only lambda events are consicered (!! The MC Information in the file is necessary !!)
//           false : all events are considered
// string save_file : file to save plots 
// string snapshot : name of snapshot to save data to
// float pT_high pT_low: values to cut on pT

// save used constants
constexpr double kK0sMass{0.497614};
constexpr double kPiMass{0.13957061};
constexpr double kProtonMass{0.9382720813};

// help functions to calculate the variables of the AP Plot using the momenta of the daughter particles
// calculates beta of the decaying particle using the toal momentum
float betaCalc(float ptotal){
    return std::sqrt(1./(std::pow((kK0sMass/ptotal),2.)+1.));
}
// total momentum of the mother is the sum of the momenta of the daughters
float ptotal(float posx, float posy, float posz, float negx, float negy, float negz){
    float massPos;
    float massNeg;

    //define vectors for momenta and calclate the length 
    TVector3 pPos(posx,posy,posz);
    TVector3 pNeg(negx,negy,negz);
    TVector3 momentum=pPos+pNeg;
    float momentum_norm=std::sqrt(momentum.Dot(momentum));

    return momentum_norm;
}

// takes the momenta of the negative and positive particle as argument
float alpha(float posx, float posy, float posz, float negx, float negy, float negz){
    // create momentum vector for both particles
    TVector3 pPos(posx,posy,posz);
    TVector3 pNeg(negx,negy,negz);
    // the momentum of V0 is the sum of the other two
    TVector3 pV0 =pPos+pNeg;
    // normalize pV0
    float pV0_norm=std::sqrt(pV0.Dot(pV0));
    // calculate longitudinal part of pPos and pNeg
    float pPlL = pPos.Dot(pV0)/pV0_norm;
    float pNegL = pNeg.Dot(pV0)/pV0_norm;
    // calculate alpha
    float alpha = (pPlL-pNegL)/((pPlL+pNegL));

    return alpha;
}

 // calculate pT of the V0 particle
float pt_v0(float posx, float posy, float posz, float negx, float negy, float negz){
    TVector3 pPos(posx,posy,posz);
    TVector3 pNeg(negx,negy,negz);
    TVector3 pV0 = pPos+pNeg;
    float pT = std::hypot(pV0[0],pV0[1]);

    return pT;
}

float qT_AP(float posx, float posy, float posz, float negx, float negy, float negz){
    // create momentum vector for both particles
    TVector3 pPos(posx,posy,posz);
    TVector3 pNeg(negx,negy,negz);
    // the momentum of V0 is the sum of the other two
    TVector3 pV0 = pPos+pNeg;
    // normalize pV0
    float pV0_norm = std::sqrt(pV0.Dot(pV0));
    // cross product of pPos and pV0
    TVector3 cross = pPos.Cross(pV0);
    // calculate qT
    float qT = (std::sqrt(cross.Dot(cross)))/pV0_norm;

    return qT;
}

// calculation of the invariant mass of lambdas
float invariantmass_lambda(float alpha, float posx, float posy, float posz, float negx, float negy, float negz){
    float massPos;
    float massNeg;

    //define vectors for momenta and calclate the length 
    TVector3 pPos(posx,posy,posz);
    TVector3 pNeg(negx,negy,negz);
    TVector3 momentum=pPos+pNeg;
    float pPos_norm=std::sqrt(pPos.Dot(pPos));
    float pNeg_norm=std::sqrt(pNeg.Dot(pNeg));
    float momentum_norm=std::sqrt(momentum.Dot(momentum));


    // use alpha to know the decay and the resulting particles
    // alpha>0 lambda -> p +pi-
    // alpha<0 pi+
    // assign the mass 
    if (alpha>0){
        massPos=kProtonMass;
        massNeg=kPiMass;
    }
    else{
        massPos=kPiMass;
        massNeg=kProtonMass;
    }

    // use mass and momentum to calculate the energy of the two particles
    float EPos=std::sqrt(std::pow(pPos_norm,2.)+std::pow((massPos),2.));
    float ENeg=std::sqrt(std::pow(pNeg_norm,2.)+std::pow((massNeg),2.));

    // total energy 
    float Energy = EPos + ENeg;

    // calculate the invariant mass using the fourmomentum vector
    float invariantmass = std::sqrt((std::pow(Energy,2.)-std::pow(momentum_norm,2.)));

    return invariantmass;
}

//calculation of the invariant mass of K0s
float invariantmass_K0(float posx, float posy, float posz, float negx, float negy, float negz){
    float massPos;
    float massNeg;

    //define vectors for momenta and calclate the length 
    TVector3 pPos(posx,posy,posz);
    TVector3 pNeg(negx,negy,negz);
    TVector3 momentum=pPos+pNeg;
    float pPos_norm=std::sqrt(pPos.Dot(pPos));
    float pNeg_norm=std::sqrt(pNeg.Dot(pNeg));
    float momentum_norm=std::sqrt(momentum.Dot(momentum));


    // decays in pi+ and pi-, masses are equal
    massPos=kPiMass;
    massNeg=kPiMass;

    // use mass and momentum to calculate the energy of the two particles
    float EPos=std::sqrt(std::pow(pPos_norm,2.)+std::pow((massPos),2.));
    float ENeg=std::sqrt(std::pow(pNeg_norm,2.)+std::pow((massNeg),2.));

    // total energy 
    float Energy = EPos + ENeg;

    // calculate the invariant mass using the fourmomentum vector
    float invariantmass = std::sqrt((std::pow(Energy,2.)-std::pow(momentum_norm,2.)));

    return invariantmass;
}



void armenterosPlot(TString file="Snapshot.root", TString tree="NewVariables") {

    // saving plots in a root file
    //std::unique_ptr<TFile> myFile( TFile::Open(file_save, "RECREATE") );

    // define a data frame to use 
    ROOT::RDataFrame complete_df(tree, file);
    
    auto hpT = complete_df.Histo1D({"p_V0","total momentum V0; p_{V0}; counts",1000,0.,8.},"ptotal");
    //hpT->Write();
    // plot invariant mass after cut
    auto hK0_new = complete_df.Histo1D({"MinvK0","Invariant mass K0; M_{inv}(K_{0}^{s})",1000,0.4,0.55},"Minv_K0");
    //hK0_new->Write();

    // define and plot a histogram which is the resulting armenteros plot after appliying the cuts
    TH2D *h = new TH2D("AP_2D","Armenteros-Podolanski Plot;#alpha;p_{T}",100,-1.,1.,100,0.,0.25);
    
    auto histo = complete_df.Histo2D(*h,"alpha","qT");

    histo->GetXaxis()->SetTitle("#alpha");
    histo->GetYaxis()->SetTitle("q_{T}");

    //histo->Write();

    //myFile->Close();

}