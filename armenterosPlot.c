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
constexpr double kK0sMass{0.497164};
constexpr double kPiMass{0.13957061};
constexpr double kProtonMass{0.9382720813};

// help functions to calculate the variables of the AP Plot using the momenta of the daughter particles
// calculates beta of the decaying particle using the toal momentum
float beta(float ptotal){
    return sqrt(1./(pow((kK0sMass/ptotal),2.)+1.));
}
// total momentum of the mother is the sum of the momenta of the daughters
float ptotal(float posx, float posy, float posz, float negx, float negy, float negz){
    float massPos;
    float massNeg;

    //define vectors for momenta and calclate the length 
    TVector3 pPos(posx,posy,posz);
    TVector3 pNeg(negx,negy,negz);
    TVector3 momentum=pPos+pNeg;
    float momentum_norm=sqrt(momentum.Dot(momentum));

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
    float pV0_norm=sqrt(pV0.Dot(pV0));
    // calculate longitudinal part of pPos and pNeg
    float pPlL = pPos.Dot(pV0)/pV0_norm;
    float pNegL = pNeg.Dot(pV0)/pV0_norm;
    // calculate alpha
    float alpha= (pPlL-pNegL)/((pPlL+pNegL));
    float beta_calc=beta(ptotal(posx, posy, posz, negx, negy, negz));

    return alpha/beta_calc;
}

 // calculate pT of the V0 particle
float pt_v0(float posx, float posy, float posz, float negx, float negy, float negz){
    TVector3 pPos(posx,posy,posz);
    TVector3 pNeg(negx,negy,negz);
    TVector3 pV0 =pPos+pNeg;
    float pT=std::hypot(pV0[0],pV0[1]);

    return pT;
}

float qT(float posx, float posy, float posz, float negx, float negy, float negz){
    // create momentum vector for both particles
    TVector3 pPos(posx,posy,posz);
    TVector3 pNeg(negx,negy,negz);
    // the momentum of V0 is the sum of the other two
    TVector3 pV0 =pPos+pNeg;
    // normalize pV0
    float pV0_norm=sqrt(pV0.Dot(pV0));
    // cross product of pPos and pV0
    TVector3 cross=pPos.Cross(pV0);
    // calculate qT
    float qT = (sqrt(cross.Dot(cross)))/pV0_norm;

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
    float pPos_norm=sqrt(pPos.Dot(pPos));
    float pNeg_norm=sqrt(pNeg.Dot(pNeg));
    float momentum_norm=sqrt(momentum.Dot(momentum));


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
    float EPos=sqrt(pow(pPos_norm,2.)+pow((massPos),2.));
    float ENeg=sqrt(pow(pNeg_norm,2.)+pow((massNeg),2.));

    // total energy 
    float Energy = EPos + ENeg;

    // calculate the invariant mass using the fourmomentum vector
    float invariantmass = sqrt((pow(Energy,2.)-pow(momentum_norm,2.)));

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
    float pPos_norm=sqrt(pPos.Dot(pPos));
    float pNeg_norm=sqrt(pNeg.Dot(pNeg));
    float momentum_norm=sqrt(momentum.Dot(momentum));


    // decays in pi+ and pi-, masses are equal
    massPos=kPiMass;
    massNeg=kPiMass;

    // use mass and momentum to calculate the energy of the two particles
    float EPos=sqrt(pow(pPos_norm,2.)+pow((massPos),2.));
    float ENeg=sqrt(pow(pNeg_norm,2.)+pow((massNeg),2.));

    // total energy 
    float Energy = EPos + ENeg;

    // calculate the invariant mass using the fourmomentum vector
    float invariantmass = sqrt((pow(Energy,2.)-pow(momentum_norm,2.)));

    return invariantmass;
}



void armenterosPlot(TString file="data/AnalysisResults_treesAP_data_LHC22o_apass6_small.root", TString tree="O2v0tableap",bool MC=false, TString file_save="APPlot_lowPt.root", TString snapshot="APPlot_df_lowPt_09.root",float pT_low=0.,float pT_high=10.) {

    // saving plots in a root file
    std::unique_ptr<TFile> myFile( TFile::Open(file_save, "RECREATE") );

    // define a data frame to use 
    ROOT::RDataFrame df_MC(tree, file);
    // filter it
    //     for MC==true only the recombined (isreco==true) lambdas (PDGCODE=3122) or K0s(310) are considered
    //     else there is a dummy cut to use the same dataframe
    //auto df_cos = df_MC.Filter(MC? "(fPDGCode == 3122 || fPDGCode== -3122) && fIsReco" : "fLen>0");
    auto df_cos = df_MC.Filter(MC? "(fPDGCode == 310) && fIsReco" : "fLen>0");
    //auto df_cos = df_MC.Filter(MC? "fIsReco" : "fLen>0");
    std::cout<<"filtered PDG"<<std::endl;

    //apply cut on cosPA
    auto df = df_cos.Filter("fCosPA>0.999");

    // define new columns with the two new variables 
    auto new_column_alpha = df.Define("alpha", alpha, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    std::cout<<"Calculated alpha"<<std::endl;
    // get the final data frame
    auto new_column_qT = new_column_alpha.Define("qT", qT, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    std::cout<<"calculated qT"<<std::endl;

    auto new_column_pT = new_column_qT.Define("pTV0", pt_v0, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    std::cout<<"calculated pT V0"<<std::endl;
    
    // add both masses as new column
    auto new_minv_l = new_column_pT.Define("Minv_lambda", invariantmass_lambda, {"alpha", "fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    std::cout<<"Calculated inv mass lambda"<<std::endl;
    auto new_minv_k0 = new_minv_l.Define("Minv_K0", invariantmass_K0, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    std::cout<<"Calculated invariant mass K0"<<std::endl;

    // create histograms to filter invariant mass and one to plot pT
    auto hL = new_minv_k0.Histo1D({"MinvL_all","Invariant mass Lambda; M_{inv}(#Lambda)",1000,1.07,1.35},"Minv_lambda");
    auto hK0 = new_minv_k0.Histo1D({"MinvK0_all","Invariant mass K0; M_{inv}(K_{0}^{s})",1000,0.275,0.6},"Minv_K0");
    auto hpTV0 = new_minv_k0.Histo1D({"pTV0_before","pT V0; pT [GeV/c]",1000,0.,5.},"pTV0");

    // plot invariant mass
    hL->Write();
    hK0->Write();
    // and pT
    hpTV0->Write();
    
    
    // add filter on calculated invariant mass
    float lambda_low=1.1;
    float lambda_high=1.13;
    float K0_low=0.49;
    float K0_high=0.502;
    // use lambda function
    auto insidePeakK =[&](float Minv){
        bool isInside;
        float low=K0_low; 
        float high=K0_high;
        if (Minv>low && Minv < high){
            isInside=true;
        }else{
            isInside=false;
        }
        return isInside;
    };
    auto insidePeakL =[&](float Minv){
        bool isInside;
        float low=lambda_low; 
        float high=lambda_high;
        if (Minv>low && Minv < high){
            isInside=true;
        }else{
            isInside=false;
        }
        return isInside;
    };
    auto outsidePeakL =[&](float Minv){
        return !insidePeakL(Minv);
    };
    auto insidePtRange =[&](float pT){
        bool isInside;
        if (pT>pT_low && pT < pT_high){
            isInside=true;
        }else{
            isInside=false;
        }
        return isInside;
    };
    
    // one filter just to ckeck on lambdas, not needed for the result in the end, just to plot the filter on invariant mass of lambda 
    auto new_filter_l = new_minv_k0.Filter(insidePeakL, {"Minv_lambda"});
    std::cout<<"Filter 1"<<std::endl;
    // filter on K0s
    auto new_filter_k = new_minv_k0.Filter(insidePeakK, {"Minv_K0"});
    std::cout<<"Filter 2"<<std::endl;
    // filter on particles which are in both peaks to reduce possible background 
    auto new_filter_both = new_filter_k.Filter(outsidePeakL, {"Minv_lambda"});
    std::cout<<"Filter 3"<<std::endl;
    // Filter on pT of V0
    auto new_filter_pT = new_filter_both.Filter(insidePtRange,{"pTV0"});
    std::cout<<"Filter 4"<<std::endl;
    // add beta and ptotal as variables
    auto new_ptotal = new_filter_pT.Define("ptotal", ptotal, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    std::cout<<"Calculated total momentum"<<std::endl;
    auto complete_df = new_ptotal.Define("beta", beta, {"ptotal"});
    std::cout<<"Calculated beta using the total momentum of V0s"<<std::endl;
    // plot pTV0 of all particles after the mass cut and the ones in the chosen range
    //auto hpTV0_a = complete_df.Histo1D({"pTV0_after","pT V0; pT [GeV/c]",1000,0.,4.},"pTV0");
    //hpTV0_a->Write();
    //auto hpTV0_all = new_filter_both.Histo1D({"pTV0_massCut","pT V0; pT [GeV/c]",1000,0.,5.},"pTV0");
    //hpTV0_all->Write();*/
    auto hpT = complete_df.Histo1D({"p_V0","total momentum V0; p_{V0}; counts",1000,0.,8.},"ptotal");
    hpT->Write();


    // define and plot a histogram which is the resulting armenteros plot after appliying the cuts
    TH2D *h = new TH2D("AP_2D","Armenteros-Podolanski Plot;#alpha;p_{T}",100,-1.,1.,100,0.,0.25);
    
    auto histo = complete_df.Histo2D(*h,"alpha","qT");

    histo->GetXaxis()->SetTitle("#alpha");
    histo->GetYaxis()->SetTitle("q_{T}");

    histo->Write();
    


    // save snapshot to use tree in next steps
    // complete_df.Snapshot("NewVariables", snapshot ,{"alpha","qT","Minv_lambda","Minv_K0","pTV0","ptotal","beta"});

    // save AP plot before cuts 
    TH2D *h1 = new TH2D("AP_2D_before","Armenteros-Podolanski Plot before cuts;#alpha;q_{T}",100,-1.,1.,100,0.,0.25);
    
    auto histo1 = new_filter_both.Histo2D(*h1,"alpha","qT");

    histo1->GetXaxis()->SetTitle("#alpha");
    histo1->GetYaxis()->SetTitle("q_{T}");

    histo1->Write();

    myFile->Close();

}

// possibility to plot all the variables 
/*

    auto cos = above.Histo1D({"CosPA","CosPA; CosPA",100,0.9988,1.0002},"fCosPA");
    cos->Write();
    auto radius = above.Histo1D({"Radius","radius; Radius",100,0.,31.},"fRadius");
    radius->Write();
    auto dcaV0 = above.Histo1D({"dcaV0PV","dcaV0PV; dcaV0PV",100,0.,31.},"fDcaV0PV");
    dcaV0->Write();
    auto dcaV0T = above.Histo1D({"fDCAV0Tracks","fDCAV0Tracks; fDCAV0Tracks",100,0.,1.1},"fDcaV0Tracks");
    dcaV0T->Write();
    auto dcaV0P = above.Histo1D({"fDCANegPV","fDCANegPV; fDCANegPV",100,0.,7.5},"fDcaNegPV");
    dcaV0P->Write();
    auto dcaV0N = above.Histo1D({"fDCAV0Tracks","fDCAV0Tracks; fDCAV0Tracks",100,0.,1.1},"fDcaV0Tracks");
    dcaV0N->Write();
    auto pxp = above.Histo1D({"fPxPos","fPxPos; fPxPos",100,-0.6,0.6},"fPxPos");
    pxp->Write();
    auto pyp = above.Histo1D({"fPyPos","fPyPos; fPyPos",100,-0.6,0.6},"fPyPos");
    pyp->Write();
    auto pzp = above.Histo1D({"fPzPos","fPzPos; fPzPos",100,-0.6,0.6},"fPzPos");
    pzp->Write();
    auto pxn = above.Histo1D({"fPxNeg","fPxNeg; fPxNeg",100,-0.6,0.6},"fPxNeg");
    pxn->Write();
    auto pyn = above.Histo1D({"fPyNeg","fPyNeg; fPyNeg",100,-0.6,0.6},"fPyNeg");
    pyn->Write();
    auto pzn = above.Histo1D({"fPzNeg","fPzNeg; fPzNeg",100,-0.6,0.6},"fPzNeg");
    pzn->Write();
    auto len = above.Histo1D({"fLen","fLen; fLen",100,-1,35.},"fLen");
    len->Write();
    auto eta = above.Histo1D({"fEta","fEta; fEta",100,-1.,1.},"fEta");
    eta->Write();
 */   