#include <ROOT/RDataFrame.hxx>

// calculating the armenteros-podolanski plot for given data and saving new variables in root file 
// string file : name of file
// string tree : name of tree
// bool MC : true : only lambda events are consicered (!! The MC Information in the file is necessary !!)
//           false : all events are considered
// string save_file : file to save plots 

void armenterosPlot(TString file="data/AnalysisResults_treesAP_data_LHC22o_apass6_small.root", TString tree="O2v0tableap",bool MC=false, TString file_save="APPlot.root") {

    // saving plots in a root file
    std::unique_ptr<TFile> myFile( TFile::Open(file_save, "RECREATE") );

    // define used masses
    float massProton=0.9382720813;
    float massPion= 0.13957061;

    // define a data frame to use 
    ROOT::RDataFrame df_MC(tree, file);
    // filter it
    //     for MC==true only the recombined (isreco==true) lambdas (PDGCODE=3122) are considered
    //     else there is a dummy cut to use the same dataframe
    auto df_cos = df_MC.Filter(MC? "(fPDGCode == 3122 || fPDGCode== -3122) && fIsReco" : "fLen>0");
    std::cout<<"filtered PDG"<<std::endl;

    //apply cut on cosPA
    auto df = df_cos.Filter("fCosPA>0.999");

    // lambda functions to calculate the variables of the AP Plot: alpha and pT
    // takes the momenta of the negative and positive particle as argument
    auto alpha = [&](float posx, float posy, float posz, float negx, float negy, float negz){
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
        float alpha= (pPlL-pNegL)/(pPlL+pNegL);

        return alpha;
    };

    auto pT = [&](float posx, float posy, float posz, float negx, float negy, float negz){
        // create momentum vector for both particles
        TVector3 pPos(posx,posy,posz);
        TVector3 pNeg(negx,negy,negz);
        // the momentum of V0 is the sum of the other two
        TVector3 pV0 =pPos+pNeg;
        // normalize pV0
        float pV0_norm=sqrt(pV0.Dot(pV0));
        // cross product of pPos and pV0
        TVector3 cross=pPos.Cross(pV0);
        // calculate pT
        float pT = (sqrt(cross.Dot(cross)))/pV0_norm;

        return pT;
    };


    // define new columns with the two new variables 
    auto new_column_alpha = df.Define("alpha", alpha, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    std::cout<<"Calculated alpha"<<std::endl;
    // get the final data frame
    auto new_column_pT = new_column_alpha.Define("pT", pT, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    std::cout<<"calculated pT"<<std::endl;

    // calculation of the invariant mass of lambdas
    auto invariantmass_lambda = [&](float alpha, float posx, float posy, float posz, float negx, float negy, float negz){
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
            massPos=massProton;
            massNeg=massPion;
        }
        else{
            massPos=massPion;
            massNeg=massProton;
        }

        // use mass and momentum to calculate the energy of the two particles
        float EPos=sqrt(pow(pPos_norm,2.)+pow((massPos),2.));
        float ENeg=sqrt(pow(pNeg_norm,2.)+pow((massNeg),2.));

        // total energy 
        float Energy = EPos + ENeg;

        // calculate the invariant mass using the fourmomentum vector
        float invariantmass = sqrt((pow(Energy,2.)-pow(momentum_norm,2.)));

        return invariantmass;
    };

    //calculation of the invariant mass of K0s
    auto invariantmass_K0 = [&](float posx, float posy, float posz, float negx, float negy, float negz){
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
        massPos=massPion;
        massNeg=massPion;

        // use mass and momentum to calculate the energy of the two particles
        float EPos=sqrt(pow(pPos_norm,2.)+pow((massPos),2.));
        float ENeg=sqrt(pow(pNeg_norm,2.)+pow((massNeg),2.));

        // total energy 
        float Energy = EPos + ENeg;

        // calculate the invariant mass using the fourmomentum vector
        float invariantmass = sqrt((pow(Energy,2.)-pow(momentum_norm,2.)));

        return invariantmass;
    };

    // add both masses as new column
    auto new_minv_l = new_column_pT.Define("Minv_lambda", invariantmass_lambda, {"alpha", "fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    std::cout<<"Calculated inv mass lambda"<<std::endl;
    auto new_minv_k0 = new_minv_l.Define("Minv_K0", invariantmass_K0, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    std::cout<<"Calculated invariant mass K0"<<std::endl;

    // create histograms to filter invariant mass
    auto hL = new_minv_k0.Histo1D({"MinvL","Invariant mass Lambda; M_{inv}(#Lambda)",1000,1.07,1.35},"Minv_lambda");
    auto hK0 = new_minv_k0.Histo1D({"MinvK0","Invariant mass K0; M_{inv}(K_{0}^{s})",1000,0.275,0.6},"Minv_K0");

    // plot invariant mass
    hL->Write();
    hK0->Write();


    // add filter
    float lambda_low=1.113;
    float lambda_high=1.118;
    float K0_low=0.488;
    float K0_high=0.505;
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
    // one filter just to ckeck on lambdas, not needed for the result in the end
    auto new_filter_l = new_minv_k0.Filter(insidePeakL, {"Minv_lambda"});
    std::cout<<"Filter 1"<<std::endl;
    // filter on K0s
    auto new_filter_k = new_minv_k0.Filter(insidePeakK, {"Minv_K0"});
    std::cout<<"Filter 2"<<std::endl;
    // filter on particles which are in both peaks to reduce possible background 
    auto complete_df = new_filter_k.Filter(outsidePeakL, {"Minv_lambda"});
    std::cout<<"Filter 3"<<std::endl;

    //plot invariant mass of lambda after filtering 
    auto hL_after = new_filter_l.Histo1D({"MinvL_cut","Invariant mass #Lambda cut; M_{inv}(#Lambda)",100,1.11,1.12},"Minv_lambda");
    hL_after->Write();


    // define and plot a histogram which is the armenteros plot
    TH2D *h = new TH2D("AP_2D","Armenteros-Podolanski Plot;#alpha;p_{T}",100,-1.,1.,100,0.,0.25);
    
    auto histo = complete_df.Histo2D(*h,"alpha","pT");

    histo->GetXaxis()->SetTitle("#alpha");
    histo->GetYaxis()->SetTitle("q_{T}");

    histo->Write();


    complete_df.Snapshot("NewVariables", "APPlot_df.root",{"alpha","pT","Minv_lambda","Minv_K0"});

    myFile->Close();

}