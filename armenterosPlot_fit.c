#include <ROOT/RDataFrame.hxx>

//calculating the armenteros-podolanski plot for given data 
// string file : name of file
// string tree : name of tree
// bool MC : true : only lambda events are consicered (!! The MC Information in the file is necessary !!)
//           false : all events are considered

void armenterosPlot_fit(TString file="data/AnalysisResults_treesAP_data.root", TString tree="O2v0tableap",bool MC=false) {

    float massProton=0.9382720813;
    float massPion= 0.13957061;

    // define a data frame to use 
    ROOT::RDataFrame df_MC(tree, file);
    // filter it
    //     for MC==true only the recombined (isreco==true) lambdas (PDGCODE=3122) are considered
    //     else there is a dummy cut to use the same dataframe
    auto df = df_MC.Filter(MC? "(fPDGCode == 3122 || fPDGCode== -3122) && fIsReco" : "fLen>0");
    std::cout<<"filtered PDG"<<std::endl;

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

    // use mass to filter the data
    // define values to cut on both peaks 
    // lambda peak should be      1.1115683
    // K0 short peak should be    0.497164
    // 1.107 1.116 0.493 0.5012
    // 1.108 1.115 0.494 0.5011
    // 1.101 1.120 0.4925 0.5017 
    // 1.101 1.120 0.485 0.505
    // 1.112 1.120 0.485 0.505
    // 1.101 1.119 0.488 0.505
    // 1.112 1.119 0.488 0.505
    float lambda_low=1.112;
    float lambda_high=1.119;
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


    // plot the different variables
    auto h1 = complete_df.Histo1D("alpha");
    TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
    h1->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    auto h2 = complete_df.Histo1D("pT");
    TCanvas * c2 = new TCanvas("c2", "c2", 600, 600);
    h2->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    auto h3 = complete_df.Histo1D("Minv_K0");
    TCanvas * c3 = new TCanvas("c3", "c3", 600, 600);
    h3->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    auto h4 = complete_df.Histo1D("Minv_lambda");
    TCanvas * c4 = new TCanvas("c4", "c4", 600, 600);
    h4->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    auto h5 = new_filter_l.Histo1D("Minv_lambda");
    TCanvas * c5 = new TCanvas("c5", "c5", 600, 600);
    h5->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    
    // define and plot a histogram which is the armenteros plot
    std::cout<<"start AP plot"<<std::endl;
    TH2D *h = new TH2D("h","Armenteros-Podolanski Plot",100,-1.,1.,40,0.,0.25);
    
    auto histo = complete_df.Histo2D(*h,"alpha","pT");
    std::cout<<"plotted histogramm"<<std::endl;

    histo->GetXaxis()->SetTitle("alpha");
    histo->GetYaxis()->SetTitle("pT");

    TCanvas * c = new TCanvas("c", "c", 900, 600);
    histo->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    // this one should only contain the lambdas
    std::cout<<"start AP plot 2"<<std::endl;
    TH2D *h0 = new TH2D("h","Armenteros-Podolanski Plot",100,-1.,1.,40,0.,0.2);
    
    auto histo1 = new_filter_l.Histo2D(*h0,"alpha","pT");
    std::cout<<"Plotted"<<std::endl;

    histo->GetXaxis()->SetTitle("alpha");
    histo->GetYaxis()->SetTitle("pT");

    TCanvas * c0 = new TCanvas("c0", "c0", 900, 600);
    histo1->DrawClone();
    gROOT->GetListOfCanvases()->Draw();
    
    // fit gaussian to bins of alpha

    // number of bins to fit
    const Int_t n = 20;
    // and resulting stepsize
    float stepsize=2./n;

    double x[n];
    double xerr[n];

    double plotvalues[2][n];

    std::cout<<"Start slicing"<<std::endl;
    for (int i = 0; i < n; i+=1) {
        // calculate x
        // use part of histogram the histogram which corresponds to one bin
        std::cout<<"Start plotting"<<std::endl;
        TH2D *h_temp = new TH2D("h","Armenteros-Podolanski Plot",5,-1+i*0.1,-1+(i+1)*0.1,40,0.,0.25);
        auto histo_temp = complete_df.Histo2D(*h_temp,"alpha","pT");
        std::cout<<"created histogramm"<<std::endl;

        // project to y axis
        TH1D *histo_bin = histo_temp->ProjectionY("projection",1,5);
        std::cout<<"Created projection"<<std::endl;

        // use gaus to fit
        TF1 *f2 = new TF1("f2", "gaus");
        histo_bin->Fit(f2);
        std::cout<<"Fitted"<<std::endl;
        // use mean-2sigma as a new fitting range
        double mean=f2->GetParameter(1);
        double sigma = f2->GetParameter(2);
        // second fit
        histo_bin->Fit(f2, "R", "", mean-2*sigma, mean+2*sigma);
        std::cout<<"Fitted 2"<<std::endl;
        // save mean and corresponding error for the plot
        plotvalues[0][i]=f2->GetParameter(1);
        plotvalues[1][i]=f2->GetParError(1);

        // calculate values to plot for alpha 
        // middle of the bin
        x[i]=-1+i*0.1+0.05;
        // error is the distance to next point 
        xerr[i]=0.05; 
    }

    // plot the values
    auto c20 = new TCanvas("c20","AP Plot",900,600); 
    auto gr = new TGraphErrors(20,x,plotvalues[0],xerr,plotvalues[1]);
    gr->SetTitle("Armenteros Podolanski plot from fitting to bins");
    gr->Draw("AP");

}