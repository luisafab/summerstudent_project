#include <ROOT/RDataFrame.hxx>

//calculating the armenteros-podolanski plot for given data 
// string file : name of file
// string tree : name of tree
// bool MC : true : only lambda events are consicered (!! The MC Information in the file is necessary !!)
//           false : all events are considered

void armenterosPlot_fit(TString file="data/AnalysisResults_treesAP_data.root", TString tree="O2v0tableap",bool MC=false) {


    // define a data frame to use 
    ROOT::RDataFrame df_MC(tree, file);
    // filter it
    //     for MC==true only the recombined (isreco==true) lambdas (PDGCODE=3122) are considered
    //     else there is a dummy cut to use the same dataframe
    auto df = df_MC.Filter(MC? "(fPDGCode == 3122 || fPDGCode== -3122) && fIsReco" : "fLen>0");

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
    auto new_column = df.Define("alpha", alpha, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    // get the final data frame
    auto complete_df = new_column.Define("pT", pT, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});

    // plot the two new variables
    auto h1 = complete_df.Histo1D("alpha");
    TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
    h1->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    auto h2 = complete_df.Histo1D("pT");
    TCanvas * c2 = new TCanvas("c2", "c2", 600, 600);
    h2->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    // define and plot a histogram which is the armenteros plot
    TH2D *h = new TH2D("h","Armenteros-Podolanski Plot",100,-1.,1.,40,0.,0.2);
    
    auto histo = complete_df.Histo2D(*h,"alpha","pT");

    histo->GetXaxis()->SetTitle("alpha");
    histo->GetYaxis()->SetTitle("pT");

    TCanvas * c3 = new TCanvas("c3", "c3", 900, 600);
    histo->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

}