#include <ROOT/RDataFrame.hxx>
#include "ellipse_fit.c"
#include "armenterosPlot.c"
#include "plot_ellipse.c"
#include "beta_calculation.c"

// do the calculation of an armenterosPlot for different pT 
void differentPt(TString file="data/AnalysisResults_treesAP_data_LHC22o_apass6_small.root", TString tree="O2v0tableap") {
    double m=0.13957061;
    // variables to save data in files
    TString save_df;
    TString save_file;
    TString save_ellipse;
    // variables to save used pT
    double pT;
    double low_pT;
    double high_pT;
    // variables to name plots
    TString label;
    int number;
    
    // save results in file
    std::unique_ptr<TFile> myFile(TFile::Open("DifferentPt_finerBinning.root","RECREATE") );
    std::cout<<"opened file"<<std::endl;

    // set values for the iterations/ binning of pT
    // number of iterations with smaller/larger bins
    // n=n_small+n_large
    int n_small=24; //24
    int n_large=12; //12
    // total number of iterations
    int n=n_small+n_large;
    // set values of pT
    double pT_start=0.3;
    double pT_middle=1.5; 
    double pT_end=4.0;
    // calculate stepsize and startvalue for different binnings
    double stepsize_small=(pT_middle-pT_start)/n_small;
    double stepsize_large=(pT_end-pT_middle)/n_large;
    double startvalue_small=pT_start+stepsize_small/2.;
    double startvalue_large=pT_middle+stepsize_large/2.;

    // variables in loop
    double stepsize;
    double start_pT;
    double factor;
    double M_temp;

    // variables to save results
    double position;
    double positions[n];
    double pTs[n];
    double beta_pT;

    // do the calculation for two different binnings (low pT has more change in beta)
    for (int i=0; i<n; i++){
        // fix stepsize and startvalue according to binning
        if (i<n_small){
            stepsize=stepsize_small;
            start_pT=startvalue_small;
            factor=i;
        }else{
            stepsize=stepsize_large;
            start_pT=startvalue_large;
            factor=i-n_small;
        }
        // calculate values of pT
        pT=start_pT+stepsize*factor;
        low_pT=pT-stepsize/2;
        high_pT=pT+stepsize/2;

        // outputs to check
        std::cout<<"calculation for pT bin from "<<low_pT<<" to "<<high_pT<<std::endl;

        // create name to save file
        number = 1000*pT;
        label=Form("pT_%i",number);
        std::cout<<label<<std::endl;
        // create directon for this pT value
        myFile->cd();
        gDirectory->mkdir("DifferentPt/"+label);
        // and one subdirectory to save different ellipses
        myFile->cd();
        gDirectory->mkdir("DifferentPt/"+label+"/differentEllipse");

        // files to save 
        save_df="DifferentPt/APPlot_df_"+label+".root";
        save_file="DifferentPt/APPlot_"+label+".root";
        save_ellipse="DifferentPt/ellipse_fit_"+label+".root";

        // call armenterosPlot to save df for corresponding pT values
        myFile->cd();
        myFile->cd("DifferentPt/"+label);
        armenterosPlot(file, tree, false,save_file,save_df,low_pT,high_pT);
        // get data frame 
        ROOT::RDataFrame df("NewVariables", save_df);
        // calculate beta for corresponding pT
        myFile->cd();
        myFile->cd("DifferentPt/"+label);
        beta_pT = beta_calculation(df, low_pT, high_pT);
        // calculate ellipse for cooresponding pT
        myFile->cd();
        myFile->cd("DifferentPt/"+label);
        position=ellipse_fit(save_df, "NewVariables",save_ellipse, beta_pT);
        // plot different ellipses around the miimum to compare
        myFile->cd();
        myFile->cd("DifferentPt/"+label+"/differentEllipse");
        for (int i=0;i<17;i++){
            M_temp=position/1000.-0.008+i*0.001;
            plot_ellipse(M_temp,m,beta_pT);
        }
        // save values of pT and M (position of minimum)
        pTs[i]=pT;
        positions[i]=position;
        // outputs to check
        std::cout<<"pT of bin: "<<pT<<std::endl;
        std::cout<<"beta of pT bin: "<<beta_pT<<std::endl;
        std::cout<<"position of minimum: "<<position<<std::endl;

    }
    

    // plot results 
    auto gr = new TGraph(n,pTs,positions);
    gr->SetTitle("Position of minima");
    gr->GetXaxis()->SetTitle("pT [GeV]");
    gr->GetYaxis()->SetTitle("M [MeV]");
    gr->SetLineWidth(0);
    gr->SetMarkerSize(1.5);
    gr->SetMarkerStyle(70);
    
    myFile->cd();
    gr->Write("Result");
    myFile->Close();
}