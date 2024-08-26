#include <ROOT/RDataFrame.hxx>


#include "armenterosPlot.c"

// do the masscalculation for different pT sampling in one bin to get an uncertainty
void differentPtMasses(TString file="SnapshotMC_genDecay.root", TString tree="NewVariables") {
    ROOT::EnableImplicitMT(4);
    // variables to save used pT
    double pT;
    double low_pT;
    double high_pT;
    // variables to name plots
    TString label;
    int number;
    
    // save results in file
    std::unique_ptr<TFile> myFile(TFile::Open("differentPt_Masses_genDecay.root","RECREATE") );
    std::cout<<"opened file"<<std::endl;

    // set values for the iterations/ binning of pT
    // number of iterations with smaller/larger bins
    int n_small=26; //24
    int n_large=24; //12
    // total number of iterations
    int n=n_small+n_large;
    // set values of pT
    double pT_start=0.3;
    double pT_middle=1.6; 
    double pT_end=4.0;
    // calculate stepsize and startvalue for different binnings
    double stepsize_small=(pT_middle-pT_start)/n_small;
    double stepsize_large=(pT_end-pT_middle)/n_large;
    std::cout<<stepsize_large<<std::endl;
    double startvalue_small=pT_start+stepsize_small/2.;
    double startvalue_large=pT_middle+stepsize_large/2.;
    std::cout<<startvalue_large<<std::endl;

    // variables in loop
    double stepsize;
    double start_pT;
    double factor;
    double M_temp;

    // variables to save results
    double position;
    double positions[n];
    double pTs[n];
    double error[n];
    // get df from snapshot 
    ROOT::RDataFrame df(tree, file);
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
        pT=start_pT+stepsize/2+factor*stepsize;
        low_pT=pT-stepsize/2;
        high_pT=pT+stepsize/2;
        std::cout<<pT<<std::endl;

        // outputs to check
        std::cout<<"calculation for pT bin from "<<low_pT<<" to "<<high_pT<<std::endl;

        // create name to save file
        number = 1000*pT;
        label=Form("pT_%i",number);
        std::cout<<label<<std::endl;
        // create directon for this pT value
        myFile->cd();
        gDirectory->mkdir("DifferentPt/"+label);

        auto insidePtRange =[&](float pT){
        bool isInside;
        if (pT>low_pT && pT < high_pT){
            isInside=true;
        }else{
            isInside=false;
        }
        return isInside;
        };
        myFile->cd();
        myFile->cd("DifferentPt/"+label);
        auto temp_df = df.Filter(insidePtRange,{"pTV0"});
        double mu = *temp_df.Mean("Minv_K0");
        // plot histogram 
        auto hK0 = temp_df.Histo1D({"MinvK0","Invariant mass K0; M_{inv}(K_{0}^{s})",100,0.497,0.49766},"Minv_K0");
        hK0->Write();
        // save values of pT and M (position of minimum)
        pTs[i]=pT;
        positions[i]=mu;
        // outputs to check
        std::cout<<"pT of bin: "<<pT<<std::endl;
        //std::cout<<"beta of pT bin: "<<beta_pT<<std::endl;

    }
    

    // plot results 
    auto gr = new TGraph(n,pTs,positions);
    gr->SetTitle("Mean of M_{inv}(K^{0})");
    gr->GetXaxis()->SetTitle("pT [GeV]");
    gr->GetYaxis()->SetTitle("M [MeV]");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    
    myFile->cd();
    gr->Write("Result");
    myFile->Close();
}