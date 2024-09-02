#include <ROOT/RDataFrame.hxx>
#include "armenterosPlot.c"
#include "fitMasses.c"

// do the masscalculation for different pT sampling in one bin to get an uncertainty
void differentPtMasses(TString file="data/SnapshotAP_noMassCut.root", TString tree="NewVariables") {
    ROOT::EnableImplicitMT(10);
    // variables to save used pT
    double pT;
    double low_pT;
    double high_pT;
    // variables to name plots
    TString label;
    int number;
    
    // save results in file
    std::unique_ptr<TFile> myFile(TFile::Open("differentPt_Masses_AP_final.root","RECREATE") );
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
        //auto hK0 = temp_df.Histo1D({"MinvK0","Invariant mass K0; M_{inv}(K_{0}^{s})",200,0.497605,0.497625},"Minv_K0");
        //0.497613,0.497615
        //hK0->Write();
        // plot histogram 
        auto mass = temp_df.Histo1D({"MinvK0","Invariant mass K0; M_{inv}(K_{S}^{0})",125,0.47,0.52},"Minv_K0");
        // use gaus to fit
        /*TF1 *f2 = new TF1("f2", "gaus");
        Int_t b_max = hK0->GetMaximumBin();
        Double_t x_max = hK0->GetBinCenter(b_max);
        Double_t y_max = hK0->GetBinContent(b_max);
        f2->SetParameter(0,y_max);
        f2->SetParameter(1,hK0->GetMean());
        f2->SetParameter(2,hK0->GetStdDev());
        hK0->Fit(f2);
        std::cout<<"Fitted"<<std::endl;
        hK0->Write();*/

        // fit crystalball to data 
        
        auto values=fitMasses(*mass);
        double mean = std::get<0>(values);
        double err = std::get<1>(values);

        positions[i]=mean*1000;
        error[i]=err*1000;
        pTs[i]=pT;
        
        

        /*
        //hK0->Write();
        // use mean-2sigma as a new fitting range
        double mean= f2->GetParameter(1);
        double sigma = f2->GetParameter(2);
        // second fit
        hK0->Fit(f2, "R", "", mean-1.*sigma, mean+1.*sigma);
        std::cout<<"Fitted 2"<<std::endl;
        hK0->Write();
        // save mean and corresponding error for the plot
        // save values of pT and M (position of minimum)
        pTs[i]=pT;
        positions[i]=f2->GetParameter(1)*1000;
        error[i]=f2->GetParError(1)*1000;
        pTs[i]=pT;
        positions[i]=mu*1000;
        // outputs to check
        std::cout<<"pT of bin: "<<pT<<std::endl;
        //std::cout<<"beta of pT bin: "<<beta_pT<<std::endl;
        */

    }
    

    // plot results 
    auto gr = new TGraphErrors(n,pTs,positions,nullptr,error);
    gr->SetTitle("Mean of cb fit");
    gr->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr->GetYaxis()->SetTitle("M [MeV]");
    gr->SetMarkerStyle(kOpenCircle);
    gr->SetMarkerColor(kRed);
    gr->SetLineColor(kRed);
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    
    myFile->cd();
    gr->Write("Result");
    myFile->Close();
}