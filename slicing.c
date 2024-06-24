#include <ROOT/RDataFrame.hxx>
// function to fit
double pT_fit (double* alpha, double* p) {
    // fix mass of k0 short or pi
    double M=0.497164;
    double m=0.13957061;
    double pT2=pow(M,2.)/4*(1-pow(alpha[0],2.))-pow(p[0],2.);
    //double pT2=pow(p[0],2.)/4*(1-pow(alpha[0],2.))-pow(m,2.);
    return sqrt(pT2);
};
// calculates pT in every bin of alpha
// file, tree: name of file and tree with variables calculated for the AP-plot
void slicing(TString file="APPlot_df.root", TString tree="NewVariables") {

    // saving plots in a root file
    std::unique_ptr<TFile> myFile( TFile::Open("results_slicing_fit.root", "RECREATE") );

    // get DataFrame out of file 
    ROOT::RDataFrame df(tree, file);
    // define histogram which is the armenteros plot
    TH2D *h = new TH2D("AP_2D","Armenteros-Podolanski Plot; #alpha; q_{T}; z",100,-1.,1.,1000,0.,0.25);
    auto histo = df.Histo2D(*h,"alpha","pT");
    // save number of bins
    int n=histo->GetNbinsX();
    // arrays to save values to plot
    double x[n];
    double xerr[n];
    double plotvalues[2][n];
    // Fit to bins 
    std::cout<<"Start slicing"<<std::endl;
    // create subdirectory to save files 
    myFile->cd();
    gDirectory->mkdir("Projections");
    myFile->cd();
    myFile->cd("Projections");
    for (int i = 1; i <=n; i+=1) {
        // calculate x
        // use part of histogram the histogram which corresponds to one bin
        std::cout<<"Iteration "<<i<<std::endl;

        double bin_low=histo->GetXaxis()->GetBinLowEdge(i);
        double bin_width=histo->GetXaxis()->GetBinWidth(i);
        double bin_up=bin_low+bin_width;

        std::cout<<"Lower edge: "<<bin_low<<"; Upper edge: "<<bin_low<<"; Width: "<<bin_width<<std::endl;

        // calculate values to plot for alpha 
        // middle of the bin
        x[i-1]=bin_low+bin_width/2;
        std::cout<<x[i-1]<<std::endl;
        // error is the distance to next point 
        xerr[i-1]=bin_width/2;

        // project to y axis
        TH1D *histo_bin = histo->ProjectionY(Form("Projection_%.2lf_%lf.2",bin_low,bin_up),i,i);
        std::cout<<"Created projection"<<std::endl;

        histo_bin->SetTitle(Form("#alpha=%.2lf",x[i-1]));
        histo_bin->GetXaxis()->SetTitle("q_{T} [GeV/c]");
        histo_bin->GetYaxis()->SetTitle("N");

        // use gaus to fit
        TF1 *f2 = new TF1("f2", "gaus");
        histo_bin->Fit(f2);
        std::cout<<"Fitted"<<std::endl;
        histo_bin->Write();

        /*histo_bin->Write();*/
        // use mean-2sigma as a new fitting range
        double mean= f2->GetParameter(1);
        double sigma = f2->GetParameter(2);
        // second fit
        histo_bin->Fit(f2, "R", "", mean-2*sigma, mean+2*sigma);
        std::cout<<"Fitted 2"<<std::endl;
        histo_bin->Write();
        // save mean and corresponding error for the plot
        plotvalues[0][i-1]=f2->GetParameter(1);
        plotvalues[1][i-1]=f2->GetParError(1);

    }
    // plot the values
    myFile->cd();
    auto gr = new TGraphErrors(n,x,plotvalues[0],xerr,plotvalues[1]);
    gr->SetTitle("Armenteros Podolanski plot from fitting to bins");
    gr->GetXaxis()->SetTitle("#alpha");
    gr->GetYaxis()->SetTitle("q_{T} [GeV/c]");
    gr->Write("Calculated values");


    // add fitting to ellipse
    // use mass of pion/ k0s to initialize parameter
    double m= 0.13957061;
    double M=0.497164;
    auto *func = new TF1("func",pT_fit,-1,1,1);
    func->SetParameters(m);
    // fitting: does not work for larger range
    gr->Fit(func,"","",-0.8,0.8);
    // Plot the graph
    gr->SetTitle("Ellipse fitted to calculated values");
    gr->GetYaxis()->SetRangeUser(0.,0.21);
    gr->Write("Fitted graph");

    myFile->Close();
    

}

