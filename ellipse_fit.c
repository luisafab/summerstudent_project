#include <ROOT/RDataFrame.hxx>
#include <Math/RootFinder.h>
#include <Math/Functor.h>
// function to calculate ellipse qT for given masses M and m 
double qT(double alpha, double M, double m, double beta){
    // including a factor beta to compensate a lower beta at lower pT
    alpha=alpha*beta;
    return sqrt(pow(M,2.)/4*(1-pow(alpha,2.))-pow(m,2.));
};
// distance of given point x,y to ellipse at point alpha
// euclidean distance: sqrt((x1-x2)^2+(y1-y2)^2)
double distance(double alpha, double M, double m, double x, double y, double beta){
    double qT_alpha=qT(alpha,M,m,beta);
    return sqrt(pow((qT_alpha-y),2.)+pow((alpha-x),2.));
}
// calculate beta for different pT (pT=m*gamma*beta)
double beta(double pT){
    double M=0.497164;
    return sqrt(1./(pow((M/pT),2.)+1.));
}

// function to calculate the distance to an ellipse with given m for APPlot saved in data frame with number of bins nx and ny
double fit_ellipse (ROOT::RDataFrame df, double M, double nx, double ny, double beta_pT){
    //double M=0.497164;
    // pion mass is fixed for calculation
    double m=0.13957061;
    //plot ellipse the data is fitted to
    int n=180;
    double alpha_plot[n];
    double qT_plot[n];
    for (int l=0; l<n; l++){
        double alpha_temp=(-0.9+l*0.01);
        alpha_plot[l]=alpha_temp;
        qT_plot[l]=qT(alpha_temp,M,m,beta_pT);
    }
    auto graph = new TGraph(n,alpha_plot,qT_plot);
    graph->SetTitle("Ellipse fitted to data");
    graph->GetXaxis()->SetTitle("#alpha");
    graph->GetYaxis()->SetTitle("q_{T}");
    graph->SetLineWidth(1);
    TAxis *axis = graph->GetXaxis();
    axis->SetLimits(-1.,1.);
    graph->GetHistogram()->SetMinimum(-0.02);
    graph->GetHistogram()->SetMaximum(0.25);
    graph->Write(Form("Ellipse_%.6lf",M));

    // Histogram for AP Plot
    // use given binning
    TH2D *h = new TH2D("AP_2D","Armenteros-Podolanski Plot; #alpha; q_{T}; z",nx,-1.,1.,ny,0.,0.25);
    auto histo = df.Histo2D(*h,"alpha","qT");
    // histograms to save results
    // calculated distance plus distance weighted by entries for every bin
    TH2D *gr = new TH2D(Form("Distance_%.6lf",M),"Distance to Ellipse; #alpha; q_{T}; z",nx,-1.,1.,ny,0.,0.25);
    TH2D *wd = new TH2D(Form("Weighted Distance_%.6lf",M),"Weighted Distance to Ellipse; #alpha; q_{T}; z",nx,-1.,1.,ny,0.,0.25);
    // variable to save the sum of distances
    double sumOfDistance=0.;
    // iterate over every bin to calculate smallest distance
    for (int i = 0; i <nx; i++) {
        for (int j = 0; j <ny; j++) {
            double content=histo->GetBinContent(i,j);
            // skip if bin is empty
            if (content==0) continue;
            // calculate position and width of bin
            double bin_x_low=histo->GetXaxis()->GetBinLowEdge(i);
            double bin_x_width=histo->GetXaxis()->GetBinWidth(i);
            double x2=bin_x_low+bin_x_width/2;
            double bin_y_low=histo->GetYaxis()->GetBinLowEdge(j);
            double bin_y_width=histo->GetYaxis()->GetBinWidth(j);
            double y2=bin_y_low+bin_y_width/2;
            
            // variables to save smallest distance and corresponding alpha and qt
            double smallestdistance=100;
            double alpha_s;
            double qT_s;
            double alpha_test=0.;
            // go throug points in alpha 
            // for points with alpha (x2) < 0 only on the left side of the ellipse 
            if (x2<0){
                alpha_test=-0.9;
            }
            
            double distance_test;
            // First rough iteration to find region of point
            for (int k = 0; k < 10; k++) {
                alpha_test+=0.1*k;
                distance_test=distance(alpha_test,M,m,x2,y2,beta_pT);
                if (distance_test<smallestdistance){
                    smallestdistance=distance_test;
                    qT_s=qT(alpha_test,M,m,beta_pT);
                    alpha_s=alpha_test;
                }
            }
            // second iteration around the chosen point
            alpha_test=alpha_s-0.2;
            for (int k = 0; k < 10000; k++) {
                alpha_test+=0.00004*k;
                distance_test=distance(alpha_test,M,m,x2,y2,beta_pT);
                if (distance_test<smallestdistance){
                    smallestdistance=distance_test;
                    qT_s=qT(alpha_test,M,m,beta_pT);
                    alpha_s=alpha_test;
                }
            }
            // save distance
            // use |xi-xellipse|^2*Ni
            double smallestsquared=pow(smallestdistance,2.)*content;
            // save results to histograms
            gr->SetBinContent(i,j,smallestdistance);
            wd->SetBinContent(i,j,smallestsquared);
            // add up
            sumOfDistance+=smallestsquared;
        }
    }
    // plot results
    gr->Write();
    wd->Write();
    //std::cout<<sumOfDistance<<std::endl;
    return sumOfDistance;
}

// file, tree: name of file and tree with variables calculated for the AP-plot
void ellipse_fit(TString file="APPlot_df_lowPt.root", TString tree="NewVariables",TString output="results_ellipse_fit_finerBinning_alpha09.root", double pT=4.) {
    //double m=0.13957061;
    double M=0.497164;
    // calculate beta for given pT to adapt in qT calculation
    double beta_pT=beta(pT);
    std::cout<<"calculate qT using beta (pT="<<pT<<") = "<<beta_pT<<std::endl;
    // stop time of the iteration
    TStopwatch watch;
    
    // saving plots in a root file
    std::unique_ptr<TFile> myFile( TFile::Open(output,"RECREATE") );
    // create subdirectory to save files 
    myFile->cd();
    gDirectory->mkdir("DifferentMasses");

    // get DataFrame out of file 
    ROOT::RDataFrame df(tree, file);
    // define histogram which is the armenteros plot
    // with number of bins
    int nx=1500;
    int ny=1000;
    TH2D *h = new TH2D("AP_2D","Armenteros-Podolanski Plot; #alpha; q_{T}; z",nx,-1.,1.,ny,0.,0.25);
    auto histo = df.Histo2D(*h,"alpha","qT");
    // plot
    histo->Write();
    // iterate over n different masses around the k0s mass M
    int n=100;
    double dist=0.2/n;
    // save different masses and calculated smallest distance
    double error[n];
    double masses[n];
    for (int i=0; i<n; i++){
        myFile->cd();
        myFile->cd("DifferentMasses");
        M=0.497164-0.08+i*dist;
        error[i]=fit_ellipse(df,M,nx,n,beta_pT);
        masses[i]=M*1000; // MeV
    }
    // plot for different masses 
    myFile->cd();
    auto gr = new TGraph(n,masses,error);
    gr->SetTitle("Sum of distances for different masses of K^{0}_{s}");
    gr->GetXaxis()->SetTitle("M [MeV]");
    gr->GetYaxis()->SetTitle("Weighted Distance");
    gr->SetLineWidth(0);
    gr->SetMarkerSize(1.5);
    gr->SetMarkerStyle(70);
    gr->GetHistogram()->SetMinimum(-0.02);
    gr->Write("Result");
    
    myFile->Close();    // close file

    std::cout<<watch.RealTime()<<std::endl;     // print time 
    

}