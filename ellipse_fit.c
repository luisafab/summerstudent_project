#include <ROOT/RDataFrame.hxx>
#include <Math/RootFinder.h>
#include <Math/Functor.h>
// function to calculate ellipse qT for given masses M and m 
double qT(double alpha, double M, double m, double beta_pT){
    // including a factor beta to compensate a lower beta at lower pT
    alpha=alpha*beta_pT;
    return sqrt(pow(M,2.)/4*(1-pow(alpha,2.))-pow(m,2.));
};
// distance of given point x,y to ellipse at point alpha
// euclidean distance: sqrt((x1-x2)^2+(y1-y2)^2)
double distance(double alpha, double M, double m, double x, double y, double beta_pT){
    double qT_alpha=qT(alpha,M,m,beta_pT);
    return sqrt(pow((qT_alpha-y),2.)+pow((alpha-x),2.));
}


// function to calculate the distance to an ellipse with given m for APPlot saved in data frame with number of bins nx and ny
double fit_ellipse (ROOT::RDataFrame df, double M, double nx, double ny, double beta_pT){
    //double M=0.497164;
    // pion mass is fixed for calculation
    double m=0.13957061;
    // plot ellipse the data is fitted to
    // use if wanted
    /*int n=180;
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
    graph->GetHistogram()->SetMaximum(0.25);*/
    //graph->Write(Form("Ellipse_%.6lf",M));

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
    // return result
    return sumOfDistance;
}

// file, tree: name of file and tree with variables calculated for the AP-plot
double ellipse_fit(TString file="APPlot_df_lowPt.root", TString tree="NewVariables",TString output="results_ellipse_fit_finerBinning_alpha09.root", double beta_pT=0.9) {
    double m=0.13957061;
    double M=0.497164;
    
    std::cout<<"calculate qT using beta  = "<<beta_pT<<std::endl;
    // stop time of the iteration
    TStopwatch watch;

    // get DataFrame out of file 
    ROOT::RDataFrame df(tree, file);
    // define histogram which is the armenteros plot
    // with number of bins
    int nx=1500;
    int ny=1000;
    TH2D *h = new TH2D("AP_2D","Armenteros-Podolanski Plot; #alpha; q_{T}; z",nx,-1.,1.,ny,0.,0.25);
    auto histo = df.Histo2D(*h,"alpha","qT");
    //plot
    histo->Write();

    // iterate over n different masses around the k0s mass M
    int n=50;
    std::cout<<"calculate error for "<<n<<" different masses"<<std::endl;
    double dist=0.1/n;
    // save different masses and calculated smallest distance
    double error[n];
    double masses[n];
    for (int i=0; i<n; i++){
        M=0.497164-0.04+i*dist;
        error[i]=fit_ellipse(df,M,nx,n,beta_pT);
        masses[i]=M*1000; // MeV
    }
    // plot for different masses 
    auto gr = new TGraph(n,masses,error);
    gr->SetTitle("Sum of distances for different masses of K^{0}_{s}");
    gr->GetXaxis()->SetTitle("M [MeV]");
    gr->GetYaxis()->SetTitle("Weighted Distance");
    gr->SetLineWidth(0);
    gr->SetMarkerSize(1.5);
    gr->SetMarkerStyle(70);
    gr->GetHistogram()->SetMinimum(-0.02);
    

    std::cout<<"Fit pol2 to area around the minimum"<<std::endl;
    // try to fit pol2 to mass plot
    TF1 *f1 = new TF1("f1", "pol2");
    gr->Fit(f1,"","",470.,520.);
    double p0 = f1->GetParameter(0);
    double p1 = f1->GetParameter(1);
    double p2 = f1->GetParameter(2);
    gr->Write("Result");
    double center=-p1/(2*p2);

    // plot matched ellipse
    // calculate values for corresponding ellipse
    int n2=160;
    double alpha_plot[n2];
    double qT_plot[n2];
    for (int l=0; l<n2; l++){
        double alpha_temp=(-0.8+l*0.01);
        alpha_plot[l]=alpha_temp;
        qT_plot[l]=qT(alpha_temp,center/1000,m,beta_pT);
    }
    auto ellipse = new TGraph(n2,alpha_plot,qT_plot);
    ellipse->SetTitle(Form("Ellipse with #beta = %.3lf",beta_pT));
    ellipse->GetXaxis()->SetTitle("#alpha");
    ellipse->GetYaxis()->SetTitle("q_{T}");
    ellipse->SetLineWidth(3);
    ellipse->SetLineColor(2);
    ellipse->Write("Ellipse");

    std::cout<<"Fitted ellipse to AP-plot with resulting mass M="<<center<<" MeV"<<std::endl;
    
    std::cout<<"Time needed: "<<watch.RealTime()<<std::endl;     // print time 
    
    return center;

}