#include <ROOT/RDataFrame.hxx>
#include "ellipse_fit.c"
#include "armenterosPlot.c"
#include "plot_ellipse.c"
#include "beta_calculation.c"
#include <tuple>


double fit_ellipse_sample (TH3D* histo, double M, int nx, int ny, int bin_z, double beta_pT){
    //double M=0.497164;
    // pion mass is fixed for calculation
    double m=0.13957061;

    // use given binning
    // variable to save the sum of distances
    double sumOfDistance=0;
    // iterate over the twenty samples
    // iterate over every bin to calculate smallest distance
    for (int i = 0; i <nx; i++) {   
        for (int j = 0; j <ny; j++) {
            double content=histo->GetBinContent(i,j,bin_z);
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
            //gr->SetBinContent(i,j,smallestdistance);
            //wd->SetBinContent(i,j,smallestsquared);
            // add up
            sumOfDistance+=smallestsquared;
        }

    }
    // return result
    return sumOfDistance;
}




std::tuple<double,double> samples(ROOT::RDataFrame df_all, double low_pT, double high_pT) {
    //std::unique_ptr<TFile> myFile( TFile::Open("test.root", "RECREATE") );
    // create df 
    //ROOT::RDataFrame df(tree, file);
    auto insidePtRange =[&](float pT){
        bool isInside;
        if (pT>low_pT && pT < high_pT){
            isInside=true;
        }else{
            isInside=false;
        }
        return isInside;
    };
    auto df = df_all.Filter(insidePtRange,{"pTV0"});
    auto beta_pT_ptr=df.Mean("beta");
    auto beta_pT=*beta_pT_ptr;
    // want to loop over samples and pT and calculate a AP plot for every bin
    int Nsamples = 50;
    int nx=200;
    int ny=200;
    int Nmasses=20;
    // use a TH3D
    TH3D *histogram = new TH3D("AP_3D", "Armenteros-Podolanski Plot;#alpha;q_{T};rnd",nx,-1.,1.,ny,0.,0.25,Nsamples,0.,1.);
    auto histo = df.Histo3D(*histogram,"alpha","qT","rnd");

    histo->GetXaxis()->SetTitle("#alpha");
    histo->GetYaxis()->SetTitle("q_{T}");
    histo->GetZaxis()->SetTitle("rnd");

    histo->Write();

    TH2D *hAP = new TH2D("AP_2D", "Armenteros-Podolanski Plot;#alpha;q_{T}",nx,-1.,1.,ny,0.,0.25);
    auto histAP = df.Histo2D(*hAP,"alpha","qT");

    histAP->GetXaxis()->SetTitle("#alpha");
    histAP->GetYaxis()->SetTitle("q_{T}");

    histAP->Write();

    // histograms to save distance and distance squared ( to check on)
    // histograms to save results
    // calculated distance plus distance weighted by entries for every bin
    //TH3D *gr = new TH2D(Form("Distance_%.6lf",M),"Distance to Ellipse; #alpha; q_{T}; sample",nx,-1.,1.,ny,0.,0.25,Nsamples,0.,1.);
    //TH3D *wd = new TH2D(Form("Weighted Distance_%.6lf",M),"Weighted Distance to Ellipse; #alpha; q_{T}; sample",nx,-1.,1.,ny,0.,0.25,Nsamples,0.,1.);
    
    double calculatedMass[Nsamples];
    double centers[Nsamples];
    // create histogram to plot masses 
    TH1D *h1 = new TH1D("h1", "Masses for different samples; M [MeV]", 400, 496, 498);
    // iterate over samples
    for (int s=1; s<=Nsamples; s++){
        // iterate over different masses 
        double error[Nmasses];
        double masses[Nmasses];
        double dist=0.06/Nmasses;
        for (int i=0; i<Nmasses; i++){
            float M=kK0sMass-0.03+i*dist;
            error[i]=fit_ellipse_sample(histo.ROOT::RDF::RResultPtr<TH3D>::GetPtr(),M,nx,ny,s,beta_pT);
            masses[i]=M*1000; // MeV
        }
        // plot for different masses 
        auto gr = new TGraph(Nmasses,masses,error);
        gr->SetTitle("Sum of distances for different masses of K^{0}_{s}");
        gr->GetXaxis()->SetTitle("M [MeV]");
        gr->GetYaxis()->SetTitle("#Chi^{2}");
        gr->SetLineWidth(0);
        gr->SetMarkerSize(1.5);
        gr->SetMarkerStyle(70);
        gr->GetHistogram()->SetMinimum(-0.02);
        std::cout<<"Fit pol2 to area around the minimum"<<std::endl;
        // try to fit pol2 to mass plot
        TF1 *f1 = new TF1("f1", "pol2");
        gr->Fit(f1,"","",467.,524.5);
        double p0 = f1->GetParameter(0);
        double p1 = f1->GetParameter(1);
        double p2 = f1->GetParameter(2);
        gr->Write("Result");
        double center=-p1/(2*p2);
        std::cout<<"Fitted ellipse to AP-plot with resulting mass M="<<center<<" MeV"<<std::endl;
        centers[s-1]=center;
        h1->Fill(center);
        
    }
    // calculate mu and sigma for the number of samples
    auto gr = new TGraph(Nsamples,centers);
    
    gr->GetXaxis()->SetTitle("Sample");
    gr->GetYaxis()->SetTitle("M [MeV]");
    gr->SetLineWidth(0);
    gr->SetMarkerSize(1.5);
    gr->SetMarkerStyle(70);
    double musum=0;
    for (int i=0; i<Nsamples;i++){
        musum+=centers[i];
    }
    double mu = musum/Nsamples;
    double sigmasum=0;
    for (int i=0; i<Nsamples;i++){
        sigmasum+=std::abs(centers[i]-mu);
    }
    double sigma = std::sqrt(sigmasum/(Nsamples*(Nsamples-1)));
    std::cout<<(Nsamples*(Nsamples))<<std::endl;
    std::cout<<mu<<" "<<sigma<<std::endl;
    gr->SetTitle(Form("M=%.3lf#pm%.3lf",mu,sigma));
    gr->Write("mu");
    h1->SetTitle(Form("M=%.3lf#pm%.3lf",mu,sigma));
    h1->Write();
    // plot resulting ellipse and summed AP Plot
    int n2=160;
    double alpha_plot[n2];
    double qT_plot[n2];
    for (int l=0; l<n2; l++){
        double alpha_temp=(-0.8+l*0.01);
        alpha_plot[l]=alpha_temp;
        qT_plot[l]=qT(alpha_temp,mu/1000,kPiMass,beta_pT);
    }
    auto ellipse = new TGraph(n2,alpha_plot,qT_plot);
    ellipse->SetTitle(Form("Ellipse with #beta = %.3lf",beta_pT));
    ellipse->GetXaxis()->SetTitle("#alpha");
    ellipse->GetYaxis()->SetTitle("q_{T}");
    ellipse->SetLineWidth(3);
    ellipse->SetLineColor(2);
    ellipse->Write("Ellipse");
    

    //myFile->Close();

    return std::make_tuple(mu, sigma);

}