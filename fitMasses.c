#include <ROOT/RDataFrame.hxx>

//double alpha1, double alpha2, double n1, double n2, double sigma1, double sigma2, double mean
Double_t crystalball_function(double *x, double *par) {

    // evaluate the crystal ball function
    Double_t cb = 0;
    if (par[4] < 0.) return 0.;

    if (par[5] < 0.) return 0.;

    // calculate z for both sides

    double z1 = (x[0] - par[5])/par[4];
    double z2 = (x[0] - par[5])/par[4];

    if (par[0] < 0) z1 = -z1;
    if (par[1] < 0) z2 = -z2;

    double abs_alpha1 = std::abs(par[0]);
    double abs_alpha2 = std::abs(par[1]);

    // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);

    // double D = std::sqrt(M_PI/2.)*(1.+ROOT::math::erf(abs_alpha/std::sqrt(2.)));

    // double N = 1./(sigma*(C+D));

    if (z1 < - abs_alpha1){

        double nDivAlpha = par[2]/abs_alpha1;

        double AA = std::exp(-0.5*abs_alpha1*abs_alpha1);

        double B = nDivAlpha -abs_alpha1;

        double arg = nDivAlpha/(B-z1);

        cb=  par[6]*(AA * std::pow(arg,par[2]));
    }
    else if (z1<=abs_alpha2){
        cb=  par[6]*std::exp(- 0.5 * z1 * z1);
    }
    else if (z2<=0){
        cb=  par[6]*std::exp(- 0.5 * z1 * z1);
    }
    else {

        //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);

        double nDivAlpha = par[3]/abs_alpha2;

        double AA = std::exp(-0.5*abs_alpha2*abs_alpha2);

        double B = nDivAlpha -abs_alpha2;

        double arg = nDivAlpha/(B+z2);

        cb=  par[6]*(AA * std::pow(arg,par[3]));

    }
    return cb + par[7]*TMath::Exp(par[8]*x[0]);

}

Double_t exponential(double *x, double *par){
    return par[0]*TMath::Exp(par[1]*x[0]);
}

std::tuple<double,double> fitMasses(TH1D mass){
    ROOT::EnableImplicitMT(3);
    //TH1D *h = new TH1D("MinvK0", "Minv", 125, 0.45, 0.52);
    //auto mass = df.Histo1D(*h,"Minv_K0");
    // define functions
    TF1 *f_cb = new TF1("crystalball", crystalball_function, 0.47, 0.52, 9 );
    TF1 *f_exp = new TF1("exp", exponential, 0.47, 0.52, 2);

    // swt Parameters of the functions
    //f_exp->SetParameters(10, -1000.);
    f_cb->SetParLimits(0,1.0,2.5);
    f_cb->SetParLimits(1,1.0,2.5);
    f_cb->SetParLimits(2,0.1,20.);
    f_cb->SetParLimits(3,0.1,30.);
    f_cb->SetParLimits(4,0.0025,0.007);
    f_cb->SetParLimits(5,0.49,0.51);
    f_cb->SetParLimits(6,1.,1e7);
    f_cb->SetParLimits(7,1., 1e6);
    f_cb->SetParLimits(8,-15.,10.);

   
    mass.Fit("crystalball");
    f_exp->SetParameter(0,f_cb->GetParameter(7));
    f_exp->SetParameter(1,f_cb->GetParameter(8));
    f_exp->Write();
    mass.Write();
    
    double mean=f_cb->GetParameter(5);
    double error=f_cb->GetParError(5);

    return std::make_tuple(mean, error);

    

}