#include <ROOT/RDataFrame.hxx>
#include "armenterosPlot.c"

float deltapT(float M, float pT1, float pT2, float E1, float E2, float p1squared, float p2squared, float p1p2){
    float deltaM= M-kK0sMass;
    return deltaM*2*M*std::sqrt(E1*E2 - p1p2)/((p1squared*E2/(pT1*E1))-(p1p2)/pT1+(p2squared*E1/(pT2*E2))-(p1p2)/pT2);
}

/*float deltapTone(float M, float pT1, float E1, float E2, float p1squared, float p1p2){
    float deltaM= M-kK0sMass;
    return deltaM*M/((2*p1squared*E2/(pT1*E1))-(p1p2)/pT1);
}*/

float deltapTone(float M, float pT1, float E1, float E2, float p1squared, float p1p2){
    float deltaM= M-kK0sMass;
    return deltaM*2*M*std::sqrt(E1*E2 - p1p2)/((p1squared*E2/(pT1*E1))-(p1p2)/pT1);
}

float deltapTonePDG(float M, float pT1, float E1, float E2, float p1squared, float p1p2){
    float deltaM= M-kK0sMass;
    return deltaM*kK0sMass/((2*p1squared*E2/(pT1*E1))-(p1p2)/pT1);
}

float LHCb(float M, float E1, float E2){
    return (M-kK0sMass)*2*kK0sMass/(std::pow(kK0sMass,2.)-2*kPiMass*(1+E2/E1));
}

float pT(float px, float py){

    return std::hypot(px,py);
}

float psquare(float px, float py, float pz){
    return std::pow(px,2.) + std::pow(py,2.) + std::pow(pz,2.);
}

float p1p2(float px1, float py1, float pz1, float px2, float py2, float pz2){
    return px1*px2 + py1*py2 + pz1*pz2;
}

float energy(float psquare){
    return std::sqrt(psquare+std::pow(kPiMass,2.));
}

float eta(float pz, float psquared){
    float pnorm=std::sqrt(psquared);
    return 0.5*std::log((pnorm+pz)/(pnorm-pz));
}

std::tuple<double,double,double,double> fitGaus(TH1D* histo) {
    // use gaus to fit
    TF1 *f2 = new TF1("f2", "gaus");
    Int_t b_max = histo->GetMaximumBin();
    Double_t x_max = histo->GetBinCenter(b_max);
    Double_t y_max = histo->GetBinContent(b_max);
    f2->SetParameter(0,y_max);
    f2->SetParameter(1,histo->GetMean());
    f2->SetParameter(2,histo->GetStdDev());
    histo->Fit(f2);
    std::cout<<"Fitted"<<std::endl;
    histo->Write();

    //hK0->Write();
    // use mean-2sigma as a new fitting range
    double mean= f2->GetParameter(1);
    double sigma = f2->GetParameter(2);
    // second fit
    histo->Fit(f2, "R", "", mean-1.5*sigma, mean+1.5*sigma);
    std::cout<<"Fitted 2"<<std::endl;
    histo->Write();
    // save mean, sigma, and corresponding errors for the plot
    mean= f2->GetParameter(1);
    sigma = f2->GetParameter(2);
    double mean_error = f2->GetParError(1);
    double sigma_error = f2->GetParError(1);
    return std::make_tuple(mean, sigma, mean_error, sigma_error);
}

void calculateShift(TString file="data/SnapshotAP.root", TString tree="NewVariables"){
    ROOT::EnableImplicitMT(10);
    ROOT::RDataFrame df_data(tree, file);
    ROOT::RDataFrame df_MC("NewVariables","data/SnapshotMC_symmetricMassCut.root");

    auto df=df_data.Filter("alpha>-0.05 && alpha<0.05").Define("pT1", pT, {"fPxPos","fPyPos"}).Define("pT2", pT, {"fPxNeg","fPyNeg"})//.Filter("pT2 > 1 && pT2 < 1.1")
                .Define("pTavg", "(pT1+pT2)/2")
                .Define("p1square", psquare, {"fPxPos","fPyPos","fPzPos"}).Define("p2square", psquare, {"fPxNeg","fPyNeg","fPzNeg"})
                .Define("eta1",eta,{"fPzPos","p1square"}).Define("eta2",eta,{"fPzNeg","p2square"}).Filter("abs(eta1) < 0.3").Filter("abs(eta2) < 0.3")
                .Define("E1", energy, {"p1square"}).Define("E2", energy, {"p2square"}).Define("p1p2", p1p2, {"fPxPos","fPyPos","fPzPos","fPxNeg","fPyNeg","fPzNeg"})
                .Define("shiftplus", deltapTone, {"Minv_K0", "pT1", "E1", "E2", "p1square", "p1p2"})
                .Define("shiftminus", deltapTone, {"Minv_K0", "pT2", "E2", "E1", "p2square", "p1p2"})
                .Define("shift", deltapT, {"Minv_K0", "pT1", "pT2", "E1", "E2", "p1square", "p2square", "p1p2"})
                .Define("shiftplusPDG", deltapTonePDG, {"Minv_K0", "pT1", "E1", "E2", "p1square", "p1p2"})
                .Define("shiftminusPDG", deltapTonePDG, {"Minv_K0", "pT2", "E2", "E1", "p2square", "p1p2"})
                .Define("shiftplusLHCb", LHCb, {"Minv_K0", "E1", "E2"})
                .Define("shiftminusLHCb", LHCb, {"Minv_K0", "E2", "E1"});

    // define columns for MC data
    auto df_MC_complete = df_MC.Filter("alpha>-0.05 && alpha<0.05").Define("pT1", pT, {"fPxPos","fPyPos"}).Define("pT2", pT, {"fPxNeg","fPyNeg"})//.Filter("pT2 > 1 && pT2 < 1.1")
                .Define("pT1gen", pT, {"fPxPosMC","fPyPosMC"}).Define("pT2gen", pT, {"fPxNegMC","fPyNegMC"})
                .Define("pTavg", "(pT1+pT2)/2").Define("pTavgGen", "(pT1gen+pT2gen)/2")
                .Define("difference1", "pT1-pT1gen").Define("difference2", "pT2-pT2gen")
                .Define("difference", "pTavg-pTavgGen")
                .Define("p1square", psquare, {"fPxPos","fPyPos","fPzPos"}).Define("p2square", psquare, {"fPxNeg","fPyNeg","fPzNeg"})
                .Define("eta1",eta,{"fPzPos","p1square"}).Define("eta2",eta,{"fPzNeg","p2square"}).Filter("abs(eta1) < 0.3").Filter("abs(eta2) < 0.3")
                .Define("E1", energy, {"p1square"}).Define("E2", energy, {"p2square"}).Define("p1p2", p1p2, {"fPxPos","fPyPos","fPzPos","fPxNeg","fPyNeg","fPzNeg"})
                .Define("shiftplus", deltapTone, {"Minv_K0", "pT1", "E1", "E2", "p1square", "p1p2"})
                .Define("shiftminus", deltapTone, {"Minv_K0", "pT2", "E2", "E1", "p2square", "p1p2"})
                .Define("shift", deltapT, {"Minv_K0", "pT1", "pT2", "E1", "E2", "p1square", "p2square", "p1p2"})
                .Define("p1squaregen", psquare, {"fPxPosMC","fPyPosMC","fPzPosMC"}).Define("p2squaregen", psquare, {"fPxNegMC","fPyNegMC","fPzNegMC"})
                .Define("E1gen", energy, {"p1squaregen"}).Define("E2gen", energy, {"p2squaregen"}).Define("p1p2gen", p1p2, {"fPxPosMC","fPyPosMC","fPzPosMC","fPxNegMC","fPyNegMC","fPzNegMC"})
                .Define("shiftplusPDG", deltapTonePDG, {"Minv_K0", "pT1gen", "E1gen", "E2gen", "p1squaregen", "p1p2gen"})
                .Define("shiftminusPDG", deltapTonePDG, {"Minv_K0", "pT2gen", "E2gen", "E1gen", "p2squaregen", "p1p2gen"})
                .Define("shiftplusLHCb", LHCb, {"Minv_K0", "E1", "E2"})
                .Define("shiftminusLHCb", LHCb, {"Minv_K0", "E2", "E1"})
                .Define("shiftplusLHCbgen", LHCb, {"Minv_K0", "E1gen", "E2gen"})
                .Define("shiftminusLHCbgen", LHCb, {"Minv_K0", "E2gen", "E1gen"});




    std::unique_ptr<TFile> myFile( TFile::Open("pT_shift_new.root", "RECREATE") );
    int n=25;

    // data

    TH2D *h = new TH2D("shift_2D_pos","Shift for positive particles;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_plus = df.Histo2D(*h,"pT1","shiftplus");

    histo_data_plus->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    histo_data_plus->GetYaxis()->SetTitle("#delta p_{T}^{+} [GeV]");

    histo_data_plus->Write();

    TH2D *h1 = new TH2D("shift_2D_neg","Shift for negative particles;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_minus = df.Histo2D(*h1,"pT2","shiftminus");

    histo_data_minus->GetXaxis()->SetTitle("p_{T}^{-} [GeV]");
    histo_data_minus->GetYaxis()->SetTitle("#delta p_{T}^{-} [GeV]");

    histo_data_minus->Write();

    // negative and positive together

    TH2D *h0 = new TH2D("shift_2D","Shift for negative particles;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data = df.Histo2D(*h0,"pTavg","shift");

    histo_data->GetXaxis()->SetTitle("#bar{p_{T}^{#pm}} [GeV]");
    histo_data->GetYaxis()->SetTitle("#delta p_{T}^{#pm} [GeV]");

    histo_data->Write();

    // data PDG Mass

    TH2D *h20 = new TH2D("shift_2D_pos_PDGMass","Shift for positive particles using PDG Mass;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_plus_PDG = df.Histo2D(*h20,"pT1","shiftplusPDG");

    histo_data_plus_PDG->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    histo_data_plus_PDG->GetYaxis()->SetTitle("#delta p_{T}^{+} [GeV]");

    histo_data_plus_PDG->Write();

    TH2D *h21 = new TH2D("shift_2D_neg_PDGMass","Shift for negative particles using PDG Mass;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_minus_PDG = df.Histo2D(*h21,"pT2","shiftminusPDG");

    histo_data_minus_PDG->GetXaxis()->SetTitle("p_{T}^{-} [GeV]");
    histo_data_minus_PDG->GetYaxis()->SetTitle("#delta p_{T}^{-} [GeV]");

    histo_data_minus_PDG->Write();

    // data LHCb formula

    TH2D *h22 = new TH2D("shift_2D_pos_LHCb","Shift for positive particles using PDG Mass;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_plus_LHCb = df.Histo2D(*h22,"pT1","shiftplusLHCb");

    histo_data_plus_LHCb->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    histo_data_plus_LHCb->GetYaxis()->SetTitle("#delta p_{T}^{+} [GeV]");

    histo_data_plus_LHCb->Write();

    TH2D *h23 = new TH2D("shift_2D_neg_LHCb","Shift for negative particles using PDG Mass;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_minus_LHCb = df.Histo2D(*h23,"pT2","shiftminusLHCb");

    histo_data_minus_LHCb->GetXaxis()->SetTitle("p_{T}^{-} [GeV]");
    histo_data_minus_LHCb->GetYaxis()->SetTitle("#delta p_{T}^{-} [GeV]");

    histo_data_minus_LHCb->Write();

    // momenta
    TH2D *h24 = new TH2D("PTs","PTs",400,0.,2.5,400,0.,2.5);
    
    auto histo_momenta = df.Histo2D(*h24,"pT1","pT2");

    histo_momenta->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    histo_momenta->GetYaxis()->SetTitle("p_{T}^{-} [GeV]");

    histo_momenta->Write();

    TH2D *h25 = new TH2D("PTs_MC_avg_pos","PTs",400,0.,2.5,400,0.,2.5);
    
    auto histo_momenta_data_avg_pos = df.Histo2D(*h25,"pTavg","pT1");

    histo_momenta_data_avg_pos->GetXaxis()->SetTitle("p_{T}^{avg} [GeV]");
    histo_momenta_data_avg_pos->GetYaxis()->SetTitle("p_{T}^{+} [GeV]");

    histo_momenta_data_avg_pos->Write();

    TH2D *h26 = new TH2D("PTs_MC_avg_neg","PTs",400,0.,2.5,400,0.,2.5);
    
    auto histo_momenta_data_avg_neg = df.Histo2D(*h26,"pTavg","pT2");

    histo_momenta_data_avg_neg->GetXaxis()->SetTitle("p_{T}^{avg} [GeV]");
    histo_momenta_data_avg_neg->GetYaxis()->SetTitle("p_{T}^{-} [GeV]");

    histo_momenta_data_avg_neg->Write();

    // MC

    TH2D *h7 = new TH2D("shift_2D_pos_MC","Shift for positive particles MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_MC_plus = df_MC_complete.Histo2D(*h7,"pT1","shiftplus");

    histo_MC_plus->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    histo_MC_plus->GetYaxis()->SetTitle("#delta p_{T}^{+} [GeV]");

    histo_MC_plus->Write();

    TH2D *h8 = new TH2D("shift_2D_neg_MC","Shift for negative particles MC;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_MC_minus = df_MC_complete.Histo2D(*h8,"pT2","shiftminus");

    histo_MC_minus->GetXaxis()->SetTitle("p_{T}^{-} [GeV]");
    histo_MC_minus->GetYaxis()->SetTitle("#delta p_{T}^{-} [GeV]");

    histo_MC_minus->Write();

    // negative and positive together

    TH2D *h10 = new TH2D("shift_2D_MC","Shift;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_MC = df_MC_complete.Histo2D(*h10,"pTavg","shift");

    histo_MC->GetXaxis()->SetTitle("#bar{p_{T}^{#pm}} [GeV]");
    histo_MC->GetYaxis()->SetTitle("#delta p_{T}^{#pm} [GeV]");

    histo_MC->Write();

    // difference of pTgen pTreco
    TH2D *h50 = new TH2D("difference_MC","Shift MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_MC_check = df_MC_complete.Histo2D(*h50,"pTavg","difference");

    histo_MC_check->GetXaxis()->SetTitle("p_{T}^{avg}^{reco} [GeV]");
    histo_MC_check->GetYaxis()->SetTitle("p_{T}^{avg}^{reco}-p_{T}^{avg}^{gen} [GeV]");

    histo_MC_check->Write();

    TH2D *h5 = new TH2D("difference_MC_pos","Shift for positive particles MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_MC_plus_check = df_MC_complete.Histo2D(*h5,"pT1","difference1");

    histo_MC_plus_check->GetXaxis()->SetTitle("p_{T}^{+}^{reco} [GeV]");
    histo_MC_plus_check->GetYaxis()->SetTitle("p_{T}^{+}^{reco}-p_{T}^{+}^{gen} [GeV]");

    histo_MC_plus_check->Write();

    TH2D *h6 = new TH2D("difference_MC_neg","Shift for positive particles MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_MC_minus_check = df_MC_complete.Histo2D(*h6,"pT2","difference2");

    histo_MC_minus_check->GetXaxis()->SetTitle("p_{T}^{-}^{reco} [GeV]");
    histo_MC_minus_check->GetYaxis()->SetTitle("p_{T}^{-}^{reco}-p_{T}^{-}^{gen} [GeV]");

    histo_MC_minus_check->Write();

    // PDG Mass and generated

    TH2D *h40 = new TH2D("shift_2D_pos_PDG_MC","Shift for positive particles using PDG Mass MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_plus_PDG_MC = df_MC_complete.Histo2D(*h40,"pT1","shiftplusPDG");

    histo_data_plus_PDG_MC->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    histo_data_plus_PDG_MC->GetYaxis()->SetTitle("#delta p_{T}^{+} [GeV]");

    histo_data_plus_PDG_MC->Write();

    TH2D *h41 = new TH2D("shift_2D_neg_PDG_MC","Shift for negative particles using PDG Mass;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_minus_PDG_MC = df_MC_complete.Histo2D(*h41,"pT2","shiftminusPDG");

    histo_data_minus_PDG_MC->GetXaxis()->SetTitle("p_{T}^{-} [GeV]");
    histo_data_minus_PDG_MC->GetYaxis()->SetTitle("#delta p_{T}^{-} [GeV]");

    histo_data_minus_PDG_MC->Write();

    // LHCb formula

    TH2D *h30 = new TH2D("shift_2D_pos_LHCb_MC","Shift for positive particles using PDG Mass MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_plus_LHCb_MC = df_MC_complete.Histo2D(*h30,"pT1","shiftplusLHCb");

    histo_data_plus_LHCb_MC->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    histo_data_plus_LHCb_MC->GetYaxis()->SetTitle("#delta p_{T}^{+} [GeV]");

    histo_data_plus_LHCb_MC->Write();

    TH2D *h31 = new TH2D("shift_2D_neg_LHCb_MC","Shift for negative particles using PDG Mass;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_minus_LHCb_MC = df_MC_complete.Histo2D(*h31,"pT2","shiftminusLHCb");

    histo_data_minus_LHCb_MC->GetXaxis()->SetTitle("p_{T}^{-} [GeV]");
    histo_data_minus_LHCb_MC->GetYaxis()->SetTitle("#delta p_{T}^{-} [GeV]");

    histo_data_minus_LHCb_MC->Write();

    // LHCb formula generated energy

    TH2D *h32 = new TH2D("shift_2D_pos_LHCb_MC_gen","Shift for positive particles using PDG Mass MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_plus_LHCb_MC_gen = df_MC_complete.Histo2D(*h32,"pT1","shiftplusLHCbgen");

    histo_data_plus_LHCb_MC_gen->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    histo_data_plus_LHCb_MC_gen->GetYaxis()->SetTitle("#delta p_{T}^{+} [GeV]");

    histo_data_plus_LHCb_MC_gen->Write();

    TH2D *h33 = new TH2D("shift_2D_neg_LHCb_MC_gen","Shift for negative particles using PDG Mass;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_minus_LHCb_MC_gen = df_MC_complete.Histo2D(*h33,"pT2","shiftminusLHCbgen");

    histo_data_minus_LHCb_MC_gen->GetXaxis()->SetTitle("p_{T}^{-} [GeV]");
    histo_data_minus_LHCb_MC_gen->GetYaxis()->SetTitle("#delta p_{T}^{-} [GeV]");

    histo_data_minus_LHCb_MC_gen->Write();

    // momenta
    TH2D *h34 = new TH2D("PTs_MC","PTs",400,0.,2.5,400,0.,2.5);
    
    auto histo_momenta_MC = df_MC_complete.Histo2D(*h34,"pT1","pT2");

    histo_momenta_MC->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    histo_momenta_MC->GetYaxis()->SetTitle("p_{T}^{-} [GeV]");

    histo_momenta_MC->Write();

    TH2D *h35 = new TH2D("PTs_MC_avg_pos_MC","PTs",400,0.,2.5,400,0.,2.5);
    
    auto histo_momenta_MC_avg_pos = df_MC_complete.Histo2D(*h35,"pTavg","pT1");

    histo_momenta_MC_avg_pos->GetXaxis()->SetTitle("p_{T}^{avg} [GeV]");
    histo_momenta_MC_avg_pos->GetYaxis()->SetTitle("p_{T}^{+} [GeV]");

    histo_momenta_MC_avg_pos->Write();

    TH2D *h36 = new TH2D("PTs_MC_avg_neg_MC","PTs",400,0.,2.5,400,0.,2.5);
    
    auto histo_momenta_MC_avg_neg = df_MC_complete.Histo2D(*h36,"pTavg","pT2");

    histo_momenta_MC_avg_neg->GetXaxis()->SetTitle("p_{T}^{avg} [GeV]");
    histo_momenta_MC_avg_neg->GetYaxis()->SetTitle("p_{T}^{-} [GeV]");

    histo_momenta_MC_avg_neg->Write();


    TH1D *h11 = new TH1D("Minv_MC","Minv_MC",500,0.48,0.51);
    
    auto histo11 = df_MC_complete.Histo1D(*h11,"Minv_K0");

    histo11->Write();

    TH1D *h12 = new TH1D("etaPos","eta",500,-1.,1.);
    
    auto histo12 = df.Histo1D(*h12,"eta1");

    histo12->Write();

    TH1D *h13 = new TH1D("etaNeg","eta",500,-1.,1.);
    
    auto histo13 = df.Histo1D(*h13,"eta2");

    histo13->Write();


    TH1D *h14 = new TH1D("etaPos_MC","eta",500,-1.,1.);
    
    auto histo14 = df_MC_complete.Histo1D(*h14,"eta1");

    histo14->Write();

    TH1D *h15 = new TH1D("etaNeg_MC","eta",500,-1.,1.);
    
    auto histo15 = df_MC_complete.Histo1D(*h15,"eta2");

    histo15->Write();

    // calculate projections 
    // project on xaxis
    // calculate values to plot for alpha 
    int n_iterate=15;
    double x[n_iterate];

    double mu_data[n_iterate];
    double sigma_data[n_iterate];
    double mu_error_data[n_iterate];
    double sigma_error_data[n_iterate];
    double mu_MC[n_iterate];
    double sigma_MC[n_iterate];
    double mu_error_MC[n_iterate];
    double sigma_error_MC[n_iterate];
    double mu_MC_diff[n_iterate];
    double sigma_MC_diff[n_iterate];
    double mu_error_MC_diff[n_iterate];
    double sigma_error_MC_diff[n_iterate];

    double mu_data_pos[n_iterate];
    double sigma_data_pos[n_iterate];
    double mu_error_data_pos[n_iterate];
    double sigma_error_data_pos[n_iterate];
    double mu_MC_pos[n_iterate];
    double sigma_MC_pos[n_iterate];
    double mu_error_MC_pos[n_iterate];
    double sigma_error_MC_pos[n_iterate];
    double mu_MC_diff_pos[n_iterate];
    double sigma_MC_diff_pos[n_iterate];
    double mu_error_MC_diff_pos[n_iterate];
    double sigma_error_MC_diff_pos[n_iterate];

    double mu_data_neg[n_iterate];
    double sigma_data_neg[n_iterate];
    double mu_error_data_neg[n_iterate];
    double sigma_error_data_neg[n_iterate];
    double mu_MC_neg[n_iterate];
    double sigma_MC_neg[n_iterate];
    double mu_error_MC_neg[n_iterate];
    double sigma_error_MC_neg[n_iterate];
    double mu_MC_diff_neg[n_iterate];
    double sigma_MC_diff_neg[n_iterate];
    double mu_error_MC_diff_neg[n_iterate];
    double sigma_error_MC_diff_neg[n_iterate];
    for (int i = 1; i <=15; i+=1) {
        // calculate x
        // use part of the histogram which corresponds to one bin
        std::cout<<"Iteration "<<i<<std::endl;
        int lower_bin;
        int upper_bin;
        if (i<=8){
            lower_bin=i;
            upper_bin=i;
        }else{
            lower_bin=9+2*(i-9);
            upper_bin=lower_bin+1;
        }
        std::cout<<i<<" "<<lower_bin<<" "<<upper_bin<<std::endl;
        double bin_low=histo_MC_minus_check->GetXaxis()->GetBinLowEdge(lower_bin);
        double bin_width=histo_MC_minus_check->GetXaxis()->GetBinWidth(lower_bin);
        if (i>8){
            bin_width=2*bin_width;
        }
        double bin_up=bin_low+bin_width;

        std::cout<<"Lower edge: "<<bin_low<<"; Upper edge: "<<bin_low<<"; Width: "<<bin_width<<std::endl;
        // middle of the bin
        x[i-1]=bin_low+bin_width/2;
        std::cout<<x[i-1]<<std::endl;
        // error is the distance to next point 
        //xerr[i-1]=bin_width/2;
        // project to y axis
        // Difference MC 
        
        TH1D *histo_bin_MC_minus_check = histo_MC_minus_check->ProjectionY("MC_neg_difference",lower_bin,upper_bin);
        TH1D *histo_bin_MC_plus_check = histo_MC_plus_check->ProjectionY("MC_pos_difference",lower_bin,upper_bin);
        TH1D *histo_bin_MC_check = histo_MC_check->ProjectionY("MC_difference",lower_bin,upper_bin);
        // calculated MC 
        TH1D *histo_bin_MC_minus = histo_MC_minus->ProjectionY("MC_neg",lower_bin,upper_bin);
        TH1D *histo_bin_MC_plus = histo_MC_plus->ProjectionY("MC_pos",lower_bin,upper_bin);
        TH1D *histo_bin_MC = histo_MC->ProjectionY("MC",lower_bin,upper_bin);
        // calculated data 
        TH1D *histo_bin_minus = histo_data_minus->ProjectionY("data_neg",lower_bin,upper_bin);
        TH1D *histo_bin_plus = histo_data_plus->ProjectionY("data_pos",lower_bin,upper_bin);
        TH1D *histo_bin = histo_data->ProjectionY("data",lower_bin,upper_bin);
        std::cout<<"Created projections"<<std::endl;
        // Use gaussian fit to calculate mu and sigma
        // use function
        myFile->cd();
        gDirectory->mkdir(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));
        // data 
        //both
        myFile->cd();
        myFile->cd(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));
        auto values_data=fitGaus(histo_bin);
        mu_data[i-1] = std::get<0>(values_data)/x[i-1];
        sigma_data[i-1] = std::get<1>(values_data)/x[i-1];
        mu_error_data[i-1] = std::get<2>(values_data)/x[i-1];
        sigma_error_data[i-1] = std::get<3>(values_data)/x[i-1];
        // plus
        myFile->cd();
        myFile->cd(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));
        auto values_data_pos=fitGaus(histo_bin_plus);
        mu_data_pos[i-1] = std::get<0>(values_data_pos)/x[i-1];
        sigma_data_pos[i-1] = std::get<1>(values_data_pos)/x[i-1];
        mu_error_data_pos[i-1] = std::get<2>(values_data_pos)/x[i-1];
        sigma_error_data_pos[i-1] = std::get<3>(values_data_pos)/x[i-1];
        std::cout<<mu_data_pos[i-1]<<" "<<sigma_data_pos[i-1]<<" "<<std::endl;
        // minus
        myFile->cd();
        myFile->cd(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));
        auto values_data_neg=fitGaus(histo_bin_minus);
        mu_data_neg[i-1] = std::get<0>(values_data_neg)/x[i-1];
        sigma_data_neg[i-1] = std::get<1>(values_data_neg)/x[i-1];
        mu_error_data_neg[i-1] = std::get<2>(values_data_neg)/x[i-1];
        sigma_error_data_neg[i-1] = std::get<3>(values_data_neg)/x[i-1];
        std::cout<<bin_low<<" "<<bin_up<<" "<<mu_data_neg[i-1]<<" "<<sigma_data_neg[i-1]<<" "<<std::endl;
        // MC 
        // both
        myFile->cd();
        myFile->cd(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));
        auto values_MC=fitGaus(histo_bin_MC);
        mu_MC[i-1] = std::get<0>(values_MC)/x[i-1];
        sigma_MC[i-1] = std::get<1>(values_MC)/x[i-1];
        mu_error_MC[i-1] = std::get<2>(values_MC)/x[i-1];
        sigma_error_MC[i-1] = std::get<3>(values_MC)/x[i-1];
        //plus
        myFile->cd();
        myFile->cd(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));
        auto values_MC_pos=fitGaus(histo_bin_MC_plus);
        mu_MC_pos[i-1] = std::get<0>(values_MC_pos)/x[i-1];
        sigma_MC_pos[i-1] = std::get<1>(values_MC_pos)/x[i-1];
        mu_error_MC_pos[i-1] = std::get<2>(values_MC_pos)/x[i-1];
        sigma_error_MC_pos[i-1] = std::get<3>(values_MC_pos)/x[i-1];
        //minus
        myFile->cd();
        myFile->cd(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));
        auto values_MC_neg=fitGaus(histo_bin_MC_minus);
        mu_MC_neg[i-1] = std::get<0>(values_MC_neg)/x[i-1];
        sigma_MC_neg[i-1] = std::get<1>(values_MC_neg)/x[i-1];
        mu_error_MC_neg[i-1] = std::get<2>(values_MC_neg)/x[i-1];
        sigma_error_MC_neg[i-1] = std::get<3>(values_MC_neg)/x[i-1];
        //difference
        //both
        myFile->cd();
        myFile->cd(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));
        auto values_MC_diff=fitGaus(histo_bin_MC_check);
        mu_MC_diff[i-1] = std::get<0>(values_MC_diff)/x[i-1];
        sigma_MC_diff[i-1] = std::get<1>(values_MC_diff)/x[i-1];
        mu_error_MC_diff[i-1] = std::get<2>(values_MC_diff)/x[i-1];
        sigma_error_MC_diff[i-1] = std::get<3>(values_MC_diff)/x[i-1];
        //plus
        myFile->cd();
        myFile->cd(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));
        auto values_MC_diff_pos=fitGaus(histo_bin_MC_plus_check);
        mu_MC_diff_pos[i-1] = std::get<0>(values_MC_diff_pos)/x[i-1];
        sigma_MC_diff_pos[i-1] = std::get<1>(values_MC_diff_pos)/x[i-1];
        mu_error_MC_diff_pos[i-1] = std::get<2>(values_MC_diff_pos)/x[i-1];
        sigma_error_MC_diff_pos[i-1] = std::get<3>(values_MC_diff_pos)/x[i-1];
        // minus
        myFile->cd();
        myFile->cd(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));
        auto values_MC_diff_neg=fitGaus(histo_bin_MC_minus_check);
        mu_MC_diff_neg[i-1] = std::get<0>(values_MC_diff_neg)/x[i-1];
        sigma_MC_diff_neg[i-1] = std::get<1>(values_MC_diff_neg)/x[i-1];
        mu_error_MC_diff_neg[i-1] = std::get<2>(values_MC_diff_neg)/x[i-1];
        sigma_error_MC_diff_neg[i-1] = std::get<3>(values_MC_diff_neg)/x[i-1];

        myFile->cd();
        myFile->cd(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));
        histo_bin_MC_plus_check->Write();
        histo_bin_MC_plus->Write();
        histo_bin_plus->Write();
        histo_bin_MC_minus_check->Write();
        histo_bin_MC_minus->Write();
        histo_bin_minus->Write();

    }

    auto gr100 = new TGraphErrors(n_iterate,x,sigma_data,nullptr,sigma_error_data);
    gr100->SetTitle("Data all tracks");
    gr100->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr100->GetYaxis()->SetTitle("#frac{#sigma(p_{T})}{p_{T}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr100->Write("sigma_data");

    auto gr = new TGraphErrors(n_iterate,x,sigma_data_pos,nullptr,sigma_error_data_pos);
    gr->SetTitle("Data positive tracks");
    gr->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    gr->GetYaxis()->SetTitle("#frac{#sigma(p_{T}^{+})}{p_{T}^{+}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr->Write("sigma_data_pos");

    auto gr1 = new TGraphErrors(n_iterate,x,sigma_data_neg,nullptr,sigma_error_data_neg);
    gr1->SetTitle("Data negative tracks");
    gr1->GetXaxis()->SetTitle("p_{T}^{-} [GeV]");
    gr1->GetYaxis()->SetTitle("#frac{#sigma(p_{T}^{-})}p_{T}^{-}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr1->Write("sigma_data_neg");

    auto gr101 = new TGraphErrors(n_iterate,x,mu_data,nullptr,mu_error_data);
    gr101->SetTitle("Data all tracks");
    gr101->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr101->GetYaxis()->SetTitle("#frac{#mu(p_{T})}{p_{T}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr101->Write("mu_data");

    auto gr2 = new TGraphErrors(n_iterate,x,mu_data_pos,nullptr,mu_error_data_pos);
    gr2->SetTitle("Data positive");
    gr2->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    gr2->GetYaxis()->SetTitle("#frac{#mu(p_{T}^{+})}p_{T}^{+}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr2->Write("mu_data_pos");

    auto gr3 = new TGraphErrors(n_iterate,x,mu_data_neg,nullptr,mu_error_data_neg);
    gr3->SetTitle("Data negative tracks");
    gr3->GetXaxis()->SetTitle("p_{T}^{-}[GeV]");
    gr3->GetYaxis()->SetTitle("#frac{#mu(p_{T}^{-})}{p_{T}^{-}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr3->Write("mu_data_neg");

    auto gr102 = new TGraphErrors(n_iterate,x,sigma_MC,nullptr,sigma_error_MC);
    gr102->SetTitle("MC");
    gr102->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr102->GetYaxis()->SetTitle("#frac{#sigma(p_{T})}{p_{T}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr102->Write("sigma_MC");

    auto gr4 = new TGraphErrors(n_iterate,x,sigma_MC_pos,nullptr,sigma_error_MC_pos);
    gr4->SetTitle("MC positive tracks");
    gr4->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    gr4->GetYaxis()->SetTitle("#frac{#sigma(p_{T}^{+})}{p_{T}^{+}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr4->Write("sigma_MC_pos");

    auto gr5 = new TGraphErrors(n_iterate,x,sigma_MC_neg,nullptr,sigma_error_MC_neg);
    gr5->SetTitle("MC negative tracks");
    gr5->GetXaxis()->SetTitle("p_{T}^{-} [GeV]");
    gr5->GetYaxis()->SetTitle("#frac{#sigma(p_{T}^{-})}{p_{T}^{-}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr5->Write("sigma_MC_neg");

    auto gr103 = new TGraphErrors(n_iterate,x,mu_MC,nullptr,mu_error_MC);
    gr103->SetTitle("MC");
    gr103->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr103->GetYaxis()->SetTitle("#frac{#mu(p_{T})}{p_T}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr103->Write("mu_MC");

    auto gr6 = new TGraphErrors(n_iterate,x,mu_MC_pos,nullptr,mu_error_MC_pos);
    gr6->SetTitle("MC positive");
    gr6->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    gr6->GetYaxis()->SetTitle("#frac{#mu(p_{T}^{+})}{p_{T}^{+}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr6->Write("mu_MC_pos");

    auto gr7 = new TGraphErrors(n_iterate,x,mu_MC_neg,nullptr,mu_error_MC_neg);
    gr7->SetTitle("MC negative tracks");
    gr7->GetXaxis()->SetTitle("p_{T}^{-} [GeV]");
    gr7->GetYaxis()->SetTitle("#frac{#mu(p_{T}^{-})}{p_{T}^{-}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr7->Write("mu_MC_neg");

    auto gr104 = new TGraphErrors(n_iterate,x,sigma_MC_diff,nullptr,sigma_error_MC_diff);
    gr104->SetTitle("MC difference");
    gr104->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr104->GetYaxis()->SetTitle("#frac{#sigma(p_{T})}{p_{T}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr104->Write("sigma_MC_diff");

    auto gr8 = new TGraphErrors(n_iterate,x,sigma_MC_diff_pos,nullptr,sigma_error_MC_diff_pos);
    gr8->SetTitle("MC difference positive tracks");
    gr8->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    gr8->GetYaxis()->SetTitle("#frac{#sigma(p_{T}^{+})}{p_{T}^{+}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr8->Write("sigma_MC_diff_pos");

    auto gr9 = new TGraphErrors(n_iterate,x,sigma_MC_diff_neg,nullptr,sigma_error_MC_diff_neg);
    gr9->SetTitle("MC difference negative tracks");
    gr9->GetXaxis()->SetTitle("p_{T}^{-} [GeV]");
    gr9->GetYaxis()->SetTitle("#frac{#sigma(p_{T}^{-})}{p_{T}^{-}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr9->Write("sigma_MC_diff_neg");

    auto gr105 = new TGraphErrors(n_iterate,x,mu_MC_diff,nullptr,mu_error_MC_diff);
    gr105->SetTitle("MC difference all tracks");
    gr105->GetXaxis()->SetTitle("p_{T} [GeV]");
    gr105->GetYaxis()->SetTitle("#frac{#mu(p_{T})}{p_T}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr105->Write("mu_MC_diff");

    auto gr10 = new TGraphErrors(n_iterate,x,mu_MC_diff_pos,nullptr,mu_error_MC_diff_pos);
    gr10->SetTitle("MC differencce positive tracks");
    gr10->GetXaxis()->SetTitle("p_{T}^{+} [GeV]");
    gr10->GetYaxis()->SetTitle("#frac{#mu(p_{T}^{+})}{p_{T}^{+}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr10->Write("mu_MC_diff_pos");

    auto gr11 = new TGraphErrors(n_iterate,x,mu_MC_diff_neg,nullptr,mu_error_MC_diff_neg);
    gr11->SetTitle("MC difference negative tracks");
    gr11->GetXaxis()->SetTitle("p_{T}^{-} [GeV]");
    gr11->GetYaxis()->SetTitle("#frac{#mu(p_{T}^{-})}{p_{T}^{-}}");
    //gr->SetLineWidth(0);
    //gr->SetMarkerSize(1.5);
    //gr->SetMarkerStyle(70);
    myFile->cd();
    gr11->Write("mu_MC_diff_neg");

    myFile->Close();

}