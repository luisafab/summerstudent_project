#include <ROOT/RDataFrame.hxx>
#include "armenterosPlot.c"

float deltapT(float M, float pT1, float pT2, float E1, float E2, float p1squared, float p2squared, float p1p2){
    float deltaM= M-kK0sMass;
    return deltaM*M/((2*p1squared*E2/(pT1*E1))+(2*p2squared*E1/(pT2*E2))-(p1p2)*(1/pT1+1/pT2));
}

float deltapTone(float M, float pT1, float E1, float E2, float p1squared, float p1p2){
    float deltaM= M-kK0sMass;
    return deltaM*M/((2*p1squared*E2/(pT1*E1))-(p1p2)/pT1);
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


void calculateShift(TString file="data/SnapshotAP.root", TString tree="NewVariables"){
    ROOT::EnableImplicitMT(10);
    ROOT::RDataFrame df_data(tree, file);
    ROOT::RDataFrame df_MC("NewVariables","SnapshotMC.root");

    auto df=df_data.Filter("alpha>-0.05 && alpha<0.05").Define("pT1", pT, {"fPxPos","fPyPos"}).Define("pT2", pT, {"fPxNeg","fPyNeg"})//.Filter("pT2 > 1 && pT2 < 1.1")
                .Define("p1square", psquare, {"fPxPos","fPyPos","fPzPos"}).Define("p2square", psquare, {"fPxNeg","fPyNeg","fPzNeg"})
                .Define("eta1",eta,{"fPzPos","p1square"}).Define("eta2",eta,{"fPzNeg","p2square"}).Filter("abs(eta1) < 0.3").Filter("abs(eta2) < 0.3")
                .Define("E1", energy, {"p1square"}).Define("E2", energy, {"p2square"}).Define("p1p2", p1p2, {"fPxPos","fPyPos","fPzPos","fPxNeg","fPyNeg","fPzNeg"})
                .Define("shiftplus", deltapTone, {"Minv_K0", "pT1", "E1", "E2", "p1square", "p1p2"})
                .Define("shiftminus", deltapTone, {"Minv_K0", "pT2", "E2", "E1", "p2square", "p1p2"})
                .Define("shiftplusPDG", deltapTonePDG, {"Minv_K0", "pT1", "E1", "E2", "p1square", "p1p2"})
                .Define("shiftminusPDG", deltapTonePDG, {"Minv_K0", "pT2", "E2", "E1", "p2square", "p1p2"})
                .Define("shiftplusLHCb", LHCb, {"Minv_K0", "E1", "E2"})
                .Define("shiftminusLHCb", LHCb, {"Minv_K0", "E2", "E1"});

    // define columns for MC data
    auto df_MC_complete = df_MC.Filter("alpha>-0.05 && alpha<0.05").Define("pT1", pT, {"fPxPos","fPyPos"}).Define("pT2", pT, {"fPxNeg","fPyNeg"})//.Filter("pT2 > 1 && pT2 < 1.1")
                .Define("pT1gen", pT, {"fPxPosMC","fPyPosMC"}).Define("pT2gen", pT, {"fPxNegMC","fPyNegMC"})
                .Define("difference1", "pT1-pT1gen").Define("difference2", "pT2-pT2gen")
                .Define("p1square", psquare, {"fPxPos","fPyPos","fPzPos"}).Define("p2square", psquare, {"fPxNeg","fPyNeg","fPzNeg"})
                .Define("eta1",eta,{"fPzPos","p1square"}).Define("eta2",eta,{"fPzNeg","p2square"}).Filter("abs(eta1) < 0.3").Filter("abs(eta2) < 0.3")
                .Define("E1", energy, {"p1square"}).Define("E2", energy, {"p2square"}).Define("p1p2", p1p2, {"fPxPos","fPyPos","fPzPos","fPxNeg","fPyNeg","fPzNeg"})
                .Define("shiftplus", deltapTone, {"Minv_K0", "pT1", "E1", "E2", "p1square", "p1p2"})
                .Define("shiftminus", deltapTone, {"Minv_K0", "pT2", "E2", "E1", "p2square", "p1p2"})
                .Define("p1squaregen", psquare, {"fPxPosMC","fPyPosMC","fPzPosMC"}).Define("p2squaregen", psquare, {"fPxNegMC","fPyNegMC","fPzNegMC"})
                .Define("E1gen", energy, {"p1squaregen"}).Define("E2gen", energy, {"p2squaregen"}).Define("p1p2gen", p1p2, {"fPxPosMC","fPyPosMC","fPzPosMC","fPxNegMC","fPyNegMC","fPzNegMC"})
                .Define("shiftplusPDG", deltapTonePDG, {"Minv_K0", "pT1gen", "E1gen", "E2gen", "p1squaregen", "p1p2gen"})
                .Define("shiftminusPDG", deltapTonePDG, {"Minv_K0", "pT2gen", "E2gen", "E1gen", "p2squaregen", "p1p2gen"})
                .Define("shiftplusLHCb", LHCb, {"Minv_K0", "E1", "E2"})
                .Define("shiftminusLHCb", LHCb, {"Minv_K0", "E2", "E1"})
                .Define("shiftplusLHCbgen", LHCb, {"Minv_K0", "E1gen", "E2gen"})
                .Define("shiftminusLHCbgen", LHCb, {"Minv_K0", "E2gen", "E1gen"});




    std::unique_ptr<TFile> myFile( TFile::Open("pT_shift.root", "RECREATE") );
    int n=50;

    // data

    TH2D *h = new TH2D("shift_2D_pos","Shift for positive particles;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.01,0.01);
    
    auto histo_data_plus = df.Histo2D(*h,"pT1","shiftplus");

    histo_data_plus->GetXaxis()->SetTitle("p_{T}^{+}");
    histo_data_plus->GetYaxis()->SetTitle("#delta p_{T}^{+}");

    histo_data_plus->Write();

    TH2D *h1 = new TH2D("shift_2D_neg","Shift for negative particles;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.01,0.01);
    
    auto histo_data_minus = df.Histo2D(*h1,"pT2","shiftminus");

    histo_data_minus->GetXaxis()->SetTitle("p_{T}^{-}");
    histo_data_minus->GetYaxis()->SetTitle("#delta p_{T}^{-}");

    histo_data_minus->Write();

    // data PDG Mass

    TH2D *h20 = new TH2D("shift_2D_pos_PDGMass","Shift for positive particles using PDG Mass;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.01,0.01);
    
    auto histo_data_plus_PDG = df.Histo2D(*h20,"pT1","shiftplusPDG");

    histo_data_plus_PDG->GetXaxis()->SetTitle("p_{T}^{+}");
    histo_data_plus_PDG->GetYaxis()->SetTitle("#delta p_{T}^{+}");

    histo_data_plus_PDG->Write();

    TH2D *h21 = new TH2D("shift_2D_neg_PDGMass","Shift for negative particles using PDG Mass;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.01,0.01);
    
    auto histo_data_minus_PDG = df.Histo2D(*h21,"pT2","shiftminusPDG");

    histo_data_minus_PDG->GetXaxis()->SetTitle("p_{T}^{-}");
    histo_data_minus_PDG->GetYaxis()->SetTitle("#delta p_{T}^{-}");

    histo_data_minus_PDG->Write();

    // data LHCb formula

    TH2D *h22 = new TH2D("shift_2D_pos_LHCb","Shift for positive particles using PDG Mass;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_plus_LHCb = df.Histo2D(*h22,"pT1","shiftplusLHCb");

    histo_data_plus_LHCb->GetXaxis()->SetTitle("p_{T}^{+}");
    histo_data_plus_LHCb->GetYaxis()->SetTitle("#delta p_{T}^{+}");

    histo_data_plus_LHCb->Write();

    TH2D *h23 = new TH2D("shift_2D_neg_LHCb","Shift for negative particles using PDG Mass;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_minus_LHCb = df.Histo2D(*h23,"pT2","shiftminusLHCb");

    histo_data_minus_LHCb->GetXaxis()->SetTitle("p_{T}^{-}");
    histo_data_minus_LHCb->GetYaxis()->SetTitle("#delta p_{T}^{-}");

    histo_data_minus_LHCb->Write();

    // momenta
    TH2D *h24 = new TH2D("PTs","PTs",n,0.,4.,200,0.,4.);
    
    auto histo_momenta = df.Histo2D(*h24,"pT1","pT2");

    histo_momenta->GetXaxis()->SetTitle("p_{T}^{+}");
    histo_momenta->GetYaxis()->SetTitle("p_{T}^{-}");

    histo_momenta->Write();

    // MC

    TH2D *h7 = new TH2D("shift_2D_pos_MC","Shift for positive particles MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.01,0.01);
    
    auto histo_MC_plus = df_MC_complete.Histo2D(*h7,"pT1","shiftplus");

    histo_MC_plus->GetXaxis()->SetTitle("p_{T}^{+}");
    histo_MC_plus->GetYaxis()->SetTitle("#delta p_{T}^{+}");

    histo_MC_plus->Write();

    TH2D *h8 = new TH2D("shift_2D_neg_MC","Shift for negative particles MC;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.01,0.01);
    
    auto histo_MC_minus = df_MC_complete.Histo2D(*h8,"pT2","shiftminus");

    histo_MC_minus->GetXaxis()->SetTitle("p_{T}^{-}");
    histo_MC_minus->GetYaxis()->SetTitle("#delta p_{T}^{-}");

    histo_MC_minus->Write();

    // difference of pTgen pTreco

    TH2D *h5 = new TH2D("difference_MC_pos","Shift for positive particles MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.05,0.05);
    
    auto histo_MC_plus_check = df_MC_complete.Histo2D(*h5,"pT1","difference1");

    histo_MC_plus_check->GetXaxis()->SetTitle("p_{T}^{+}^{reco}");
    histo_MC_plus_check->GetYaxis()->SetTitle("p_{T}^{+}^{reco}-p_{T}^{+}^{gen}");

    histo_MC_plus_check->Write();

    TH2D *h6 = new TH2D("difference_MC_neg","Shift for positive particles MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.05,0.05);
    
    auto histo_MC_minus_check = df_MC_complete.Histo2D(*h6,"pT2","difference2");

    histo_MC_minus_check->GetXaxis()->SetTitle("p_{T}^{-}^{reco}");
    histo_MC_minus_check->GetYaxis()->SetTitle("p_{T}^{-}^{reco}-p_{T}^{-}^{gen}");

    histo_MC_minus_check->Write();

    // PDG Mass and generated

    TH2D *h40 = new TH2D("shift_2D_pos_PDG_MC","Shift for positive particles using PDG Mass MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_plus_PDG_MC = df_MC_complete.Histo2D(*h40,"pT1","shiftplusPDG");

    histo_data_plus_PDG_MC->GetXaxis()->SetTitle("p_{T}^{+}");
    histo_data_plus_PDG_MC->GetYaxis()->SetTitle("#delta p_{T}^{+}");

    histo_data_plus_PDG_MC->Write();

    TH2D *h41 = new TH2D("shift_2D_neg_PDG_MC","Shift for negative particles using PDG Mass;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_minus_PDG_MC = df_MC_complete.Histo2D(*h41,"pT2","shiftminusPDG");

    histo_data_minus_PDG_MC->GetXaxis()->SetTitle("p_{T}^{-}");
    histo_data_minus_PDG_MC->GetYaxis()->SetTitle("#delta p_{T}^{-}");

    histo_data_minus_PDG_MC->Write();

    // LHCb formula

    TH2D *h30 = new TH2D("shift_2D_pos_LHCb_MC","Shift for positive particles using PDG Mass MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_plus_LHCb_MC = df_MC_complete.Histo2D(*h30,"pT1","shiftplusLHCb");

    histo_data_plus_LHCb_MC->GetXaxis()->SetTitle("p_{T}^{+}");
    histo_data_plus_LHCb_MC->GetYaxis()->SetTitle("#delta p_{T}^{+}");

    histo_data_plus_LHCb_MC->Write();

    TH2D *h31 = new TH2D("shift_2D_neg_LHCb_MC","Shift for negative particles using PDG Mass;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_minus_LHCb_MC = df_MC_complete.Histo2D(*h31,"pT2","shiftminusLHCb");

    histo_data_minus_LHCb_MC->GetXaxis()->SetTitle("p_{T}^{-}");
    histo_data_minus_LHCb_MC->GetYaxis()->SetTitle("#delta p_{T}^{-}");

    histo_data_minus_LHCb_MC->Write();

    // LHCb formula generated energy

    TH2D *h32 = new TH2D("shift_2D_pos_LHCb_MC_gen","Shift for positive particles using PDG Mass MC;p_{T}^{+} #delta p^{+}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_plus_LHCb_MC_gen = df_MC_complete.Histo2D(*h32,"pT1","shiftplusLHCbgen");

    histo_data_plus_LHCb_MC_gen->GetXaxis()->SetTitle("p_{T}^{+}");
    histo_data_plus_LHCb_MC_gen->GetYaxis()->SetTitle("#delta p_{T}^{+}");

    histo_data_plus_LHCb_MC_gen->Write();

    TH2D *h33 = new TH2D("shift_2D_neg_LHCb_MC_gen","Shift for negative particles using PDG Mass;p_{T}^{-} #delta p^{-}",n,0.,2.5,200,-0.1,0.1);
    
    auto histo_data_minus_LHCb_MC_gen = df_MC_complete.Histo2D(*h33,"pT2","shiftminusLHCbgen");

    histo_data_minus_LHCb_MC_gen->GetXaxis()->SetTitle("p_{T}^{-}");
    histo_data_minus_LHCb_MC_gen->GetYaxis()->SetTitle("#delta p_{T}^{-}");

    histo_data_minus_LHCb_MC_gen->Write();

    // momenta
    TH2D *h34 = new TH2D("PTs_MC","PTs",n,0.,4.,200,0.,4.);
    
    auto histo_momenta_MC = df_MC_complete.Histo2D(*h34,"pT1","pT2");

    histo_momenta_MC->GetXaxis()->SetTitle("p_{T}^{+}");
    histo_momenta_MC->GetYaxis()->SetTitle("p_{T}^{-}");

    histo_momenta_MC->Write();


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
    double x[n];
    for (int i = 1; i <=n; i+=1) {
        // calculate x
        // use part of the histogram which corresponds to one bin
        std::cout<<"Iteration "<<i<<std::endl;

        double bin_low=histo_MC_minus_check->GetXaxis()->GetBinLowEdge(i);
        double bin_width=histo_MC_minus_check->GetXaxis()->GetBinWidth(i);
        double bin_up=bin_low+bin_width;

        std::cout<<"Lower edge: "<<bin_low<<"; Upper edge: "<<bin_low<<"; Width: "<<bin_width<<std::endl;
        // middle of the bin
        x[i-1]=bin_low+bin_width/2;
        std::cout<<x[i-1]<<std::endl;
        // error is the distance to next point 
        //xerr[i-1]=bin_width/2;
        // project to y axis
        TH1D *histo_bin_MC_minus_check = histo_MC_minus_check->ProjectionY("MC_neg_difference",i,i);
        TH1D *histo_bin_MC_plus_check = histo_MC_plus_check->ProjectionY("MC_pos_difference",i,i);
        TH1D *histo_bin_MC_minus = histo_MC_minus->ProjectionY("MC_neg",i,i);
        TH1D *histo_bin_MC_plus = histo_MC_plus->ProjectionY("MC_pos",i,i);
        TH1D *histo_bin_minus = histo_data_minus->ProjectionY("data_neg",i,i);
        TH1D *histo_bin_plus = histo_data_plus->ProjectionY("data_pos",i,i);
        std::cout<<"Created projections"<<std::endl;

        myFile->cd();
        gDirectory->mkdir(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));

        myFile->cd();
        myFile->cd(Form("Projection_%.2lf_%.2lf",bin_low,bin_up));
        histo_bin_MC_plus_check->Write();
        histo_bin_MC_plus->Write();
        histo_bin_plus->Write();
        histo_bin_MC_minus_check->Write();
        histo_bin_MC_minus->Write();
        histo_bin_minus->Write();

    }

    myFile->Close();

}