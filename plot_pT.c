#include <ROOT/RDataFrame.hxx>

// function to calculate difference in 1/pT for MC and data 
// Snapshot needs to be one of MC data, containing pT for data and MC
void plot_pT(TString file="SnapshotMC.root",TString tree="NewVariables"){

    std::unique_ptr<TFile> myFile( TFile::Open("pTplot.root", "RECREATE") );

    ROOT::RDataFrame df(tree, file);

    auto oneoverpT =[&](float pT){
        return 1./pT;
    };


    auto new_df = df.Define("oneoverpT","1./pTreco");
    auto new_column = new_df.Define("oneoverpTgen","1./pTV0");
    auto complete_df = new_column.Define("difference","oneoverpT-oneoverpTgen");

    // save Snapshot to use for pT binning
    complete_df.Snapshot("pTs","Snapshot_pTs.root",{"difference","oneoverpT","oneoverpTgen"});

    
    TH2D *h = new TH2D("pTPlot","Compare generated and calculated data;1/pT;1/pT-1/pT",500,0,3.5,500,-0.23,0.23);
    
    auto histo = complete_df.Histo2D(*h,"oneoverpT","difference");

    histo->GetXaxis()->SetTitle("#frac{1}{p_{T}_{reco}}");
    histo->GetYaxis()->SetTitle("#frac{1}{p_{T}_{reco}}-#frac{1}{p_{T}_{gen}}");

    histo->Write();

    TH2D *h5 = new TH2D("pTPlot2","Compare generated and calculated data;1/pT;1/pT-1/pT",500,0,3.5,500,-0.23,0.23);
    
    auto histo5 = complete_df.Histo2D(*h5,"oneoverpTgen","difference");

    histo5->GetXaxis()->SetTitle("#frac{1}{p_{T}_{gen}}");
    histo5->GetYaxis()->SetTitle("#frac{1}{p_{T}_{reco}}-#frac{1}{p_{T}_{gen}}");

    histo5->Write();


    TH1D *h1 = new TH1D("oneoverpT","1/pT;1/pT;1/pT-1/pT",100,0,3.5);
    
    auto histo1 = complete_df.Histo1D(*h1,"oneoverpT");

    histo1->GetXaxis()->SetTitle("#frac{1}{p_{T}_{reco}}");

    histo1->Write();

    TH1D *h2 = new TH1D("oneoverpTgen","1/pT;1/pT;1/pT-1/pT",100,0,3.5);
    
    auto histo2 = complete_df.Histo1D(*h2,"oneoverpTgen");

    histo2->GetXaxis()->SetTitle("#frac{1}{p_{T}_{gen}}");

    histo2->Write();

    TH1D *h3 = new TH1D("pTreco","pTV0;1/pT;1/pT-1/pT",100,0,3.5);
    
    auto histo3 = complete_df.Histo1D(*h3,"pTreco");

    histo3->GetXaxis()->SetTitle("p_{T}_{reco}");

    histo3->Write();

    TH1D *h4 = new TH1D("pTgen","pTgen;1/pT;1/pT-1/pT",100,0,3.5);
    
    auto histo4 = complete_df.Histo1D(*h4,"pTV0");

    histo4->GetXaxis()->SetTitle("p_{T}_{gen}");

    histo4->Write();

    myFile->Close();
}