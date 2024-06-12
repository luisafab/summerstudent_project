#include <ROOT/RDataFrame.hxx>

int test() {
    // Erstellen eines RDataFrame-Objekts ohne explizite Spaltennamen
    ROOT::RDataFrame df("O2v0tableap", "data/AnalysisResults_treesAP_data.root");

    // Definieren einer neuen Spalte mit einer Lambda-Funktion
    // auto lambda_function = [&](int PDGCode, float Len) {
    auto alpha = [&](float posx, float posy, float posz, float negx, float negy, float negz){
        TVector3 pPos(posx,posy,posz);
        TVector3 pNeg(negx,negy,negz);
        TVector3 pV0 =pPos+pNeg;
        float pV0_norm=sqrt(pV0.Dot(pV0));
        float pPlL = pPos.Dot(pV0)/pV0_norm;
        float pNegL = pNeg.Dot(pV0)/pV0_norm;
        float alpha= (pPlL-pNegL)/(pPlL+pNegL);

        return alpha;

    };

    auto pT = [&](float posx, float posy, float posz, float negx, float negy, float negz){
        TVector3 pPos(posx,posy,posz);
        TVector3 pNeg(negx,negy,negz);
        TVector3 pV0 =pPos+pNeg;
        float pV0_norm=sqrt(pV0.Dot(pV0));
        float pPlL = pPos.Dot(pV0)/pV0_norm;
        float pNegL = pNeg.Dot(pV0)/pV0_norm;
        TVector3 cross=pPos.Cross(pV0);
        float pT = (sqrt(cross.Dot(cross)))/pV0_norm;

        return pT;

    };
    auto new_column = df.Define("alpha", alpha, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    auto complete_column = new_column.Define("pT", pT, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});

    // Ausgabe der neu erstellten Spalte
    //new_column.Snapshot("output.root", {"new_column"});

    auto h1 = complete_column.Histo1D("alpha");
    TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
    h1->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    auto h2 = complete_column.Histo1D("pT");
    TCanvas * c2 = new TCanvas("c2", "c2", 600, 600);
    h2->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    // define and plot a histogram which is the armenteros plot
    TH2D *h = new TH2D("h","Armenteros-Podolanski Plot",100,-1.,1.,40,0.,0.2);
    
    //TH1D *h = new TH1D("h","h",10,-1.,1.);
    auto histo = complete_column.Histo2D(*h,"alpha","pT");

    histo->GetXaxis()->SetTitle("alpha");
    histo->GetYaxis()->SetTitle("pT");

    TCanvas * c3 = new TCanvas("c3", "c3", 900, 600);
    histo->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    //auto h3 = complete_column.Histo2D("alpha","pT");
    //TCanvas * c3 = new TCanvas("c3", "c3", 600, 600);
    //h3->DrawClone();
    //gROOT->GetListOfCanvases()->Draw();


    return 0;
}