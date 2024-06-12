using namespace ROOT; // RDataFrame's namespace
void dataframe_test() {
    // ein dataframe definieren (name des trees und der Datei)
    //RDataFrame d1("O2v0tableap", "data/AnalysisResults_treesAP_data.root");
    RDataFrame d1("O2mcv0tableap", "data/AnalysisResults_treesAP.root");
    
    
    
    // Cuts definieren
    // Beispielcut für fLen
    auto cutPDG = [](int PDG) { return (PDG == 3122 || PDG== -3122); };
    auto cutIsRecoPDG = [](bool IsReco, int PDG) { return ((PDG == 3122 || PDG== -3122) && IsReco); };
    // Beispielcut für fEta
    //auto cutb1b2 = [](int b2, double b1) { return b2 % 2 && b1 < 4.; };

    auto df_MCcuts = d1.Filter(cutIsRecoPDG,{"fIsReco","fPDGCode"}); // <- no column name specified here!
                      //.Filter(cutb1b2, {"b2", "b1"})
    //               .Count();

    // Histogramme mit den Variablen
    //auto h = d1.Histo1D("fPDGCode");
    auto h1 = df_MCcuts.Histo1D("fPDGCode");
    //auto ppos=df_MCcuts.Take("fPDGCode");
    // Histogramme zeichnen
    // es muss zuerst ein Canvas definiert werden und es braucht einen work around für den pointer
    //TCanvas * c1 = new TCanvas("c", "c", 600, 600);
    //h->DrawClone();
    //gROOT->GetListOfCanvases()->Draw();
    TCanvas * c2 = new TCanvas("c2", "c2", 600, 600);
    h1->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    //std::cout << fPDGCode << " entries passed all filters" << std::endl;
}