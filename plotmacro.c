void plotmacro() {
    TF1 efunc("efunc","exp([0]+[1]*x)",0.,5.);
    efunc.SetParameter(0,1);
    efunc.SetParameter(1,-1);
    TH1F h("h","example histogram",100,0.,5.);
    for (int i=0;i<1000;i++) {
        h.Fill(efunc.GetRandom());
    }
    h.Draw();
}