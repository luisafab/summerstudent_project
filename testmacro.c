void testmacro() {
    TCanvas * c1 = new TCanvas("c", "c", 600, 600);
    TF1 *f1 = new TF1("f1","sin(x)",0.,10.);
    f1->Draw();
}