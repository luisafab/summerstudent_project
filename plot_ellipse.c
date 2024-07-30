// calculate values for corresponding ellipse
#include <ROOT/RDataFrame.hxx>


// function to plot ellipse for given m,M and beta
void plot_ellipse(double M, double m, double beta) {
    int n=160;
    double alpha_plot[n];
    double qT_plot[n];
    for (int i=0; i<n; i++){
        double alpha_temp=(-0.8+i*0.01);
        alpha_plot[i]=alpha_temp;
        qT_plot[i]=qT(alpha_temp,M,m,beta);
    }
    auto ellipse = new TGraph(n,alpha_plot,qT_plot);
    ellipse->SetTitle(Form("Ellipse with #beta = %.3lf and M=%.4lf",beta, M));
    ellipse->GetXaxis()->SetTitle("#alpha");
    ellipse->GetYaxis()->SetTitle("q_{T}");
    ellipse->SetLineWidth(3);
    ellipse->SetLineColor(1);
    ellipse->Write(Form("Ellipse_M%.4lf",M));

}