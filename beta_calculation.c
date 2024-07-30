#include <ROOT/RDataFrame.hxx>

// calculate beta for different pT (pT=m*gamma*beta)
double beta(double pT){
    double M=0.497164;
    return sqrt(1./(pow((M/pT),2.)+1.));
}
// calculate weighted beta in one pT bin from pT_low to pT_high
// use 10 bins in pT range and weight beta in every bin by entries
// df needs to contain pT 
double beta_calculation(ROOT::RDataFrame df, double low_pT, double high_pT){
    TH1D *histo= new TH1D("pT_range","used pT range; p_{T} [GeV]; counts",10,low_pT,high_pT);
    auto pTV0 = df.Histo1D(*histo,"pTV0");
    std::cout<<"calculate beta for pT bin from "<<low_pT<<" to "<<high_pT<<std::endl;
    // loop over every bin
    // calculate beta in every pT bin and weight by content
    double sum_beta = 0.;
    int totalcounts=0;
    for (int i=1; i<11; i++){
        // save values for evry bin
        int content=pTV0->GetBinContent(i);
        double pT_bin=pTV0->GetBinCenter(i);
        double beta_bin=beta(pT_bin);
        // sum up counts and weighted beta
        totalcounts+=content;
        sum_beta+=beta_bin*content;
        std::cout<<"pT bin "<<i<<" with center "<<pT_bin<<" beta: "<<beta_bin<<" and entries: "<<content<<std::endl;
        std::cout<<"new values of weighted sum: "<<sum_beta<<" and total number of counts "<<totalcounts<<std::endl;
    }
    // divide by number of counts to get new beta
    double beta_pT=sum_beta/totalcounts;
    // calculate the "old" beta (beta for middle of the bin to compare)
    double old_beta=beta(low_pT+(high_pT-low_pT)/2)
    // write new and old beta in title
    pTV0->SetTitle(Form("#beta_{old}=%.3lf and #beta_{new}=%.3lf",old_beta,beta_pT));
    // plot PT
    pTV0->Write();
    std::cout<<"calculated weighted beta: beta_new="<<beta_pT<<" (beta_old="<<old_beta<<")"<<std::endl;
    return beta_pT;
}