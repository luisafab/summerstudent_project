#include <TF1.h>
#include <TFile.h>
#include <TGenPhaseSpace.h>
#include <TH1D.h>
#include <TLorentzVector.h> // ugly, but TGenPhaseSpace uses this.
#include <TMath.h>
#include <TRandom.h>
#include <ROOT/RDataFrame.hxx>
#include "armenterosPlot.c"

constexpr int kNtrials{100000000};

void genDecay() {
  gRandom->SetSeed(1234);
  //define exponential function to describe distribution of pT
  TF1 mtExpo("mtExpo","[0]*x*std::exp(-std::hypot([2], x)/[1])", 0, 10);
  mtExpo.SetParameter(0, 1.);
  mtExpo.SetParameter(1, 0.5);
  mtExpo.SetParameter(2, kK0sMass);

  // plot exponential
  int n=160;
  double pT_plot[n];
  double exp_plot[n];
  for (int i=0; i<n; i++){
    double pT_temp=i*0.04;
    pT_plot[i]=pT_temp;
    exp_plot[i]=mtExpo(pT_temp);
  }
  auto gr = new TGraph(n,pT_plot,exp_plot);
  gr->SetTitle("Used exponential function");
  gr->GetXaxis()->SetTitle("p_{T}");
  gr->GetYaxis()->SetTitle("f(p_{T})");
  gr->SetLineWidth(2);
  gr->SetLineColor(1);
  
  // create histograms 
  TH1D hV0Pt("hV0Pt", ";#it{p}_{T} (GeV/#it{c});Entries; Calculated pT", 100, 0.5, 1.5);
  TH1D hV0Pt_1("hV0Pt_before", ";#it{p}_{T} (GeV/#it{c});Entries; Calculated pT", 100, 0.5, 1.5);
  TH1D hptotal_1("hPtotal_before", ";#it{p}_{T} (GeV/#it{c});Entries; Input to decay generator", 100, 0.5, 2.);
  TH1D hptotal("hPtotal", ";#it{p}_{Total} (GeV/#it{c});Entries;Calculated total momentum", 100, 0.5, 2.);

  TH1D hPiPPt("hPipPt", ";#it{p}_{T} (GeV/#it{c});Entries", 100, 0, 2);
  TH1D hPiMPt("hpimPt", ";#it{p}_{T} (GeV/#it{c});Entries", 100, 0, 2);
  TH1D hPiPPt_1("hPipPt_before", ";#it{p}_{T} (GeV/#it{c});Entries", 100, 0, 2);
  TH1D hPiMPt_1("hpimPt_before", ";#it{p}_{T} (GeV/#it{c});Entries", 100, 0, 2);

  TH1D hPiPPx("hPiPPx", ";#it{p}_{x} (GeV/#it{c});Entries", 100, -3, 3);
  TH1D hPiMPx("hPiMPx", ";#it{p}_{x} (GeV/#it{c});Entries", 100, -3, 3);
  TH1D hPiPPy("hPiPPy", ";#it{p}_{y} (GeV/#it{c});Entries", 100, -3, 3);
  TH1D hPiMPy("hPiMPy", ";#it{p}_{y} (GeV/#it{c});Entries", 100, -3, 3);
  TH1D hPiPPz("hPiPPz", ";#it{p}_{z} (GeV/#it{c});Entries", 100, -3, 3);
  TH1D hPiMPz("hPiMPz", ";#it{p}_{z} (GeV/#it{c});Entries", 100, -3, 3);
  TH1D hPiPPx_1("hPiPPx_before", ";#it{p}_{x} (GeV/#it{c});Entries;shifted", 100, -3, 3);
  TH1D hPiMPx_1("hPiMPx_before", ";#it{p}_{x} (GeV/#it{c});Entries;shifted", 100, -3, 3);
  TH1D hPiPPy_1("hPiPPy_before", ";#it{p}_{y} (GeV/#it{c});Entries;shifted", 100, -3, 3);
  TH1D hPiMPy_1("hPiMPy_before", ";#it{p}_{y} (GeV/#it{c});Entries;shifted", 100, -3, 3);
  TH1D hPiPPz_1("hPiPPz_before", ";#it{p}_{z} (GeV/#it{c});Entries;shifted", 100, -3, 3);
  TH1D hPiMPz_1("hPiMPz_before", ";#it{p}_{z} (GeV/#it{c});Entries;shifted", 100, -3, 3);

  TH1D halpha("alpha", ";#alpha ;Entries", 100, -1.1, 1.1);
  TH1D hqT("qT", ";#it{q}_{T} (GeV/#it{c});Entries", 100, 0., 0.3);
  TH2D hAP("AP", ";#alpha; #it{q}_{T} (GeV/#it{c});Entries", 150, -1.1, 1.1, 150, 0., 0.3);
  TH1D hMK0("MK0", ";M(K^{0}_{S});Entries", 1000, 0.3, 0.7);

  TH1D halpha_1("alpha_before", ";#alpha ;Entries", 100, -1.1, 1.1);
  TH1D hqT_1("qT_before", ";#it{q}_{T} (GeV/#it{c});Entries", 100, 0., 0.3);
  TH2D hAP_1("AP_before", ";#alpha; #it{q}_{T} (GeV/#it{c});Entries", 150, -1.1, 1.1, 150, 0., 0.3);
  TH1D hMK0_1("MK0_before", ";M(K^{0}_{S});Entries", 1000, 0.3, 0.7);

  
  TH1D hPiPPL("hpipPL", ";#it{p}_{L}^{+} (GeV/#it{c});Entries", 100, 0, 2);
  TH1D hPiNPL("hpimPL", ";#it{p}_{L}^{-} (GeV/#it{c});Entries", 100, 0, 2);


  // vectors to save variables for RDataFrame
  //auto pxp_a=new TObjArray();
  /*float pxp_a[kNtrials];
  float pxn_a[kNtrials];
  float pyp_a[kNtrials];
  float pyn_a[kNtrials];
  float pzp_a[kNtrials];
  float pzn_a[kNtrials];*/

  
  //TFile f("trees_MC_genDecay.root","recreate");
  //TTree t1("mctable","generated momenta");
  TFile* f = new TFile("trees_MC_genDecay.root","recreate");
  TTree* t1 = new TTree("mctable","generated momenta");
  Float_t pxp, pyp, pzp, pxn, pyn, pzn;
  t1->Branch("fPxPosMC",&pxp,"pxp/F");
  t1->Branch("fPyPosMC",&pyp,"pyp/F");
  t1->Branch("fPzPosMC",&pzp,"pzp/F");
  t1->Branch("fPxNegMC",&pxn,"pxn/F");
  t1->Branch("fPyNegMC",&pyn,"pyn/F");
  t1->Branch("fPzNegMC",&pzn,"pzn/F");

  // fill the tree
   // save the Tree heade; the file will be automatically closed
   // when going out of the function scope
  
  
  // lorentz vectors to save particles
  TLorentzVector mother, pip, pim;
  // GenPhaseSpace to generate decay
  TGenPhaseSpace genPi;
  TGenPhaseSpace genPi_test;
  //masses of resulting particles
  const double massesPi[2]{kPiMass, kPiMass};
  
  // simulate decay
  // loop over number of trials
  for (int i=0; i < kNtrials; ++i) {
    // initialise variables 
    // shift and smirring of 1/pT
    float pshift=0.;
    float psmirr=0.;
    // pT eta and phi of mother
    float pT_original=mtExpo.GetRandom();
    float eta=gRandom->Uniform(-1, 1);
    float phi=gRandom->Uniform(0, TMath::TwoPi());

    float pT_=1./pT_original;
    float pT_new=gRandom->Gaus(pT_+pshift,pT_*psmirr);
    float pT_shift=1/pT_new;

    // generate mother and decay
    mother.SetPtEtaPhiM(pT_shift, eta, phi, kK0sMass);
    genPi.SetDecay(mother, 2, massesPi);
    genPi.Generate();

    // p+ and p- of daughters
    float pxp_1=genPi.GetDecay(0)->Px();
    float pxn_1=genPi.GetDecay(1)->Px();
    float pyp_1=genPi.GetDecay(0)->Py();
    float pyn_1=genPi.GetDecay(1)->Py();
    float pzp_1=genPi.GetDecay(0)->Pz();
    float pzn_1=genPi.GetDecay(1)->Pz();
    
    // get eta and phi for the daughters
    float etaPlus=genPi.GetDecay(0)->Eta();
    float etaMinus=genPi.GetDecay(1)->Eta();
    float phiPlus=genPi.GetDecay(0)->Phi();
    float phiMinus=genPi.GetDecay(1)->Phi();
    // get pT of the particles 
    float pT_pos=1./genPi.GetDecay(0)->Pt();
    float pT_neg=1./genPi.GetDecay(1)->Pt();
    float pT_pos_new=gRandom->Gaus(pT_pos+pshift,pT_pos*psmirr);
    float pT_neg_new=gRandom->Gaus(pT_neg+pshift,pT_neg*psmirr);
    // use the new 1/pT, eta and phi 
    // assign lorentz vector to generate components 
    TLorentzVector pip_new, pim_new;
    pip_new.SetPtEtaPhiM(1/pT_pos_new, etaPlus, phiPlus, kPiMass);
    pim_new.SetPtEtaPhiM(1/pT_neg_new, etaMinus, phiMinus, kPiMass);
    pxp=pip_new.Px();
    pxn=pim_new.Px();
    pyp=pip_new.Py();
    pyn=pim_new.Py();
    pzp=pip_new.Pz();
    pzn=pim_new.Pz();
    // cut on ptotal to make it comparable
    float ptotal_calc=ptotal(pxp,pyp,pzp,pxn,pyn,pzn);
    /*if (ptotal_calc<0.95 || ptotal_calc>1.05){
      continue;
    }*/
    // calculate pt of V0 particle 
    float ptV0_calc=pt_v0(pxp,pyp,pzp,pxn,pyn,pzn);
    float alpha_calc=alpha(pxp,pyp,pzp,pxn,pyn,pzn);
    float qT_calc=qT_AP(pxp,pyp,pzp,pxn,pyn,pzn);

    float ptV0_calc_before=pt_v0(pxp_1,pyp_1,pzp_1,pxn_1,pyn_1,pzn_1);
    float alpha_calc_before=alpha(pxp_1,pyp_1,pzp_1,pxn_1,pyn_1,pzn_1);
    float qT_calc_before=qT_AP(pxp_1,pyp_1,pzp_1,pxn_1,pyn_1,pzn_1);
    
    float mass_calc=invariantmass_K0(pxp,pyp,pzp,pxn,pyn,pzn);
    float mass_calc_before=invariantmass_K0(pxp_1,pyp_1,pzp_1,pxn_1,pyn_1,pzn_1);
    float ptotal_calc_before=ptotal(pxp_1,pyp_1,pzp_1,pxn_1,pyn_1,pzn_1);
  

    hV0Pt.Fill(ptV0_calc);
    hV0Pt_1.Fill(ptV0_calc_before);
    hptotal.Fill(ptotal_calc);
    hptotal_1.Fill(ptotal_calc_before);

    hPiPPt.Fill(pip_new.Pt());
    hPiMPt.Fill(pim_new.Pt());
    hPiPPt_1.Fill(genPi.GetDecay(0)->Pt());
    hPiMPt_1.Fill(genPi.GetDecay(1)->Pt());

    hPiPPx_1.Fill(pxp_1);
    hPiMPx_1.Fill(pxn_1);
    hPiPPy_1.Fill(pyp_1);
    hPiMPy_1.Fill(pyn_1);
    hPiPPz_1.Fill(pzp_1);
    hPiMPz_1.Fill(pzn_1);
    hPiPPx.Fill(pxp);
    hPiMPx.Fill(pxn);
    hPiPPy.Fill(pyp);
    hPiMPy.Fill(pyn);
    hPiPPz.Fill(pzp);
    hPiMPz.Fill(pzn);

    hAP.Fill(alpha_calc,qT_calc);
    hAP_1.Fill(alpha_calc_before,qT_calc_before);
    hqT.Fill(qT_calc);
    halpha.Fill(alpha_calc);
    hMK0.Fill(mass_calc);
    hqT_1.Fill(qT_calc_before);
    halpha_1.Fill(alpha_calc_before);
    hMK0_1.Fill(mass_calc_before);
    t1->Fill();
  }
  t1->Write();
  
  std::unique_ptr<TFile> myFile(TFile::Open("genDecay_test.root","RECREATE") );
  // create directon for this pT value
  myFile->cd();
  gDirectory->mkdir("components");
  myFile->cd();
  gDirectory->mkdir("components_before");
  myFile->cd();
  gDirectory->mkdir("V0");
  myFile->cd();
  gDirectory->mkdir("pT_daughters");
  myFile->cd();
  gDirectory->mkdir("APComponents");
  // plot used exponential
  gr->Write("exponential");
  // Pt and ptotal of V0
  myFile->cd();
  myFile->cd("V0");

  hV0Pt.Write();
  hptotal.Write();
  hV0Pt_1.Write();
  hptotal_1.Write();

  myFile->cd();
  // AP before and after
  hAP.Write();
  hAP_1.Write();
  
  // cariables for p+ and p-
  // pt before and after
  myFile->cd();
  myFile->cd("pT_daughters");
  hPiPPt.Write();
  hPiMPt.Write();
  hPiPPt_1.Write();
  hPiMPt_1.Write();

  // x, y, z coponents before and after
  myFile->cd();
  myFile->cd("components");
  hPiPPx.Write();
  hPiMPx.Write();
  hPiPPy.Write();
  hPiMPy.Write();
  hPiPPz.Write();
  hPiMPz.Write();

  myFile->cd();
  myFile->cd("components_before");
  hPiPPx_1.Write();
  hPiMPx_1.Write();
  hPiPPy_1.Write();
  hPiMPy_1.Write();
  hPiPPz_1.Write();
  hPiMPz_1.Write();
  // variables used for final AP plot
  myFile->cd();
  hMK0.Write();
  hMK0_1.Write();
  myFile->cd();
  myFile->cd("APComponents");
  halpha.Write();
  hqT.Write();
  halpha_1.Write();
  hqT_1.Write();

  
  myFile->Close();
}