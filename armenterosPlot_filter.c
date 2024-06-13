#include <ROOT/RDataFrame.hxx>

//calculating the armenteros-podolanski plot for given data 
// string file : name of file
// string tree : name of tree
// bool MC : true : only lambda events are consicered (!! The MC Information in the file is necessary !!)
//           false : all events are considered

void armenterosPlot_filter(TString file="data/AnalysisResults_treesAP_data.root", TString tree="O2v0tableap",bool MC=false) {

    float massProton=0.9382720813;
    float massPion= 0.13957061;

    // define a data frame to use 
    ROOT::RDataFrame df_MC(tree, file);
    // filter it
    //     for MC==true only the recombined (isreco==true) lambdas (PDGCODE=3122) are considered
    //     else there is a dummy cut to use the same dataframe
    auto df = df_MC.Filter(MC? "(fPDGCode == 3122 || fPDGCode== -3122) && fIsReco" : "fLen>0");

    // lambda functions to calculate the variables of the AP Plot: alpha and pT
    // takes the momenta of the negative and positive particle as argument
    auto alpha = [&](float posx, float posy, float posz, float negx, float negy, float negz){
        // create momentum vector for both particles
        TVector3 pPos(posx,posy,posz);
        TVector3 pNeg(negx,negy,negz);
        // the momentum of V0 is the sum of the other two
        TVector3 pV0 =pPos+pNeg;
        // normalize pV0
        float pV0_norm=sqrt(pV0.Dot(pV0));
        // calculate longitudinal part of pPos and pNeg
        float pPlL = pPos.Dot(pV0)/pV0_norm;
        float pNegL = pNeg.Dot(pV0)/pV0_norm;
        // calculate alpha
        float alpha= (pPlL-pNegL)/(pPlL+pNegL);

        return alpha;
    };

    auto pT = [&](float posx, float posy, float posz, float negx, float negy, float negz){
        // create momentum vector for both particles
        TVector3 pPos(posx,posy,posz);
        TVector3 pNeg(negx,negy,negz);
        // the momentum of V0 is the sum of the other two
        TVector3 pV0 =pPos+pNeg;
        // normalize pV0
        float pV0_norm=sqrt(pV0.Dot(pV0));
        // cross product of pPos and pV0
        TVector3 cross=pPos.Cross(pV0);
        // calculate pT
        float pT = (sqrt(cross.Dot(cross)))/pV0_norm;

        return pT;
    };

    auto invariantmass = [&](float alpha, float posx, float posy, float posz, float negx, float negy, float negz){
        float massPos;
        float massNeg;

        //define vectors for momenta and calclate the length 
        TVector3 pPos(posx,posy,posz);
        TVector3 pNeg(negx,negy,negz);
        TVector3 momentum=pPos+pNeg;
        float pPos_norm=sqrt(pPos.Dot(pPos));
        float pNeg_norm=sqrt(pNeg.Dot(pNeg));
        float momentum_norm=sqrt(momentum.Dot(momentum));


        // use alpha to know the decay and the resulting particles
        // alpha>0 lambda -> p +pi-
        // alpha<0 pi+
        // assign the mass 
        if (alpha>0){
            massPos=massProton;
            massNeg=massPion;
        }
        else{
            massPos=massPion;
            massNeg=massProton;
        }

        // use mass and momentum to calculate the energy of the two particles
        float EPos=sqrt(pow(pPos_norm,2.)+pow((massPos),2.));
        float ENeg=sqrt(pow(pNeg_norm,2.)+pow((massNeg),2.));

        // total energy 
        float Energy = EPos + ENeg;

        // calculate the invariant mass using the fourmomentum vector
        float invariantmass = sqrt((pow(Energy,2.)-pow(momentum_norm,2.)));

        return invariantmass;
    };

    // define new columns with the new variables 
    auto new_alpha = df.Define("alpha", alpha, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    auto new_pT = new_alpha.Define("pT", pT, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
    auto all_df = new_pT.Define("Minv", invariantmass, {"alpha","fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});

    // filter on invariant mass
    auto complete_df = all_df.Filter("(Minv>1) && (Minv<2)");


    // plot the two new variables
    auto h1 = complete_df.Histo1D("alpha");
    TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
    h1->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    auto h2 = complete_df.Histo1D("pT");
    TCanvas * c2 = new TCanvas("c2", "c2", 600, 600);
    h2->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    auto h3 = complete_df.Histo1D("Minv");
    TCanvas * c0 = new TCanvas("c2", "c2", 600, 600);
    h3->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

    // define and plot a histogram which is the armenteros plot
    TH2D *h = new TH2D("h","Armenteros-Podolanski Plot",100,-1.,1.,40,0.,0.2);
    
    auto histo = complete_df.Histo2D(*h,"alpha","pT");

    histo->GetXaxis()->SetTitle("#alpha");
    histo->GetYaxis()->SetTitle("p_{T}");


    TCanvas * c3 = new TCanvas("c3", "c3", 900, 600);
    histo->DrawClone();
    gROOT->GetListOfCanvases()->Draw();

}