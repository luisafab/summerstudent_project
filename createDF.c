#include <ROOT/RDataFrame.hxx>
#include "armenterosPlot.c"


// this function creates a dataframe with the calculated values and saves them in a dataframe to use 
// the tree is different for MC data, data, and generated data with the function gendecay, also the resulting dataframe is different
// two booleans to indicate which tree is used

//void createDF(TString file="data/AnalysisResults_treesAP_data_LHC22o_apass6_small.root",TString tree="O2v0tableap", bool MC = false){
//void createDF(TString file="data/trees_test_k0_newMC.root",TString tree="DF_2263624207437933/O2mcv0tableap", int mode=0){
void createDF(TString file="data/trees_MC.root",TString tree="DF_2261906121493398/O2mcv0tableap", bool MC=true, bool genDecay=false){
//void createDF(TString file="trees_MC_genDecay.root",TString tree="mctable", bool MC=true){
    
    // define a data frame to use 
    ROOT::RDataFrame df_MC(tree, file);
    // use different momenta for MC and not MC
    vector<basic_string<char> > labels1 = {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"};
    vector<basic_string<char> > labels2 = {"alpha","fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"};
    if (MC){
        labels1 = {"fPxPosMC", "fPyPosMC","fPzPosMC","fPxNegMC", "fPyNegMC","fPzNegMC"};
        labels2 = {"alpha","fPxPosMC", "fPyPosMC","fPzPosMC","fPxNegMC", "fPyNegMC","fPzNegMC"};
    }
    // filter it
    //     for MC==true only the recombined (isreco==true) lambdas (PDGCODE=3122) or K0s(310) are considered
    //     else there is a dummy cut to use the same dataframe
    //auto df_cos = df_MC.Filter(MC? "(fPDGCode == 3122 || fPDGCode== -3122) && fIsReco" : "fLen>0");
    auto df_cos = df_MC.Filter(MC? "(fPDGCode == 310) && fIsReco" : "fLen>0");
    std::cout<<"filtered PDG"<<std::endl;
    // apply cut on cosPA and dummy cut for gendecay
    auto df = df_MC.Filter(genDecay? "fPxPosMC>-1000" : "fCosPA>0.999");

    // calculate values
    // variables of AP plot: alpha and qT
    auto new_column_alpha = df.Define("alpha", alpha, labels1);
    std::cout<<"Calculated alpha"<<std::endl;
    auto new_column_qT = new_column_alpha.Define("qT", qT_AP, labels1);
    std::cout<<"calculated qT"<<std::endl;

    // calculate invariant mass to cut on 
    // add both masses as new column to cut on 
    auto new_minv_l = new_column_qT.Define("Minv_lambda", invariantmass_lambda, labels2);
    std::cout<<"Calculated inv mass lambda"<<std::endl;
    auto new_minv_k0 = new_minv_l.Define("Minv_K0", invariantmass_K0, labels1);
    std::cout<<"Calculated invariant mass K0"<<std::endl;
    // add filter on calculated invariant mass
    float lambda_low=1.1;
    float lambda_high=1.13;
    float K0_low=0.49;
    float K0_high=0.502;
    // use lambda function
    auto insidePeakK =[&](float Minv){
        bool isInside;
        float low=K0_low; 
        float high=K0_high;
        if (Minv>low && Minv < high){
            isInside=true;
        }else{
            isInside=false;
        }
        return isInside;
    };
    auto insidePeakL =[&](float Minv){
        bool isInside;
        float low=lambda_low; 
        float high=lambda_high;
        if (Minv>low && Minv < high){
            isInside=true;
        }else{
            isInside=false;
        }
        return isInside;
    };
    auto outsidePeakL =[&](float Minv){
        return !insidePeakL(Minv);
    };
    
    // filter on K0s
    auto new_filter_k = new_minv_k0.Filter(insidePeakK, {"Minv_K0"});
    std::cout<<"Filter 2"<<std::endl;
    // filter on particles which are in both peaks to reduce possible background 
    auto new_filter_both = new_filter_k.Filter(outsidePeakL, {"Minv_lambda"});
    std::cout<<"Filter 3"<<std::endl;

    // add beta and ptotal as variables
    auto new_ptotal = new_filter_both.Define("ptotal", ptotal, labels1);
    std::cout<<"Calculated total momentum"<<std::endl;
    auto complete_df = new_ptotal.Define("beta", beta, {"ptotal"});
    std::cout<<"Calculated beta using the total momentum of V0s"<<std::endl;

    // generate a random number for to sample data later
    auto new_column_random = complete_df.Define("rnd", "gRandom->Rndm()");
    std::cout<<"Calculated random"<<std::endl;

    // calculate pT of V0
    auto new_column_pT = new_column_random.Define("pTV0", pt_v0, labels1);
    std::cout<<"calculated pT V0"<<std::endl;
    // save dataframe
    if (MC && genDecay) {
        new_column_pT.Snapshot("NewVariables","SnapshotGenDecay.root");
    }else if(MC){
        // calculate also pTreco for MCParticles to compare
        auto new_column_pTreco = new_column_pT.Define("pTreco", pt_v0, {"fPxPos", "fPyPos","fPzPos","fPxNeg", "fPyNeg","fPzNeg"});
        new_column_pTreco.Snapshot("NewVariables","SnapshotMC.root");
    }else{
        new_column_pT.Snapshot("NewVariables", "test.root");
    }

}