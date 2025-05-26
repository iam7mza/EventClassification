#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        cerr << "Usage: " << argv[0] << " <root_file> [sample_size]" << endl;
        cerr << "  root_file: Input ROOT file" << endl;
        cerr << "  sample_size: Number of events to randomly sample (optional, default: all events)" << endl;
        return 1;
    }
    
    string filename = argv[1];
    Long64_t sampleSize = -1; // -1 means use all events
    
    if (argc == 3) {
        sampleSize = stoll(argv[2]);
        if (sampleSize <= 0) {
            cerr << "Error: Sample size must be positive" << endl;
            return 1;
        }
    }
    
    cout << "Opening file: " << filename << endl;
    
    TFile *file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return 1;
    }
    
    TTree *tree = (TTree*)file->Get("mini");
    if (!tree) {
        cerr << "Error: Cannot find tree in file" << endl;
        cout << "Available keys in file:" << endl;
        file->ls();
        file->Close();
        return 1;
    }
    
    // Print available branches
    cout << "Available branches:" << endl;
    TObjArray *branches = tree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); i++) {
        cout << branches->At(i)->GetName() << endl;
    }
    
    // Variables to store all data
    // Basic event info
    int runNumber = 0;
    int eventNumber = 0;
    int channelNumber = 0;
    float mcWeight = 0;
    float XSection = 0;
    float SumWeights = 0;
    
    // Scale factors
    float scaleFactor_PILEUP = 0;
    float scaleFactor_ELE = 0;
    float scaleFactor_MUON = 0;
    float scaleFactor_PHOTON = 0;
    float scaleFactor_TAU = 0;
    float scaleFactor_BTAG = 0;
    float scaleFactor_LepTRIGGER = 0;
    float scaleFactor_PhotonTRIGGER = 0;
    
    // Triggers
    bool trigE = false;
    bool trigM = false;
    bool trigP = false;
    
    // Leptons
    unsigned int lep_n = 0;
    vector<bool> *lep_truthMatched = nullptr;
    vector<bool> *lep_trigMatched = nullptr;
    vector<float> *lep_pt = nullptr;
    vector<float> *lep_eta = nullptr;
    vector<float> *lep_phi = nullptr;
    vector<float> *lep_E = nullptr;
    vector<float> *lep_z0 = nullptr;
    vector<int> *lep_charge = nullptr;
    vector<int> *lep_type = nullptr;
    vector<bool> *lep_isTightID = nullptr;
    vector<float> *lep_ptcone30 = nullptr;
    vector<float> *lep_etcone20 = nullptr;
    vector<float> *lep_trackd0pvunbiased = nullptr;
    vector<float> *lep_tracksigd0pvunbiased = nullptr;
    
    // Missing ET
    float met_et = 0;
    float met_phi = 0;
    
    // Jets
    unsigned int jet_n = 0;
    vector<float> *jet_pt = nullptr;
    vector<float> *jet_eta = nullptr;
    vector<float> *jet_phi = nullptr;
    vector<float> *jet_E = nullptr;
    vector<float> *jet_jvt = nullptr;
    vector<int> *jet_trueflav = nullptr;
    vector<bool> *jet_truthMatched = nullptr;
    vector<float> *jet_MV2c10 = nullptr;
    
    // Photons
    unsigned int photon_n = 0;
    vector<bool> *photon_truthMatched = nullptr;
    vector<bool> *photon_trigMatched = nullptr;
    vector<float> *photon_pt = nullptr;
    vector<float> *photon_eta = nullptr;
    vector<float> *photon_phi = nullptr;
    vector<float> *photon_E = nullptr;
    vector<bool> *photon_isTightID = nullptr;
    vector<float> *photon_ptcone30 = nullptr;
    vector<float> *photon_etcone20 = nullptr;
    vector<int> *photon_convType = nullptr;
    
    // Large R jets
    unsigned int largeRjet_n = 0;
    vector<float> *largeRjet_pt = nullptr;
    vector<float> *largeRjet_eta = nullptr;
    vector<float> *largeRjet_phi = nullptr;
    vector<float> *largeRjet_E = nullptr;
    vector<float> *largeRjet_m = nullptr;
    vector<int> *largeRjet_truthMatched = nullptr;
    vector<float> *largeRjet_D2 = nullptr;
    vector<float> *largeRjet_tau32 = nullptr;
    
    // Taus
    unsigned int tau_n = 0;
    vector<float> *tau_pt = nullptr;
    vector<float> *tau_eta = nullptr;
    vector<float> *tau_phi = nullptr;
    vector<float> *tau_E = nullptr;
    vector<int> *tau_charge = nullptr;
    vector<bool> *tau_isTightID = nullptr;
    vector<bool> *tau_truthMatched = nullptr;
    vector<bool> *tau_trigMatched = nullptr;
    vector<int> *tau_nTracks = nullptr;
    vector<float> *tau_BDTid = nullptr;
    
    // Di-tau
    float ditau_m = 0;
    
    // Systematic uncertainties
    vector<float> *lep_pt_syst = nullptr;
    float met_et_syst = 0;
    vector<float> *jet_pt_syst = nullptr;
    vector<float> *photon_pt_syst = nullptr;
    vector<float> *largeRjet_pt_syst = nullptr;
    vector<float> *tau_pt_syst = nullptr;
    
    // Set branch addresses - Basic event info
    tree->SetBranchAddress("runNumber", &runNumber);
    tree->SetBranchAddress("eventNumber", &eventNumber);
    tree->SetBranchAddress("channelNumber", &channelNumber);
    tree->SetBranchAddress("mcWeight", &mcWeight);
    tree->SetBranchAddress("XSection", &XSection);
    tree->SetBranchAddress("SumWeights", &SumWeights);
    
    // Scale factors
    tree->SetBranchAddress("scaleFactor_PILEUP", &scaleFactor_PILEUP);
    tree->SetBranchAddress("scaleFactor_ELE", &scaleFactor_ELE);
    tree->SetBranchAddress("scaleFactor_MUON", &scaleFactor_MUON);
    tree->SetBranchAddress("scaleFactor_PHOTON", &scaleFactor_PHOTON);
    tree->SetBranchAddress("scaleFactor_TAU", &scaleFactor_TAU);
    tree->SetBranchAddress("scaleFactor_BTAG", &scaleFactor_BTAG);
    tree->SetBranchAddress("scaleFactor_LepTRIGGER", &scaleFactor_LepTRIGGER);
    tree->SetBranchAddress("scaleFactor_PhotonTRIGGER", &scaleFactor_PhotonTRIGGER);
    
    // Triggers
    tree->SetBranchAddress("trigE", &trigE);
    tree->SetBranchAddress("trigM", &trigM);
    tree->SetBranchAddress("trigP", &trigP);
    
    // Leptons
    tree->SetBranchAddress("lep_n", &lep_n);
    tree->SetBranchAddress("lep_truthMatched", &lep_truthMatched);
    tree->SetBranchAddress("lep_trigMatched", &lep_trigMatched);
    tree->SetBranchAddress("lep_pt", &lep_pt);
    tree->SetBranchAddress("lep_eta", &lep_eta);
    tree->SetBranchAddress("lep_phi", &lep_phi);
    tree->SetBranchAddress("lep_E", &lep_E);
    tree->SetBranchAddress("lep_z0", &lep_z0);
    tree->SetBranchAddress("lep_charge", &lep_charge);
    tree->SetBranchAddress("lep_type", &lep_type);
    tree->SetBranchAddress("lep_isTightID", &lep_isTightID);
    tree->SetBranchAddress("lep_ptcone30", &lep_ptcone30);
    tree->SetBranchAddress("lep_etcone20", &lep_etcone20);
    tree->SetBranchAddress("lep_trackd0pvunbiased", &lep_trackd0pvunbiased);
    tree->SetBranchAddress("lep_tracksigd0pvunbiased", &lep_tracksigd0pvunbiased);
    
    // Missing ET
    tree->SetBranchAddress("met_et", &met_et);
    tree->SetBranchAddress("met_phi", &met_phi);
    
    // Jets
    tree->SetBranchAddress("jet_n", &jet_n);
    tree->SetBranchAddress("jet_pt", &jet_pt);
    tree->SetBranchAddress("jet_eta", &jet_eta);
    tree->SetBranchAddress("jet_phi", &jet_phi);
    tree->SetBranchAddress("jet_E", &jet_E);
    tree->SetBranchAddress("jet_jvt", &jet_jvt);
    tree->SetBranchAddress("jet_trueflav", &jet_trueflav);
    tree->SetBranchAddress("jet_truthMatched", &jet_truthMatched);
    tree->SetBranchAddress("jet_MV2c10", &jet_MV2c10);
    
    // Photons
    tree->SetBranchAddress("photon_n", &photon_n);
    tree->SetBranchAddress("photon_truthMatched", &photon_truthMatched);
    tree->SetBranchAddress("photon_trigMatched", &photon_trigMatched);
    tree->SetBranchAddress("photon_pt", &photon_pt);
    tree->SetBranchAddress("photon_eta", &photon_eta);
    tree->SetBranchAddress("photon_phi", &photon_phi);
    tree->SetBranchAddress("photon_E", &photon_E);
    tree->SetBranchAddress("photon_isTightID", &photon_isTightID);
    tree->SetBranchAddress("photon_ptcone30", &photon_ptcone30);
    tree->SetBranchAddress("photon_etcone20", &photon_etcone20);
    tree->SetBranchAddress("photon_convType", &photon_convType);
    
    // Large R jets
    tree->SetBranchAddress("largeRjet_n", &largeRjet_n);
    tree->SetBranchAddress("largeRjet_pt", &largeRjet_pt);
    tree->SetBranchAddress("largeRjet_eta", &largeRjet_eta);
    tree->SetBranchAddress("largeRjet_phi", &largeRjet_phi);
    tree->SetBranchAddress("largeRjet_E", &largeRjet_E);
    tree->SetBranchAddress("largeRjet_m", &largeRjet_m);
    tree->SetBranchAddress("largeRjet_truthMatched", &largeRjet_truthMatched);
    tree->SetBranchAddress("largeRjet_D2", &largeRjet_D2);
    tree->SetBranchAddress("largeRjet_tau32", &largeRjet_tau32);
    
    // Taus
    tree->SetBranchAddress("tau_n", &tau_n);
    tree->SetBranchAddress("tau_pt", &tau_pt);
    tree->SetBranchAddress("tau_eta", &tau_eta);
    tree->SetBranchAddress("tau_phi", &tau_phi);
    tree->SetBranchAddress("tau_E", &tau_E);
    tree->SetBranchAddress("tau_charge", &tau_charge);
    tree->SetBranchAddress("tau_isTightID", &tau_isTightID);
    tree->SetBranchAddress("tau_truthMatched", &tau_truthMatched);
    tree->SetBranchAddress("tau_trigMatched", &tau_trigMatched);
    tree->SetBranchAddress("tau_nTracks", &tau_nTracks);
    tree->SetBranchAddress("tau_BDTid", &tau_BDTid);
    
    // Di-tau
    tree->SetBranchAddress("ditau_m", &ditau_m);
    
    // Systematic uncertainties
    tree->SetBranchAddress("lep_pt_syst", &lep_pt_syst);
    tree->SetBranchAddress("met_et_syst", &met_et_syst);
    tree->SetBranchAddress("jet_pt_syst", &jet_pt_syst);
    tree->SetBranchAddress("photon_pt_syst", &photon_pt_syst);
    tree->SetBranchAddress("largeRjet_pt_syst", &largeRjet_pt_syst);
    tree->SetBranchAddress("tau_pt_syst", &tau_pt_syst);
    
    // Get total number of entries
    Long64_t totalEntries = tree->GetEntries();
    cout << "Total events in file: " << totalEntries << endl;
    
    // Determine actual sample size and create event indices
    Long64_t actualSampleSize = (sampleSize == -1 || sampleSize > totalEntries) ? totalEntries : sampleSize;
    cout << "Sampling " << actualSampleSize << " events..." << endl;
    
    // Create vector of event indices
    vector<Long64_t> eventIndices;
    for (Long64_t i = 0; i < totalEntries; i++) {
        eventIndices.push_back(i);
    }
    
    // Random sampling if needed
    if (actualSampleSize < totalEntries) {
        random_device rd;
        mt19937 gen(rd());
        shuffle(eventIndices.begin(), eventIndices.end(), gen);
        eventIndices.resize(actualSampleSize);
        sort(eventIndices.begin(), eventIndices.end()); // Sort for efficient tree access
    }
    
    // Create output CSV filename
    string csvFilename = filename.substr(0, filename.find_last_of("."));
    if (actualSampleSize < totalEntries) {
        csvFilename += "_sample" + to_string(actualSampleSize);
    }
    csvFilename += ".csv";
    cout << "Output CSV file: " << csvFilename << endl;
    
    // Open CSV file for writing
    ofstream csvFile(csvFilename);
    if (!csvFile.is_open()) {
        cerr << "Error: Cannot create CSV file " << csvFilename << endl;
        file->Close();
        return 1;
    }
    
    // Helper function to write vector as list
    auto writeVectorList = [&csvFile](const auto* vec, unsigned int size) {
        csvFile << ",\"[";
        for (unsigned int i = 0; i < size; i++) {
            csvFile << (*vec)[i];
            if (i < size - 1) csvFile << ",";
        }
        csvFile << "]\"";
    };
    
    // Write CSV header
    csvFile << "event,runNumber,eventNumber,channelNumber,mcWeight,XSection,SumWeights,"
            << "scaleFactor_PILEUP,scaleFactor_ELE,scaleFactor_MUON,scaleFactor_PHOTON,scaleFactor_TAU,scaleFactor_BTAG,scaleFactor_LepTRIGGER,scaleFactor_PhotonTRIGGER,"
            << "trigE,trigM,trigP,"
            << "lep_n,lep_truthMatched,lep_trigMatched,lep_pt,lep_eta,lep_phi,lep_E,lep_z0,lep_charge,lep_type,lep_isTightID,lep_ptcone30,lep_etcone20,lep_trackd0pvunbiased,lep_tracksigd0pvunbiased,"
            << "met_et,met_phi,"
            << "jet_n,jet_pt,jet_eta,jet_phi,jet_E,jet_jvt,jet_trueflav,jet_truthMatched,jet_MV2c10,"
            << "photon_n,photon_truthMatched,photon_trigMatched,photon_pt,photon_eta,photon_phi,photon_E,photon_isTightID,photon_ptcone30,photon_etcone20,photon_convType,"
            << "largeRjet_n,largeRjet_pt,largeRjet_eta,largeRjet_phi,largeRjet_E,largeRjet_m,largeRjet_truthMatched,largeRjet_D2,largeRjet_tau32,"
            << "tau_n,tau_pt,tau_eta,tau_phi,tau_E,tau_charge,tau_isTightID,tau_truthMatched,tau_trigMatched,tau_nTracks,tau_BDTid,"
            << "ditau_m,lep_pt_syst,met_et_syst,jet_pt_syst,photon_pt_syst,largeRjet_pt_syst,tau_pt_syst" << endl;
    
    // Process selected events
    for (Long64_t idx = 0; idx < actualSampleSize; idx++) {
        Long64_t eventNum = eventIndices[idx];
        tree->GetEntry(eventNum);
        
        if (idx % 10000 == 0) {
            cout << "Processing event " << idx + 1 << "/" << actualSampleSize << " (original event " << eventNum << ")" << endl;
        }
        
        // Write basic event info
        csvFile << eventNum << "," << runNumber << "," << eventNumber << "," << channelNumber << "," 
                << mcWeight << "," << XSection << "," << SumWeights << ","
                << scaleFactor_PILEUP << "," << scaleFactor_ELE << "," << scaleFactor_MUON << "," 
                << scaleFactor_PHOTON << "," << scaleFactor_TAU << "," << scaleFactor_BTAG << "," 
                << scaleFactor_LepTRIGGER << "," << scaleFactor_PhotonTRIGGER << ","
                << trigE << "," << trigM << "," << trigP << ","
                << lep_n;
        
        // Write lepton vectors
        writeVectorList(lep_truthMatched, lep_n);
        writeVectorList(lep_trigMatched, lep_n);
        writeVectorList(lep_pt, lep_n);
        writeVectorList(lep_eta, lep_n);
        writeVectorList(lep_phi, lep_n);
        writeVectorList(lep_E, lep_n);
        writeVectorList(lep_z0, lep_n);
        writeVectorList(lep_charge, lep_n);
        writeVectorList(lep_type, lep_n);
        writeVectorList(lep_isTightID, lep_n);
        writeVectorList(lep_ptcone30, lep_n);
        writeVectorList(lep_etcone20, lep_n);
        writeVectorList(lep_trackd0pvunbiased, lep_n);
        writeVectorList(lep_tracksigd0pvunbiased, lep_n);
        
        // Write missing ET
        csvFile << "," << met_et << "," << met_phi << "," << jet_n;
        
        // Write jet vectors
        writeVectorList(jet_pt, jet_n);
        writeVectorList(jet_eta, jet_n);
        writeVectorList(jet_phi, jet_n);
        writeVectorList(jet_E, jet_n);
        writeVectorList(jet_jvt, jet_n);
        writeVectorList(jet_trueflav, jet_n);
        writeVectorList(jet_truthMatched, jet_n);
        writeVectorList(jet_MV2c10, jet_n);
        
        // Write photon info
        csvFile << "," << photon_n;
        writeVectorList(photon_truthMatched, photon_n);
        writeVectorList(photon_trigMatched, photon_n);
        writeVectorList(photon_pt, photon_n);
        writeVectorList(photon_eta, photon_n);
        writeVectorList(photon_phi, photon_n);
        writeVectorList(photon_E, photon_n);
        writeVectorList(photon_isTightID, photon_n);
        writeVectorList(photon_ptcone30, photon_n);
        writeVectorList(photon_etcone20, photon_n);
        writeVectorList(photon_convType, photon_n);
        
        // Write large R jet info
        csvFile << "," << largeRjet_n;
        writeVectorList(largeRjet_pt, largeRjet_n);
        writeVectorList(largeRjet_eta, largeRjet_n);
        writeVectorList(largeRjet_phi, largeRjet_n);
        writeVectorList(largeRjet_E, largeRjet_n);
        writeVectorList(largeRjet_m, largeRjet_n);
        writeVectorList(largeRjet_truthMatched, largeRjet_n);
        writeVectorList(largeRjet_D2, largeRjet_n);
        writeVectorList(largeRjet_tau32, largeRjet_n);
        
        // Write tau info
        csvFile << "," << tau_n;
        writeVectorList(tau_pt, tau_n);
        writeVectorList(tau_eta, tau_n);
        writeVectorList(tau_phi, tau_n);
        writeVectorList(tau_E, tau_n);
        writeVectorList(tau_charge, tau_n);
        writeVectorList(tau_isTightID, tau_n);
        writeVectorList(tau_truthMatched, tau_n);
        writeVectorList(tau_trigMatched, tau_n);
        writeVectorList(tau_nTracks, tau_n);
        writeVectorList(tau_BDTid, tau_n);
        
        // Write di-tau and systematic uncertainties
        csvFile << "," << ditau_m;
        writeVectorList(lep_pt_syst, lep_n);
        csvFile << "," << met_et_syst;
        writeVectorList(jet_pt_syst, jet_n);
        writeVectorList(photon_pt_syst, photon_n);
        writeVectorList(largeRjet_pt_syst, largeRjet_n);
        writeVectorList(tau_pt_syst, tau_n);
        
        csvFile << endl;
    }
    
    csvFile.close();
    file->Close();
    
    cout << "Successfully exported " << actualSampleSize << " events to " << csvFilename << endl;
    return 0;
}