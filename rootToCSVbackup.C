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
    
    // Variables to store jet data
    unsigned int jet_n = 0;
    vector<float> *jet_pt = nullptr;
    vector<float> *jet_eta = nullptr;
    vector<float> *jet_E = nullptr;
    vector<int> *jet_trueflav = nullptr;
    
    // Set branch addresses
    int status = tree->SetBranchAddress("jet_n", &jet_n);
    if (status < 0) {
        cerr << "Error: Could not set branch address for jet_n" << endl;
        file->Close();
        return 1;
    }
    
    tree->SetBranchAddress("jet_pt", &jet_pt);
    tree->SetBranchAddress("jet_eta", &jet_eta);
    tree->SetBranchAddress("jet_E", &jet_E);
    tree->SetBranchAddress("jet_trueflav", &jet_trueflav);
    
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
    
    // Write CSV header
    csvFile << "event,jet_n,jet_pt_list,jet_eta_list,jet_E_list,jet_trueflav_list" << endl;
    
    // Process selected events
    for (Long64_t idx = 0; idx < actualSampleSize; idx++) {
        Long64_t eventNum = eventIndices[idx];
        tree->GetEntry(eventNum);
        
        if (idx % 10000 == 0) {
            cout << "Processing event " << idx + 1 << "/" << actualSampleSize << " (original event " << eventNum << ")" << endl;
        }
        
        // Write event number and jet count
        csvFile << eventNum << "," << jet_n;
        
        // Write jet_pt list
        csvFile << ",\"[";
        for (unsigned int j = 0; j < jet_n; j++) {
            csvFile << (*jet_pt)[j];
            if (j < jet_n - 1) csvFile << ",";
        }
        csvFile << "]\"";
        
        // Write jet_eta list
        csvFile << ",\"[";
        for (unsigned int j = 0; j < jet_n; j++) {
            csvFile << (*jet_eta)[j];
            if (j < jet_n - 1) csvFile << ",";
        }
        csvFile << "]\"";
        
        // Write jet_E list
        csvFile << ",\"[";
        for (unsigned int j = 0; j < jet_n; j++) {
            csvFile << (*jet_E)[j];
            if (j < jet_n - 1) csvFile << ",";
        }
        csvFile << "]\"";
        
        // Write jet_trueflav list
        csvFile << ",\"[";
        for (unsigned int j = 0; j < jet_n; j++) {
            csvFile << (*jet_trueflav)[j];
            if (j < jet_n - 1) csvFile << ",";
        }
        csvFile << "]\"" << endl;
    }
    
    csvFile.close();
    file->Close();
    
    cout << "Successfully exported " << actualSampleSize << " events to " << csvFilename << endl;
    return 0;
}