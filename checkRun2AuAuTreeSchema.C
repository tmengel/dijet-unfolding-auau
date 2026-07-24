#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"

// Compare the run-2 Au+Au jet tree schema with makeUnfoldingHists.C.
// Run with:
// root -l -b -q 'checkRun2AuAuTreeSchema.C("/path/to/input.root")'
void checkRun2AuAuTreeSchema(
    const char* filename = "/home/tmengel/PPG14/rootfiles/run2auau_ana509_2024p022_v001_r03_jets.root")
{
  TFile input(filename, "READ");
  if (input.IsZombie()) {
    std::cerr << "Could not open " << filename << std::endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(input.Get("T"));
  if (!tree) tree = dynamic_cast<TTree*>(input.Get("ttree"));
  if (!tree) {
    std::cerr << "Neither T nor ttree is present in " << filename << std::endl;
    return;
  }

  struct Mapping {
    const char* expected;
    const char* available;
  };

  // cone_size = 3, the default in makeUnfoldingHists.C
  const std::vector<Mapping> mappings = {
      {"centrality", "cent"},
      {"mbd_vertex_z", "zvrtx"},
      {"jet_pt_3_sub", "jet_pT"},
      {"jet_pt_unsub_3_sub", "jet_unsub_pT"},
      {"jet_e_3_sub", "jet_E"},
      {"jet_e_unsub_3_sub", "jet_unsub_E"},
      {"jet_eta_3_sub", "jet_eta"},
      {"jet_phi_3_sub", "jet_phi"},
      {"gl1_scaled", nullptr},
      {"minbias", nullptr},
      {"mbd_time_zero", nullptr},
  };

  std::cout << "Input tree: " << tree->GetName()
            << " (makeUnfoldingHists.C expects: ttree)\n\n";

  for (const auto& mapping : mappings) {
    const char* source = mapping.available ? mapping.available : mapping.expected;
    const bool present = tree->GetBranch(source) != nullptr;
    std::cout << (present ? "MATCH   " : "MISSING ")
              << mapping.expected << " <- "
              << (mapping.available ? mapping.available : "no corresponding branch")
              << '\n';
  }

  std::cout << "\nUse T (not ttree) and the listed branch mappings in an adapter. "
            << "The three MISSING fields require a documented replacement or selection change; "
            << "they cannot be recovered by renaming branches." << std::endl;
}
