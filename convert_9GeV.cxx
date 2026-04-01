#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "TH2F.h"
#include "TFile.h"
#include "TString.h"
#include "TEfficiency.h"

using std::map;
using std::vector;
using std::string;

struct SystematicCuts {
    string tag;
    double dca;
    int nHitsFit;
    double nSigmaProton;
    double mass2Low;
    double mass2High;
};

void convert() {
    // change those configuration in need
    const char* fNameIn1 = "9aGeV.hadd.root";
    const char* fNameIn2 = "9bGeV.hadd.root";
    const char* fNameOut = "CIF.TPC.EFF.9.root";
    string energy = "9";

    // don't need to do anything else below
    std::cout << "[LOG] Input file Region 1: " << fNameIn1 << std::endl;
    std::cout << "[LOG] Input file Region 2: " << fNameIn2 << std::endl;
    std::cout << "[LOG] Outpurfile: " << fNameOut << std::endl;
    std::cout << "[LOG] Energy tag " << energy << std::endl;

    
    auto fIn1 = TFile::Open(fNameIn1, "read");
    auto fIn2 = TFile::Open(fNameIn2, "read");

    map<string, 
    map<int,
    map<int, TEfficiency*>>> effPro1;
    map<string, 
    map<int,
    map<int, TEfficiency*>>> effPbar1;
    map<string, 
    map<int,
    map<int, TH2F*>>> hPro1;
    map<string, 
    map<int,
    map<int, TH2F*>>> hPbar1;
    map<string, 
    map<int,
    map<int, TEfficiency*>>> effPro2;
    map<string, 
    map<int,
    map<int, TEfficiency*>>> effPbar2;
    map<string, 
    map<int,
    map<int, TH2F*>>> hPro2;
    map<string, 
    map<int,
    map<int, TH2F*>>> hPbar2;

    int nvz;
    int nCent = 24;
    vector<SystematicCuts> sysCuts;

    if (energy == "27") {
		sysCuts = {
			{"default", 1.0, 15, 2, 0.6, 1.2},
			{"dca0p8", 0.8, 15, 2, 0.6, 1.2},
			{"dca0p9", 0.9, 15, 2, 0.6, 1.2},
			{"dca1p1", 1.1, 15, 2, 0.6, 1.2},
			{"dca1p2", 1.2, 15, 2, 0.6, 1.2},
			{"nhit12", 1.0, 12, 2, 0.6, 1.2},
			{"nhit17", 1.0, 17, 2, 0.6, 1.2}
		};
		nvz = 3;
	} else {
		sysCuts = {
			{"default", 1.0, 20, 2, 0.6, 1.2},
			{"dca0p8", 0.8, 20, 2, 0.6, 1.2},
			{"dca0p9", 0.9, 20, 2, 0.6, 1.2},
			{"dca1p1", 1.1, 20, 2, 0.6, 1.2},
			{"dca1p2", 1.2, 20, 2, 0.6, 1.2},
			{"nhit15", 1.0, 15, 2, 0.6, 1.2},
			{"nhit18", 1.0, 18, 2, 0.6, 1.2},
			{"nhit22", 1.0, 22, 2, 0.6, 1.2},
			{"nhit25", 1.0, 25, 2, 0.6, 1.2}
		};
		nvz = 5;
	}

    std::cout << "[LOG] Number of centrality bins: " << nCent << std::endl;
    std::cout << "[LOG] Number of Vz bins: " << nvz << std::endl;
    std::cout << "[LOG] Systematic tags:\n";
    for (const auto& cut: sysCuts) {
        std::cout << "\t" << cut.tag << "\n";
    }
    std::cout << "[LOG] Now converting..." << std::endl;
	for (const auto& cut : sysCuts) {
		for (int icent=0; icent<nCent; icent++) {
			for (int ivz=0; ivz<nvz; ivz++) {
                effPro1[cut.tag][icent][ivz] = (TEfficiency*)fIn1->Get(Form("EffPro_cent%d_vz%d_%s", icent, ivz, cut.tag.c_str()));
                effPbar1[cut.tag][icent][ivz] = (TEfficiency*)fIn1->Get(Form("EffPbar_cent%d_vz%d_%s", icent, ivz, cut.tag.c_str()));
				hPro1[cut.tag][icent][ivz] = (TH2F*)(effPro1[cut.tag][icent][ivz]->CreateHistogram());
                hPro1[cut.tag][icent][ivz]->SetName(Form("TpcEff_cent%d_vz%d_Pro_%s_0", icent, ivz, cut.tag.c_str()));
				hPbar1[cut.tag][icent][ivz] = (TH2F*)(effPbar1[cut.tag][icent][ivz]->CreateHistogram());
                hPbar1[cut.tag][icent][ivz]->SetName(Form("TpcEff_cent%d_vz%d_Pbar_%s_0", icent, ivz, cut.tag.c_str()));

                effPro2[cut.tag][icent][ivz] = (TEfficiency*)fIn2->Get(Form("EffPro_cent%d_vz%d_%s", icent, ivz, cut.tag.c_str()));
                effPbar2[cut.tag][icent][ivz] = (TEfficiency*)fIn2->Get(Form("EffPbar_cent%d_vz%d_%s", icent, ivz, cut.tag.c_str()));
				hPro2[cut.tag][icent][ivz] = (TH2F*)(effPro2[cut.tag][icent][ivz]->CreateHistogram());
                hPro2[cut.tag][icent][ivz]->SetName(Form("TpcEff_cent%d_vz%d_Pro_%s_1", icent, ivz, cut.tag.c_str()));
				hPbar2[cut.tag][icent][ivz] = (TH2F*)(effPbar2[cut.tag][icent][ivz]->CreateHistogram());
                hPbar2[cut.tag][icent][ivz]->SetName(Form("TpcEff_cent%d_vz%d_Pbar_%s_1", icent, ivz, cut.tag.c_str()));
			}
		}
	}

    auto fOut = TFile::Open(fNameOut, "recreate");
    fOut->cd();
    std::cout << "[LOG] Now saving..." << std::endl;

	for (const auto& cut : sysCuts) {
		for (int icent=0; icent<nCent; icent++) {
			for (int ivz=0; ivz<nvz; ivz++) {
				hPro1[cut.tag][icent][ivz]->Write();
				hPbar1[cut.tag][icent][ivz]->Write();
				hPro2[cut.tag][icent][ivz]->Write();
				hPbar2[cut.tag][icent][ivz]->Write();
			}
		}
	}

    fOut->Close();
    fIn1->Close();
    fIn2->Close();
    std::cout << "[LOG] All done!" << std::endl;

}