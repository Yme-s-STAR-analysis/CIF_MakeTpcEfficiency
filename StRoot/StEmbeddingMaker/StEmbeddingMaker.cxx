#include "StEmbeddingMaker.h"

#include <TMath.h>

#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StPicoEvent/StPicoMcTrack.h"
#include "StPicoEvent/StPicoMcVertex.h"
#include "StThreeVectorF.hh"
#include "StLorentzVector.hh"
#include "Stiostream.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include "phys_constants.h"

#include "StRoot/CentCorrTool/CentCorrTool.h"
#include "StRoot/MeanDcaTool/MeanDcaTool.h"
#include "StRoot/TpcShiftTool/TpcShiftTool.h"
#include "StRoot/TriggerTool/TriggerTool.h"
#include "StRoot/StCFMult/StCFMult.h"
#include "StRoot/VtxShiftTool/VtxShiftTool.h"

StEmbeddingMaker::StEmbeddingMaker(
	const char* name, 
	StPicoDstMaker* picoMaker,
    const char* outName,
	const std::string& energy,
	int targetID
) : StMaker(name) {
	mOutputName = outName;
	mPicoDstMaker = picoMaker;
	mPicoDst = 0;
	this->energy = energy;
	this->targetID = targetID;
	nvz = 5;
	nCent = 24;
}

StEmbeddingMaker::~StEmbeddingMaker() {}

Int_t StEmbeddingMaker::Init() {
  	mFileOut = new TFile(mOutputName, "recreate");

	// patch 3.0: 
	// - 1) Add a flag 'match', 1 for RC matched, 0 for NOT matched
	// -> thus, we don't need MC information:
	// -> When a track is matched, we will use RC pt, y, eta...
	// -> When a track is NOT matched, we will use MC pt, y, eta...
	// -> And if it is NOT matched, the pt, y and eta of RC tracks would be ones from MC (indeed, there is no RC track for this case)
	// - 2) Some quantities will not be used, like nHitsPoss, dedx, dcaxy/z
	//	Previous quantities:
    //	"cent:vz:"
    //  "pTMc:etaMc:yMc:"
    //  "pTRc:etaRc:yRc:"
    // 	"nHitsFit:nHitsPoss:nHitsRatio:nHitsDedx:dedx:dca:dcaXY:dcaZ";

	// initialize costume modules

	mtDca = new BES2Processing::MeanDcaTool();
	mtCent = new BES2Processing::CentCorrTool();
	mtMult = new BES2Processing::StCFMult(energy == "27");
	mtShift = new BES2Processing::TpcShiftTool(energy);
	mtTrg = new BES2Processing::TriggerTool(energy);
	mtVtx = new BES2Processing::VtxShiftTool(energy);

	// mean dca tool setup
	mtDca->ReadParams(energy);

	// centrality tool
	mtCent->InitParams(energy);

	// Multiplicity and shift tool
	mtShift->Init();
	mtMult->ImportShiftTool(mtShift);

	// tag, dca, nHitsFit, nSigmaProton, mass2Low, mass2High
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
		vzCut = 27;
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
		vzCut = 50;
		nvz = 5;
	}

	// set vz cut for extended rapidity (y0p6)
	std::map<std::string, double> vzMapExt = {
		{"7", 20},
		{"9", 30},
		{"11", 30},
		{"14", 40},
		{"17", 40},
		{"19", 40},
		{"27", 27}
	};
	auto it = vzMapExt.find(energy);
	if (it != vzMapExt.end()) {
		vzCutExt = it->second;
	} else {
		vzCutExt = vzCut;
		std::cout << "[WARNING] - DstMaker: Energy tag [" << energy << "] not found!" << std::endl;
	}

	// init TEfficiency here
	std::string partName;
	if (targetID == 14) {
		partName = "EffPro";
	} else if (targetID == 15) {
		partName = "EffPbar";
	} else {
		std::cout << "[ERROR] - DstMaker: Invalid particle type (" << targetID << ")" << std::endl;
		return kStFATAL;
	}
	for (const auto& cut : sysCuts) {
		for (int icent=0; icent<nCent; icent++) {
			for (int ivz=0; ivz<nvz; ivz++) {
				std::string effName = partName + 
				"_cent" + std::to_string(icent) + 
				"_vz" + std::to_string(ivz) + 
				"_" + cut.tag;
				efficiency[cut.tag][icent][ivz] = new TEfficiency(
					effName.c_str(),
					";y;p_{T} (GeV/c);Efficiency",
					30, -1.5, 1.5, 
					21, 0.0, 2.1
				);
			}
		}
	}

	return kStOK;
}

//---------------------------------------------------------
Int_t StEmbeddingMaker::Finish() {
	std::cout << "[LOG] Number of events: " << nEvents << ".\n";
	mFileOut->cd();
	for (const auto& cut : sysCuts) {
		for (int icent=0; icent<nCent; icent++) {
			for (int ivz=0; ivz<nvz; ivz++) {
				efficiency[cut.tag][icent][ivz]->Write();
			}
		}
	}
	mFileOut->Close();
	std::cout << "[LOG] This is the end of this job.\n";
	return kStOK;
}

void StEmbeddingMaker::Clear(Option_t* opt) {}

//---------------------------------------------------------------
Int_t StEmbeddingMaker::Make() {
	if (!mPicoDstMaker) {
		LOG_WARN << " No PicoDstMaker! Skip! " << endm;
		return kStWarn;
	}

	mPicoDst = mPicoDstMaker->picoDst();
	if (!mPicoDst) {
		LOG_WARN << " No PicoDst! Skip! " << endm;
		return kStWarn;
	}

	if (!mPicoDst) {
		return kStOK;
	}

	// Load event
	event = (StPicoEvent*)mPicoDst->event();
	if (!event) {
		cerr << "Error opening picoDst Event, skip!" << endl;
		return kStOK;
	}
	nEvents += 1;

	// Vertex Cut
	TVector3 pVtx = event->primaryVertex();
	Float_t vx = event->primaryVertex().X();
	Float_t vy = event->primaryVertex().Y();
	Float_t vz = event->primaryVertex().Z();

	if (fabs(vx) < 1.e-5 && 
		fabs(vy) < 1.e-5 &&
		fabs(vz) < 1.e-5) {
		return kStOK;
	}

	// wide vertex cut -> tight cut now
	if (fabs(vz) > vzCut) { return kStOK; }
	auto vr = mtVtx->GetShiftedVr(vx, vy);
	auto vzBin = vz_split(vz);
	if (vr > 1) { return kStOK; }

	// check trigger ID
	Int_t trgid = mtTrg->GetTriggerID(event);
	if (trgid < 0) { return kStOK; }

	// centrality
	mtMult->make(mPicoDst);
	Int_t refMult = mtMult->mRefMult;
	Int_t tofMult = mtMult->mTofMult;
	Int_t nTofMatch = mtMult->mNTofMatch;
	Int_t nTofBeta = mtMult->mNTofBeta;
	Int_t refMult4 = mtMult->mRefMult4;

	Int_t refMult3X = mtMult->mRefMult3X;// this will be RefMult3 for 27 GeV (implemented in StCFMult)
	refMult3X = mtCent->GetCorrectedRefMult3(
		refMult3X, refMult, tofMult, nTofMatch, nTofBeta,
		vz, refMult4
	);

	if (refMult3X < 0) { return kStOK; } // pile up events: -1

	int centX = mtCent->GetCentralityClass24(refMult3X);
	if (centX < 0) { return kStOk; }

	// check DCA
	if (!mtDca->Make(mPicoDst)) { return kStOK; }
	if (mtDca->IsBadMeanDcaZEvent(mPicoDst) || mtDca->IsBadMeanDcaXYEvent(mPicoDst)) {
		return kStOK;
	}

	Int_t numberOfMcTracks = mPicoDst->numberOfMcTracks();
	Int_t numberOfRcTracks = mPicoDst->numberOfTracks();
	Int_t numberOfMcVertices = mPicoDst->numberOfMcVertices();
	if (!numberOfMcVertices || !numberOfMcTracks) { 
		// this event has no MC information, skip it
		return kStOK;
	}
	
	const Float_t mField = event->bField();
	// Reconstructed track loop
	// to construct the map of RcTrack ID -> McTrack ID
	std::map<Int_t, Int_t> mMc2Rc;
	for (Int_t iRcTrk=0; iRcTrk<numberOfRcTracks; iRcTrk++) {
		rcTrack = (StPicoTrack*)mPicoDst->track(iRcTrk);
		if (!rcTrack) {
			continue;
		}
		if (!rcTrack->isPrimary()) {
			continue;
		}

		Int_t idTruth = rcTrack->idTruth(); // index of corresponding MC track
		if (idTruth <= 0 || idTruth > 10000) {
			continue;
		}

		// note that, here idTruth - 1 will be the quantity iMcTrk in later codes
		// not just idTruth. VERY IMPORTANT!
		mMc2Rc.insert(std::pair<Int_t, Int_t>(idTruth - 1, iRcTrk)); 
	}

	// MC track loop
	// to record the MC track and its RC track
	bool flag = false;
	for (Int_t iMcTrk=0; iMcTrk<numberOfMcTracks; iMcTrk++){
		mcTrack = (StPicoMcTrack*)mPicoDst->mcTrack(iMcTrk);
		if (!is_McTrack_from_PV()) {
			continue;
		}
		if (!is_target_particle()) {
			// only fill Geant ID histogram itself when the track is not proton
			continue;
		}
		if (mcTrack->idVtxStart() != 1) {
			// make sure this track is from PV
			continue;
		}

		// don't need this, the flag will be set later
		// if (mMc2Rc.count(iMcTrk) == 0) { // this mc track does not have rc track
		// 	flag = false;
		// }

		// find the reconstruct track
		// Indeed, mcTrack->id() is iMcTrk + 1;
		// cannot just use iRcTrk = mMc2Rc.find(iMcTrk)->second to get the RC track ID
		pair<map<Int_t, Int_t>::iterator, map<Int_t, Int_t>::iterator> ret;
		ret = mMc2Rc.equal_range(iMcTrk);
		map<Int_t, Int_t>::iterator iter;
		Int_t count = 0;
		Int_t iRcTrk = -1;
		Int_t iRcTrk_best = -1;
		Int_t qaTruth_best = -1;
		for (iter=ret.first; iter!=ret.second; iter++, count++) { // loop over the possible reconstructed tracks and find the best one with greatest qaTruth
			iRcTrk = iter->second;
			rcTrack = (StPicoTrack*)mPicoDst->track(iRcTrk);
			if (!rcTrack) {
				continue;
			}
			if (!rcTrack->isPrimary()){
				continue;
			}
			if (rcTrack->qaTruth() > qaTruth_best){
				qaTruth_best = rcTrack->qaTruth();
				iRcTrk_best = iRcTrk;
			}
		} 
		if (count > 0) { // indicates that rc track was found, and set the best one
			rcTrack = (StPicoTrack*)mPicoDst->track(iRcTrk_best);
			if (!rcTrack || !rcTrack->isPrimary()){
				flag = false;
			} else {
				flag = true;
			}
		} else {
			rcTrack = 0;
			flag = false;
		}

		TVector3 rcMom;
		StThreeVector<Float_t> rcMom3;
		StLorentzVector<Float_t> rcMom4;
		Double_t dca = -999.0;
		if (flag) {
			rcMom = rcTrack->pMom();
			rcMom3 = StThreeVector<Float_t>(rcMom.X(), rcMom.Y(), rcMom.Z());
			// calculate rapidity
			Float_t MP = 0.938272;
			Float_t EP = sqrt(rcMom3.mag2() + MP*MP);
			rcMom4 = StLorentzVector<Float_t>(rcMom3, EP);
			StPicoPhysicalHelix helix = rcTrack->helix(mField);
			dca = fabs(helix.geometricSignedDistance(pVtx));
		}

		// fill the efficiency map
		auto dpt = flag ? (Float_t)rcMom.Pt() : (Float_t)mcTrack->pt();
		auto dy = flag ? (Float_t)rcMom4.rapidity() : (Float_t)mcTrack->rapidity();
		auto ddca = flag ? dca : 999;
		auto dnhitsfit = flag ? rcTrack->nHitsFit() * 1.0 : -999;
		auto dnhitsratio = flag ? (rcTrack->nHitsFit() * 1.0) / (rcTrack->nHitsPoss() * 1.0) : -999;
		auto dnhitsdedx = flag ? rcTrack->nHitsDedx() * 1.0 : -999;
		for (const auto& cut : sysCuts) {
			bool acc = true;
			if (ddca > cut.dca) { acc = false; }
			if (dnhitsratio < 0.52) { acc = false; }
			if (dnhitsdedx <= 5) { acc = false; }
			if (dnhitsfit <= cut.nHitsFit) { acc =false; }
			efficiency[cut.tag][centX][vzBin]->Fill(acc, dy, dpt);
		}
	}

	return kStOK;
}

bool StEmbeddingMaker::is_target_particle() {
	return (mcTrack->geantId() == targetID);
}

bool StEmbeddingMaker::is_McTrack_from_PV() {
	if (!mcTrack) {
		return false;
	}
	Int_t idMcVx = mcTrack->idVtxStart();
	while (idMcVx != 1) {
		StPicoMcVertex* mcVertex = (StPicoMcVertex*)mPicoDst->mcVertex(idMcVx - 1);
		Int_t idMcTrack = mcVertex->idOfParentTrack();
		if (!idMcTrack) {
			break;
		}
		StPicoMcTrack* mcTrackP = (StPicoMcTrack*)mPicoDst->mcTrack(idMcTrack - 1);
		idMcVx = mcTrackP->idVtxStart();
		if (!idMcVx) {
			break;
		}
	}
	if (idMcVx != 1) {
		return false;
	}
	return true;
}

int StEmbeddingMaker::vz_split(double val) {
	if (energy == "27") {
		if (-27 < val && val < -10) {
			return 0;
		} else if (-10 < val && val < 10) {
			return 1;
		} else if (10 < val && val < 27) {
			return 2;
		} else {
			return -1;
		}
	} else {
		if (-30 < val && val < -10) {
			return 0;
		} else if (-10 < val && val < 10) {
			return 1;
		} else if (10 < val && val < 30) {
			return 2;
		} else if (-50 < val && val < -30) {
			return 3;
		} else if (30 < val && val < 50) {
			return 4;
		} else {
			return -1;
		}
	}
}