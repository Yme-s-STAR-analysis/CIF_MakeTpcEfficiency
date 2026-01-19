#ifndef _StEmbeddingMaker_head
#define _StEmbeddingMaker_head
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "TString.h"
#include "TVector3.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVectorD.hh"
#include "TNtuple.h"

#include "StRoot/CentCorrTool/CentCorrTool.h"
#include "StRoot/MeanDcaTool/MeanDcaTool.h"
#include "StRoot/TpcShiftTool/TpcShiftTool.h"
#include "StRoot/StCFMult/StCFMult.h"
#include "StRoot/TriggerTool/TriggerTool.h"
#include "StRoot/VtxShiftTool/VtxShiftTool.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoMcTrack;
class StPicoDstMaker;
class TH1F;
class TH2F;
class TEfficiency;
class TProfile;
class TTree;
class TH2D;
class TNtuple;

struct SystematicCuts {
	std::string tag;
	double dca;
	int nHitsFit;
	double nSigmaProton;
	double mass2Low;
	double mass2High;
};

class StEmbeddingMaker : public StMaker {
	public:
		StEmbeddingMaker(const char *name, 
			StPicoDstMaker *picoMaker, 
			const char *outName, 
			const std::string& energy,
			int targetID
		);
		virtual ~StEmbeddingMaker();

		virtual Int_t Init();
		virtual Int_t Make();
		virtual void  Clear(Option_t *opt="");
		virtual Int_t Finish();

		bool is_target_particle();
		bool is_McTrack_from_PV();

		int vz_split(double val);

	private:
		StPicoDstMaker *mPicoDstMaker;
		StPicoDst      *mPicoDst;
		StPicoEvent	   *event;
		StPicoTrack    *rcTrack;
		StPicoMcTrack  *mcTrack;

		Int_t targetID;

		BES2Processing::StCFMult* mtMult;
		BES2Processing::MeanDcaTool* mtDca;
		BES2Processing::CentCorrTool* mtCent;
		BES2Processing::TpcShiftTool* mtShift;
		BES2Processing::TriggerTool* mtTrg;
		BES2Processing::VtxShiftTool* mtVtx;

		std::map<std::string, 
			std::map<int, // centrality bin
			std::map<int, // vz bin
			TEfficiency*>>
		> efficiency;

		std::vector<SystematicCuts> sysCuts;
		double vzCut;
		double vzCutExt;
		int nCent;
		int nvz;

		TString mOutputName;
		TFile* mFileOut;

		std::string energy;
		Int_t nEvents;

		ClassDef(StEmbeddingMaker, 1)
};


ClassImp(StEmbeddingMaker)

#endif
