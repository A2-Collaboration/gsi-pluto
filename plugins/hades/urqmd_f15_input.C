#include "PHUrReader.h"
#include "PReaction.h"
#include "PPlutoBulkDecay.h"
#include "PChannel.h"
#include "PParticle.h"
#include "PData.h"


#include "TString.h"
#include "TClonesArray.h"


#include "hphysicsconstants.h"

#include <iostream>

#include <iostream>
using namespace std;

Int_t UserAnalysis(PParticle** pA,Int_t size){
    

    TClonesArray* addon = PHUrReader::getAddon();
    TClonesArray* event = PHUrReader::getEvent();
    if(!addon) return 1;

    Int_t sizeAdd = addon->GetEntries();
    cout<<"UserAnalysis : nparticle ="<<size<<" ninput ="<<sizeAdd<<endl;
    for(Int_t i = 0; i < size; i ++){
	PParticle* p = pA[i];
	if(0&&i<sizeAdd) {
	    PParticle* p1 = (PParticle*)event->At(i);
	    cout<<i<<" "<<p->ID()<<" "<<p1->ID()<<endl;
	}
	if(i >= sizeAdd) {
	    PParticle* p1=0;
	    Int_t parentIndex = p->GetParentIndex() ;
	    cout<<"daughter "<<i<<" id "<<HPhysicsConstants::pid(p->ID())<<flush;
	    while(parentIndex!=-1){
		p1=0;
		if(parentIndex < sizeAdd ){
		    p1 = (PParticle*) event->At(parentIndex);
		} else {
		    p1 = pA[parentIndex];
		}
		if(p1) {
		    cout<<" -> "<<parentIndex<<" id "<<HPhysicsConstants::pid(p1->ID())<<flush;
		    parentIndex = p1->GetParentIndex();
		} else cout<<"p1=0, should not happen"<<endl;
	    }
            cout<<endl;
	}
    }
    return 1;
}

void urqmd_f15_input(TString infile = "/hera/hades/dstsim/apr12/urqmd/bmax9/f15/Au_Au_1230MeV_1000evts_1.f15",
		     TString outfile ="out-dilep_pluto",
                     Int_t nEvents=10
		    )
{
    // Read in urqmd f15 file (collision history) and put the particles
    // on the PLUTO bulk stack. Either dileptons (D0,D+,w,rho0,phi,eta,eta' and pi0)
    // can be decayed via shining (translated for original UrQMD fortran routines,
    // authors: Katharina Schmidt, Sascha Vogel, Christian Sturm 2007-2009) into
    // e+/e- or by PLUTO bulk decay. Optional stable non dileptons can be
    // written out. PHUrReader adds a TClonesArray "Addon" to the PLUTO Tree
    // that contains original UrQMD infos about creation and absorbtion times , densities etc.
    //

    Int_t maxnum                      = 19999;
    Int_t writeNonLeptons             = 2;   // default : 0=no nonleptons, 1=stable non leptons, 2=all nonlepton
    Bool_t runUrmdShining             = kFALSE;   // default : kTRUE  = run urqmd dilep shining code
    Bool_t writeFreezeOutLeptonsOnly  = kFALSE;   // default : kFALSE = store intermediate dileptons to output
    Bool_t doPlutoDecay               = kFALSE;   // default : kFALSE = do not decay dileptons with PLUTO
    TString outputUrmdDilepton        = "test_dilep.dat";

    if(doPlutoDecay && runUrmdShining){
	cout<<"WRONG CONFIG: RUNNING UrQMD SHINING and PLUTO DECAY AT THE SAME TIME MAKES NO SENSE! PLUTO SELECTED! "<<endl;
        runUrmdShining = kFALSE;
    }

    makeStaticData()->AddDecay(-1,"eta' -> dilepton, gamma","eta'","dilepton,gamma",0.0009); // missing channel: BR from PDG upper limit
    *(makeStaticData()->GetBatchValue("_system_particle_stacksize")) = maxnum;   // we have a lot of particles from urqmd

    //------------------------------------------------------
    // setup Urqmd Reader
    PHUrReader *input = new PHUrReader();
    input->Input                      (infile);
    input->Output                     (outputUrmdDilepton);
    input->SetOutputNonDileptons      (writeNonLeptons);
    input->SetOutputLeptons           (runUrmdShining );
    input->SetOutputFreezeoutDileptons(writeFreezeOutLeptonsOnly);
    input->SetMaxNumParticles(maxnum);
    //------------------------------------------------------
    // setup PLUTO

    //######################################################
    // Database work for selected channels of bulk decay (needed for PLUTO decay only)
    //------------------------------------------------------
    // switch all channels off
    makeStaticData()->DisableAllChannelBR("D0");
    makeStaticData()->DisableAllChannelBR("D+");
    makeStaticData()->DisableAllChannelBR("w");
    makeStaticData()->DisableAllChannelBR("rho0");
    makeStaticData()->DisableAllChannelBR("phi");
    makeStaticData()->DisableAllChannelBR("eta");
    makeStaticData()->DisableAllChannelBR("eta'");
    makeStaticData()->DisableAllChannelBR("pi0");


    //------------------------------------------------------
    // switch all wanted channels on
    makeStaticData()->SetEnhanceChannelBR("D0"  , "dilepton, n"    ,1.);
    makeStaticData()->SetEnhanceChannelBR("D+"  , "dilepton, p"    ,1.);
    makeStaticData()->SetEnhanceChannelBR("w"   , "dilepton, pi0"  ,0.000595/(0.000595+0.000071));   // both dal +dir  : Branching ratio=0.000595
    makeStaticData()->SetEnhanceChannelBR("eta" , "dilepton, gamma",1.);
    makeStaticData()->SetEnhanceChannelBR("eta'", "dilepton, gamma",1.); // missing channel!!
    makeStaticData()->SetEnhanceChannelBR("pi0" , "dilepton, gamma",1.);

    makeStaticData()->SetEnhanceChannelBR("phi" , "e+, e-"         ,1.);
    makeStaticData()->SetEnhanceChannelBR("w"   , "e+, e-"         ,0.000071/(0.000595+0.000071));   // both dal +dir   : Branching ratio=0.000071
    makeStaticData()->SetEnhanceChannelBR("rho0", "e+, e-"         ,1.);
    //######################################################

    PPlutoBulkDecay *pl = new PPlutoBulkDecay(); // needed to decay Dileptons with PLUTO
    pl->SetRecursiveMode(1);                     // Let also the products decay
    pl->SetTauMax(0.001);                        // maxTau in ns


    PReaction my_reaction((Char_t*)outfile.Data());
    //my_reaction.SetUserAnalysis(UserAnalysis);
    my_reaction.AddBulk(input);
    if(doPlutoDecay){
	my_reaction.AddBulk(pl);
    }

    my_reaction.Print();

    cout << my_reaction.Loop(nEvents) << " events converted" << endl;


}
