////////////////////////////////////////////////////////
//  
//This class reads UniGen-files
//see
//http://www.gsi.de/forschung/helmholtz-vi/unigen/
//for further information
//
//Parts of the source code have been taken from 
//the UniGen project.
//
//
//                    Author:  Ingo Froehlich
//                    Written: 01/02/2011
//
////////////////////////////////////////////////////////

#include "PUniGenInput.h"
#include "PChannel.h"
#include "PDistributionManager.h"

#include "UParticle.h"

PUniGenInput::PUniGenInput() {

    fPriority = FILEINPUT_PRIORITY;
    unknown_pdg_pointer = 0;
    reset_mass = 0;
}

PUniGenInput::PUniGenInput(char *filename) {

    fPriority = FILEINPUT_PRIORITY;
    unknown_pdg_pointer = 0;
    reset_mass = 0;

    Input(filename);
}

Bool_t PUniGenInput::Input(char *filename) {
    //Copied from UManager
    fInFile = new TFile(filename);
    if(NULL == fInFile) {
	Fatal("Input", "Input file <%s> does not exist",filename);
    }

    fRun = (URun*) fInFile->Get("run");
    
    fInTree = (TTree*) fInFile->Get("events");
    if(NULL == fInTree) {
	Fatal("SetInputFile", "Input file <%s> has no events tree", filename);
    }

    Info("Input", "Input file has %lli events", fInTree->GetEntries());
    
    nentries = fInTree->GetEntries();
    centry = 0;

    if(NULL == fEvent) {
	fEvent = new UEvent();
    }

    fInTree->SetBranchAddress("event", &fEvent);

    return kTRUE;

};

Bool_t PUniGenInput::Modify(PParticle **mstack, int *, int *num, int stacksize) {

    if (centry == 0) { 
	//init before 1st event
	makeDistributionManager()->Exec("pdg:init");
	pdg_param    = makeDataBase()->GetParamInt("pdg");
	if (pdg_param < 0) 
	    return kFALSE;
	pid_param    = makeDataBase()->GetParamInt("pid");
    }

    if (centry == nentries) {
	Info("Modify", "Unigen file: number of events reached");
	return kFALSE;
    }
    
    int cur = *num;
    PParticle dummy("dummy");

    fInTree->GetEntry(centry);
    centry++;

    Int_t npart     = fEvent->GetNpa();
#if 0
    //not yet used....
    Int_t evnr      = fEvent->GetEventNr();
    Int_t stepnr    = fEvent->GetStepNr();
#endif

    for (int i=0; i<npart; i++) {

	if (cur == stacksize) {
	    Warning("Modify", "Stack size too small, increase '_system_particle_stacksize'");
	    return kTRUE;
	}

	UParticle *part = fEvent->GetParticle(i);


	if (!makeDataBase()->GetParamInt (pdg_param, part->GetPdg(), pid_param, &i_result)) {
	    bool found = kFALSE;
	    for (int i=0; i<unknown_pdg_pointer; i++) {
		if (unknown_pdg[i] == part->GetPdg()) 
		    found=kTRUE;
	    }
	    if (!found) {
		Warning("Modify",
			"Particle with PDG code %i is not defined in Pluto, further particles will get the 'dummy' pid",
			part->GetPdg());
		if (unknown_pdg_pointer != MAX_UNIGENUNKNOWNPDG) {
		    unknown_pdg[unknown_pdg_pointer] = part->GetPdg();
		    unknown_pdg_pointer++;
		}	    
	    }
	    
	} else {
	    dummy.SetID(*i_result);
	}
	
	dummy.SetVect4(part->GetMomentum());
	Double_t x = part->GetPosition().X();
	Double_t y = part->GetPosition().Y();
	Double_t z = part->GetPosition().Z();
	Double_t t = part->GetPosition().T();
	dummy.SetVertex(x,y,z,t);
	dummy.SetIndex(part->GetIndex());
	dummy.SetParentIndex(part->GetParent());
	dummy.SetSiblingIndex(part->GetMate());
	dummy.SetSourceId(part->GetParentDecay());
	dummy.SetW(part->GetWeight());
	
	if (dummy.M() < PData::LMass(dummy.ID()) || dummy.M() > PData::UMass(dummy.ID())) {
	    if (reset_mass == 1) {
		dummy.SetM(makeStaticData()->GetParticleMass(dummy.ID()));
	    } else if (reset_mass == 0) {
		Warning("Modify", "Found a particle mass which is out of bounds, it is recommended to use ResetInvalidMass()");
		Warning("Modify", "(following warnings will be skipped)");
		reset_mass=-1;
	    } 
	}
	
	*(mstack[cur]) = dummy;
			
	cur++;	
    }

    *num = cur;

    return kTRUE;
}

ClassImp(PUniGenInput)
