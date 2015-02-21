/////////////////////////////////////////////////////////////////////
//
// Generic distribution to change any kind of distribution 
// by using the batch script language to define equations
//
// This is a rejection method. It works as an afterburner
// for the normal distribution(s).
//
// Usage:
// 1.) Define a template, e.g.:
// PAnyDistribution* decay = 
//        new PAnyDistribution("name","description");
// decay->Add("....,     parent");
// decay->Add("....,     daughter");
//             .....
// 2.) If a non-uniform distribution should be defined (otherwise go to step 4),
// one has to add a "cache" which stores the undistorted
// distribution. The result is calculated by comparing
// the undistorted cache (as filled in 3) with the equation in step 4
// The size and binning of the cache must
// be chosen such, that during runtime (or even better: during Preheating) the statistic
// is sufficiently filled
// TH1F * cache  = new TH1F ("cache","......",...);
// (2- or 3-dim histograms are possible as well - if one has enough time!) 
//
// 3.) Add the equation to fill the cache:
// decay->AddEquation(cache,"_x = ...(calculate something here)...;");
// Here it is important to know that the daughters are in 
// their rest frame (i.e. the parent).
// But the parent itself is in lab frame.
// The parent is indicated by "_parent"
// Moreover, all particles of the decay (daughters and parent)
// can be accessed via the usual '[parname]' identifier.
// N.B.: "x,y,z,t" are reserved in TFormula, do NOT use such variable names.
//
// 4.) Add the usual equation:
// The distribution (the probability function)
// must be stored in "_f"
// decay->AddEquation("_f = ...(calculate something here)...");
//
// 5.) Remember, AnyDistribution is a rejection method. Therefore
// it can happen, that parts of the phase space, where _f has a 
// large probability, is not well populated by the generated events.
// In this case, the event loop will run forever, as Pluto tries
// to match the shape defined by _f.
// The following factor is the maximum enhancement factor to avoid such
// deadlocks.
// N.B.: It directly scales with the computing time!!!
// decay->SetMaxEnhancementFactor(10);
//
// 6.) Add this model to the Pluto data base:
// makeDistributionManager()->Add(decay);
// 
// An example can be found in macros/demo_PAnyDistribution.C
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PAnyDistribution.h"
#include "PDynamicData.h"

ClassImp(PAnyDistribution)

PAnyDistribution::PAnyDistribution() {
    Fatal("PAnyDistribution","Wrong ctor called");
} ;

PAnyDistribution::PAnyDistribution(const Char_t *id, const Char_t *de) :
    PDistribution(id, de) {

    parent    = NULL;
    projector = NULL;
    cache1    = NULL;
    cache2    = NULL;
    cache3    = NULL;
    max = 0;
    max_enhance_factor = 10;

    daughter_pos=0;
    for (int i=0; i< ANY_DISTRIBUTION_MAX_DAUGHTERS;i++) {
        daughters[i]= NULL;
    }
} ;

PDistribution* PAnyDistribution::Clone(const char*delme) const {
    return new PAnyDistribution((const PAnyDistribution &)* this);
};

void PAnyDistribution::MakeVars() {
    if (!projector) projector = new PProjector();

    vparent   = makeDynamicData()->GetBatchParticle("_parent");  
    vf        = makeStaticData()->GetBatchValue("_f"); 

    x = makeStaticData()->GetBatchValue("_x"); 
    y = makeStaticData()->GetBatchValue("_y"); 
    z = makeStaticData()->GetBatchValue("_z"); 


};

Bool_t PAnyDistribution::AddEquation(char * command) {
    MakeVars();
    return projector->AddCommand(command);
};

Bool_t PAnyDistribution::AddEquation(TH1 * histo, char * command) {
    MakeVars();
    cache1 = histo;
    return projector->AddHistogram(histo,command,2);
};

Bool_t PAnyDistribution::AddEquation(TH2 * histo, char * command) {
    MakeVars();
    cache2 = histo;
    return projector->AddHistogram(histo,command,2);
};

Bool_t PAnyDistribution::AddEquation(TH3 * histo, char * command) {
    MakeVars();
    cache3 = histo;
    return projector->AddHistogram(histo,command,2);
};

// Bool_t PAnyDistribution::AddHistogram(TH2 * histo, const char * command) {
//     MakeVars();
//     return  projector->AddHistogram(histo,command,0);
// };


Bool_t PAnyDistribution::Init(void) {

    daughter_pos=0; //clear stuff because the Attach function makes a clone

    //get the parent
    for (int i=0; i<position; i++) {
	if (particle_flag[i] == PARTICLE_LIST_PARENT)
	    parent=particle[i];
    }
    if (!parent) {
	Warning("Init","Parent not found");
	return kFALSE;
    }

    
    //Now get all daughters
    for (int i=0; i < ANY_DISTRIBUTION_MAX_DAUGHTERS;i++) {
        daughters[i]= GetParticle("daughter");
        if (daughters[i]) {         
            daughter_pos++;
        }
    }
 
    return kTRUE;    
};

Bool_t PAnyDistribution::Prepare(void) {
    return kTRUE;
};

Bool_t PAnyDistribution::Finalize(void) {
    return kTRUE;
};

Bool_t PAnyDistribution::CheckAbort(void) {
    return kFALSE;
};

Bool_t PAnyDistribution::IsValid(void) {

    double factor = 1.;

    if (projector) {	

	*vparent = parent;
	//vparent->Boost(-parent->BoostVector()); 

	stack[0]=vparent;
	for (int i=0; i < ANY_DISTRIBUTION_MAX_DAUGHTERS;i++) {
	    stack[i+1] = daughters[i];
	}
	//projector->Execute();
	Int_t num = daughter_pos + 1;
	projector->Modify((PParticle **)&stack,(int *)&decay_done, &num, 
			  ANY_DISTRIBUTION_MAX_DAUGHTERS);
	factor =  *vf;
    } else {
	Fatal("IsValid","No equation set");
    }
    
    //cout << "Factor: " << factor <<  endl;


    Double_t myres=0.;
    if (cache1) {
	int bin = cache1->FindBin(*x);
	myres = cache1->GetBinContent(bin);
	//dw= hist3D->GetBinError(bin);
	if (cache1->GetEntries()) myres /= cache1->GetEntries();
	myres *=cache1-> GetNbinsX();
    } else if (cache2) {
	int bin = cache2->FindBin(*x,*y);
	myres = cache2->GetBinContent(bin);
	if (cache2->GetEntries()) myres /= cache2->GetEntries();
	myres *= cache1-> GetNbinsX() * cache1-> GetNbinsY();
    } else if (cache3) {
	int bin = cache3->FindBin(*x,*y,*z);
	myres = cache3->GetBinContent(bin);
	if (cache3->GetEntries()) myres /= cache3->GetEntries();
	///gsffsfsfsdfsdf
    } 

    
    Double_t reduction_factor = 1.;
    if (myres) {
	reduction_factor /= myres;
    }

    if (reduction_factor > max_enhance_factor) reduction_factor = max_enhance_factor;

    factor *= reduction_factor;

    if (factor>max) {
	max = factor * 1.1;
	Warning("IsValid","Factor > max, increased (%f)",max);
    }
 
    if ((factor/max)>PUtils::sampleFlat()) return kTRUE; // sample now 
    
    return kFALSE;

};


