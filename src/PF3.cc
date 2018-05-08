////////////////////////////////////////////////////////
//  
// Adapted TF1
//
////////////////////////////////////////////////////////

#include <math.h>
#include <iostream>
#include "PUtils.h"
using namespace std;

#ifdef WIN32
#pragma optimize("",off)
#endif

#include "PF3.h"

//Standard ROOT ctors

PF3::PF3() :
    TF3() {
    projector = NULL;    
};

PF3::PF3(const char *name, const char *formula, Double_t xmin, Double_t xmax, Double_t, Double_t, Double_t, Double_t) :
    TF3(name, formula, xmin, xmax) {
    projector = NULL;
};

#if 0
PF3::PF3(const char *name, void *fcn, Double_t xmin, Double_t xmax, Double_t, Double_t, Double_t, Double_t, Int_t npar) :
    TF3(name, fcn, xmin, xmax, npar) {
    projector = NULL;
};
#endif

PF3::PF3(const char *name, Double_t (*fcn)(Double_t *, Double_t *), 
	 Double_t xmin, Double_t xmax, Double_t, Double_t, Double_t, Double_t, Int_t npar) :
    TF3(name, fcn, xmin, xmax, npar) {
    projector = NULL;
};

PF3::PF3(const char *name, ROOT::Math::ParamFunctor f, Double_t xmin, Double_t xmax, Double_t, Double_t, Double_t, Double_t, Int_t npar) :
    TF3(name, f, xmin, xmax, npar) {
    projector = NULL;
};  

PF3::PF3(const PF3 &f3) :
    TF3(f3) {
    projector = NULL;
};

//Pluto ctor
PF3::PF3(const char *name, Double_t xmin, Double_t xmax, Double_t, Double_t, Double_t, Double_t) :
    TF3(name,"x+y+z", xmin, xmax) {
    SetTitle(name);
    projector = new PProjector();
    vf = makeStaticData()->GetBatchValue("_f"); 
    vx = makeStaticData()->GetBatchValue("_x"); 
    vy = makeStaticData()->GetBatchValue("_y"); 
    vz = makeStaticData()->GetBatchValue("_z"); 
}

Bool_t PF3::Add(char *command) {
    if (!projector) {
	Warning("Add", "No PProjector, wrong ctor called?");
	return kFALSE;
    }
    return projector->AddCommand(command);
};

Bool_t PF3::Add(TH1 *histo, const char *command) {
    if (!projector) {
	Warning("Add", "No PProjector, wrong ctor called?");
	return kFALSE;
    }
    return  projector->AddHistogram(histo, command, 0);
};

Bool_t PF3::Add(TH2 *histo, const char *command) {
    if (!projector) {
	Warning("Add", "No PProjector, wrong ctor called?");
	return kFALSE;
    }
    return  projector->AddHistogram(histo, command, 0);
};

Bool_t PF3::Add(TH3 *histo, const char *command) {
    if (!projector) {
	Warning("Add", "No PProjector, wrong ctor called?");
	return kFALSE;
    }
    return  projector->AddHistogram(histo, command, 0);
};

Bool_t PF3::Add(TGraph *graph, const char *command) {
    if (!projector) {
	Warning("Add", "No PProjector, wrong ctor called?");
	return kFALSE;
    }
    return  projector->AddTGraph(graph, command);
};

Bool_t PF3::Add(TGraph2D *graph, const char *command) {
    if (!projector) {
	Warning("Add", "No PProjector, wrong ctor called?");
	return kFALSE;
    }
    return  projector->AddTGraph2D(graph, command);
};

Double_t PF3::EvalPar(const Double_t *x, const Double_t *params) {
    if (!projector) return TF3::EvalPar(x, params);

    //Extension via batch script:
    *vx = x[0];
    *vy = x[1];
    *vz = x[2];

    Int_t num=0;
    if (projector->Modify(NULL, NULL, &num, 0))
	return  *vf;
    return 0;
}



ClassImp(PF3)
