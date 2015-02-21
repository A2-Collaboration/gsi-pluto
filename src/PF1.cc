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

#include "PF1.h"

//Standard ROOT ctors

PF1::PF1() :
    TF1() {
    projector = NULL;    
};

PF1::PF1(const char *name, const char *formula, Double_t xmin, Double_t xmax) :
    TF1(name, formula, xmin, xmax) {
    projector = NULL;
};

PF1::PF1(const char *name, void *fcn, Double_t xmin, Double_t xmax, Int_t npar) :
    TF1(name, fcn, xmin, xmax, npar) {
    projector = NULL;
};


PF1::PF1(const char *name, Double_t (*fcn)(Double_t *, Double_t *), 
	 Double_t xmin, Double_t xmax, Int_t npar) :
    TF1(name,fcn,xmin,xmax,npar) {
    projector = NULL;
};

PF1::PF1(const char *name, ROOT::Math::ParamFunctor f, Double_t xmin, Double_t xmax, Int_t npar) :
    TF1(name, f, xmin, xmax,npar) {
    projector = NULL;
};  

PF1::PF1(const PF1 &f1) :
    TF1(f1) {
    projector = NULL;
};

//Pluto ctor
PF1::PF1(const char *name, Double_t xmin, Double_t xmax) :
    TF1(name,"x", xmin, xmax) {
    SetTitle(name);
    projector = new PProjector();
    vf = makeStaticData()->GetBatchValue("_f"); 
    vx = makeStaticData()->GetBatchValue("_x"); 
}

Bool_t PF1::Add(char * command) {
    if (!projector) {
	Warning("Add","No PProjector, wrong ctor called?");
	return kFALSE;
    }
    return projector->AddCommand(command);
};

Bool_t PF1::Add(TH1 * histo, const char * command) {
    if (!projector) {
	Warning("Add","No PProjector, wrong ctor called?");
	return kFALSE;
    }
    return  projector->AddHistogram(histo,command,0);
};

Bool_t PF1::Add(TH2 * histo, const char * command) {
    if (!projector) {
	Warning("Add","No PProjector, wrong ctor called?");
	return kFALSE;
    }
    return  projector->AddHistogram(histo,command,0);
};

Bool_t PF1::Add(TH3 * histo, const char * command) {
    if (!projector) {
	Warning("Add","No PProjector, wrong ctor called?");
	return kFALSE;
    }
    return  projector->AddHistogram(histo,command,0);
};

Bool_t PF1::Add(TGraph * graph, const char * command) {
    if (!projector) {
	Warning("Add","No PProjector, wrong ctor called?");
	return kFALSE;
    }
    return  projector->AddTGraph(graph,command);
};

Bool_t PF1::Add(TGraph2D * graph, const char * command) {
    if (!projector) {
	Warning("Add","No PProjector, wrong ctor called?");
	return kFALSE;
    }
    return  projector->AddTGraph2D(graph,command);
};

Double_t PF1::EvalPar(const Double_t *x, const Double_t *params) {
    if (!projector) return TF1::EvalPar(x,params);

    //Extension via batch script:
    *vx = x[0];

    Int_t num=0;
    if (projector->Modify(NULL, NULL, &num, 0))
	return  *vf;
    return 0;
}



ClassImp(PF1)
