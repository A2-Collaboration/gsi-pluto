// Author: I. Froehlich
// Written: 14.9.2012
// Revised: 

#ifndef _PF3_H_
#define _PF3_H_

#include "TF3.h"
#include "PProjector.h"


class PF3: public TF3 {

 public:

    
    //Standard ROOT ctors
    PF3();
    PF3(const char *name, const char *formula, Double_t xmin=0, Double_t xmax=1, Double_t ymin=0, Double_t ymax=1, Double_t zmin=0, Double_t zmax=1);
    //PF3(const char *name, void *fcn, Double_t xmin=0, Double_t xmax=1, Double_t ymin=0, Double_t ymax=1, Double_t zmin=0, Double_t zmax=1, Int_t npar=0);

    PF3(const char *name, Double_t (*fcn)(Double_t *, Double_t *), Double_t xmin=0, Double_t xmax=1, 
	Double_t ymin=0, Double_t ymax=1, Double_t zmin=0, Double_t zmax=1, Int_t npar=0);
    PF3(const char *name, Double_t (*fcn)(const Double_t *, const Double_t *), Double_t xmin=0, Double_t xmax=1, 
	Double_t ymin=0, Double_t ymax=1, Double_t zmin=0, Double_t zmax=1, Int_t npar=0);
    PF3(const char *name, ROOT::Math::ParamFunctor f, Double_t xmin = 0, Double_t xmax = 1,
	Double_t ymin=0, Double_t ymax=1, Double_t zmin=0, Double_t zmax=1, Int_t npar = 0);  
    PF3(const PF3 &f3);
  
    //Pluto ctor
    PF3(const char *name, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, Double_t zmin, Double_t zmax);

    virtual Double_t EvalPar(const Double_t *x, const Double_t *params);
    //TF wrapper

  
    Bool_t Add(char *command);
    Bool_t Add(TH1  *histo, const char *command = "");
    Bool_t Add(TH2  *histo, const char *command = "");
    Bool_t Add(TH3  *histo, const char *command = "");
    Bool_t Add(TGraph *graph, const char *command = "");
    Bool_t Add(TGraph2D *graph, const char *command = "");


 protected:

    PProjector *projector;
    Double_t   *vf, *vx, *vy, *vz;


    ClassDef(PF3, 0)  // Pluto TF3 extension

};

#endif
