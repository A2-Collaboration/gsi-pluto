// Author: I. Froehlich
// Written: 14.9.2012
// Revised: 

#ifndef _PF1_H_
#define _PF1_H_

#include "TF1.h"
#include "PProjector.h"


class PF1: public TF1 {

 public:

    
    //Standard ROOT ctors
    PF1();
    PF1(const char *name, const char *formula, Double_t xmin=0, Double_t xmax=1);
    //PF1(const char *name, void *fcn, Double_t xmin=0, Double_t xmax=1, Int_t npar=0);

    PF1(const char *name, Double_t (*fcn)(Double_t *, Double_t *), Double_t xmin=0, Double_t xmax=1, Int_t npar=0);
    PF1(const char *name, Double_t (*fcn)(const Double_t *, const Double_t *), Double_t xmin=0, Double_t xmax=1, Int_t npar=0);
    PF1(const char *name, ROOT::Math::ParamFunctor f, Double_t xmin = 0, Double_t xmax = 1, Int_t npar = 0);  
    PF1(const PF1 &f1);
  
    //Pluto ctor
    PF1(const char *name, Double_t xmin, Double_t xmax);

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
    Double_t   *vf, *vx, *vy;


    ClassDef(PF1, 0)  // Pluto TF1 extension

};

#endif
