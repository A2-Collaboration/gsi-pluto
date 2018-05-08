// Author: I. Froehlich
// Written: 18.8.2010
// Revised: 

#ifndef _PF2_H_
#define _PF2_H_

#include "TF2.h"
#include "PProjector.h"


class PF2: public TF2 {

 public:

    
    //Standard ROOT ctors
    PF2();
    PF2(const char *name, const char *formula, Double_t xmin=0, Double_t xmax=1, Double_t ymin=0, Double_t ymax=1);
    //PF2(const char *name, void *fcn, Double_t xmin=0, Double_t xmax=1, Double_t ymin=0, Double_t ymax=1, Int_t npar=0);

    PF2(const char *name, Double_t (*fcn)(Double_t *, Double_t *), Double_t xmin=0, Double_t xmax=1, Double_t ymin=0, Double_t ymax=1, Int_t npar=0);
    PF2(const char *name, Double_t (*fcn)(const Double_t *, const Double_t *), Double_t xmin=0, Double_t xmax=1, Double_t ymin=0, Double_t ymax=1, Int_t npar=0);
    PF2(const char *name, ROOT::Math::ParamFunctor f, Double_t xmin = 0, Double_t xmax = 1, Double_t ymin = 0, Double_t ymax = 1, Int_t npar = 0);  
    PF2(const PF2 &f2);
  
    //Pluto ctor
    PF2(const char *name, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax);

    virtual Double_t    EvalPar(const Double_t *x, const Double_t *params);
    //TF wrapper

  
    Bool_t Add(const char *command);
    Bool_t Add(TH1 *histo, const char *command = "");
    Bool_t Add(TH2 *histo, const char *command = "");
    Bool_t Add(TH3 *histo, const char *command = "");
    Bool_t Add(TGraph *graph, const char *command = "");
    Bool_t Add(TGraph2D *graph, const char *command = "");

    void SetEpsilon(Double_t e) {
	epsilon = e;
    };

    Bool_t MakeIntegral(void);
    void GetRandom2(Double_t &xrandom, Double_t &yrandom);


 protected:

    using TF2::Integral;
    Double_t Integral(Double_t ax, Double_t bx, Double_t ay, Double_t by, Double_t epsilon);

    Double_t epsilon;

    PProjector *projector;
    Double_t   *vf, *vx, *vy;


    ClassDef(PF2, 0)  // Pluto TF2 extension

};

#endif
