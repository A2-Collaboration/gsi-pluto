// Author: I. Froehlich
// Written: 14.8.2012
// Revised: 

#ifndef _PSCATTERCROSSSECTION_H_
#define _PSCATTERCROSSSECTION_H_


#include "TF1.h"
#include "TF2.h"
#include "PAngularDistribution.h"
#include "PProjector.h"
#include "PF2.h"


class PScatterCrossSection : public PAngularDistribution  {
  
 public:
    PScatterCrossSection();
    PScatterCrossSection(const Char_t *id, const Char_t *de);
    PDistribution *Clone(const char *delme=NULL) const;

    Bool_t Init(void);
    Bool_t Prepare(void);

    Bool_t AddEquation(char *command);
    Bool_t AddHistogram(TH2 *histo, char *command);
    Bool_t AddHistogram(TH1 *histo, char *command);

    using TF1::SetRange;
    void SetRange(Double_t mymin, Double_t mymax) {
	TF1::SetRange(qmin, qmax);
	qmax = mymax;
	qmin = mymin;
    };

    void SetNpx(Int_t my_npx);
    void SetNpy(Int_t my_npy);
  
    PF2 *GetFunction(void) {return pf2;};

 protected:
    
    double SamplePolarAngle(double);

    Double_t qmin,qmax;
    Double_t costheta, q;

    PF2 *pf2;
    
    PParticle *vprimary;
    Double_t  *vf;
    Bool_t MakeVars(void);

    Int_t npx, npy;

    ClassDef(PScatterCrossSection, 0)  //Cross section model for the reaction a+b -> c+d

};

#endif


