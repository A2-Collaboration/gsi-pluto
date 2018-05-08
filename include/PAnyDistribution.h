// Author: I. Froehlich
// Written: 12.8.2011

#ifndef _PANYDISTRIBUTION_H_
#define _PANYDISTRIBUTION_H_

#define ANY_DISTRIBUTION_MAX_DAUGHTERS 8

#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "PDistribution.h"
#include "PProjector.h"


class PAnyDistribution : public PDistribution  {
  
 public:
    PAnyDistribution();
    PAnyDistribution(const Char_t *id, const Char_t *de);
    PDistribution *Clone(const char *delme=NULL) const;

    Bool_t Init(void);
    Bool_t Prepare(void);
    Bool_t Finalize(void);
    Bool_t IsNotRejected(void);
    Bool_t CheckAbort(void);

    Bool_t AddEquation(TH1 *histo, const char *command);  //adds an equation + cache for non-uniform distributions
    Bool_t AddEquation(TH2 *histo, const char *command);  //adds an equation + cache for non-uniform distributions
    Bool_t AddEquation(TH3 *histo, const char *command);  //adds an equation + cache for non-uniform distributions
    Bool_t AddEquation(const char *command);
    
    void SetMaxEnhancementFactor(Double_t _max_enhance_factor) {
	max_enhance_factor = _max_enhance_factor;
    };

 private:
    PParticle *daughters[ANY_DISTRIBUTION_MAX_DAUGHTERS]; //up to 8 daughters
    PParticle *stack[ANY_DISTRIBUTION_MAX_DAUGHTERS+1];
    int decay_done[ANY_DISTRIBUTION_MAX_DAUGHTERS+1];
    int daughter_pos;
    PParticle *parent;

    PProjector *projector;

    void MakeVars(void);

    TH1 *cache1;
    TH2 *cache2;
    TH3 *cache3;

    PParticle *vparent;
    Double_t  *vf;
    Double_t *x, *y, *z;
    Double_t max, max_enhance_factor;

    ClassDef(PAnyDistribution, 0)  //Any possible distribution for a->b+c+d 

};

#endif


