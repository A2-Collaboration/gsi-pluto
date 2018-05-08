// Author: I. Froehlich
// Written: 3.7.2006
// Revised: 

#ifndef _PDALITZDISTRIBUTION_H_
#define _PDALITZDISTRIBUTION_H_

#include "TF1.h"
#include "TF2.h"
#include "PDistribution.h"
#include "PProjector.h"


class PDalitzDistribution : public PDistribution  {
  
 public:
    PDalitzDistribution();
    PDalitzDistribution(const Char_t *id, const Char_t *de);
    PDistribution* Clone(const char *delme=NULL) const;

    Bool_t Init(void);
    Bool_t Prepare(void);
    Bool_t Finalize(void);
    Bool_t IsNotRejected(void);
    Bool_t CheckAbort(void);

    void SetSlopes(Double_t myslope1, Double_t myslope2) {slope1=myslope1;slope2=myslope2;};
    Bool_t AddEquation(const char *command);
    Bool_t AddHistogram(TH2 *histo, const char *command = "");
    Bool_t Do(char *command) {
	return AddEquation(command);
    };
    Bool_t Do(TH2 *histo, const char *command = "") {
	return AddHistogram(histo, command);
    };

    void SetMax(Double_t mymax) {
	max = mymax;
    };

 private:
    PParticle *side_particle[2]; //2 additional particles
    PParticle *primary;
    PParticle *parent;

    PParticle *vprimary, *vs1, *vs2;
    Double_t  *vf;

    Double_t slope1, slope2,max;
    PProjector *projector;

    void MakeVars(void);

    ClassDef(PDalitzDistribution, 0)  //Dalitz plane decay slopes for a->b+c+d 
};

#endif


