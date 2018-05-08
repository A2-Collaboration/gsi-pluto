// Author: I. Froehlich
// Written: 18.7.2007
// Revised: 

#ifndef _PBEAMSMEARING_H_
#define _PBEAMSMEARING_H_


#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "PDistribution.h"


class PBeamSmearing : public PDistribution  {
  
 public:
    PBeamSmearing();
    ~PBeamSmearing();
    PBeamSmearing(const Char_t *id, const Char_t *de);
    
    PDistribution *Clone(const char *delme=NULL) const;

    Bool_t Init(void);
    Bool_t Prepare(void);
    
    void Print(const Option_t *delme=NULL) const;  //Debug info
    
    void SetMomentumSmearing(TF1 *f) {
	//A function stretching the beam momentum (x=1 means no modification)
	mom_smearing = f;
    };

    void SetMomentumSmearing(TH1 *f) {
	//A function stretching the beam momentum (x=1 means no modification)
	mom_smearing_histo = f;
    };
    
    void SetMomentumFunction(TF1 *f) {
	//A function which desribes the absolute momentum distribution
	mom_function = f;
    };
    
    void SetMomentumFunction(TH1 *f) {
	//A function which desribes the absolute momentum distribution
	mom_function_histo = f;
    };
    
    void SetBeamParameters(Double_t th, Double_t ph, Double_t sig) {
	//N.B. these are in units DEGREE because this is more convenient!
	thetaBeam = th * TMath::Pi() / 180.;
	phiBeam   = ph * TMath::Pi() / 180.;
	sigmaBeam = sig* TMath::Pi() / 180.;
    };

    void SetAngularSmearing(TF1 *f) {
	//smearing distribution in degree
	ang_smearing=f;
    };
    
    void SetReaction(const char *reaction);

 private:

    PParticle *beam,   *mybeam;
    PParticle *target, *mytarget;
    PParticle *parent;

    TF1     *mom_smearing;
    TF1     *mom_function;
    TF1     *ang_smearing;
    TH1     *mom_smearing_histo;
    TH1     *mom_function_histo;
    Double_t thetaBeam;           // polar angle of beam
    Double_t phiBeam;             // azimuthal angle of beam
    Double_t sigmaBeam;           // smearing of beam
    ClassDef(PBeamSmearing, 0)    // General purpose beam smearing models
};

#endif


