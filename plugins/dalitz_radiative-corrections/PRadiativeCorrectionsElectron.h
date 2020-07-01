#ifndef _PRADIATIVECORRECTIONSELECTRON_H_
#define _PRADIATIVECORRECTIONSELECTRON_H_

#include "PRadiativeCorrections.h"

#include "dilepton_radiative_corrections.h"

class PRadiativeCorrectionsElectron : public PRadiativeCorrections {

  public:
    PRadiativeCorrectionsElectron(){}
    PRadiativeCorrectionsElectron(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const override;

    Double_t GetWeight() override;

  private:
    void SetLimits() override;

    TGraph2D *corrections_pi0, *corrections_eta, *corrections_etap;
    TH2D *h_corrections_eta, *h_corrections_etap;
    bool pi0, eta, etap;

    ClassDef(PRadiativeCorrectionsElectron,0)
};

#endif
