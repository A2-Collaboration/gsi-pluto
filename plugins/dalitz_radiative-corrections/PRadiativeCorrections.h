#ifndef _PRADIATIVECORRECTIONS_H_
#define _PRADIATIVECORRECTIONS_H_

#include "PChannelModel.h"
#include "PDynamicData.h"
#include "PKinematics.h"

#include "TGraph2D.h"

class PRadiativeCorrections : public PChannelModel {

 public:
    PRadiativeCorrections();
    PRadiativeCorrections(const Char_t *id, const Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;

    Bool_t Init();

    using PDistribution::GetWeight;
    using PChannelModel::GetWeight;

    Double_t GetWeight(void);

    Bool_t IsValid();

 protected:

    PParticle *parent, *meson, *lp, *lm;

    PChannelModel *rad_corrections;

 private:
    TGraph2D *corrections_pi0, *corrections_eta_ee, *corrections_etap_ee, *corrections_eta_mumu, *corrections_etap_mumu;

    bool pi0, eta, etap;
    bool muon;

    ClassDef(PRadiativeCorrections,0)
};

#endif
