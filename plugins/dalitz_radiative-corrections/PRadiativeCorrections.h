#ifndef _PRADIATIVECORRECTIONS_H_
#define _PRADIATIVECORRECTIONS_H_

#include "PChannelModel.h"

#include "TGraph2D.h"

class PRadiativeCorrections : public PChannelModel {

 public:
    PRadiativeCorrections();
    PRadiativeCorrections(const Char_t *id, const Char_t *de, Int_t key);

    Bool_t Init();

    using PDistribution::GetWeight;
    using PChannelModel::GetWeight;

    virtual Double_t GetWeight() = 0;

    Bool_t IsValid();

 protected:
    PParticle *parent, *meson, *lp, *lm;

    PChannelModel *rad_corrections;

    virtual void SetMaximumWeight() = 0;
    bool weight_max_set = false;

    Double_t weight_max;


    ClassDef(PRadiativeCorrections,0)
};

#endif
