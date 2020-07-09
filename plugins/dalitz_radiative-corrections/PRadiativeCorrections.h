#ifndef _PRADIATIVECORRECTIONS_H_
#define _PRADIATIVECORRECTIONS_H_

#include "PChannelModel.h"

#include "TGraph2D.h"
#include "TH2D.h"

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

    virtual void SetLimits() = 0;
    virtual void SetXminmax(const double);
    bool limits_set;

    Double_t weight_max;
    Double_t y_max, x_min, x_max;
    std::vector<double>* y_vals = nullptr;
    std::vector<double>* x_tuples = nullptr;


    ClassDef(PRadiativeCorrections,0)
};

#endif
