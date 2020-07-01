#include "PRadiativeCorrections.h"

#include <iostream>

using namespace std;


ClassImp(PRadiativeCorrections)

PRadiativeCorrections::PRadiativeCorrections() {
}

PRadiativeCorrections::PRadiativeCorrections(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de, key) {

    if (is_channel<0)
        Warning("PRadiativeCorrections","This model should be bound to CHANNELS only");

    parent = nullptr;
    meson = nullptr;
    lp = nullptr;
    lm = nullptr;

    limits_set = false;
    weight_max = 1.;
    x_min = 0.;
}

Bool_t PRadiativeCorrections::Init() {
    //Init function called once for each PChannel

    //looking for parent. This is mandatory
    parent = GetParticle("parent");
    if (!parent) {
        Warning("Init","Parent not found");
        return kFALSE;
    }
    bool muon;
    if (parent->Is("dilepton"))
        muon = false;
    else if (parent->Is("dimuon"))
        muon = true;
    else {
        Warning("Init","Channel is not a dilepton/dimuon decay");
        return kFALSE;
    }

    lp = muon ? GetParticle("mu+") : GetParticle("e+");
    if (!lp) {
        Warning("Init","Positron / mu+ not found");
        return kFALSE;
    }
    lm = muon ? GetParticle("mu-") : GetParticle("e-");
    if (!lm) {
        Warning("Init","Electron / mu- not found");
        return kFALSE;
    }

    rad_corrections = GetSecondaryModel("corrections");

    if (!rad_corrections) {
        Warning("Init","Model for radiative corrections not found");
        return kFALSE;
    }

    return kTRUE;
}

Bool_t PRadiativeCorrections::IsValid() {
    //Use rejection mode...

    if (GetVersionFlag() & VERSION_WEIGHTING) return kTRUE;
    //...but not if weighting enabled.

    meson = parent->GetParent();
    double q2 = parent->M2();
    double im2 = meson->M2();
    double x = q2/im2;
    double y = meson->Vect4().Dot(lp->Vect4()-lm->Vect4());
    y = 2*abs(y)/meson->M2()/(1-x);
    double nu2 = 4*lm->M2()/im2;
    double beta = sqrt(1-nu2/x);
    if ((x < nu2) || (x > 1.))
        return kFALSE;
    if ((y < 0.) || (y > beta))
        return kFALSE;

    // valid x and y values, do rejection now
    double weight = GetWeight();
    // weight 1 is returned if no value has been found in the correction tables
    // return true in this case since weight_max can be smaller than 1
    // (if only negative correction values are present)
    if (weight == 1.)
        return kTRUE;

    if ((weight/weight_max) > PUtils::sampleFlat())
        return kFALSE;

    return kTRUE;
}
