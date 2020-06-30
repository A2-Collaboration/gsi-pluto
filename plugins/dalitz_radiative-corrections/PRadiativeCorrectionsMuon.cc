#include "PRadiativeCorrectionsMuon.h"

using namespace std;


ClassImp(PRadiativeCorrectionsMuon)

PRadiativeCorrectionsMuon::PRadiativeCorrectionsMuon(const Char_t *id, const Char_t *de, Int_t key) :
    PRadiativeCorrections(id, de, key) {

    eta = etap = false;

    //TODO: include pi0 corrections

    corrections_eta = new TGraph2D("eta_mumu_corrections", "Radiative Corrections #eta #to #mu^{+} #mu^{-} #gamma",
                                   eta_mumu::corr.size(), &eta_mumu::x[0], &eta_mumu::y[0], &eta_mumu::corr[0]);
    corrections_etap = new TGraph2D("etap_mumu_corrections", "Radiative Corrections #eta' #to #mu^{+} #mu^{-} #gamma",
                                    etap_mumu::corr.size(), &etap_mumu::x[0], &etap_mumu::y[0], &etap_mumu::corr[0]);

    // calling GetHistogram() method avoids Interpolate() returning 0 if an exact point is being interpolated
    corrections_eta->GetHistogram();
    corrections_etap->GetHistogram();
}

PDistribution* PRadiativeCorrectionsMuon::Clone(const char*delme) const {
    return new PRadiativeCorrectionsMuon((const PRadiativeCorrectionsMuon &)* this);
}

void PRadiativeCorrectionsMuon::SetMaximumWeight()
{
    double correction = 0.;

    if (parent->GetParent()->Is("eta"))
        correction = corrections_eta->GetZmax();
    else if (parent->GetParent()->Is("eta'"))
        correction = corrections_etap->GetZmax();

    weight_max += correction/100.;
    weight_max_set = true;
}

Double_t PRadiativeCorrectionsMuon::GetWeight()
{
    if (!weight_max_set)
        SetMaximumWeight();

    meson = parent->GetParent();
    eta = etap = false;
    if (!meson->IsMeson()) {
        Warning("GetWeight","Grandparent is not a meson");
        return 1.;
    }
    if (meson->Is("eta"))
        eta = true;
    else if (meson->Is("eta'"))
        etap = true;
    else {
        Warning("GetWeight","No pseudo-scalar meson found");
        return 1.;
    }

    double x = parent->M2()/meson->M2();
    double y = meson->Vect4().Dot(lp->Vect4()-lm->Vect4());
    // use absolute value for y since y should be symmetric
    y = 2*abs(y)/meson->M2()/(1-x);

    double weight = 1.;
    double correction = 0.;

    if (eta)
        correction = corrections_eta->Interpolate(x, y);
    else if (etap)
        correction = corrections_etap->Interpolate(x, y);
    else
        Warning("GetWeight","No correction table found for this channel");

    weight += correction/100.;

    return weight;
}
