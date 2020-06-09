#include "PRadiativeCorrectionsMuon.h"

using namespace std;


ClassImp(PRadiativeCorrectionsMuon)

PRadiativeCorrectionsMuon::PRadiativeCorrectionsMuon(const Char_t *id, const Char_t *de, Int_t key) :
    PRadiativeCorrections(id, de, key) {

    eta = etap = false;

    //TODO: include pi0 corrections

    const string pwd = get_base(__FILE__);
    corrections_eta = new TGraph2D((pwd+"/eta_mumu_corrections").c_str());
    corrections_etap = new TGraph2D((pwd+"/etap_mumu_corrections").c_str());

    // calling GetHistogram() method avoids Interpolate() returning 0 if an exact point is being interpolated
    corrections_eta->GetHistogram();
    corrections_etap->GetHistogram();
}

PDistribution* PRadiativeCorrectionsMuon::Clone(const char*delme) const {
    return new PRadiativeCorrectionsMuon((const PRadiativeCorrectionsMuon &)* this);
}

Double_t PRadiativeCorrectionsMuon::GetWeight()
{
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

    double q2 = parent->M2();
    double im2 = meson->M2();
    //double x = parent->M2()/meson->M2();
    double x = q2/im2;
    //double y = meson->Vect4()*(lp-lm);
    double y = meson->Vect4().Dot(lp->Vect4()-lm->Vect4());
    double nu2 = 4*lm->M2()/im2;
    double beta = sqrt(1-nu2/x);
    // use absolute value for y since correction values are just provided for positive y values
    // y should be symmetric so under this assumption everything should be fine
    y = 2*abs(y)/meson->M2()/(1-x);

    if ((x < nu2) || (x > 1.))
        cerr << "x value outside of kinematical bounds: x = " << x << " not in [" << nu2 << " , 1]" << endl;
    if ((y < 0.) || (y > beta))
        cerr << "y value outside of kinematical bounds: y = " << y << " not in [0 , " << beta << "]" << endl;

    double weight = 1.;
    double correction = 0.;

    if (eta)
        correction = corrections_eta->Interpolate(x, y);
    else if (etap)
        correction = corrections_etap->Interpolate(x, y);
    else
        Warning("GetWeight","No correction table found for this channel");

    weight += correction/100.;
    cout << "[DEBUG]   x = " << x << "  y = " << y << "   apply weight: " << weight << endl;

    return weight;
}
