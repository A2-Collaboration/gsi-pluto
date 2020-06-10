#include "PRadiativeCorrectionsElectron.h"

using namespace std;


ClassImp(PRadiativeCorrectionsElectron)

PRadiativeCorrectionsElectron::PRadiativeCorrectionsElectron(const Char_t *id, const Char_t *de, Int_t key) :
    PRadiativeCorrections(id, de, key) {

    pi0 = eta = etap = false;

    corrections_eta = new TGraph2D(eta_ee_corrections.size());
    corrections_etap = new TGraph2D(etap_ee_corrections.size());

    int i = 0;
    for (const auto c : eta_ee_corrections)
        corrections_eta->SetPoint(i++, c.x, c.y, c.c);
    i = 0;
    for (const auto c : etap_ee_corrections)
        corrections_etap->SetPoint(i++, c.x, c.y, c.c);

    // calling GetHistogram() method avoids Interpolate() returning 0 if an exact point is being interpolated
    corrections_eta->GetHistogram();
    corrections_etap->GetHistogram();
}

PDistribution* PRadiativeCorrectionsElectron::Clone(const char*delme) const {
    return new PRadiativeCorrectionsElectron((const PRadiativeCorrectionsElectron &)* this);
}

Double_t PRadiativeCorrectionsElectron::GetWeight()
{
    meson = parent->GetParent();
    pi0 = eta = etap = false;
    if (!meson->IsMeson()) {
        Warning("GetWeight","Grandparent is not a meson");
        return 1.;
    }
    if (meson->Is("pi0"))
        pi0 = true;
    else if (meson->Is("eta"))
        eta = true;
    else if (meson->Is("eta'"))
        etap = true;
    else {
        Warning("GetWeight","No pseudo-scalar meson found");
        return 1.;
    }

    //meson->Print();

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
//    cout << "[DEBUG]   x elem [" << nu2 << " , 1]" << endl;
//    cout << "[DEBUG]   y elem [0 , " << beta << "]" << endl;

    if ((x < nu2) || (x > 1.))
        cerr << "x value outside of kinematical bounds: x = " << x << " not in [" << nu2 << " , 1]" << endl;
    if ((y < 0.) || (y > beta))
        cerr << "y value outside of kinematical bounds: y = " << y << " not in [0 , " << beta << "]" << endl;

    double weight = 1.;
    double correction = 0.;

    if (pi0)
        correction = corrections_pi0->Interpolate(x, y);
    else if (eta)
        correction = corrections_eta->Interpolate(x, y);
    else if (etap)
        correction = corrections_etap->Interpolate(x, y);
    else
        Warning("GetWeight","No correction table found for this channel");

    weight += correction/100.;
    cout << "[DEBUG]   x = " << x << "  y = " << y << "   apply weight: " << weight << endl;

    return weight;
}
