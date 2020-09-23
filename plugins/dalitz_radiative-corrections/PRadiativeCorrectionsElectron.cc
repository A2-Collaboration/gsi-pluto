#include "PRadiativeCorrectionsElectron.h"

using namespace std;


ClassImp(PRadiativeCorrectionsElectron)

PRadiativeCorrectionsElectron::PRadiativeCorrectionsElectron(const Char_t *id, const Char_t *de, Int_t key) :
    PRadiativeCorrections(id, de, key) {

    pi0 = eta = etap = false;

    //TODO: include pi0 corrections

    corrections_eta = new TGraph2D("eta_ee_corrections", "Radiative Corrections #eta #to e^{+} e^{-} #gamma",
                                   eta_ee::corr.size(), &eta_ee::x[0], &eta_ee::y[0], &eta_ee::corr[0]);
    corrections_etap = new TGraph2D("etap_ee_corrections", "Radiative Corrections #eta' #to e^{+} e^{-} #gamma",
                                    etap_ee::corr.size(), &etap_ee::x[0], &etap_ee::y[0], &etap_ee::corr[0]);

    // calling GetHistogram() method avoids Interpolate() returning 0 if an exact point is being interpolated
    h_corrections_eta = corrections_eta->GetHistogram();
    h_corrections_etap = corrections_etap->GetHistogram();
}

PDistribution* PRadiativeCorrectionsElectron::Clone(const char*delme) const {
    return new PRadiativeCorrectionsElectron((const PRadiativeCorrectionsElectron &)* this);
}

void PRadiativeCorrectionsElectron::SetLimits()
{
    TGraph2D* corr;

    if (parent->GetParent()->Is("eta")) {
        corr = corrections_eta;
        y_vals = &eta_ee::y_vals;
        x_tuples = &eta_ee::x_tuples;
        x_vec = &eta_ee::x;
        y_vec = &eta_ee::y;
        z_vec = &eta_ee::corr;
    } else if (parent->GetParent()->Is("eta'")) {
        corr = corrections_etap;
        y_vals = &etap_ee::y_vals;
        x_tuples = &etap_ee::x_tuples;
        x_vec = &etap_ee::x;
        y_vec = &etap_ee::y;
        z_vec = &etap_ee::corr;
    }

    weight_max += corr->GetZmax()/100.;

    limits_set = true;
}

Double_t PRadiativeCorrectionsElectron::GetWeight()
{
    // Set limits here since the grandparent (meson) is not accesible during Init()
    if (!limits_set)
        SetLimits();

    meson = parent->GetParent();
    //pi0 = eta = etap = false;
    if (!meson->IsMeson()) {
        Warning("GetWeight","Grandparent is not a meson");
        return 1.;
    }
    /*if (meson->Is("pi0"))
        pi0 = true;
    else if (meson->Is("eta"))
        eta = true;
    else if (meson->Is("eta'"))
        etap = true;
    else {
        Warning("GetWeight","No pseudo-scalar meson found");
        return 1.;
    }*/

    double q2 = parent->M2();
    double im2 = meson->M2();
    double x = q2/im2;
    double y = meson->Vect4().Dot(lp->Vect4()-lm->Vect4());
    //double nu2 = 4*lm->M2()/im2;
    //double beta = sqrt(1-nu2/x);
    // use absolute value for y since correction values are just provided for positive y values
    // y should be symmetric so under this assumption everything should be fine
    y = 2*abs(y)/meson->M2()/(1-x);

    //if ((x < nu2) || (x > 1.))
    //    cerr << "x value outside of kinematical bounds: x = " << x << " not in [" << nu2 << " , 1]" << endl;
    //if ((y < 0.) || (y > beta))
    //    cerr << "y value outside of kinematical bounds: y = " << y << " not in [0 , " << beta << "]" << endl;

    double weight = 1.;
    //double correction = 0.;

    /*if (pi0)
        correction = corrections_pi0->Interpolate(x, y);
    else if (eta)
        correction = corrections_eta->Interpolate(x, y);
    else if (etap)
        correction = corrections_etap->Interpolate(x, y);
    else
        Warning("GetWeight","No correction table found for this channel");*/

    weight += ApproximateValue(x, y)/100.;

    return weight;
}
