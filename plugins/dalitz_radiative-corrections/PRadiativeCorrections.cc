#include "PRadiativeCorrections.h"

#include <iostream>

#include "radiative_corrections_etap.h"

using namespace std;


ClassImp(PRadiativeCorrections)

PRadiativeCorrections::PRadiativeCorrections() {
}

PRadiativeCorrections::PRadiativeCorrections(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de, key) {

    if (is_channel<0)
        Warning("PRadiativeCorrections","This model should be bound to CHANNELS only");

    pi0 = eta = etap = false;

    parent = NULL;
    meson = NULL;
    lp = NULL;
    lm = NULL;

    corrections_etap_ee = new TGraph2D(etap::corr_values.size());
    size_t i = 0;
    for (auto x : etap::x_values)
        for (auto y : etap::y_values)
            corrections_etap_ee->SetPoint(i, x, y, etap::corr_values.at(i++));

    cout << "[PRadiativeCorrections] Test some correction points" << endl;
    cout << "x = 0.1, y = 0.2: " << corrections_etap_ee->Interpolate(.1, .2) << endl;
    cout << "x = 0.34, y = 0.61: " << corrections_etap_ee->Interpolate(.34, .61) << endl;
    cout << "x = 0.621, y = 0.97: " << corrections_etap_ee->Interpolate(.621, .97) << endl;
    cout << "x = 0, y = 1: " << corrections_etap_ee->Interpolate(0, 1) << endl;
    cout << "End tests" << endl;

} ;

PDistribution* PRadiativeCorrections::Clone(const char*delme) const {
    return new PRadiativeCorrections((const PRadiativeCorrections &)* this);
};

Bool_t PRadiativeCorrections::Init() {
    //Init function called once for each PChannel

    //looking for parent. This is mandatory
    parent = GetParticle("parent");
    if (!parent) {
        Warning("Init","Parent not found");
        return kFALSE;
    }
    if (parent->Is("dilepton"))
        muon = false;
    else if (parent->Is("dimuon"))
        muon = true;
    else {
        Warning("Init","Channel is not a dilepton/dimuon decay");
        return kFALSE;
    }
//    meson = parent->GetParent();
//    if (!meson->IsMeson()) {
//        Warning("Init","Grandparent is not a meson");
//        return kFALSE;
//    }
//    if (meson->Is("pi0"))
//        pi0 = true;
//    else if (meson->Is("eta"))
//        eta = true;
//    else if (meson->Is("eta'"))
//        etap = true;
//    else {
//        Warning("Init","No pseudo-scalar meson found");
//        return kFALSE;
//    }
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


    if (rad_corrections)
        Info("Init","Found radiative corrections");

    return kTRUE;
}


Double_t PRadiativeCorrections::GetWeight(void) {

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
    if (pi0 && muon) {
        Warning("GetWeight","pi0 Dalitz decay in mu+ mu- found, impossible");
        return 1.;
    }

    //meson->Print();

    double q2 = parent->M2();
    double im2 = meson->M2();
//    cout << "[DEBUG]   q2: " << q2 << " GeV^2; IM(eta'): " << meson->M() << " GeV" << endl;
//    cout << "[DEBUG]   angle between e+ and e-: " << TMath::RadToDeg()*lm->Angle(lp->Vect()) << " deg" << endl;
//    cout << "[DEBUG]   angle between e+ and eta': " << TMath::RadToDeg()*meson->Angle(lp->Vect()) << " deg" << endl;
    //double x = parent->M2()/meson->M2();
    double x = q2/im2;
    //double y = meson->Vect4()*(lp-lm);
    double y = meson->Vect4().Dot(lp->Vect4()-lm->Vect4());
    double nu2 = 4*lm->M2()/im2;
    double beta = sqrt(1-nu2/x);
//    cout << "[DEBUG]   nu = " << sqrt(nu2) << endl;
//    cout << "[DEBUG]   beta = " << beta << endl;
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

    //TODO: split muon and electron into separate classes
    if (pi0 && !muon)
        correction = corrections_pi0->Interpolate(x, y);
    else if (eta && !muon)
        correction = corrections_eta_ee->Interpolate(x, y);
    else if (eta && muon)
        correction = corrections_eta_mumu->Interpolate(x, y);
    else if (etap && !muon)
        correction = corrections_etap_ee->Interpolate(x, y);
    else if (etap && muon)
        correction = corrections_etap_mumu->Interpolate(x, y);
    else
        Warning("GetWeight","No correction table found for this channel");

    weight += correction/100.;
    cout << "[DEBUG]   x = " << x << "  y = " << y << "   apply weight: " << weight << endl;

    return weight;
}

Bool_t PRadiativeCorrections::IsValid() {
    //Use rejection mode...

    if (GetVersionFlag() & VERSION_WEIGHTING) return kTRUE;
    //...but not if weighting enabled.

    return kTRUE;
}
