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


    if (rad_corrections) {
        Info("Init","Found radiative corrections");
        rad_corrections->Print();
    }

    return kTRUE;
}

Bool_t PRadiativeCorrections::IsValid() {
    //Use rejection mode...

/*    if (GetVersionFlag() & VERSION_WEIGHTING) return kTRUE;
    //...but not if weighting enabled.
*/
    cout << "+++ Calculate x and y in IsValid +++" << endl;
    meson = parent->GetParent();
    double q2 = parent->M2();
    double im2 = meson->M2();
    double x = q2/im2;
    double y = meson->Vect4().Dot(lp->Vect4()-lm->Vect4());
    y = 2*abs(y)/meson->M2()/(1-x);
    double nu2 = 4*lm->M2()/im2;
    double beta = sqrt(1-nu2/x);
    cout << "[DEBUG]   x = " << x << "  y = " << y << endl;
    if ((x < nu2) || (x > 1.))
        return kFALSE;
    if ((y < 0.) || (y > beta))
        return kFALSE;

    return kTRUE;
}

string PRadiativeCorrections::get_base(const std::string& path)
{
    size_t pos = path.find_last_of("/");
    return (std::string::npos == pos) ? "" : path.substr(0, pos);
}
