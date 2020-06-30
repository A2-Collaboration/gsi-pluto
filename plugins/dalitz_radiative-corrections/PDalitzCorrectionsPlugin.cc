#include "PDataBase.h"
#include "PDalitzCorrectionsPlugin.h"


PDalitzCorrectionsPlugin::PDalitzCorrectionsPlugin() {
}

PDalitzCorrectionsPlugin::PDalitzCorrectionsPlugin(const Char_t *id, const Char_t *de):
    PDistributionCollection(id, de) {

    rad_corrections_ee = nullptr;
    rad_corrections_mumu = nullptr;

}

PDalitzCorrectionsPlugin::~PDalitzCorrectionsPlugin() {
}

Bool_t PDalitzCorrectionsPlugin::Activate(void) {

    return kTRUE;
};


Bool_t PDalitzCorrectionsPlugin::ExecCommand(const char * command, Double_t value) {

    PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();

    //called the 1st time?
    if (!rad_corrections_ee) {
        rad_corrections_ee = new PRadiativeCorrectionsElectron("radiative_corrections_dilepton@dilepton_to_e+_e-/corrections", "Radiative corrections Dalitz decay", -1);
        pdmutil->Add(rad_corrections_ee);
    }
    if (!rad_corrections_mumu) {
        rad_corrections_mumu = new PRadiativeCorrectionsMuon("radiative_corrections_dimuon@dimuon_to_mu+_mu-/corrections", "Radiative corrections Dalitz decay", -1);
        pdmutil->Add(rad_corrections_mumu);
    }

    return kTRUE;
}

ClassImp(PDalitzCorrectionsPlugin)
