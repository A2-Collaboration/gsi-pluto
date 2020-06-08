#include "PDataBase.h"
#include "PDalitzCorrectionsPlugin.h"


PDalitzCorrectionsPlugin::PDalitzCorrectionsPlugin() {
}


PDalitzCorrectionsPlugin::PDalitzCorrectionsPlugin(const Char_t *id, const Char_t *de):
    PDistributionCollection(id, de) {

    rad_corrections = NULL;

}

PDalitzCorrectionsPlugin::~PDalitzCorrectionsPlugin() {
}

Bool_t PDalitzCorrectionsPlugin::Activate(void) {

    return kTRUE;

};


Bool_t PDalitzCorrectionsPlugin::ExecCommand(const char * command, Double_t value) {

    PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();

    //called the 1st time?
    if (!rad_corrections) {
        rad_corrections = new PRadiativeCorrections("radiative_corrections@dilepton_to_e+_e-/corrections", "Radiative corrections Dalitz decay", -1);
        rad_corrections->EnableWeighting();
        pdmutil->Add(rad_corrections);
    }

    return kTRUE;
}

ClassImp(PDalitzCorrectionsPlugin)

