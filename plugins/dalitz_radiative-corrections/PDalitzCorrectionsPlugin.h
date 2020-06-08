#ifndef _PDALITZCORRECTIONSPLUGIN_H_
#define _PDALITZCORRECTIONSPLUGIN_H_

#include "PDistributionManagerUtil.h"
#include "PDistributionCollection.h"

#include "PRadiativeCorrections.h"


using namespace std;

class PDalitzCorrectionsPlugin : public PDistributionCollection {

 public:

    PDalitzCorrectionsPlugin();
    PDalitzCorrectionsPlugin(const Char_t *id, const Char_t *de);
    ~PDalitzCorrectionsPlugin();

    Bool_t ExecCommand(const char * command, Double_t value);

    Bool_t Activate(void);

 private:

    PRadiativeCorrections *rad_corrections;


    ClassDef(PDalitzCorrectionsPlugin,0)
};

#endif // _PDALITZCORRECTIONSPLUGIN_H_
