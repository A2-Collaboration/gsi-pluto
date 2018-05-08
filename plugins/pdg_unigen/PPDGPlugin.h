// Author: I. Froehlich
// Written: 02.02.2011
// Revised: 
// 


#ifndef _PPDG_PLUGIN_H_
#define _PPDG_PLUGIN_H_

#include "TROOT.h"

#include "PChannelModel.h"
#include "PDistributionManagerUtil.h"
#include "PDistributionCollection.h"

using namespace std;

class PPDGPlugin : public PDistributionCollection {
    
 public:

    //constructor
    PPDGPlugin(const Char_t *id, const Char_t *de);
    //destructor
    ~PPDGPlugin();

    Bool_t ExecCommand(const char *command, Double_t value); 

    Bool_t Activate(void);

 private:

    Int_t is_initialized, is_extend_resonances;

    ClassDef(PPDGPlugin, 0) //Plugin to add pdg code for UniGen
};

#endif //_PPDG_PLUGIN_H_







