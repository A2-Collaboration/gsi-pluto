// Author: Ingo Froehlich
// Written: 01/02/2011


#ifndef _PUNIGENINPUT_H_
#define _PUNIGENINPUT_H_

#include "PBulkInterface.h"
#include "UEvent.h"
#include "URun.h"
#include "TTree.h"
#include "TFile.h"

#define MAX_UNIGENPARTICLES 500
#define MAX_UNIGENUNKNOWNPDG 10

class PUniGenInput: public PBulkInterface {

 private:

    //PParticle *local      [MAX_UNIGENPARTICLES];

    //From batch system
    //Double_t *vertex_x,*vertex_y,*vertex_z;

    TFile  *fInFile;               // Input file
    TTree  *fInTree;               // Input tree
    URun   *fRun;                  // Run object
    UEvent *fEvent;                // Event
    
    Int_t centry, nentries, reset_mass;

    Int_t pdg_param, pid_param, *i_result;

    Int_t unknown_pdg[MAX_UNIGENUNKNOWNPDG];
    Int_t unknown_pdg_pointer;

 protected:
    
    
 public:
    
    PUniGenInput();
    PUniGenInput(char *filename);
    
    Bool_t Input(char *filename);

    Bool_t Modify(PParticle **stack, int *decay_done, int *num, int stacksize);  //bulk interface
    
    void ResetInvalidMass(void) {reset_mass=1;};
    
    ClassDef(PUniGenInput, 0) // Adds particles from an UniGen file into a PReaction
};
#endif 

















