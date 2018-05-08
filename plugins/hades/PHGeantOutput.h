// Author: Jochen Markert
// Written: 05/12/2007
// Modified: 
// PHGeantOutput Class Header

#ifndef _PHGEANTOUTPUT_H_
#define _PHGEANTOUTPUT_H_

#include "PFileOutput.h"
#include "PParticle.h"
#include <stdio.h>

class PHGeantOutput: public PFileOutput {

 private:

    Double_t *event_impact_param, *seqnr;

 protected:

     FILE  *asciiFile;
     Int_t  ctEvt;
     Int_t  ctParticlePerEvt;
     Bool_t writeSEQNUMBER;
 public:

    PHGeantOutput();

    bool OpenFile(const char *_filename);      //filename
    bool CloseFile(void);                //
    bool WriteEventHeader(void);          
    bool WriteParticle(PParticle *par);  //write one particle
    void SetWriteSeqNumber(Bool_t write) { writeSEQNUMBER = write; }

    ClassDef(PHGeantOutput, 0) // Pluto ascii file output for HGEANT
};
#endif 

















