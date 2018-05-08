////////////////////////////////////////////////////////
//  File input interface implementation file
//
//  This class serves as template for file output classes
//
//                    Author:  Ingo Froehlich
//                    Written: 14/05/2007
//                    Revised:
//
////////////////////////////////////////////////////////

#include "PFileOutput.h"


PFileOutput::PFileOutput() {
    cnt       = -1; 
    getVERTEX = 0;
    filename  = "<user file>";
}

bool PFileOutput::OpenFile(const char *) {
    return kFALSE; //pure virtual not allowed
}

bool PFileOutput::CloseFile(void) {
    return kFALSE; //pure virtual not allowed
}

bool PFileOutput::WriteEvent(void) {
    //next event
    return kFALSE; //pure virtual not allowed
}

bool PFileOutput::WriteEventHeader(void) {
    return kFALSE; //pure virtual not allowed
}

bool PFileOutput::WriteBranchHeader(void) {
    return kFALSE; //pure virtual not allowed
}

bool PFileOutput::WriteParticle(PParticle *) {
    //write one particle
    return kFALSE; //pure virtual not allowed
}


ClassImp(PFileOutput) 
