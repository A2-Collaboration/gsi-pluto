// Author: I. Froehlich
// Written: 19.2.2009
// Revised: 
// PFormula Class Header

#ifndef _PFORMULA_H_
#define _PFORMULA_H_

#include "v5/TFormula.h"


class PFormula: public ROOT::v5::TFormula {

 public:

    PFormula();
    PFormula(const char *name,const char *expression);
    
    TString error_string;
    Int_t   error_code;
    TString chaine;

    void Analyze(const char *schain, Int_t &err, Int_t offset);

 protected:

    ClassDef(PFormula,0)  // Adapted TFormula class

};

#endif
