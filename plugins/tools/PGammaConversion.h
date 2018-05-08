// Author:
// Written: 12/01/2014
// Modified: 
// PGammaConversion Class Header

#ifndef _PGAMMACONVERSION_H_
#define _PGAMMACONVERSION_H_

#include "PBulkInterface.h"

#include "TH2F.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TString.h"

#include <vector>

using namespace std;

#define MAX_GAMMACONVERSION_DILEPTONS 20
#define NEFUNCS 250


class PGammaConversion: public PBulkInterface {

 private:

    PParticle ep, em;
    Float_t Z;
    Double_t x1, y1, z1;
    Float_t convpar[6];
    Float_t ConvProb;

    static Double_t dSde(Double_t *x, Double_t *par); // has to be static to be used by TF1 object
    Float_t ConversionXS(Float_t E);
    TF1 *Elepton[NEFUNCS];

    TRandom3 *RGen;
    TF1 *ulep;


    //------------------------------------
    // histograms
    Int_t nVertex;
    Int_t nbinTheta;
    Int_t nbinR;
    Int_t nbinZ;
    TH3F *hdimension;  // keep one hist for axis ranges
    TH1F *hvertexgeom; // keep target positions
    vector<TH2F*> fhzr;
    Int_t findVertexBin( Float_t zVertexEvt);
    Int_t getThetaBin(Double_t thetadeg);
    TH2F *getHist(Double_t zVertexEvt, Double_t thetadeg);
    //------------------------------------

 
 protected:
    
    
 public:
    
    PGammaConversion();
    PGammaConversion(Float_t Z, Float_t p,Int_t run);
    ~PGammaConversion();
    //void SetAverageCharge(Float_t z) {Z = z;}  // set average atomic charge of converter material
    void    SetConvLookupParameters(Float_t c1, Float_t c2, Float_t c3) {convpar[0]=c1; convpar[1]=c2; convpar[2]=c3;};
    Float_t GetConversionProb(Float_t E, Float_t theta,Int_t run);
    void    SetConversionProb(Float_t p) {ConvProb = (p>0 ? p : 0);}  // fix conversion probability
    Int_t   readHist(TString inputroot);

    void    SetVertex(Double_t x, Double_t y, Double_t z) {
	x1 = x; y1 = y; z1 = z;
    };
    
    Bool_t Modify(PParticle **stack, int *decay_done, int *num, int stacksize);  //bulk interface

    ClassDef(PGammaConversion, 0) // Fast post processing gamma conversion in material
};
#endif 

















