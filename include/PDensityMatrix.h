// Author: Ingo Froehlich
// Written: 24/04/2013

#ifndef _PDENSITYMATRIX_H_
#define _PDENSITYMATRIX_H_

#include "PEmbeddedParticles.h"
#include "PProjector.h"

#define DENSITYMATRIX_MAX_MATRICES 16
#define DENSITYMATRIX_MAX          1000

class PDensityMatrix: public PEmbeddedParticles {

 private:

    Double_t *matrix[DENSITYMATRIX_MAX_MATRICES], 
	*matrix_integral[DENSITYMATRIX_MAX_MATRICES];

    Double_t *axes[3];
    Int_t     axes_size[3];
    Int_t     dimension;
    Int_t     matrixsize;
    Int_t     current_matrix;

    Bool_t    IsBorder(Int_t bin);
    
    PProjector *projector;
    Double_t *x,*y,*z;
    
 public:
    
    PDensityMatrix();
    
    Bool_t Modify(PParticle **stack, int *decay_done, int *num, int stacksize);  //bulk interface
    
    Bool_t ReadDensityMatrix(const char *filename, Int_t dim, Bool_t use_bin_width,
			     Double_t min_selection, Double_t max_selection);

    Bool_t   GetBin(Double_t *x,    Int_t *bins);
    Bool_t   GetBin(Int_t     bin,  Int_t *bins);
    Double_t GetBinWidth(Int_t dim, Int_t bin);
    Int_t    GetRandomBin(Int_t num);

    void     SetMatrix(Int_t i) {current_matrix = i;};

    Bool_t   Do(const char *command);

    ClassDef(PDensityMatrix, 0) // Add particles from a density matrix
};
#endif 

















