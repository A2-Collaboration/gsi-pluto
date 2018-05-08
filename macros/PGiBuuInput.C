////////////////////////////////////////////////////////
//  This file reads the GiBUU dipelton data files
//  and convert the values to PParticles
//  
//  The contructed PParticles are added to the particle 
//  array of the event loop
//
////////////////////////////////////////////////////////

#include "../src/PBulkInterface.h"
#include "../src/PUtils.h"

class PGiBuuInput : public PBulkInterface {

private:

    FILE * dataFile;
    PParticle eplus,eminus;

 public:
    
    PGiBuuInput(char * filename);

    Bool_t Modify(PParticle ** stack, int *decay_done, int * num, int stacksize);  //bulk interface

    ClassDef(PGiBuuInput,0) 

};

PGiBuuInput::PGiBuuInput(char * filename) {
    dataFile = fopen(filename,"r");
    if (!dataFile) Fatal("PiBuuInput","File %s not found",filename);

    eplus.SetID(2);  //id of e+
    eminus.SetID(3); //id of e-
}

bool PGiBuuInput::Modify(PParticle ** stack, int *decay_done, int * num, int stacksize) {

    double e, px, py, pz, w;
    int    id, flag;

    //Read until EOF. If the number of particles
    //in the event loop is larger, the event loop is stopped
    if (fscanf (dataFile, "+ %lf %lf %lf %lf %lf %d %d\n" , 
		&e, &px, &py, &pz, &w, &id, &flag) == EOF) return kFALSE;
    
    eplus.SetPxPyPzE(px,py,pz,e);
    eplus.SetW(w);
    eplus.SetSourceId(id);
    
    if (fscanf (dataFile, "- %lf %lf %lf %lf %lf %d %d\n" , 
		&e, &px, &py, &pz, &w, &id, &flag) == EOF) return kFALSE;
    eminus.SetPxPyPzE(px,py,pz,e);
    eminus.SetW(w);
    eminus.SetSourceId(id);
    
    //put the particles on the stack:
    *(stack[(*num)++]) = eplus;
    *(stack[(*num)++]) = eminus;

    return kTRUE;
};

















