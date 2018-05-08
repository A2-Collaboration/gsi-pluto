////////////////////////////////////////////////////////
//  This file reads the Pluto particle array
//  and writes the values to the dilepton data file
//
////////////////////////////////////////////////////////

#include "../src/PFileOutput.h"
#include "../src/PUtils.h"

class PGiBuuOutput : public PFileOutput {

private:

    FILE * dataFile;           //Pointer to data file
    PParticle *eplus,*eminus;  //Private pointers
    Int_t event_number;

 public:
    
    PGiBuuOutput(char * filename);

    Bool_t Modify(PParticle ** stack, int *decay_done, int * num, int stacksize);  //bulk interface
    bool CloseFile(void);


    ClassDef(PGiBuuOutput,0) 


};

PGiBuuOutput::PGiBuuOutput(char * filename) {
    //Constructor
    dataFile = fopen(filename,"w"); //open the data file for writing
    if (!dataFile) Fatal("PiBuuOutput","File %s cannot be opened",filename);
    fPriority =99;  //This is very important, never change this value
    event_number=0;
}

bool PGiBuuOutput::CloseFile(void) {
    if (dataFile) fclose(dataFile);
    cout << "Data file closed" << endl;
    return kTRUE;
}

bool PGiBuuOutput::Modify(PParticle ** stack, int *decay_done, int * num, int stacksize) {

    double e, px, py, pz, w;
    int    id, flag;

    fprintf(dataFile, "# %i\n",event_number);
    event_number++;

    for (int i = 0; i< (*num); i++) {
	//loop over all particles

	if (stack[i]->Is("e+")) {
	    //found e+
	    eplus = stack[i];
	    e    = eplus->E();
	    px   = eplus->Px();
	    py   = eplus->Py();
	    pz   = eplus->Pz();
	    w    = eplus->W();
	    id   = eplus->GetSourceId();
	    flag = eplus->IsActive();
	    if (flag)
		if (!fprintf(dataFile, "+ %lf %lf %lf %lf %lf %d %d\n" , 
			     e, px, py, pz, w, id, flag)) return kFALSE;
	    
	} else if (stack[i]->Is("e-")) {
	    //found e-
	    eminus = stack[i];
	    e    = eminus->E();
	    px   = eminus->Px();
	    py   = eminus->Py();
	    pz   = eminus->Pz();
	    w    = eminus->W();
	    id   = eminus->GetSourceId();
	    flag = eminus->IsActive();
	    if (flag)
		if (!fprintf(dataFile, "- %lf %lf %lf %lf %lf %d %d\n" , 
			     e, px, py, pz, w, id, flag)) return kFALSE;
	} 
 
    }
    return kTRUE;
};

















