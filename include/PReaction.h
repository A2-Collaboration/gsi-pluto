// Author: M.A. Kagarlis
// Written: 03.02.99
// Revised: 21.04.05  R.H.
// PReaction Class Header

#ifndef _PREACTION_H_
#define _PREACTION_H_

#include "TClonesArray.h"
#include "TTree.h"
#include "TFile.h"
#include "TMethodCall.h"
#include "PChannel.h"
// #include "G__ci.h"

#include "PDistributionManager.h"
#include "PStdData.h"
#include "PFileOutput.h"
#include "PBulkInterface.h"
#include "PProjector.h"
#include "PPlutoBulkDecay.h"

#define MAX_FILEOUTPUT 5


#define MIN_PARTICLE_STACKSIZE 500
#define INC_PARTICLE_STACKSIZE 100
#define MAX_PARTICLE_STACKSIZE 10000
#define MAX_REACTION_FILTERS   100


class TPythia6;
class PFilter;

Int_t select(PParticle*);

class PReaction : public TObject  {

    friend class PDecayManager;
 public:

    PReaction();

    PReaction(const char *filename);

    PReaction(PChannel **, const char *, int n=2, int f0=1, int f1=0, int f2=0, int f3=0, TTree *ttree=NULL);

    PReaction(PChannel **, int n=2, unsigned int ff=0, TTree *ttree=NULL, const char *filename=NULL);

    PReaction(Double_t momentum, const char *beam, const char *target, const char *reaction, const char *file_name=NULL,
	      Int_t f0=1, Int_t f1=0, Int_t f2=0, Int_t f3=0, TTree *ttree=NULL); 
    // build reaction from a descriptor string
    PReaction(const char *command,  const char *beam, const char *target, 
	      const char *reaction, const char *file_name=NULL,
	      Int_t f0=1, Int_t f1=0, Int_t f2=0, Int_t f3=0, TTree *ttree=NULL); 
    // build reaction from a descriptor string but uses energy parser
    PReaction(PParticle *, const char *reaction, const char *file_name=NULL,
	      Int_t f0=1, Int_t f1=0, Int_t f2=0, Int_t f3=0, TTree *ttree=NULL); 
    // build reaction from a descriptor string but uses a seed particle

    void AddReaction(const char *reaction);
    //add a reaction to the list

    ~PReaction();
    // Reaction destructor

    void InitChannels();

    int loop(int i=-1, int wf=0, int verbose=1) {
	//kept for backward compatibility
	return Loop(i, wf, verbose);
    };

    int Loop(int i=-1, int wf=0, int verbose=1);
    // number of events to simulate

    int GetNumAll() { return ndpar; }
    // returns number of all the particles in the reaction

    int *GetCIndex() { 
	return cindex.GetArray(); }
    // returns product channel index

    void VertexYes() {
	// turn on vertex calculation after initialization
	if (IsExtTree()) return;
	getVERTEX = 1;
	ropt = ropt | ff2;
    }

    void VertexNo()  {
	// turn off vertex calculation after initialization
	if (IsExtTree()) return;
	getVERTEX = 0;
	ropt=~((~ropt)|ff2);
    }


    void allParticles() {
	// put all particles on root file after initialization
	if (IsExtTree()) return;
	allPARTICLES = 1;
	ropt = ropt | ff0;
    }

    void trackedParticles() {
	// put tracked particles on root file file after initialization
	if (IsExtTree()) return;
	allPARTICLES = 0;
	ropt=~((~ropt)|ff0);
    }

    void asciiYes() {
	// enable ascii output for GEANT after initialization
	if (IsExtTree()) return;
	asciiOUTPUT = 1;
	ropt = ropt | ff3;
	SetName(filename);
    }

    void asciiNo() {
	// disable ascii output for GEANT after initialization
	if (IsExtTree()) return;
	asciiOUTPUT = 0;
	ropt=~((~ropt)|ff3);
	SetName(filename);
    }

    void SetName(const char *);
    // sets up output file names

    PParticle **GetProducts() { return particle_stack; }
    // returns product pointer

    void SetFilter(int ch_id, PFilter *filter);
    // This function is invoked by the Filter class. It directs loop to
    // impose a filter at address (*filter) following the decay represented
    // by the reaction Channel number ch_id. See also PFilter.

    PFilter **GetFilter() { return NULL; }
    // returns pointer to array of existing filters

    void Print(const Option_t *delme=NULL) const;
    void PrintReport() const; //Print a final report after sampling

    void Close();
    // close root output file

    void SetHGeant(int flag) {HGeant = (flag!=0);}  // set to 1, if PLUTO is run
    // from the HGeant prompt 
    void setHGeant(int flag) {SetHGeant(flag);}; //backward comp.

    void SetUserSelection(Int_t (*f)(PParticle*)) { // compiled user selection
	SetUserSelection((void*)f);
    }
    void SetUserSelection(void *f);                 // interpreted user selection

    void SetUserAnalysis(Int_t (*f)(PParticle**,Int_t)) { // compiled analysis
	SetUserAnalysis((void*)f);
    }
    void SetUserAnalysis(void *f);                  // interpreted user analysis

    void SetTrigCond(Int_t n) {nTrigCond = n;}      // set trigger level

    void setDecayAll(Float_t tau=1.) { 
	Warning ("setDecayAll", "This method is depreciated: use SetDecayAll()");
	SetDecayAll(tau);
    }

    void SetDecayAll(Float_t tau=1.) {              
	// decay all particles with
	// lifetime < tau (in ns)
	// choose either Pluto or Pythia, depending on availability

	//if (allPARTICLES && (asciiOUTPUT||HGeant) ) {
	//    Warning ("setDecayAll",
	//	     "\nOptions decayALL & allPARTICLES & (asciiOUTPUT||HGeant) incompatible!\n\n");
	//    return;
	//}

	//Wrapper to the new scheme
	tauMax = tau;   // go to sec
    

#ifdef USE_PYTHIA6
    
	PPythiaBulkDecay *pl = new PPythiaBulkDecay();
	pl->SetPythia(fPythia);
	pl->SetTauMax(tauMax);
	AddBulk(pl);

#else

	PPlutoBulkDecay *pl = new PPlutoBulkDecay();
	pl->SetRecursiveMode(1);
	pl->SetTauMax(tauMax);
	AddBulk(pl);

#endif
    
    }


    void SetPythia(TPythia6 *p) {fPythia=p;}        // set pointer to Pythia

    void SetMaxFileSize(Int_t bytes) {nMaxBytes = bytes;}  // set max file size

#if 0
    Int_t testPointer(void *p2f) {   // determine type of pointer-to-function
	// (see CINT reference manual for details) 
	if (p2f == NULL) return -1;
	Int_t ret = G__isinterpretedp2f(p2f);
	char *fname;
	fname = G__p2f2funcname(p2f);
	printf("Pointer to function %s is of type %d\n",fname,ret);
	return ret;
    }
#endif

    void SetWriteIndex(Bool_t flag) {
	if (flag == kTRUE) writeINDEX = 1;
	else writeINDEX = 0;
    }

    void DisableWeightReset() {weight_reset=0;};

    Int_t GetReactionId() { return reactionId; }

    PDistributionManager *GetDistributionManager(void){
	return makeDistributionManager();
    }

    Bool_t AddFileOutput(PFileOutput *file) {
	if (fileoutput_pos == MAX_FILEOUTPUT ) {
	    Warning("AddFileOutput", "MAX_FILEOUTPUT reached");
	    return kFALSE;
	}
	files[fileoutput_pos++] = file;
	return kTRUE;
    }

    Bool_t AddBulk(PBulkInterface *mybulk);
    Bool_t AddPrologueBulk(PBulkInterface *mybulk); //Bulk IO before any decay

    Bool_t Do(const char *command) {
	return GetCurrentProjector()->AddCommand(command);
    }
    Bool_t Do(TH1 *f, const char *command, Int_t flag=1) {
	return GetCurrentProjector()->AddHistogram(f,command,flag);
    }
    Bool_t Do(TH2 *f, const char *command, Int_t flag=1) {
	return GetCurrentProjector()->AddHistogram(f,command,flag);
    }
    Bool_t Do(TH3 *f, const char *command, Int_t flag=1) {
	return GetCurrentProjector()->AddHistogram(f,command,flag);
    }
    Bool_t Output(TNtuple *f, const char *command = (char *)"") {
	return GetCurrentProjector()->AddOutputTNtuple(f,command);
    }
    Bool_t Input(TNtuple *f) {
	return GetCurrentProjector()->AddInputTNtuple(f);
    }
    Bool_t Output(const char *f, const char *command="") {
	return GetCurrentProjector()->AddOutputASCII(f,command);
    }
    Bool_t Input(const char *f, const char *command="") {
	return GetCurrentProjector()->AddInputASCII(f, command);
    }
    Bool_t CloseFile(void) {
	return GetCurrentProjector()->CloseFile();
    }

    PProjector *GetCurrentProjector(void) {
	if (!bulkdecay_pos) {
	    current_projector = new PProjector();
	    AddBulk(current_projector);
	} else {
	    if ((strcmp("PProjector", bulk[bulkdecay_pos-1]->GetName()) == 0)) {
		if ((bulk[bulkdecay_pos-1]->GetPriority() > FILTER_PRIORITY) ||
		    ((PProjector *)bulk[bulkdecay_pos-1])->IsFileOpen()) {
		    current_projector = (PProjector *)bulk[bulkdecay_pos-1];
		} else {
		    current_projector = new PProjector();
		    AddBulk(current_projector);
		}
	    } else {		
		current_projector = new PProjector();
		AddBulk(current_projector);
	    }
/* 	    if ((strcmp("PProjector",bulk[bulkdecay_pos-1] -> GetName())==0)  */
/* 		&& (bulk[bulkdecay_pos-1] -> GetPriority() > FILTER_PRIORITY)) { */
/* 		current_projector = (PProjector *)bulk[bulkdecay_pos-1]; */
/* 	    } else { */
		
/* 		current_projector = new PProjector(); */
/* 		AddBulk(current_projector); */
/* 	    } */
	}
	return current_projector;
    }

    void Preheating(Int_t num) {pre_heating=num;}; //Perform num dummy decays

    void IsInline(void) {
	//pure inline PReaction, e.g. FAIRROOT
	is_inline = 1;
    }

 private:

    Int_t reactionId;                          //  reaction identifier
    Int_t (*userSelection)(PParticle*);        //! pointer to selection function 
    Int_t (*userAnalysis)(PParticle**,Int_t);  //! pointer to analysis function 

    Int_t num_filters;           //! Filters from PBatch commands
    Int_t filter_keys[MAX_REACTION_FILTERS], filter_counter[MAX_REACTION_FILTERS]; 
    Double_t *filter_values[MAX_REACTION_FILTERS]; 

    int nchan, ntpar, ndpar, nclones, loop_count, reset_count, status;
    // number of reaction channels (channels),
    // number of non-decayed (tracked) particles,
    // number of decay products,
    // number of filters
    // number of times loop function invoked
    // number of reset calls
    // error status
    static int activeCnt;  // number of clones for root branches,

    TArrayI cindex, dindex, ftrack;
    // channel index of the decay products,
    // indices of parents to channels after the 1st
    // flag indicating whether particle is tracked (1) or not (0)

    unsigned int ropt, allPARTICLES, asciiOUTPUT, writeINDEX,
	resetCHANNELS, getVERTEX, extTREE, HGeant, decayALL, inactivate_decayed_particles;
    static const unsigned int ff0, ff1, ff2, ff3, ff4;
    // option flags (see constructor)

    Float_t tauMax;   // max lifetime for decay (see decayALL option) 

    PChannel **channel;
    // address of array of pointers to reaction channels

    TString filename, file1, file2, original_filename,reaction_string;
    // output file name, filename.evt for (GEANT), filename.root (ROOT)
    void ConvertFilename(void);
    
    TTree *tree;
    // address of the event tree

    PParticle **particle_stack;
    // Particle stack for the event loop
    int decay_done[MAX_PARTICLE_STACKSIZE];
    int stacksize;

    PParticle **particle;
    //particles from the PChannel list

    static TClonesArray *evt[MAX_NUM_BRANCHES+1];
    // particle clone

    TFile *rootfile;
    // the root output file

    int size_branches;
    int key_branches[MAX_NUM_BRANCHES];
    int doonce;

    int fileoutput_pos;
    PFileOutput *files[MAX_FILEOUTPUT];
    int bulkdecay_pos, pro_bulkdecay_pos;
    
    PBulkInterface *bulk[MAX_BULKDECAY];
    PBulkInterface *pro_bulk[MAX_BULKDECAY];
    PProjector *current_projector;

  
    Int_t nTrigCond;

    TPythia6 *fPythia;                   //! pointer to Pythia object

    void SetReactionId(void);
    // calculate reaction id from list of products of 1st channel

    void SetUp(PChannel **);
    // get the channels and particles, identify the physics, set up the defaults.

    Bool_t parse_script(const char *e, const char *beam, const char *target, const char *reaction, const char *file_name=NULL,
			Int_t f0=1, Int_t f1=0, Int_t f2=0, Int_t f3=0, TTree *ttree=NULL); 

    void InitLoop();
    // resets the dynamical objects of the reaction and opens files

    Int_t IsExtTree() { return extTREE; }
    // check whether tree is external

    Int_t ParseChannel(PParticle *parent, const char *channel,
		       TList &plutoList, Int_t &numChannels);
    // parse string describing reaction
    PParticle *MakeParticle(char * name);

    Int_t nMaxBytes;

    static FILE *asciif;
    static Int_t globalEventCounter;

    Int_t weight_reset;
    Int_t is_inline;

    // beam axis skewing and smearing
    /* Double_t thetaBeam, phiBeam, sigmaBeam; */
    /* Double_t beam_energy,beam_momentum; //for consecutive reactions */
    const char *r_beam, *r_target;

    //From batch system
    Double_t *vertex_x, *vertex_y, *vertex_z;
    Double_t *event_impact_param , *event_plane;
    Double_t *weight_version;

    // ascii file

    static void passEvent(float,  float,  float,  int, int*, int*, int*, int*,
			  float*, float*, float*, float*,
			  float*, float*, float*, float*) {;}
    // dummy HGeant interface.  Comment to make HGeant routine visible!


    Int_t pre_heating;
    PReaction *sub_reaction;

    ClassDef(PReaction, 0) //Pluto Reaction Class

};

#endif // _PREACTION_H_














	
