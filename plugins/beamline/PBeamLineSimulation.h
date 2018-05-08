
#ifndef _PBEAMLINESIMULATION_H_
#define _PBEAMLINESIMULATION_H_

#include "PBeamSmearing.h"
#include "PBulkInterface.h"
#include "PProjector.h"

#define MAX_NUM_BEAM_PARTICLES 1000


class PBeamLineSimulation : public PBeamSmearing  {
  
 public:
    PBeamLineSimulation();
    ~PBeamLineSimulation();
    PBeamLineSimulation(Char_t *id, Char_t *de);
    
    PDistribution *Clone(const char *delme=NULL) const;

    Bool_t Init(void);
    Bool_t Prepare(void);

    Bool_t InitBeamLine(char *filename);
    Int_t  AddDetector (char *name, Double_t distance);

    Bool_t TargetIsElement(Int_t n);
    void   SetGlobalMomentum(Double_t p) {global_p = p;};

    Bool_t AddEquation(char *command);
    Bool_t AddHistogram(TH2 *histo, const char *command = "");
    Bool_t Do(char *command) {return AddEquation(command);};
    Bool_t Do(TH2 *histo, const char *command = "") {return AddHistogram(histo, command);};
    
    void SetBeamtubeSize(Double_t _beamtube_size_x, Double_t _beamtube_size_y) {
	beamtube_size_x = _beamtube_size_x;
	beamtube_size_y = _beamtube_size_y;
    };
    
 private:

    PParticle *beam,   *mybeam;
    PParticle *target, *mytarget;
    PParticle *parent;

    PParticle *beam_particles[MAX_NUM_BEAM_PARTICLES];
    Int_t      beam_id       [MAX_NUM_BEAM_PARTICLES];
    Int_t      beam_counter;

    Double_t  *filter;
    Double_t   beamtube_size_x, beamtube_size_y; //size in mm

    char      *detector_name[MAX_NUM_BEAM_PARTICLES];
    Double_t   detector_distance[MAX_NUM_BEAM_PARTICLES];
    Int_t      detector_counter;
    
    Int_t AddBeamParticle(char *name);
    void MakeVars(void);

    PProjector *projector;

    //taken from the original code:
    long double T11[51],T12[51],T13[51],T14[51],T16[51];
    long double T111[51],T112[51],T113[51],T114[51],T116[51];
    long double T122[51],T123[51],T124[51],T126[51];
    long double T133[51],T134[51],T136[51];
    long double T144[51],T146[51];
    long double T166[51],along[51];

    long double T21[51],T22[51],T23[51],T24[51],T26[51];      
    long double T211[51],T212[51],T213[51],T214[51],T216[51];
    long double T222[51],T223[51],T224[51],T226[51];
    long double T233[51],T234[51],T236[51];
    long double T244[51],T246[51];
    long double T266[51];      

    long double T31[51],T32[51],T33[51],T34[51],T36[51];
    long double T311[51],T312[51],T313[51],T314[51],T316[51];
    long double T322[51],T323[51],T324[51],T326[51];
    long double T333[51],T334[51],T336[51];
    long double T344[51],T346[51];
    long double T366[51];

    long double T41[51],T42[51],T43[51],T44[51],T46[51]; 
    long double T411[51],T412[51],T413[51],T414[51],T416[51];
    long double T422[51],T423[51],T424[51],T426[51];
    long double T433[51],T434[51],T436[51];
    long double T444[51],T446[51];
    long double T466[51];
           
    Int_t nelem; 

    char txt2[3]; 
    char txt5[6]; 
    
    long double z1; 
    long double z2; 
    long double z3; 
    long double z4; 
    long double z5; 
    
    long double xi; 
    long double xpi;
    long double yi; 
    long double ypi; 
    long double dp; 
    
    Int_t iok; 
    long double xf; 
    long double yf; 
    long double xpf; 
    long double ypf; 
    long double xf1; 
    long double yf1; 
    long double xpf1; 
    long double ypf1; 
    long double xf2; 
    long double xpf2; 
    long double yf2; 
    long double ypf2; 

    long double aa;
    long double sig11; 
    long double sig22; 
    long double r12; 
    long double sig33; 
    long double r13; 
    long double r23; 
    long double sig44;
    long double r14; 
    long double r24; 
    long double r34; 
    long double sig55; 
    long double r15; 
    long double r25; 
    long double r35; 
    long double r45; 
    long double sig66; 
    long double r16; 
    long double r26; 
    long double r36; 
    long double r46; 
    long double r56; 

    Int_t ii; 
    Int_t ind1; 
    Int_t ind23; 

    Int_t num_target;

    Double_t *in_xi;
    Double_t *in_yi;
    Double_t *in_dp;
    Double_t *in_xpi;
    Double_t *in_ypi;

    Double_t out_xi[51];
    Double_t out_yi[51];
    Double_t out_dp[51];
    Double_t out_xpi[51];
    Double_t out_ypi[51];
    Double_t global_p;

    ClassDef(PBeamLineSimulation, 0)  // General purpose beam line transport simulation
};

#endif


