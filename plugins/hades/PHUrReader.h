#ifndef __PHUrReader_h__
#define __PHUrReader_h__


#include "PHUrDilep.h"
#include "PHUrEventHeader.h"
#include "PHUrCollisionHeader.h"
#include "PHUrParticle.h"
#include "PHUrAddon.h"

#include "TString.h"
#include "TClonesArray.h"

#include "PBulkInterface.h"


#include <vector>
#include <map>
#include <fstream>

using namespace std;


class PHUrReader : public PBulkInterface {


private:


    PHUrDilep dilep;       // urqmd dilepton shining code

    //-------------------------------------------------
    // ascii io
    TString     fInputName;
    TString     fOutputName;
    ifstream    fInputAscii;
    FILE       *fOutputAscii;
    //-------------------------------------------------

    //-------------------------------------------------
    // urqmd structure f15
    PHUrEventHeader     evtheader;                //! object to keep eventheader info
    PHUrEventHeader     evtheaderCP;              //! object to keep eventheader info
    PHUrCollisionHeader colheader;                //! object to keep collisonheader info
    vector< vector<PHUrParticle> > collisionsIn;  //! list of in-going particles  of collisions per event (entry i (collision) in vector is a list of particles)
    vector< vector<PHUrParticle> > collisionsOut; //! list of out-going particles of collisions per event (entry i (collision) in vector is a list of particles)
    vector< PHUrCollisionHeader  > collisions;    //! list of collision header  per event (entry i (collision) in vector is a collisionheader)
    //-------------------------------------------------


    //-------------------------------------------------
    // event handling
    Int_t fNumPart;
    Int_t fNumPartAddon;
    Int_t fNumMax;
    static TClonesArray *fAdd;
    static TClonesArray *fEvent;        // should be set to PLUTO TREE
    Int_t fEvtCt;                //!
    //-------------------------------------------------

    //-------------------------------------------------
    // particle ids etc
    map    <Int_t,Int_t>   mUrqmdToPdg;
    map    <Int_t,TString> mUrqmdProcess;
    map    <Int_t,Int_t> mapIDs;
    vector <Int_t > funknownIDs;
    map<Int_t , Int_t> fmapIds[3];
    Int_t pdg_param ;
    Int_t pid_param ;
    //-------------------------------------------------


    //-------------------------------------------------
    // flags
    Bool_t printEvtHeader; //!                     // urqmd reader
    Bool_t printColHeader; //!                     // urqmd reader
    Bool_t printParticle ; //!                     // urqmd reader
    Int_t  outputNonDiLeptons; //!                 // default: 0 ,1=put all urqmd stable particles to output 2=put all urqmd particle to output
    Bool_t outputLeptons;      //!                 // default: kTRUE,  let dileptons decay in urqmddilep shining code
    Bool_t outputFreezeoutDiLeptons; //!           // default: kFALSE, take into accout only final stage dileptons
    //-------------------------------------------------

    void          ClearVector(vector<vector<PHUrParticle> >& v);
    Bool_t        ReadEvent();
    const Char_t *GetUrQMDProcess(Int_t);
public:

    PHUrReader(TString inputfile="");
    ~PHUrReader();

    void         Input (TString filename);
    void         Output(TString filename);
    void         SetOutputNonDileptons      (Int_t doit)        { outputNonDiLeptons       = doit;}
    void         SetOutputLeptons           (Bool_t doit=kTRUE) { outputLeptons            = doit;}
    void         SetOutputFreezeoutDileptons(Bool_t doit=kTRUE) { outputFreezeoutDiLeptons = doit;}
    void         SetMaxNumParticles(Int_t max) { fNumMax = max; }
    void         SetMapID(Int_t Id_in,Int_t Id_out) ;

    virtual bool Modify(PParticle **array, int *decay_done, int *num, int maxnum);  //Modify particle bulk
    void         CreateAddonArray(TTree *T, Int_t stacksize);
    PHUrAddon   *CreateAddon(Int_t& index);
    void         CreateParticleArray(TTree *T, Int_t stacksize);
    PParticle   *CreateParticle(Int_t &index);
    Int_t        GetNFreeSlots()    { return  fNumMax - fEvent->GetEntries() ;}
    static TClonesArray *getAddon() { return fAdd  ;}
    static TClonesArray *getEvent() { return fEvent;}

    static Int_t UserAnalysis(PParticle**, Int_t);

    ClassDef(PHUrReader, 0)
};


#endif
