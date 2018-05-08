// Author: Ingo Froehlich
// Written: 14/05/2007
// Modified: 
// PFileOutput Class Header

#ifndef _PFILEOUTPUT_H_
#define _PFILEOUTPUT_H_

#include "PParticle.h"
#include "PBulkInterface.h"
#include "PChannel.h"

class PFileOutput: public PBulkInterface {

 private:

 protected:

     Int_t cnt;              //! number of total particles in event
     Int_t subcnt;           //! number of particles in branch
     Int_t getVERTEX;        //! transport getVERTEX switch from PReaction
     Int_t writeINDEX;       //! transport writeINDEX switch from PReaction
     Int_t allPARTICLES;     //! transport allPARTICLES switch from PReaction
     Int_t asciiOUTPUT;      //! transport asciiOUTPUT switch from PReaction
     PChannel **channel;     //! transport channel array from PReaction
     Int_t      nChannel;    //! transport number of channels in array

     Int_t branchNum;        //! Branch number, 0:std branch, >0:additional branches
     const char *branchName; //! Branch name

     const char *filename;

  public:
    PFileOutput();

    void SetHeader(Int_t _cnt, 
		   Int_t _allPARTICLES,
		   Int_t _getVERTEX,
		   Int_t _asciiOUTPUT,
		   Int_t _writeINDEX,
		   PChannel **_channel,
                   Int_t _nCahnnel
		  ) {
	//number of expected WriteParticles:
	cnt          = _cnt;
	//PReaction flags:
	allPARTICLES = _allPARTICLES;
	getVERTEX    = _getVERTEX;
	asciiOUTPUT  = _asciiOUTPUT;
	writeINDEX   = _writeINDEX;
	channel      = _channel;
        nChannel     = _nCahnnel;
    };

    void SetBranchHeader(Int_t _subcnt, Int_t _branchNum, const char *_branchName) {
	subcnt     = _subcnt;
	branchNum  = _branchNum;
	branchName = _branchName;
    };

    virtual bool OpenFile(const char *_filename);      //filename
    virtual bool CloseFile(void);                //
    virtual bool WriteEvent(void);               //next event
    virtual bool WriteEventHeader(void);
    virtual bool WriteBranchHeader(void);
    virtual bool WriteParticle(PParticle *par);  //write single particle

    const char  *GetFilename() {
	return filename;
    };

    ClassDef(PFileOutput, 0) // Pluto file output base class
};
#endif 

















