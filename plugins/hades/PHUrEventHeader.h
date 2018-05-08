#ifndef __PHUrEventHeader_h__
#define __PHUrEventHeader_h__

#include "PHUrCollisionHeader.h"

#include "TObject.h"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;


class PHUrEventHeader : public TObject {  
public:

    Int_t    id;               // 1 : Int_t    The first integer in the event header is a -1.
    Int_t    evtnum;           // 2 : Int_t    event number,
    Int_t    mass_proj;        // 3 : Int_t    mass of projectile
    Int_t    mass_target;      // 4 : Int_t    mass of target,
    Double_t impact_par;       // 5 : Double_t impact parameter,
    Double_t E_reaction;       // 6 : Double_t two-particle c.m. energy of heavy-ion reaction,
    Double_t tot_sig;          // 7 : Double_t the total cross section of the heavy-ion reaction ,
    Double_t E_beam;           // 8 : Double_t the beam energy ,
    Double_t mom_per_particle; // 9 : Double_t momentum (per particle) in the laboratory frame

    using TObject::Read;
    Bool_t Read(ifstream& in){
	if(in.eof()) return kFALSE;
	if(!in.good()) return kFALSE;
	in >> id >> evtnum >> mass_proj >> mass_target >> impact_par >> 
	    E_reaction >> tot_sig >> E_beam >> mom_per_particle;
	if(evtnum!=-2 && !in.eof() && !in.good()) 
	    return kFALSE;
        return kTRUE;
    }

    using TObject::Print;
    void Print(Option_t *) {
	cout << setw(5)  << dec << id  << " "
	     << setw(6)  << evtnum     << " "
	     << setw(4)  << mass_proj  << " "
	     << setw(4)  << mass_target << " "
	     << setw(15) << scientific << impact_par << " "
	     << setw(15) << E_reaction << " "
	     << setw(15) << tot_sig    << " "
	     << setw(15) << E_beam     << " "
	     << setw(15) << mom_per_particle
	     << dec << endl;
    }

    using TObject::Copy;
    void Copy(PHUrCollisionHeader& colheader){
	id               = colheader.n_in;
	evtnum           = colheader.n_out;
	mass_proj        = colheader.process_id;
	mass_target      = colheader.n_collision;
	impact_par       = colheader.t_collision;
	E_reaction       = colheader.E_total_CM;
	tot_sig          = colheader.sig_total;
	E_beam           = colheader.sig_partial;
	mom_per_particle = colheader.baryon_density;
    }

    using TObject::Clear;
    void Clear(Option_t *) {
	id               = -2;
	evtnum           = -1;
	mass_proj        = -99;
	mass_target      = -99;
	impact_par       = -1.;
	E_reaction       = -1.;
	tot_sig          = -1.;
	E_beam           = -1.;
	mom_per_particle = -1.;
    }

    PHUrEventHeader(){ Clear(NULL); };
    ~PHUrEventHeader(){};

    ClassDef(PHUrEventHeader, 0)
};


#endif
