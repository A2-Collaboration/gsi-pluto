#ifndef __PHUrCollisionHeader_h__
#define __PHUrCollisionHeader_h__


#include "TObject.h"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

class PHUrCollisionHeader : public TObject {

public:
    Int_t    n_in;            //   1  Int_t n ingoing particles
    Int_t    n_out;           //   2  Int_t n outgoing particles
    Int_t    process_id;      //   3  Int_t process id
    Int_t    n_collision;     //   4  Int_t n collision
    Double_t t_collision;     //   5  Double_t collison time fm/c
    Double_t E_total_CM;      //   6  Double total CM energy GeV
    Double_t sig_total;       //   7  Double_t  total cross section mbarn
    Double_t sig_partial;     //   8  Double_t  partial cross section of exit channel mbarn
    Double_t baryon_density;  //   9  Double_t baryon density at collision point

    using TObject::Read;
    Bool_t Read(ifstream& in) {
	if(in.eof())   return kFALSE;
	if(!in.good()) return kFALSE;
	in >> n_in >> n_out >> process_id >> n_collision >> t_collision >> 
	    E_total_CM >> sig_total >> sig_partial >> baryon_density;
	if(n_in!=-2 && !in.eof() && !in.good()) return kFALSE;
        return kTRUE;
    }

    using TObject::Print;
    void Print(Option_t *) {
	cout<< setw(5)  << dec << n_in << " "
	    << setw(6)  << n_out       << " "
	    << setw(4)  << process_id  << " "
	    << setw(4)  << n_collision << " "
	    << setw(15) << scientific  << t_collision << " "
	    << setw(15) << E_total_CM  << " "
	    << setw(15) << sig_total   << " "
	    << setw(15) << sig_partial << " "
	    << setw(15) << baryon_density
	    << dec << endl;
    }

    using TObject::Clear;
    void Clear(Option_t *) {
	n_in           = -2;
	n_out          = -1;
	process_id     = -1;
	n_collision    = -1;
	t_collision    = -1;
	E_total_CM     = -1;
	sig_total      = -1;
	sig_partial    = -1;
	baryon_density = -1;
    }
    PHUrCollisionHeader(){ Clear(NULL);}
    ~PHUrCollisionHeader(){}
    Bool_t IsScattering()            { return (n_in==2 && n_out==2 ) ; }
    Bool_t IsDecay()                 { return (n_in==1 && n_out==2 ) ; }
    Bool_t IsAnnihilation()          { return (n_in==2 && n_out==1 ) ; }
    Bool_t IsPauliBlockedCollision() { return (n_in==2 && n_out==0 ) ; }
    Bool_t IsPauliBlockedDecay()     { return (n_in==1 && n_out==0 ) ; }
    Bool_t IsStringDecay5()          { return (n_in==2 && n_out==5 ) ; }

    ClassDef(PHUrCollisionHeader,0)

};


#endif
