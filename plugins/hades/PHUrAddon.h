#ifndef __PHUrAddon_h__
#define __PHUrAddon_h__


#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"


#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;


class PHUrAddon : public TObject {

public:
    Int_t    n_in;       //  n ingoing particles
    Int_t    n_out;      //  n outgoing particles
    Int_t    process_id; //  process id
    Int_t    n_collision;//  n collision
    Double_t t_cre;     //   t   : computational frame time of particle in fm/c
    Double_t t_abs;     //   t   : computational frame time of particle in fm/c
    Double_t dens_cre;  //   baryon density at collision point
    Double_t dens_abs;  //   baryon density at collision point
    TVector3 fr_cre;    //   x,y,z coordinate in fm
    TVector3 fr_abs;    //   x,y,z coordinate in fm
    Int_t    stable;    // flag for stable particles
    Int_t    instance;  // default -1 , 0 first n last appearance of the particle
    Int_t    first;     // default -1 , 0 first n last appearance of the particle
    Int_t    last;      // default -1 , 0 first n last appearance of the particle

    using TObject::Clear;
    void  Clear(Option_t *) {
	n_in           = -2;
	n_out          = -1;
	process_id     = -1;
	n_collision    = -1;
	t_cre          = -1;
	t_abs          = -1;
	dens_cre       = -1;
        dens_abs       = -1;
	first          = -1;
        last           = -1;

        fr_cre.SetXYZ(-1.,-1.,-1.);
	fr_abs.SetXYZ(-1.,-1.,-1.);
	stable = -1;
        instance=-1;
    }

    PHUrAddon(){ Clear(NULL); }
    ~PHUrAddon(){}

    Bool_t IsScattering()            { return (n_in==2 && n_out==2 ) ; }
    Bool_t IsDecay()                 { return (n_in==1 && n_out==2 ) ; }
    Bool_t IsAnnihilation()          { return (n_in==2 && n_out==1 ) ; }
    Bool_t IsPauliBlockedCollision() { return (n_in==2 && n_out==0 ) ; }
    Bool_t IsPauliBlockedDecay()     { return (n_in==1 && n_out==0 ) ; }
    Bool_t IsStringDecay5()          { return (n_in==2 && n_out==5 ) ; }

    ClassDef(PHUrAddon, 1)
};


#endif
