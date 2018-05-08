#ifndef __PHUrParticle_h__
#define __PHUrParticle_h__


#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"


#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;


class PHUrParticle : public TLorentzVector {

public:
    Int_t ind;                //   1      ind : index of particle
    Double_t t;               //   2      t   : computational frame time of particle in fm/c
    TVector3 fr;              //          x,y,z coordinate in fm
    Double_t mass;            //   10     m   : mass of particle in GeV
    Int_t id;                 //   11     ityp: particle-ID
    Int_t I3;                 //   12   2*I3  : isospin z-projection (doubled)
    Int_t chrg;               //   13     ch  : charge of particle
    Int_t ind_part;           //   14         : index of last collision partner
    Int_t ncoll;              //   15    Ncoll: number of collisions
    Int_t s;                  //   16     S   : strangeness
    Int_t parent_process;     //   17         : history information (parent process)
    Int_t pdg;                //   ==> compact id for pdg conversion

    Int_t instance;           // default -1 , 0 first n last appearance of the particle
    Int_t first;              // first apperance ==1
    Int_t last;               // last apearance

    vector<Int_t> listInCollIndex;     // list of collisions involing this particle
    vector<Int_t> listOutCollIndex;

    using TObject::Read;
    Bool_t Read(ifstream& in) {
	if(in.eof())   return kFALSE;
	if(!in.good()) return kFALSE;
        Double_t px, py, pz, E;
	Double_t rx, ry, rz;
	in >> ind >> t >> rx >> ry >> rz >> E >> px >> py >> pz >> mass >> id >> I3 >> 
	    chrg >> ind_part >> ncoll >> s >> parent_process;
	ind      -= 1;  // fortran ->C++
	ind_part -= 1;  // fortran ->C++

	if (id >= 0)
	    pdg =  1000 * (chrg + 2) + id;
	else
            pdg = -1000 * (chrg + 2) + id;

	SetPxPyPzE(px, py, pz, E);
	fr.SetXYZ(rx, ry, rz);

	if(!in.eof() && !in.good()) return kFALSE;
        return kTRUE;
    }

    Int_t GetFirstCollision() {
        Int_t ind = -1;
	if(listInCollIndex.size()  > 0) ind = listInCollIndex[0];
        if(listOutCollIndex.size() > 0 && (listOutCollIndex[0] < ind || ind < 0)) ind = listOutCollIndex[0];
	return ind;
    }

    Int_t GetFirstCollisionInput() {
        Int_t ind = -1;
	if(listInCollIndex.size() > 0) ind = listInCollIndex[0];
	return ind;
    }

    Int_t GetFirstCollisionOutput() {
        Int_t ind = -1;
        if(listOutCollIndex.size() > 0 ) ind = listOutCollIndex[0];
	return ind;
    }

    Int_t GetLastCollision() {
        Int_t ind = -1;
	if(listInCollIndex.size()  > 0) { ind = listInCollIndex[listInCollIndex.size()-1]; }
        if(listOutCollIndex.size() > 0 && listOutCollIndex[listOutCollIndex.size()-1] > ind) ind = listOutCollIndex[listOutCollIndex.size()-1];
	return ind;
    }

    Int_t GetLastCollisionInput() {
        Int_t ind = -1;
	if(listInCollIndex.size()  > 0) { ind = listInCollIndex[listInCollIndex.size()-1]; }
	return ind;
    }

    Int_t GetLastCollisionOutput() {
        Int_t ind = -1;
        if(listOutCollIndex.size() > 0 ) ind = listOutCollIndex[listOutCollIndex.size()-1];
	return ind;
    }

    Bool_t IsPrimary() {
        Int_t ind1 = -1;
        Int_t ind2 = -1;
	if(listInCollIndex.size()  > 0 ) ind1 = listInCollIndex [0];
        if(listOutCollIndex.size() > 0 ) ind2 = listOutCollIndex[0];

        if(ind1 ==-1 && ind2 ==-1)               return kTRUE;  // no collision
        if(ind1 !=-1 && ind2 !=-1 && ind1<=ind2) return kTRUE;  // first appearance in input
        if(ind1 !=-1 && ind2 ==-1 )              return kTRUE;  // only in input

        return kFALSE;
    }

    using TObject::Print;
    void Print(Option_t *) {
	cout << setw(5)  << dec << ind << " "
	     << setw(15) << scientific
	     << setw(15) << fr.X() << " "
	     << setw(15) << fr.Y() << " "
	     << setw(15) << fr.Z() << " "
	     << setw(15) << E()    << " "
	     << setw(15) << Px()   << " "
	     << setw(15) << Py()   << " "
	     << setw(15) << Pz()   << " "
	     << setw(15) << mass   << " "
	     << setw(5)  << dec << id << " "
	     << setw(5)  << I3   << " "
	     << setw(5)  << chrg << " "
	     << setw(5)  << ind_part << " "
	     << setw(5)  << ncoll << " "
	     << setw(5)  << s << " "
	     << setw(5)  << parent_process
	     << dec << endl;
    }

    using TObject::Clear;
    void  Clear(Option_t *) {
	ind            = -1;
	t              = -1;
	SetPxPyPzE(-1.,-1.,-1.,-1.);
        fr.SetXYZ(-1.,-1.,-1.);
	mass           = -1;
	id             = -1;
	I3             = -1;
	chrg           = -1;
	ind_part       = -1;
	ncoll          = -1;
	s              = -1;
	parent_process = -1;
	pdg            = -1;
	instance       = -1;
	first          = -1;
	last           = -1;

	listInCollIndex .clear();
        listOutCollIndex.clear();
    }

    PHUrParticle(){ Clear(NULL); }
    ~PHUrParticle(){}

    ClassDef(PHUrParticle, 0)
};


#endif
