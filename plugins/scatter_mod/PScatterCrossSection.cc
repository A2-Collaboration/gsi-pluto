/////////////////////////////////////////////////////////////////////
//
//
//                                  Author:  I. Froehlich
/////////////////////////////////////////////////////////////////////


using namespace std;
#include <sstream>
#include <iostream>
#include <iomanip>

#include "PScatterCrossSection.h"
#include "PDynamicData.h"

ClassImp(PScatterCrossSection)

PScatterCrossSection::PScatterCrossSection() {
    Fatal("PScatterCrossSection", "Wrong ctor called");
};

PScatterCrossSection::PScatterCrossSection(const Char_t *id, const Char_t *de) :
    PAngularDistribution(id, de) {

    beam = target = NULL;
    pf2  = NULL;
    qmax = 0.0;
    qmin = -1;
    npx = npy = -1;

    TF1 *dummy = new TF1("dummy", "1", -1, 1);
    SetAngleFunction(dummy);
 
};

PDistribution *PScatterCrossSection::Clone(const char *) const {
    return new PScatterCrossSection((const PScatterCrossSection &)* this);
};

void PScatterCrossSection::SetNpx(Int_t my_npx) {
    npx = my_npx;
    if (pf2) pf2->SetNpx(npx);
};

void PScatterCrossSection::SetNpy(Int_t my_npy) {
    npy = my_npy;
    if (pf2) pf2->SetNpy(npy);
};
  
Bool_t PScatterCrossSection::MakeVars() {
    if (qmin < 0) {
	Error("MakeVars", "Range of the c.m. energy not defined");
	return kFALSE;
    }
    if (!pf2) {
	pf2 = new PF2("", -1., 1., qmin, qmax);
	if (npx > 0) pf2->SetNpx(npx);
	if (npy > 0) pf2->SetNpy(npy);
    }

    vprimary   = makeDynamicData()->GetBatchParticle("_primary");  
    vf         = makeStaticData()->GetBatchValue("_f"); 
    return kTRUE;
};

Bool_t PScatterCrossSection::AddEquation(char *command) {
    if (!MakeVars()) return kFALSE;
    return pf2->Add(command);
};

Bool_t PScatterCrossSection::AddHistogram(TH1 *histo, char *command) {    
    if (!MakeVars()) return kFALSE;
    return pf2->Add(histo, command);
};

Bool_t PScatterCrossSection::AddHistogram(TH2 *histo, char *command) {
    if (!MakeVars()) return kFALSE;
    return pf2->Add(histo, command);
};

Bool_t PScatterCrossSection::Init(void) {
    if (!PAngularDistribution:: Init()) return kFALSE;

    if (!beam || !target) {
	Warning("Init", "beam or target not found");
	return kFALSE;
    }

    return kTRUE;    
};

Bool_t PScatterCrossSection::Prepare(void) {

    PAngularDistribution::Prepare();

    if (!pf2) {
	Error("Prepare", "Function not defined");
	return kFALSE;
    }

    pf2->GetRandom2(costheta, q);

    Double_t m_t = target->M();
    Double_t m_b = beam->M();

    //momentum fixed target:
    Double_t mom_b = sqrt( pow(  (q*q-m_b*m_b-m_t*m_t)/(2*m_t)  ,2) - m_b*m_b );

    //cout << costheta << ":" << q << ":" << mom_b << endl;

    if (std::isnan(mom_b)) {
	return kFALSE;
    }

    //return kTRUE;
    //Make a copy of beam and target:
    PParticle my_beam(beam);
    PParticle my_target(target);

    my_beam.Boost(-my_target.BoostVector());
    my_beam.SetMom(mom_b);
    my_beam.Boost(my_target.BoostVector());

    beam->SetVect(my_beam.Vect());
    beam->SetM(m_b);

    parent->Reconstruct();

    return kTRUE;
};



double PScatterCrossSection::SamplePolarAngle(double) {
    return costheta;
}
