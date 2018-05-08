#ifndef __PHUrDilep_h__
#define __PHUrDilep_h__


#include "TObject.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TTree.h"

#include "PParticle.h"
#include "PHUrParticle.h"
#include "PHUrAddon.h"
#include "PHUrEventHeader.h"
#include "PHUrCollisionHeader.h"


#include <vector>
#include <map>

using namespace std;

class PHUrReader;


class PHUrDilep : public TObject {

public:

    PHUrReader *reader; //!


    //---------------------------------------
    // ascii io
    FILE *out;
    //---------------------------------------


    //---------------------------------------
    // urqmd structure
    PHUrEventHeader      *evtheader;
    PHUrParticle         *particleIn;
    PHUrParticle         *particleOut;
    PHUrCollisionHeader  *collisionIn;
    PHUrCollisionHeader  *collisionOut;
    TLorentzVector        ep,em;
    //---------------------------------------

    //---------------------------------------
    // particle id etc
    Int_t omega    ;
    Int_t rho      ;
    Int_t phi      ;
    Int_t pion     ;
    Int_t etaprime ;
    Int_t eta      ;
    Int_t delta    ;
    Int_t omegadir ;
    Int_t omegadal ;

    map <Int_t,Int_t> *mUrqmdToPdg;
    //---------------------------------------

    //---------------------------------------
    // flags
    Bool_t outputLeptons;      //!                 // default: kTRUE,  let dileptons decay
    //---------------------------------------

private:


    //---------------------------------------
    // constants
    Double_t  alpha_em       ;
    Double_t  vacmass_pi0    ; // pi0
    Double_t  mass_electron  ; // e-
    Double_t  mass_nucleon   ; // proton
    Double_t  vacmass_eta    ; // eta
    Double_t  vacmass_etaprime; // eta'
    Double_t  vacmass_rho     ; // rho0
    Double_t  vacmass_omega   ; // omega
    Double_t  vacmass_phi     ; // phi
    Double_t  vacmass_delta   ; // D+

    Double_t  gev            ;
    Double_t  g              ;     // delta


    Double_t  lambda_omega   ;
    Double_t  gamma_omega    ;
    Double_t  gamma_photon   ;

    Double_t  br_pi0         ;
    Double_t  b_pi0          ;

    Double_t  lambda_eta     ;
    Double_t  br_eta         ;

    Double_t  lambda_etaprime ;
    Double_t  gamma_etaprime  ;
    Double_t  br_etaprime     ;
    //---------------------------------------

    void dalpi       (Double_t mx, Double_t& dgamma);
    void daleta      (Double_t mx, Double_t& dgamma);
    void daletaprime (Double_t mx, Double_t& dgamma);
    void daldelta (Double_t tau, Double_t mx, Double_t mres, Double_t& dgamma);
    void dalomega (Double_t tau, Double_t mx, Double_t &dgamma);
    void diromega (Double_t tau, Double_t mres, Int_t multi, Double_t& weight);
    void dirphi   (Double_t tau, Double_t mres, Int_t multi, Double_t& weight);
    void dirrho   (Double_t tau, Double_t mres, Int_t multi, Double_t& weight);
    void lobo_dal (Double_t p0_gstar, Double_t p0_particle,
		   Double_t m_gstar,  Double_t beta, Double_t gamma,
		   Double_t x_gstar,  Double_t y_gstar, Double_t z_gstar, Double_t t_gstar,
		   Double_t x_particle, Double_t y_particle, Double_t z_particle,
		   Double_t m_res);
    void lobo_dir(Double_t beta, Double_t gamma, Double_t m_res);

    void  dgamma_sum(Int_t ityp, Double_t tau, Int_t multi,
		     Double_t& weight,
		     Double_t& smin, Double_t& smax, Double_t& tmax);

    void dgamma_sum_delta(Double_t tau, Double_t mres, Int_t multi,
			  Double_t& weight_delta,
			  Double_t& smin_del, Double_t& smax_del, Double_t& tmax_del);

    void  gamma_star(Double_t smin, Double_t smax, Double_t tmax,
		     Double_t& s, Double_t& t,
		     Double_t& xgstar, Double_t& ygstar, Double_t& zgstar, Double_t& tgstar,
		     Double_t& xparticle, Double_t& yparticle, Double_t& zparticle, Double_t& tparticle);

    void t_delta(Double_t tau, Double_t mres, Int_t multi, Double_t& tmax_del);

    void t_omega(Double_t tau, Int_t multi, Double_t& smin_omega, Double_t& smax_omega, Double_t& tmax_omega);
    Double_t rndfunc();

    void Output(Int_t ityp, Double_t weight, Double_t mres, Int_t multi, Double_t tau);
    void OutputROOT(Int_t ityp, Double_t weight);


public:
    PHUrDilep();
    ~PHUrDilep();

    void Output(TString infilename, TString outfilename);
    void Dilep();


    ClassDef(PHUrDilep, 0)
};


#endif
