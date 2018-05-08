
#include "../src/PChannelModel.h"


//This is the class definition
//Calculation by F. Iachello
//Provided by B. Ramstein


class PDeltaDalitzFF : public PChannelModel  {
  
 public:

    using PDistribution::GetWeight;
    PDeltaDalitzFF(Char_t *id, Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);
    
    void SetQED(Int_t i){useQED=i;};
    void SetCC(Double_t a, Double_t b, Double_t c) {
	g_m2=a*a;
	g_e2=b*b;
	g_c2=c*c;
	rescale = 1;
    };
    

 private:

    Double_t facteur_rho(Double_t x);
    Double_t gk2(Double_t t);

    Double_t beta_prime;// rho meson coupling strength
    Double_t beta;// size parameter for N-Delta transition (fitted using space-like data)
    
    Double_t a2; // phase for space-like to time-like continuation (fitted to elastic time-like
    // form factor data) 

    //Masses:
    Double_t mass_pi,fm_pis ,mass_rho,gamma_rho,mass_delta,mass_n,mu_p;

    Double_t theta;

    //Coupling constants for QED version
    Double_t g_m2, g_e2, g_c2;
    Int_t    useQED, rescale;
 

    ClassDef(PDeltaDalitzFF,0)  //Iachello form factor
};

PDistribution* PDeltaDalitzFF::Clone(const char*) const {
    //clone the object
    return new PDeltaDalitzFF((const PDeltaDalitzFF &)* this);
};

PDeltaDalitzFF::PDeltaDalitzFF(Char_t *id, Char_t *de, Int_t key) : PChannelModel(id, de,key) {
    //Constructor

    beta_prime = 0.004;
    beta       = 1.2748;
    a2         = 0.29;
    theta      = 53.;

    mass_pi    = makeStaticData()->GetParticleMass("pi+");
    fm_pis     = 4*mass_pi*mass_pi;
    mass_rho   = makeStaticData()->GetParticleMass("rho0");
    gamma_rho  = makeStaticData()->GetParticleTotalWidth("rho0");
    mass_delta = makeStaticData()->GetParticleMass("D+");
    mass_n     = makeStaticData()->GetParticleMass("n");

    mu_p = 2.793;
    g_m2 = 3.2*3.2;
    g_e2 = 0.04*0.04;
    g_c2 = 0.2 *0.2;

    useQED =0;
    rescale = 0;
} ;

Double_t PDeltaDalitzFF::gk2(Double_t t) {
    // intrinsic form_factor 
    Double_t denom2=1 + a2*a2*t*t + 2.*a2*t*cos(theta*TMath::Pi()/180.);
//cout <<"t= "<<t<< " g= "<<1-1./denom2<<endl;
    return 1/denom2;
}

Double_t PDeltaDalitzFF::facteur_rho(Double_t x) {
    // rho propagator term 
    if (x < fm_pis) return mass_rho*mass_rho
			/(mass_rho*mass_rho - x);
    //scaling to avoid "step"
    //  Double_t scale=0.91;
    Double_t scale=1.;

    Double_t logpart = (sqrt(x - fm_pis ) + sqrt(x))/(2*mass_pi);
    Double_t alpha= 2/TMath::Pi()*sqrt((x - fm_pis )/x)*log(logpart);
    Double_t denom_reel = mass_rho*mass_rho - x + (fm_pis -x)*gamma_rho*alpha/mass_pi;
    Double_t qnorm=x/fm_pis;
    Double_t beta = sqrt(pow(qnorm-1.,3)/qnorm);
    double denom_im = gamma_rho*4*mass_pi*beta;
    return scale*(mass_rho*mass_rho + 
	    8*gamma_rho*mass_pi/TMath::Pi())
	/sqrt(denom_reel*denom_reel+denom_im*denom_im);
}

Double_t PDeltaDalitzFF::GetWeight(Double_t *mass, Int_t *) {

    Double_t x=mass[0]*mass[0];
    Double_t m_d2=mass[1]*mass[1]; //mass of Delta sq.

    //q-dependent QED-Formfactor
    if (useQED)
	return (g_m2 + 3 * g_e2 + (g_c2*x/(2*m_d2)));
    

    //VMD-Formfactor
    Double_t ff = 2.4637*gk2(-x)*(beta_prime + beta*facteur_rho(x));
    if (rescale) ff *= sqrt((4.76826/4.68927)*(g_m2/(3.2*3.2)));
    return ff*ff;

}


ClassImp(PDeltaDalitzFF)


