{
    PUtils::SetSeed(1282137294);
    makeStaticData()->SetParticleUMass("w",1);

//     Float_t Eb    = 1.65;     // beam energy in AGeV
//     Float_t ctr   = 1.0;     // 
//     //Float_t bmax  = 3.9;     // max. impact parameter (corresponds to 60% of all)
//     Float_t bmax  =0.0;         // Poisson sampling
//     if (bmax>0.) ctr = 1.;   // use b sampling instead of Poisson sampling
    
//     Float_t T1    = 0.052;   // temperature in GeV (for pion 2-component spectra)
//     Float_t T2    = 0.095;   // temperature in GeV (assume this for thermalized source)
//     Float_t frac  = 0.7;     // fraction of pion low-T component (from Jehad's fit to QMD)
//     Float_t blast = 0.0;     // radial expansion velocity
    
//     Float_t A2    = 0.3;     // polar distribution
//     Float_t A4    = 0.0;
//     Float_t v1    = 0.0;     // side flow
//     Float_t v2    = 0.0;     // elliptic flow
    
    //vom CBM-macro:
Float_t Eb   = 25;
//Float_t T1    = 0.130;   // temperature in GeV
Float_t T1    = 0.250;   // temperature in GeV
Float_t T2    = 0.;      // temperature in GeV
Float_t blast = 0.3;      // radial expansion velocity
Float_t w = 0.;      // elliptic flow


//PFireball *source4=new PFireball("w",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);
 PFireball *source4=new PFireball("w",Eb,T1,T2,1.,blast,0.0,0,0,0);
 //source4->SetEpsilon(0.1);
 //source4->setSigma(0.41);
 //ource4->SetW(1.0);

PChannel *c4 = source4->makeChannel(1, 0); // omega -> e+e-

source4->Print();

PParticle** ptcls;

////////////////    omega direct decay     /////////////////////////////

TH1F * histo1     = new TH1F ("histo1","w mass",100,0.1,1.2);

PParticle* elec4 = new PParticle("e-");
PParticle* posi4 = new PParticle("e+");

//   PParticle* elec4 = new PParticle("mu-");
//   PParticle* posi4 = new PParticle("mu+");

ptcls = c4->GetParticles();
PParticle* s4_1[] = {ptcls[1], elec4, posi4};
PChannel* c4_1 = new PChannel(s4_1,2,1);
// PChannel* c4_1 = NULL;

 PChannel *ccw[] = {c4,c4_1};

 PReaction *r=new PReaction(ccw,"delme",1);
r->Do(histo1,"_x= ([w])->M()");

 //r->allParticles();
 //r->SetDecayAll(1.);   // decay all particles with tau < 1 ns
r->SetHGeant(0);      // set to 1, if PLUTO is run from HGeant prompt




r->Print();

r->loop(20000);   // do n events


}

