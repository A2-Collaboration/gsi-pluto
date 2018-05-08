{{
//
// Thermal source model of Ag+Ag at 1.65 AGeV 30% central 

//makeDistributionManager()->Disable("helicity_angles");


Float_t Eb    = 1.65;     // beam energy in AGeV
Float_t ctr   = 1.0;     // 
//Float_t bmax  = 3.9;     // max. impact parameter (corresponds to 60% of all)
Float_t bmax  =0.0;         // Poisson sampling
if (bmax>0.) ctr = 1.;   // use b sampling instead of Poisson sampling

Float_t T1    = 0.052;   // temperature in GeV (for pion 2-component spectra)
Float_t T2    = 0.095;   // temperature in GeV (assume this for thermalized source)
Float_t frac  = 0.7;     // fraction of pion low-T component (from Jehad's fit to QMD)
Float_t blast = 0.0;     // radial expansion velocity

Float_t A2    = 0.3;     // polar distribution
Float_t A4    = 0.0;
Float_t v1    = 0.0;     // side flow
Float_t v2    = 0.0;     // elliptic flow

Float_t Mpi0   = 13.41*ctr;      // pi0 multiplicity in p+C 3.5 AGeV
Float_t Mpip   = 13.41*ctr;      // pi+
Float_t Mpim   = 13.07*ctr;      // pi-
Float_t Meta   = 0.36*ctr;     // eta
Float_t Momega = 0.045*ctr;     // omega
Float_t Mrho0  = 0.52*ctr;     // rho0  
Float_t Mphi   = 0.0017*ctr;    // phi from ArKcl data
Float_t MDelta = 2.*3/2.*Mpi0; // Delta0 + Delta+ (2? enhanced for in-medium   

Float_t enhance = 20.;           // enhancement factor
Float_t enhancepi=1;            //enhancement factor for pi

Mpi0 *=enhancepi;
Meta *=enhance;
Momega *= enhance;
Mrho0 *= enhance;
Mphi *= enhance;
MDelta *= enhance;

PFireball *source1=new PFireball("pi0",Eb,T1,T2,frac,blast,A2,A4,v1,v2);
PFireball *source2=new PFireball("eta",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);
PFireball *source3=new PFireball("w",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);
PFireball *source4=new PFireball("w",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);
PFireball *source5=new PFireball("rho0",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);
PFireball *source6=new PFireball("phi",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);
PFireball *source8=new PFireball("D0",Eb,T2,0.,1.,blast,0.0,A4,v1,v2);

source1->SetW(1.0/enhancepi);
source2->SetW(1.0/enhance);
source3->SetW(1.0/enhance);
source4->SetW(1.0/enhance);
source5->SetW(1.0/enhance);
source6->SetW(1.0/enhance);
source8->SetW(1.0/enhance);

PChannel *c1 = source1->makeChannel(40, Mpi0);          // pi0 decay directly

PChannel *c2 = source2->makeChannel(1, 0.006*Meta);     // eta -> gamma e+e-
PChannel *c3 = source3->makeChannel(1, 5.9e-4*Momega);  // omega -> pi0 e+e-
PChannel *c4 = source4->makeChannel(1, 7.15e-5*Momega); // omega -> e+e-
PChannel *c5 = source5->makeChannel(1, 4.48e-5*Mrho0);  // rho0 -> e+e-
PChannel *c6 = source6->makeChannel(1, 3.00e-4*Mphi);   // phi -> e+e-
PChannel *c8 = source8->makeChannel(1, 4.0e-5*MDelta);  // Delta -> N e+e-

source1->Print();
source2->Print();
source3->Print();
source4->Print();
source5->Print();
source8->Print();

PParticle** ptcls;

/// Dalitz decays

PParticle* gam2  = new PParticle("g");
PParticle* diel2 = new PParticle("dilepton");
PParticle* elec2 = new PParticle("e-");
PParticle* posi2 = new PParticle("e+");


////////////////    eta Dalitz decay       /////////////////////////////

ptcls = c2->GetParticles();
PParticle* s2_1[] = {ptcls[1], gam2, diel2};
PParticle* s2_2[] = {diel2, elec2, posi2};
PChannel* c2_1 = new PChannel(s2_1,2,1);
PChannel* c2_2 = new PChannel(s2_2,2,1);


////////////////    omega Dalitz decay     /////////////////////////////

PParticle* pi03  = new PParticle("pi0");
PParticle* diel3 = new PParticle("dilepton");
PParticle* elec3 = new PParticle("e-");
PParticle* posi3 = new PParticle("e+");

ptcls = c3->GetParticles();
PParticle* s3_1[] = {ptcls[1], pi03, diel3};
PParticle* s3_2[] = {diel3, elec3, posi3};
PChannel* c3_1 = new PChannel(s3_1,2,1);
PChannel* c3_2 = new PChannel(s3_2,2,1);

////////////////    omega direct decay     /////////////////////////////

PParticle* elec4 = new PParticle("e-");
PParticle* posi4 = new PParticle("e+");

ptcls = c4->GetParticles();
PParticle* s4_1[] = {ptcls[1], elec4, posi4};
PChannel* c4_1 = new PChannel(s4_1,2,1);

////////////////    rho0 direct decay     /////////////////////////////

PParticle* elec5 = new PParticle("e-");
PParticle* posi5 = new PParticle("e+");

ptcls = c5->GetParticles();
PParticle* s5_1[] = {ptcls[1], elec5, posi5};
PChannel* c5_1 = new PChannel(s5_1,2,1);


////////////////    phi direct decay      /////////////////////////////

PParticle* elec6 = new PParticle("e-");
PParticle* posi6 = new PParticle("e+");

ptcls = c6->GetParticles();
PParticle* s6_1[] = {ptcls[1], elec6, posi6};
PChannel* c6_1 = new PChannel(s6_1,2,1);

////////////////    Delta Dalitz decay    /////////////////////////////

PParticle* neut8 = new PParticle("n");
PParticle* diel8 = new PParticle("dilepton");
PParticle* elec8 = new PParticle("e-");
PParticle* posi8 = new PParticle("e+");

ptcls = c8->GetParticles();
PParticle* s8_1[] = {ptcls[1], neut8, diel8};
PParticle* s8_2[] = {diel8, elec8, posi8};
PChannel* c8_1 = new PChannel(s8_1,2,1);
PChannel* c8_2 = new PChannel(s8_2,2,1);

PChannel *cc[] = {
    c1,                  /* pi0          */
    c2,c2_1,c2_2,       /* eta Dalitz   */
    c3,c3_1,c3_2,       /* omega Dalitz */
    c4,c4_1,            /* omega direct */
    c5,c5_1,            /* rho0         */
    c6,c6_1,            /* phi direct   */
    c8,c8_1,c8_2         /* Delta Dalitz */
};

 
PReaction *r=new PReaction(cc,"delme",14,0,0,0,0);


r->allParticles();
r->SetDecayAll(1.);   // decay all particles with tau < 1 ns
r->SetHGeant(0);      // set to 1, if PLUTO is run from HGeant prompt

 TH2F * histo2 = new TH2F ("histo2","pi0 pt vs. rapitidy",20,-1,2.,20,0,1);
 r->Do(histo2,"foreach(pi0); _x = [pi0]->Rapidity(); _y = [pi0]->Pt(); ");
 
 TH1F * histo1 = new TH1F ("histo1","dilepton mass",100,0,1.2);
 r->Do(histo1,"foreach(dilepton); _x = [dilepton]->M();");

TH1F * histo1a = new TH1F ("histo1a","ee inv. mass",100,0,1.2);
 r->Do(histo1a,"do_ee: _x = ([e+] + [e-])->M();");
 r->Do("formore(e+); formore(e-); goto do_ee");


TH1F * histo3 = new TH1F ("histo3","photon energy",20,0,1.);
r->Do(histo3,"foreach(g); _x = [g]->E();");

 TH2F * histo2a = new TH2F ("histo2a","dilepton pt vs. rapitidy",20,-1,2.,20,0,1);
 r->Do(histo2a,"foreach(dilepton); _x = [dilepton]->Rapidity(); _y = [dilepton]->Pt(); ");

r->Print();

r->loop(100000);   // do n events

 TCanvas *can1 = new 
     TCanvas ("can1","pi0 Dalitz decay ee mass");
 can1->SetLogy(1);
 
 histo1a->Draw();
 histo1->Draw("e1same");

 TCanvas *can2 = new 
     TCanvas ("can2","pi0");
 histo2->Draw("");

 TCanvas *can2a = new 
     TCanvas ("can2a","dilepton");
  histo2a->Draw("");
 
 TCanvas *can3 = new 
     TCanvas ("can3","gamma E");
 histo3->Draw();

  

return 0;
}

