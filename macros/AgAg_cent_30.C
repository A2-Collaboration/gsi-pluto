{
//
// Thermal source model of Ag+Ag at 1.65 AGeV 30% central 
//
// To save space, only leptons branches are generated and written to file.
// allows for different seeds each time the macro is started

    PUtils::SetSeed(123);

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

//Float_t enhance = 500.;           // enhancement factor for mesons
Float_t enhance = 100.;           // enhancement factor
//Float_t enhancepi=4.;            //enhancement factor for pi
Float_t enhancepi=1;            //enhancement factor for pi

Mpi0 *=enhancepi;
//Mpip *=enhancepi;
Meta *=enhance;
Momega *= enhance;
Mrho0 *= enhance;
Mphi *= enhance;
MDelta *= enhance;

PFireball *source1=new PFireball("pi0",Eb,T1,T2,frac,blast,A2,A4,v1,v2);
PFireball *source13=new PFireball("pi0",Eb,T1,T2,frac,blast,A2,A4,v1,v2);


source1->SetW(1.0/enhancepi);
source13->SetW(1.0/enhancepi);


PChannel *c1 = source1->makeChannel(20, Mpi0);      // pi0 decay directly
//PChannel *c1 = source1->makeChannel(1, Mpi0);


source1->Print();


PParticle** ptcls;

/// Dalitz decays


//  PParticle* gam1  = new PParticle("g");
//  PParticle* diel1 = new PParticle("dilepton");
//  PParticle* elec1 = new PParticle("e-");
//  PParticle* posi1 = new PParticle("e+");

////////////////    pi0 Dalitz decay       /////////////////////////////

//  ptcls = c1->GetParticles();
//  PParticle* s1_1[] = {ptcls[1], gam1, diel1};
//  PParticle* s1_2[] = {diel1, elec1, posi1};
 //PChannel* c1_1 = new PChannel(s1_1,2,1);
 //PChannel* c1_2 = new PChannel(s1_2,2,1);




PChannel *cc[] = {
    c1,
    //    c1_1, c1_2
                 };
//change  *cc[]if you add want to add new channels and modify param after outFile
 
PReaction *r=new PReaction(cc,"delme",1,0,0,0,0);
//PReaction *r=new PReaction(cc,"delme",1,0,0,0,0);
 r->allParticles();

 r->SetDecayAll(1.);   // decay all particles with tau < 1 ns
r->SetHGeant(0);      // set to 1, if PLUTO is run from HGeant prompt

//r->Do("echo ------------------");
//r->Do("foreach(*); pidx = [*]->GetParentIndex(); idx = [*]->GetIndex(); id = [*]->ID(); echo $id,$idx,$pidx"); 
 //r->Do("foreach(*); [*]->Print()");

 TH2F * histo2 = new TH2F ("histo2","pi0 pt vs. rapitidy",20,-1,2.,20,0,1);
 r->Do(histo2,"foreach(pi0); _x = [pi0]->Rapidity(); _y = [pi0]->Pt(); ");

TH1F * histo1 = new TH1F ("histo1","dilepton mass",20,0,0.2);
r->Do(histo1,"foreach(dilepton); _x = [dilepton]->M();");

TH1F * histo1a = new TH1F ("histo1a","ee inv. mass",20,0,0.2);
 r->Do(histo1a,"do_ee: _x = ([e+] + [e-])->M();");
 r->Do("formore(e+); formore(e-); goto do_ee");


TH1F * histo3 = new TH1F ("histo3","photon energy",20,0,1.);
r->Do(histo3,"foreach(g); _x = [g]->E();");

 TH2F * histo2a = new TH2F ("histo2a","dilepton pt vs. rapitidy",20,-1,2.,20,0,1);
 r->Do(histo2a,"foreach(dilepton); _x = [dilepton]->Rapidity(); _y = [dilepton]->Pt(); ");

r->Print();

r->loop(20000);   // do n events

 TCanvas *can1 = new 
     TCanvas ("can1","pi0 Dalitz decay ee mass");
 can1->SetLogy(1);
 
 histo1->Draw();
 histo1a->Draw("e1same");

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

