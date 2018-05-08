{

	
Float_t Eb   = 25;
//Float_t T1    = 0.130;   // temperature in GeV
Float_t T1    = 0.250;   // temperature in GeV
Float_t T2    = 0.;      // temperature in GeV
Float_t blast = 0.3;      // radial expansion velocity
Float_t w = 0.;      // elliptic flow

Int_t const num_of_react=1000;

//original:
  PFireball *source_P=new PFireball("w",Eb,T1,T2,1.,blast,0.,0.,0.,0.);
//PFireball *source_P=new PFireball("w",Eb,T1,T2,1.,blast,w,w,w,w);
//source_P->setSigma(0.41);
source_P->Print();

PParticle *P = new PParticle("w");
 PParticle *mumP = new PParticle("mu-");
 PParticle *mupP = new PParticle("mu+");
 //PParticle *mumP = new PParticle("e-");
 //PParticle *mupP = new PParticle("e+");

PParticle* s_P[]    = {source_P,P};
//PChannel* c_sP   = new PChannel(s_P,1,1);
PChannel* c_sP   = source_P->makeChannel(1, 0);
PParticle** ptcls;
ptcls = c_sP->GetParticles();


//PParticle *s_Pdimu[]  ={P,mumP,mupP};
PParticle *s_Pdimu[]  ={ptcls[1],mumP,mupP};
PChannel  *c_Pdimu  = new PChannel(s_Pdimu,2,1);
PChannel  *cc_P[]   = {c_sP,c_Pdimu};

PReaction *r_P = new PReaction(cc_P,"omega_25gev_0009",sizeof(cc_P)/sizeof(cc_P[0]),0,0,0,0);
 TH2F * histo2 = new TH2F ("histo2","Rap. vs. Pt",50,-1.5,3.5,50,0,2.0);
 r_P->Do(histo2,"_x = ([e+] + [e-])->Rapidity(); _y=([e+] + [e-])->Pt(); ");


//r_P->Print();
r_P->setHGeant(0);
r_P->loop(num_of_react);

 histo2->Draw("colz");

}


