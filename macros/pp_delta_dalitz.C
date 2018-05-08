// p + p -> p + Delta+ -> p + p + gamma* -> p + p + e- + e+
// This macro tests:
//   1. The mass-dependent width Breit-Wigner distribution of the Delta 
//   2. The anisotropic (s+p wave) production angle for the pp->pDelta channel
//   3. The dilepton mass for Delta Dalitz decay

{
    gROOT->Reset();
    
    // use to get pure QED form factors:
    ((PDalitzDecay *) makeDistributionManager()->GetDistribution("D+_dalitz"))->SetUseQED(1);
    
    PParticle *p1 = new PParticle("p",2.2);
    PParticle *p2 = new PParticle("p");
    PParticle *p3 = new PParticle("p"); 
    PParticle *p4 = new PParticle("p");
    PParticle *delta = new PParticle("D+");
    PParticle *g  = new PParticle("dilepton");
    PParticle *em = new PParticle("e-");
    PParticle *ep = new PParticle("e+");
    PParticle *q  = new PParticle(*p1 + *p2);
    PParticle *s1[] = {q,p3,delta};
    PParticle *s2[] = {delta,p4,g};
    PParticle *s3[] = {g,em,ep};
    PChannel *c1 = new PChannel(s1);
    PChannel *c2 = new PChannel(s2); 
    PChannel *c3 = new PChannel(s3); 
    PChannel *cc[] = {c1,c2,c3};

    PReaction *r = new PReaction(cc, "pp_delta_dalitz", 3, 1);

    //filtered dilepton pairs:
    r->Do("theta_ep = ([e+]->Theta() * 180.)/TMath::Pi()");
    r->Do("theta_em = ([e-]->Theta() * 180.)/TMath::Pi()");
    r->Do("opang = ([e+]->Angle([e-]) * 180.)/TMath::Pi()");
    r->Do("filter=1; if ((theta_ep<18 || theta_ep>85 || theta_em<18 || theta_em>85) || opang < 9) filter=0");
    r->Do("if (filter); [dilepton]->Push(Branch(Accepted)); [e+]->Push(Branch(Accepted)); [e-]->Push(Branch(Accepted))");
    
    makeDistributionManager()->GetDistribution("D+_dalitz")->Print();
    
    r->Print();
    r->loop(100000);
}
