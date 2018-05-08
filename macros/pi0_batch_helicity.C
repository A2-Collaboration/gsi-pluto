Double_t hel(Double_t *x, Double_t *par) {
   return par[0]*(1 + par[1]*x[0]*x[0]);
}

void pi0_batch_helicity() {


    //gSystem->CompileMacro( "../macros/PHadesAcc.C");


TH1F * histo2 = new TH1F ("theta","cos theta",20,0.,1.);
TH1F * histo1 = new TH1F ("histo1","opang",50,0,3.14);
TH1F * histo3 = new TH1F ("histo3","phi",50,-3.14,3.14);
histo2->Sumw2();
histo1->Sumw2();

//makeDistributionManager()->Disable("pp_delta+_angle");
//makeDistributionManager()->Disable("pp_delta_waves1");
//makeDistributionManager()->Disable("helicity_angles");

//PReaction my_reaction(3.13,"p","p","p p pi0 [dilepton [e+ e-] g]");
PReaction my_reaction(3.13,"p","p","p D+ [p pi0 [dilepton [e+ e-] g]]");

//my_reaction.AddBulk(new PHadesAcc());

PProjector *m1 = new PProjector(); 

//m1->AddHistogram(histo1,"if(_opang);_ep=[e+]; _em=[e-]; _x=_ep->Angle(_em)");

//m1->AddHistogram(histo2,"if(_hadacc); if(_opang);_pi0=[pi0]; _ep=[e+]; _em=[e-]; _pi0->Boost([p + p]); _em->Boost([p + p]); _ep->Boost([p + p]); _ep->Rot(_pi0); _em->Rot(_pi0); _pi0->Rot(_pi0) ; _ep->Boost(_pi0); _em->Boost(_pi0); dil=_ep+_em; _ep->Rot(dil); _em->Rot(dil); dil->Rot(dil); _ep->Boost(dil); _em->Boost(dil) ; s1= _ep->Theta();s2=_em->Theta(); _x = fabs(cos(s1)) ");

m1->AddHistogram(histo2,"_pi0=[pi0]; _ep=[e+]; _em=[e-]; _pi0->Boost([p + p]); _em->Boost([p + p]); _ep->Boost([p + p]); _ep->Rot(_pi0); _em->Rot(_pi0); _pi0->Rot(_pi0) ; _ep->Boost(_pi0); _em->Boost(_pi0); dil=_ep+_em; _ep->Rot(dil); _em->Rot(dil); dil->Rot(dil); _ep->Boost(dil); _em->Boost(dil) ; s1= _ep->Theta();s2=_em->Theta(); _x = fabs(cos(s1)) ");

 m1->AddHistogram(histo3,"_x=[e+]->Phi()");

my_reaction.AddBulk(m1);



my_reaction.Print();
my_reaction.Loop(20000);

//histo->Draw("box");

 histo2->SetMinimum(0);
histo2->Draw("e1");


TF1 *f1 = new TF1("myfunc",hel,0,1,2);
 histo2->Fit(f1);

//histo1->Draw("e1");
}
