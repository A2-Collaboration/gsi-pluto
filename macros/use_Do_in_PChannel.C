//TITLE Using scripting in a PChannel

{ 

    PFireball *source=new PFireball("w",1.0,0.055,0.0,1.0,0.0,0.,0.,0.,0.); //1AGeV temp.
    source->setTrueThermal(kFALSE);

    PParticle *w = new PParticle("w");
    PParticle *s1[] = {source,w};
    PChannel  *c1 = new PChannel(s1,1,1);

    c1->Do("loop: [w]->SetM({w_bw}->SampleTotalMass()); if ([w]->M() > 0.85) goto loop");

    PParticle *ep = new PParticle("mu+");
    PParticle *em = new PParticle("mu-");
    PParticle *pi = new PParticle("pi0");
    PParticle *dm = new PParticle("dimuon");
    PParticle *s2[] = {w,dm,pi};
    PChannel  *c2 = new PChannel(s2,2,1);

    PChannel  *cc[] = {c1,c2};
    PReaction *r = new PReaction(cc,"thermal_w",2);
    TH1F *histo = new TH1F("histo", "mumu mass", 100, 0.1, 1.2);
    TH1F *histo2 = new TH1F("histo2", "w mass", 100, 0.1, 1.2);
    r->Do(histo, "_x = ([dimuon])->M()");
    r->Do(histo2, "_x = ([w])->M()");
    r->Print();
    r->loop(10000);

    histo2->Draw();

    //data.Draw("M()","ID()==41","");

}
