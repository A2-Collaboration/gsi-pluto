//TITLE Decay Manager

{{

    //This macro tests the pd fusion to He3

    makeStaticData()->SetBatchValue("_system_inactivate_decayed_particles",1);

    PDecayChannel *c1;
    PDecayChannel *c2;
    PDecayManager *pdm = new PDecayManager;
    pdm->SetVerbose(1);
    pdm->SetDefault("pi0");
    pdm->SetDefault("dilepton");
    PParticle *p = new PParticle("p",1.00);
    PParticle *d = new PParticle("d",0.,0.,0.,1.875613,45);
    PParticle *q = new PParticle(*p+*d);
    PParticle *eta = new PParticle("eta");
    c1 = new PDecayChannel;
    c2 = new PDecayChannel;
    c1->AddChannel(1.0,"He3","eta");
    c2->AddChannel(1.0,"pi0","pi0","pi0");
    pdm->AddChannel("eta",c2);
    pdm->InitReaction(q,c1);
    pdm->Print();

    //Debug //X
    //pdm->Do("echo -------------------------------------------");//X
    //pdm->Do("foreach(*) ; [*]->Print()");//X


    //Multiplicity for stable particles
    TH1F * histo1 = new TH1F ("histo1","part mult",21,-0.5,20.5);
    pdm->Do("counter = 0;"); 
    pdm->Do("foreach(*); counter = counter +1");
    //pdm->Do("echo $counter;"); 
    pdm->Do(histo1,"_x = counter");

    int n = pdm->loop(100000,0,"pd_eta_3pi0",0,0,0,1);
    //int n = pdm->loop(100,0,"pd_eta_3pi0",1,0,0,0); //X also unstable ones for DEBUG


    TCanvas *obl_test_dp_He3_c1 = new                       //X
	TCanvas ("obl_test_dp_He3_c1","pdm");//X
    
    obl_test_dp_He3_c1->SetBorderMode(0);    //X
    obl_test_dp_He3_c1->SetFillColor(0);     //X
    obl_test_dp_He3_c1->SetTicky(1);         //X
    obl_test_dp_He3_c1->SetLogy(1);          //X
    obl_test_dp_He3_c1->SetTickx(1);         //X
    obl_test_dp_He3_c1->SetFrameFillColor(0);//X
    
    histo1->SetXTitle("particle multiplicity");//X
    histo1->Draw();//X result should have maximum at 7

    obl_test_dp_He3_c1->Update();   //X
    obl_test_dp_He3_c1->Modified(); //X
    obl_test_dp_He3_c1->Print("obl_test_dp_He3@pdm:_Particle_multiplicity_in_dp_to_He3+X.png");//X

}}
