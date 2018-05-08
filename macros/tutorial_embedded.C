{

    PReaction my_reaction("embedded");

    //Construct the bulk container:
    PEmbeddedParticles * embedded = new PEmbeddedParticles();

    //Add an eta which we emit at a single point:
    PParticle * eta = new PParticle("eta",1.,2.,3.);
    //Just add the particle to the container: 
    embedded->AddParticle(eta);

    //Add dileptons's, which we emit in a small cone:
    PParticle * dilepton = new PParticle("dilepton");  
    embedded->AddParticle(dilepton,100); //downscaling
    embedded->SetSampling(0, 1.,   //pmin and pmax in lab frame 
			  TMath::Pi()/2., //opening angle
			  0., 0., //Theta, phi of pointing vect.
			  0.01 , 0.8 //mass range
			  ); 

    
    
    //Add our container to the reaction:
    my_reaction.AddBulk(embedded);

    TH1F * histo1 = new TH1F ("histo1","dilepton momentum",100,0.0,1.5);
    my_reaction.Do(histo1,"_x =  ([dilepton])->P()");
    TH1F * histo2 = new TH1F ("histo2","dilepton cos theta",100,-1,1);
    my_reaction.Do(histo2,"_x =  ([dilepton])->CosTheta()");
        
    PPlutoBulkDecay *pl = new PPlutoBulkDecay();
    pl->SetRecursiveMode(1);  //Let also the products decay
    pl->SetTauMax(0.001);     //maxTau in ns
    my_reaction.AddBulk(pl);

    TH1F * histo3 = new TH1F ("histo3","di-lepton mass",50,0.0,1.0);
    my_reaction.Do(histo3,"_x =  ([e+,1] + [e-,1])->M()");
    my_reaction.Do(histo3,"_x =  ([e+,2] + [e-,2])->M()");


    my_reaction.Print();

    my_reaction.Loop(100000);

    histo1->Draw();


}
