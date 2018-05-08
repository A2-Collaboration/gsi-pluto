{

    Int_t pid_lambda1405 = makeStaticData()->AddParticle(-1,"Lambda1405", 1.406);
    makeStaticData()->SetParticleTotalWidth("Lambda1405",0.05);
    makeStaticData()->SetParticleBaryon("Lambda1405",1);


    makeStaticData()->AddDecay("Lambda(1405)0 --> Sigma+ + pi-", "Lambda1405", 
			       "Sigma+,pi-", .3333 );
    makeStaticData()->AddDecay("Lambda(1405)0 --> Sigma- + pi+", "Lambda1405", 
			       "Sigma-,pi+", .3333 );
    makeStaticData()->AddDecay("Lambda(1405)0 --> Sigma0 + pi0", "Lambda1405", 
			       "Sigma0,pi0", .3333 );

    
    Int_t pid_sigma1385_p = makeStaticData()->AddParticle(-1,"Sigma1385+", 1.3828);
    makeStaticData()->SetParticleTotalWidth("Sigma1385+",0.0358);
    makeStaticData()->SetParticleBaryon("Sigma1385+",1);

    makeStaticData()->AddDecay("Sigma(1385)+ --> Lambda + pi+", "Sigma1385+",
			       "Lambda,pi+", .8815 );    
    makeStaticData()->AddDecay("Sigma(1385)+ --> Sigma+ + pi0", "Sigma1385+",
			       "Sigma+, pi0", .05925);
    makeStaticData()->AddDecay("Sigma(1385)+ --> Sigma0 + pi+", "Sigma1385+",
			       "Sigma0, pi+", .05925);
    

    listParticle("Lambda1405");
    listParticle("Sigma1385+");

    

    TH1F histo("histo", "Lambda1405 mass", 100, 1.2, 1.6);
    PReaction r("3.5", "p", "p", "p Lambda1405 [Sigma+ pi-]");
    r.Do(&histo,"_x = [Lambda1405]->M()");
    r.Print();
    r.Loop(10000);
    
    histo.Draw();

}
