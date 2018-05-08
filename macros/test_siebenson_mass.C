{

    PUtils::SetSeed(123);

    makeDistributionManager()->Exec("strangeness:init");
    
     pid_resonance = makeStaticData()->AddParticle(80,"DeltaPP1950", 1.89);        // 1.25=mass
     makeStaticData()->SetParticleTotalWidth("DeltaPP1950",0.33);
     makeStaticData()->SetParticleBaryon("DeltaPP1950",1);     
     makeStaticData()->SetParticleLMass("DeltaPP1950",1.70);  
     makeStaticData()->SetParticleUMass("DeltaPP1950",2.3);  
     decay_index = makeStaticData()->AddDecay(-1,"DeltaPP1950 --> Sigma1385+ + K+", "DeltaPP1950","Sigma1385+,K+",0.5);
     
     pid_resonance = makeStaticData()->AddParticle(81,"DeltaPP2420", 2.42);        // 1.25=mass
     makeStaticData()->SetParticleTotalWidth("DeltaPP2420",0.4);
     makeStaticData()->SetParticleBaryon("DeltaPP2420",1);
     makeStaticData()->SetParticleLMass("DeltaPP2420",1.70);  
     makeStaticData()->SetParticleUMass("DeltaPP2420",2.3);  
     decay_index = makeStaticData()->AddDecay(-1,"DeltaPP2420 --> Sigma1385+ + K+", "DeltaPP2420","Sigma1385+,K+",0.5);

     TH1F *hf1new= new TH1F("hf1new","",100,1.,.2.8);
     //PReaction my_reaction("_T1=3.5","p","p","p DeltaPP2420 [Sigma1385+ K+]");
     PReaction my_reaction("_T1=3.5","p","p","p DeltaPP1950");

     PPlutoBulkDecay *pl = new PPlutoBulkDecay();
     pl->SetRecursiveMode(1);  //Let also the products decay
     pl->SetTauMax(0.001);     //maxTau in ns
     my_reaction.AddBulk(pl);
     
 
     my_reaction.Print();
     //my_reaction.Do(hf1new,"_x = [DeltaPP2420]->M()");
     //my_reaction.Do(hf1new,"_x = [DeltaPP1950]->M()");
     my_reaction.Do(hf1new,"_x = [Sigma1385+]->M()");
     my_reaction.Loop(300000);

     hf1new->Draw();

 }
