{


    *(makeStaticData()->GetBatchValue("_system_force_m1n")) = 1;

    makeStaticData()->AddAlias("DP33+","Delta1600+");

    //This is for particles which exists already:
    PResonanceDalitz * newmodel = 
	new PResonanceDalitz("Delta1600+_dalitz@Delta1600+_to_p_dilepton",
			     "dgdm from Zetenyi/Wolf",-1);
    newmodel->setGm(0.162);
    makeDistributionManager()->Add(newmodel);



    //Here we add new particles

    makeStaticData()->AddParticle(75,"N1720+", 1.72);   //PID, name, mass
    makeStaticData()->AddAlias("N1720+","N*(1720)+");
    makeStaticData()->SetParticleTotalWidth("N1720+",0.2);
    makeStaticData()->SetParticleBaryon("N1720+",1);
    makeStaticData()->SetParticleSpin("N1720+",3);
    makeStaticData()->SetParticleCharge("N1720+",1);
    makeStaticData()->SetParticleIsospin("N1720+",1);
    makeStaticData()->SetParticleParity("N1720+",1);
    
    makeStaticData()->AddDecay("N1720+ --> n + pi+", "N1720+",
			       "n,pi+", 0.1 );    //Not the real BR!
    Int_t n1720p_ee_id =makeStaticData()->AddDecay("N1720+ --> dilepton + p", "N1720+",
						   "dilepton,p", .001 );  //Not the real BR!

    PResonanceDalitz * newmodel = 
	new PResonanceDalitz("N1720+_dalitz@N1720+_to_p_dilepton",
			     "dgdm from Zetenyi/Wolf",-1);
    newmodel->setGm(0.193);
    makeDistributionManager()->Add(newmodel);



    //N*, First constant for charge =0, second constant for charge =+
	// setGm(0 , 1.98  , 1.980, 1, 3);              // D1232
       // setGm(1 , 0.098 , 0.139, 1, 1);              // N1440
       // setGm(2 , 0.719 , 0.793,-1, 3);              // N1520
       // setGm(3 , 0.549 , 0.634,-1, 1);              // N1535
       // setGm(4 , 0.285 , 0.315,-1, 1);              // N1650
       // setGm(5 , 1.52  , 0.678,-1, 5);              // N1675
       // setGm(6 , 2.74  , 0.971, 1, 5);              // N1680
       // setGm(7 , 0.193 , 0.126,-1, 3);              // N1700
       // setGm(8 , 0.019 , 0.030, 1, 1);              // N1710
       // setGm(9 , 0.386 , 0.193, 1, 3);              // N1720
       // setGm(10, 0.162 , 0.162, 1, 3);              // D1600
       // setGm(11, 0.162 , 0.162,-1, 1);              // D1620
       // setGm(12, 0.549 , 0.549,-1, 3);              // D1700
       // setGm(13, 0.713 , 0.713, 1, 5);              // D1905
       // setGm(14, 0.066 , 0.066, 1, 1);              // D1910
       // setGm(15, 0.479 , 0.479,-1, 5);              // D1930 
    
    PReaction dummy("_T1=3.5","p","p","p p N1720+ "); //dummy, just to add the composite
    makeDistributionManager()->Print("");
    makeDistributionManager()->Print("decay_models");


    //PReaction my_reaction("_T1=3.5","p","p","p p N1720+ [pi+ n]");
    PReaction my_reaction("_T1=3.5","p","p","p p N1720+ [dilepton [e+ e-] p]");
 
    TH1F * histo1 = new TH1F ("histo1","ee invariant mass",100,0.,1.2);
    my_reaction.Do(histo1,"_x = ([e+,1] + [e-,1])->M() ");

    TH1F * histo2 = new TH1F ("histo2","n* invariant mass",100,0.9,1.5);
    my_reaction.Do(histo2,"_x = ([N1720+])->M() ");

    my_reaction.Print();
    my_reaction.Loop(10000);
    
    

#if 0
    PBreitWigner * bw = (PBreitWigner *)(makeDistributionManager()->GetDistribution("N1720+_bw"));
    bw->Draw();
    bw->SetParameter(1,-1);
    PBreitWigner * bw_ee = bw->Clone();
    bw_ee->SetParameter(1,n1720p_ee_id);
    bw_ee->Draw("same");
#endif

    histo2->Draw();
	

}
