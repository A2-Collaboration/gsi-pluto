{

    //The PParticles must be created before the PProjector:
    PParticle beam("p",0,0,1.696);
    PParticle target("d");
    PParticle q=beam + target;

    PParticle p1("pi+");
    PParticle p2("pi-");
    PParticle p3("pi0");
    
    PParticle he("He3");
    PParticle *part[]={&q,&p1,&p2,&p3,&he};
 

    //Create and compile the script in advance.
    PProjector *m1 = new PProjector();
    //m1->AddCommand("");
    m1->AddCommand("he3_miss = [p + d] - [He3]"); //the p+d composite must have been created, but this is done.....
    m1->AddCommand("#miss=0.; if he3_miss->M() > 0.5 && he3_miss < 0.6; #miss=1.");

    // m1->AddCommand("foreach(*); id = [*]->ID(); echo $id");
    // m1->AddCommand("mass=he3_miss->M();  echo $mass , $#miss");
    // m1->AddCommand("echo *** end of event ***");
    

    for (Int_t i=49;i<=200;i++){

	//the PChannels *must* be created inside loop, they cannot be re-used
	PChannel prod(part,4,1);
	PChannel *c[ ]={prod};

	PReaction *pd = new PReaction(c,Form("pd2HePipPimPi0_MMHe3cut_app250k_%i",i),1,0,0,0,0);
    
	pd.AddBulk(m1); //re-use precompiled script
	
	pd->Print();
	pd->loop(250000);
    }
	

}
