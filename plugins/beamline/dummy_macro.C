{

    //Init the reaction for c.m. sampling:
    PBeamLineSimulation *sim = new PBeamLineSimulation("beam", "Beam line simulation");
    sim->SetReaction("p + p");
    //Read the data file
    sim->InitBeamLine("pibeam_set6.data");
    //Place (virtual) detector along the beam line. Units are in mm:
    sim->AddDetector("det1", -3000.0);
    sim->AddDetector("det2", -6000.0);


    for (double dist = -30000; dist<0; dist+=300)
	sim->AddDetector("complete", dist);
    //open the ROOT file and type:
    //data.Draw("complete.fV.fX:complete.fV.fZ >>h3(100,-30000,0,100,-80,80)","","colz");


    //Choose which module should be the target:
    sim->TargetIsElement(33);
    //Select the global momentum:
    sim->SetGlobalMomentum(3.0);

    //Beam profile at the production target
    //position in mm:
    sim->Do("_beam_x = 0;  _beam_y = 0;");
    //divergence in px/pz and py/pz:
    //sim->Do("_beam_px = 0; _beam_py = 0.05;");
    sim->Do("_beam_px = 0; _beam_py = 0;");
    //momentum spread:
    sim->Do("_beam_dp = 0.16*sampleFlat() - 0.08; "); //+- 8%
    //sim->Do("_beam_dp = -0.08; ");

    //Add and enable module:
    makeDistributionManager()->Add(sim);

    PReaction my_reaction("_T1 = 2.2","p","p","p p eta [dilepton [e+ e-] g]", "beam_line");
    //PReaction my_reaction("_T1 = 2.2","p","p","p p", "beam_line");

    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(1000);


}
