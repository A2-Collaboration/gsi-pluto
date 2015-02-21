//TITLE Using a batch script to simulate pileup

{

    //first, we define in the data base:
    //1.) The mean time difference between 2 events
    //    N.B. the time in Pluto is in mm/s, so we have to convert it
    //    let us start with 100000 events/s
    makeStaticData()->SetBatchValue("beam_particle_mean_time", 10e-5 * TMath::C() / 10e-3);
    
    //2.) The time of the first event (we choose 0)
    makeStaticData()->SetBatchValue("absolute_event_time", 0);
    
    //3.) The threshold for a pile-up event (just e.g. 50% of beam_particle_mean_time)
    makeStaticData()->SetBatchValue("pileup_time_theshold", 0.5 * 10e-5 * TMath::C() / 10e-3);

    //4.) Just a flag for pileup
    makeStaticData()->SetBatchValue("pileup_flag" ,0);

    //this just constructs the temporary particles:
    makeDynamicData()->GetBatchParticle("p1");
    makeDynamicData()->GetBatchParticle("p2");
    makeDynamicData()->GetBatchParticle("eta1");

    //**************************************************************
    PReaction my_reaction("_T1=2.2","p","p","p p eta","pileup_eta");
    my_reaction.trackedParticles();  //this removes the "empty events"
    //**************************************************************
    
    //the following loop sets the absolute time for each particle
    my_reaction.Do("foreach(*); event_time = [*]->T(); [*]->SetT(event_time + absolute_event_time)");

    //set a event time difference based on a "flat" distribution
    my_reaction.Do("timediff = sampleFlat()*2*beam_particle_mean_time");
    
    //check if the last event was pileup and has been removed:
    my_reaction.Do("if !pileup_flag; goto check_pileup");
    my_reaction.Do("push(p1); push(p2); push(eta1);");
    my_reaction.Do("pileup_flag = 0");
    my_reaction.Do("goto event_ok");

    //check if we are below threshold
    my_reaction.Do("check_pileup:");
    my_reaction.Do("if timediff > pileup_time_theshold; goto event_ok");
    
    //if yes, store particles and kill event
    my_reaction.Do("p1 = [p,1]; p2 = [p,2]; eta1 = [eta];");
    my_reaction.Do("[p,1]->SetInActive(); [p,2]->SetInActive(); [eta]->SetInActive();");
    my_reaction.Do("pileup_flag = 1");
    

    my_reaction.Do("event_ok:");
    //set the new time for the next event
    my_reaction.Do("absolute_event_time = absolute_event_time + timediff");

    

    //this is just for debugging:
    //my_reaction.Do("time = [p,1]->T(); echo $time");
    //my_reaction.Do("foreach(*); [*]->Print();");
    

    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(100);
}
