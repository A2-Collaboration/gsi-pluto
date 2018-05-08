//TITLE Re-read the ASCII-file from batch_write_ntuple.C and make a bulk decay

{
    //Define an "empty" reaction
    PReaction my_reaction("eta_decays");
    
    PProjector *input = new PProjector();
    input->AddCommand("myeta = P3M(0.,0.,0.,eta.mass)");
    input->AddCommand("myeta->SetID(eta.pid)");
    input->AddInputASCII("eta_sample.txt", 
			 "readline{@px @py @pz @mass}; myeta->SetXYZM(px,py,pz,mass)");
    input->AddCommand("push(myeta)");
    my_reaction.AddPrologueBulk(input); //The "prolog" is done before the decay

    my_reaction.SetDecayAll(1.);
    
    my_reaction.Print();

    cout << my_reaction.Loop() << " events recovered" << endl;
}
