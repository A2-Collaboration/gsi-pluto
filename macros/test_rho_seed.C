{

    PParticle rho("rho0");

    // decay of the rho -> dilepton
    PParticle em("e-");
    PParticle ep("e+");
    PParticle *dileptons[]={&rho,&em,&ep};
    PChannel dilepton_decay(dileptons,2,1);
    PChannel *c[ ]={&dilepton_decay};
    PReaction r(c,"rho",1,1,0,0,0);

    PProjector *input = new PProjector();
    input->AddInputASCII("rho.txt", "readline{@rho_px @rho_py @rho_pz @rho_mass}; echo $rho_px $rho_py $rho_pz $rho_mass");
    //input->AddInputASCII("rho.txt", "readline{@rho_px @rho_py @rho_pz}; echo $rho_px");
    //input->AddInputASCII("rho.txt", "readline{@rho_px}; echo $rho_px");


    //input->AddInputASCII("rho.txt", "readline{@rho_px @x @e}; echo $rho_px");
    r->AddPrologueBulk(input);

    //r.Do("mass = [rho0]->M(); echo $mass");

    r.Print();
    r.loop(100);


}
