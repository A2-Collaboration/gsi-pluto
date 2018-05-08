{

    PParticle omega("w");

    // decay of the w -> dilepton
    PParticle em("e-");
    PParticle ep("e+");
    PParticle *dileptons[]={&omega,&em,&ep};
    PChannel dilepton_decay(dileptons,2,1);
    PChannel *c[ ]={&dilepton_decay};
    PReaction r(c,"omega",1,1,0,0,0);

    PProjector *input = new PProjector();
    input->AddInputASCII("omega_PYTHIA.txt", 
			 "readline{@px @py @pz @mass}; [w]->SetXYZM(px,py,pz,mass)");
    r->AddPrologueBulk(input); //The "prolog" is done before the decay

    // This is just for debugging:
    // r.Do("mass2 = [w]->M(); px2=[w]->Px(); py2=[w]->Py(); echo $px2,$py2,$mass2");
    // r.Do("mass3 = ([e+]+[e-])->M(); px3=([e+]+[e-])->Px(); py3=([e+]+[e-])->Py(); echo $px3,$py3,$mass3");

    r.Print();
    r.loop(100);


}
