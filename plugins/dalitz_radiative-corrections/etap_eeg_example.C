//TITLE Activate radiative corrections for the Dalitz decay eta' --> e+ e- g

{
    // activate the corrections
    makeDistributionManager()->Exec("dalitz_corrections");

    // add the eta' Dalitz decay
    makeStaticData()->AddDecay("eta' Dalitz", "eta'", "dilepton,g", 0.0009);

    // change the form factor
    PSimpleVMDFF *ff = new PSimpleVMDFF("etaprime_ff@eta'_to_dilepton_g/formfactor", "Eta prime form factor",-1);
    ff->AddEquation("_ff2 = .5776*(.5776+.01)/((.5776-_q2)*(.5776-_q2)+.5776*.01)");  //equation from pluto paper page 11
    makeDistributionManager()->Add(ff);

    // define the reaction as usual
    PReaction my_reaction(1.600,"g","p","p eta' [dilepton [e+ e-] g]","etap_e+e-g_radCorr",1,0,0,0);

    my_reaction.Print();
    my_reaction.Loop(100000);
}

