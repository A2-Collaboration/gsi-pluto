//TITLE <b>Setting up reactions:</b> Use the one-liner parser to produce eta Dalitz events

{
    PUtils::SetSeed(123); //this is to have a fixed SEED. By default, the systime is used....

    PReaction my_reaction("_T1 = 2.2","p","p","p p eta [dilepton [e+ e-] g]", "eta_dalitz",1,0,0,0);
    
    //_T1 is the beam kinetic energy
    //one can use _P1 for the beam momentum
    //and _T2, _P2 in addition for collider experiments (separated by ';')

    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(100000);
}
