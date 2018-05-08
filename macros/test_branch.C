{
    PReaction my_reaction("_T1=2.2","p","p","p p eta [dilepton [e+ e-] g]" ,"eta_dalitz");
    my_reaction.Do("bla = Branch(newparticles);echo $bla");
    my_reaction.Do("bla = Branch(newparticles2);echo $bla");
    my_reaction.Do("bla = Branch(newparticles);echo $bla");
    my_reaction.Do("bla = Branch(newparticles3);echo $bla");
    my_reaction.Do("echo ----------");
    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(100);
}
