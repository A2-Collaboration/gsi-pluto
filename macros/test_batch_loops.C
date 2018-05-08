{

    
    PReaction my_reaction("2.2","p","p","p p eta [dilepton [e+ e-] g]");

    my_reaction.Do("echo ********");

    my_reaction.Do("num_p = p.npar;  num_eta = eta.npar; echo $num_p , $num_eta");

    my_reaction.Do("foreach(p); [p]->Print();");
    my_reaction.Do("[p,1]->Print();");
    my_reaction.Do("[p,2]->Print();");

    my_reaction.Do("echo ******** 1 ********");

    my_reaction.Do("[*,1]->Print();");
    my_reaction.Do("num=1; echo ***************************************num is now $num");
    my_reaction.Do("[*,$num]->Print();");
    
    
    my_reaction.Do("echo ******** 2 ********");

    my_reaction.Do("[*,2]->Print();");
    my_reaction.Do("num=2");
    my_reaction.Do("[*,$num]->Print();");

    my_reaction.Do("echo ******** LOOP ********");
    my_reaction.Do("foreach(p); [p]->Print(); pos=p.cpos; echo $pos");

    my_reaction.Do("foreach(*); id = [*]->ID(); pos=*.cpos; echo found particle with $id at pos $pos");

    //Geht das:
    my_reaction.Do("echo ******** ???????????? ********");
    //my_reaction.Do("[p,3]->Print() ;");
    //my_reaction.Do("num=3; [p,$num]->Print() ;");

    my_reaction.Do("echo loop started; foreach(*); id = [*]->ID(); pos=*.cpos; echo found particle with $id at pos $pos");

    my_reaction.Do("echo ******** SELFMADE LOOP ********");
    // my_reaction.Do("cur = 1; if (cur < (p.npar +1.5)); echo $cur ; [p,$cur]->Print() ; cur = cur +1;");
    my_reaction.Do("cur = 0; myloop: if !(cur ~ p.npar); cur = cur +1; echo proton $cur ; [p,$cur]->Print() ;goto myloop");


    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(1);

    my_reaction.GetCurrentProjector()->Print();

}
