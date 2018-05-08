{


    PBatch x;
    x.AddCommand("my = 0.;");
    x.Execute();
    
    //First define our histograms
    TH1F * histo1 = new TH1F ("histo1","combinatorial background",100,0.,1.0);
    
    //Define the reaction as usual
    PReaction my_reaction(100,"p","p","p p pi0 [g dilepton [e+ e-]] pi0 [g dilepton [e+ e-]]");
    //PReaction my_reaction(100,"p","p","p p omega [e+ e-]");
    //PReaction my_reaction(100,"p","p","p N*(1535)+ [p eta [g dielectron [e+ e-] ]]");
      

        
    my_reaction.Do("echo ***************");

    //my_reaction.Do("label: ([p,1]->P())->Print();");
    //my_reaction.Do("xx: my = my + 0.1 ;  yy: my->Print(); ([p,2]->P())->Print();  if (my < 0.5);  goto(xx)  ");
    //my_reaction.Do("xx: my = my + 0.1 ;  yy: my->Print();   ");
    //my_reaction.Do("([p,2]->P())->Print();  if (my < 0.5);  goto(xx)  ");
    
    //my_reaction.Do("label: ([p]->P())->Print();   formore(p); goto(label)");
    //my_reaction.Do(histo1,"label2: _x =([e+] + [e-])->M() ;   ");
    //my_reaction.Do("formore(e+); goto label2 ");
    //my_reaction.Do("formore(e-); goto label2 ");

    my_reaction.Do("[dilepton]->SetM(my)");
    my_reaction.Do("my = my + 1");
    my_reaction.Do("([dilepton]->M())->Print()");
    
    //my_reaction.Do("[dilepton]->SetXYZM(0.1,0.2,0.3,my)");
    my_reaction.Do("[dilepton]->SetXYZM(0.1,0.2,0.3,my)");
    my_reaction.Do("[dilepton]->Print()");
    

    //my_reaction.Loop(5000);
    my_reaction.Print();
    my_reaction.Loop(5);
    
    // m1->Print();
    
    //histo_p1->Draw("");

}
