 {

     
     makeDistributionManager()->Print();
     
     //PReaction my_reaction(2.2,"p","p","p p eta [dilepton [e+ e-] pi+ pi-]");
     PReaction my_reaction(2.2,"p","p","p p eta [e+ e- pi+ pi-]");


     
     my_reaction.Print();

     //Test the dilepton mass distribution:
     //makeDistributionManager()->GetDistribution("eta_ee_pipi")->Draw();
     
   
     //This is for on-line debugging
     TH1F * histo1 = new TH1F ("histo1","ee invariant mass",100,0.,0.6);
     histo1->Sumw2();
     my_reaction.Do(histo1,"_x = ([e+,1] + [e-,1])->M() ");
     
     my_reaction.Loop(100000);
     
     histo1->Draw("");
     
     
   

 }
