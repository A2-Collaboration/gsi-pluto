 {


     PReaction my_reaction("2.2","p","p","p p w [dilepton [e+ e-] pi0]");
     my_reaction.Print();

   
     //This is for on-line debugging
     TH1F * histo1 = new TH1F ("histo1","ee invariant mass",100,0.,0.8);
     //histo1->Sumw2();
     my_reaction.Do(histo1,"_x = ([e+,1] + [e-,1])->M() ");
     
     my_reaction.Loop(10000);
     
     histo1->Draw("");


     PSimpleVMDFF *ff = new PSimpleVMDFF("vmd_ff_dd@w_to_dilepton_pi0/formfactor",
					 "VMD form factor",-1);
     ff->SetVectorMesonMass(0.77);
     ff->SetWeightMax(5.);
     
     //0.17918225/( (0.4225-m2)*(0.4225-m2) + 0.000676);
     //ff->AddEquation("_ff2 = 0.17918225/( (0.4225- _q2)*(0.4225- _q2) + 0.000676)");
     ff->AddEquation("_ff2 = 5.");

     makeDistributionManager()->Add(ff);
     makeDistributionManager()->LinkDB();

     PReaction my_reaction2("2.2","p","p","p p w [dilepton [e+ e-] pi0]");
     my_reaction2.Print();

   
     //This is for on-line debugging
     TH1F * histo2 = new TH1F ("histo2","ee invariant mass",100,0.,0.8);
     histo2->Sumw2();
     my_reaction2.Do(histo2,"_x = ([e+,1] + [e-,1])->M() ");
     
     my_reaction2.Loop(10000);
     
     histo2->Draw("e1same");
     
     
   

 }
