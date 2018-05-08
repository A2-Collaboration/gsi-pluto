{
    

    //old version:
    PReaction my_reaction(3.13,"p","p","p p eta [dilepton [e+ e-] g]");
    TH1F * histo1 = new TH1F ("histo1","dilepton mass",100,0.0,0.7);
    my_reaction.Do(histo1,"_x =  ([e+] + [e-])->M()");
 
    my_reaction.Print();
    my_reaction.Loop(10000);





    //Disable old eta Dalitz model and enable it for weighting:
    makeDistributionManager()->GetDistribution("eta_dalitz")->EnableWeighting();
    
    //TF1 object representing the di-lepton statistics:
    TF1 *flat = new TF1("flat","1",0,1);
    
    //The "PInclusiveModel" can be used as a generator:
    PInclusiveModel *dilepton_generator = 
	new PInclusiveModel("flat@eta_to_g_dilepton/generator","Dilepton generator",-1);
    
    //The distribution template:
    dilepton_generator->Add("eta,parent");
    dilepton_generator->Add("g,daughter");
    dilepton_generator->Add("dilepton,daughter,primary");
    dilepton_generator->SetSampleFunction(flat);

    //Enable distribution as a generator
    dilepton_generator->EnableGenerator();
    
    makeDistributionManager()->Add(dilepton_generator);
    
    PReaction my_reaction2(3.13,"p","p","p p eta [dilepton [e+ e-] g]");
    TH1F * histo2 = new TH1F ("histo2","dilepton mass",100,0.0,0.7);
    my_reaction2.Do(histo2,"_x =  ([e+] + [e-])->M()");
    

    my_reaction2.Print();

    //100 dummy events to have a good starting point for the normalization:
    my_reaction2.Preheating(100);
    my_reaction2.Loop(10000);
 
    histo2->Draw();
    histo2->SetLineColor(2);
    histo1->Draw("same");

   

    


}

