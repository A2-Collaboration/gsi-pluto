{
 
    TH1F *histo = new TH1F("histo", "ee mass", 20, 0., 0.6);

    makeDistributionManager()->Unpack("pluto_ee_filter.root");
    makeDistributionManager()->Startup("_filter_debug=1");
    //makeDistributionManager()->Startup("_filter_remove_particles=1");

    PReaction my_reaction("2.2","p","p","p p eta [dilepton [e+ e-] g]","filtered_events");
    //my_reaction.Do("#eefilter=bool_acc");
    my_reaction.Do(histo,"if (bool_acc); _x = [dilepton]->M()");
    

    //my_reaction.Do("theta_ep = ([e+]->Theta() * 180.)/TMath::Pi()");
    //my_reaction.Do("theta_em = ([e-]->Theta() * 180.)/TMath::Pi()");
    //my_reaction.Do("#filter = 1; if theta_ep<18 || theta_ep>85 || theta_em<18 || theta_em>85; #filter = 0");

    my_reaction.Print();
    my_reaction.Loop(10000);


}
