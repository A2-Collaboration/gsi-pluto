{

    TF1 *angles1=new TF1("angles1","(1+x*x)/2",-1,1);
    PAngularDistribution *ang = new PAngularDistribution("my_angle","My angular dist");
    ang->Add("q,parent,reference");
    ang->Add("p,daughter");
    ang->Add("NP11+,daughter,primary");
    ang->SetAngleFunction(angles1);
    makeDistributionManager()->Add(ang);
 
    PReaction my_reaction (3,"p","p","p NP11+ [p pi0]","angular_distribution",1,0,0,0);
    
    TH1F * histo1 = new TH1F ("histo1","cos theta of p1 ",40,-1,-1);
    my_reaction.Do(histo1,"p1 = [p,1]; p1->Boost([p+p]); _x =  p1->CosTheta()");
    TH1F * histo2 = new TH1F ("histo2","cos theta of p2 ",40,-1,-1);
    my_reaction.Do(histo2,"p2 = [p,2]; p2->Boost([p+p]); _x =  p2->CosTheta()");

    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(100000);

    histo1->Draw();
    histo2->Draw("same");

}
