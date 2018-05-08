{

    PReaction my_reaction("2.2","p","p","p p eta [dilepton [e+ e-] g]", "eta_dalitz",1,0,0,0);
//If first number in quotation marks, it is the beam energy
//If not, it is the momentum 

    TH1F * histo1 = new TH1F ("histo1","dilepton mass",100,0.0,0.7);
    my_reaction.Do(histo1,"_x =  ([e+] + [e-])->M()");

    my_reaction.Do("opang = ([e+]->Angle([e-]))* 180. /TMath::Pi() ");
    my_reaction.Do("#opangfilter = 1; if opang<9 ; #opangfilter = 0");

    TH1F * histo2 = new TH1F ("histo2","dilepton mass",100,0.0,0.7);
    my_reaction.Do(histo2,"if opang>9 ; _x =  ([e+] + [e-])->M()");

    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(10000);

    histo1->Draw();
    histo2->Draw("same");

}
