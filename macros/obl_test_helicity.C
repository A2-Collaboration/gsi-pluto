//TITLE Angles

Double_t hel(Double_t *x, Double_t *par) {//X
   return par[0]*(1 + par[1]*x[0]*x[0]);  //X
}                                         //X

void obl_test_helicity() {

    TCanvas *obl_test_helicity_c1 = new                       //X
	TCanvas ("obl_test_helicity_c1","pi0 helicity angle");//X
    
    obl_test_helicity_c1->SetBorderMode(0);    //X
    obl_test_helicity_c1->SetFillColor(0);     //X
    obl_test_helicity_c1->SetTicky(1);         //X
    obl_test_helicity_c1->SetTickx(1);         //X
    obl_test_helicity_c1->SetFrameFillColor(0);//X


    TH1F * histo1 = new TH1F ("histo1","D+ angle",100,-1,1);
    histo1->SetXTitle("cos (#theta ^{#Delta +})");//X
    TH1F * histo2 = new TH1F ("histo2","helicity angle",50,0,1);
    histo2->SetXTitle("cos (#alpha ^{ee}_{e})");//X
    histo2->Sumw2();//X

    //makeDistributionManager()->Disable("helicity_angles");
    
    PReaction my_reaction(3.13,"p","p","p D+ [p pi0 [dilepton [e+ e-] g]]");
    
    my_reaction.Do(histo2,"_pi0=[pi0]; _ep=[e+]; _em=[e-]; _pi0->Boost([p + p]); _em->Boost([p + p]); _ep->Boost([p + p]); _ep->Rot(_pi0); _em->Rot(_pi0); _pi0->Rot(_pi0) ; _ep->Boost(_pi0); _em->Boost(_pi0); dil=_ep+_em; _ep->Rot(dil); _em->Rot(dil); dil->Rot(dil); _ep->Boost(dil); _em->Boost(dil) ; s1= _ep->Theta();s2=_em->Theta(); _x = fabs(cos(s1)) ");
    
    my_reaction.Do(histo1,"mydelta = [D+]; mydelta->Boost([p + p]); _x = mydelta->CosTheta() ");
    
    my_reaction.Print();
    my_reaction.Loop(100000);
        
    histo2->SetMinimum(0);//X
    histo2->Draw("e1");   //X
    
    
    TF1 *f1 = new TF1("myfunc",hel,0,1,2);//X
    histo2->Fit(f1);                      //X

    
    TCanvas *obl_test_helicity_c2 = new                        //X
	TCanvas ("obl_test_helicity_c2","D+ production angle");//X
    obl_test_helicity_c2->SetBorderMode(0);    //X
    obl_test_helicity_c2->SetFillColor(0);     //X
    obl_test_helicity_c2->SetTicky(1);         //X
    obl_test_helicity_c2->SetTickx(1);         //X
    obl_test_helicity_c2->SetFrameFillColor(0);//X


    histo1->Draw("");//X

    obl_test_helicity_c1->Update();   //X
    obl_test_helicity_c1->Modified(); //X
    obl_test_helicity_c2->Print("obl_test_helicity@Angles:_Polar_angle_of_the_Delta+.png");//X
    obl_test_helicity_c1->Print("obl_test_helicity@Angles:_Helicity_angle_of_eta.png");    //X

}
