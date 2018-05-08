//TITLE pp elastic

{{

  //makeDistributionManager()->Disable("pp_elastic");
  
  Int_t nev = 10000;
  //nev = 1; //X
  
  TCanvas *obl_test_elastic_c1 = new        //X
    TCanvas ("obl_test_elastic_c1","pp elastic");//X
  obl_test_elastic_c1->SetBorderMode(0);    //X
  obl_test_elastic_c1->SetFillColor(0);     //X
  obl_test_elastic_c1->SetLogy(1);          //X
  obl_test_elastic_c1->SetTicky(1);         //X
  obl_test_elastic_c1->SetTickx(1);         //X
  obl_test_elastic_c1->SetFrameFillColor(0);//X

  obl_test_elastic_c1->Divide(2,2);                    //X
  obl_test_elastic_c1->GetPad(1)->SetFrameFillColor(0);//X
  obl_test_elastic_c1->GetPad(2)->SetFrameFillColor(0);//X
  obl_test_elastic_c1->GetPad(3)->SetFrameFillColor(0);//X
  obl_test_elastic_c1->GetPad(4)->SetFrameFillColor(0);//X

  TCanvas *obl_test_elastic_c2 = new        //X
    TCanvas ("obl_test_elastic_c2","pp elastic 2");//X
  obl_test_elastic_c2->SetBorderMode(0);    //X
  obl_test_elastic_c2->SetFillColor(0);     //X
  obl_test_elastic_c2->SetTicky(1);         //X
  obl_test_elastic_c2->SetTickx(1);         //X
  obl_test_elastic_c2->SetFrameFillColor(0);//X

  PReaction ela1("0.050","p","p","p p");
  ela1.Print();
  TH1F * theta1 = 
    new TH1F ("theta1","Polar angle at 50 MeV",100,-1,1);
  ela1.Do(theta1,"p1 = [p,1]; p1->Boost([p + p]); _x = p1->CosTheta()");
  ela1.Loop(nev);
  obl_test_elastic_c1->cd(1);       //X
  theta1->SetXTitle("cos (#theta)");//X
  theta1->Draw("");                 //X



  PReaction ela2("1.25","p","p","p p");
  ela2.Print();
  TH1F * theta2 = 
    new TH1F ("theta2","Polar angle at 1.25 GeV",100,-1,1);
  ela2.Do(theta2,"p1 = [p,1]; p1->Boost([p + p]); _x = p1->CosTheta()");
  ela2.Loop(nev);
  obl_test_elastic_c1->cd(2);       //X
  theta2->SetXTitle("cos (#theta)");//X
  theta2->Draw("");                 //X


    
  PReaction ela3("2.2","p","p","p p");
  ela3.Print();
  TH1F * theta3 = 
    new TH1F ("theta3","Polar angle at 2.2 GeV",100,-1,1);
  ela3.Do(theta3,"p1 = [p,1]; p1->Boost([p + p]); _x = p1->CosTheta()");
  ela3.Loop(nev);
  obl_test_elastic_c1->cd(3);       //X
  theta3->SetXTitle("cos (#theta)");//X
  theta3->Draw("");                 //X



  PReaction ela4("0.050","p","d","p p (n)");
  ela4.Print();
  TH1F * theta4 = 
    new TH1F ("theta4","Polar angle at 50 MeV (quasi-free d)",100,-1,1);
  TH1F * phi4 = 
    new TH1F ("phi4","Azimuthal angle at 50 MeV (quasi-free d)",100,-180.,180.);
  ela4.Do(theta4,"p1 = [p,1]; p1->Boost([p + p]); _x = p1->CosTheta()");
  ela4.Do(phi4,"_x = [p,1]->Phi()* 180. /TMath::Pi()");
  ela4.Do(theta4,"_x = 1.");
  ela4.Loop(nev);
  obl_test_elastic_c1->cd(4);       //X
  theta4->SetXTitle("cos (#theta)");//X
  theta4->Draw("");                 //X


    
  obl_test_elastic_c2->cd();    //X
  phi4->SetMinimum(0);          //X
  phi4->SetXTitle("cos (#phi)");//X
  phi4->Draw("");               //X

  obl_test_elastic_c1->Update();   //X
  obl_test_elastic_c1->Modified(); //X
  obl_test_elastic_c2->Print("obl_test_elastic@Elastic_pp_scattering:_azimuthal_angle_in_p+d.png");//X
  obl_test_elastic_c1->Print("obl_test_elastic@Elastic_pp_scattering:_polar_angle.png");           //X

}}
