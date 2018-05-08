//TITLE pn elastic

{{

  
  Int_t nev = 10000;
  //nev = 1; //X
  
  TCanvas *obl_test_pn_elastic_c1 = new        //X
    TCanvas ("obl_test_pn_elastic_c1","pn elastic");//X
  obl_test_pn_elastic_c1->SetBorderMode(0);    //X
  obl_test_pn_elastic_c1->SetFillColor(0);     //X
  obl_test_pn_elastic_c1->SetLogy(1);          //X
  obl_test_pn_elastic_c1->SetTicky(1);         //X
  obl_test_pn_elastic_c1->SetTickx(1);         //X
  obl_test_pn_elastic_c1->SetFrameFillColor(0);//X

  obl_test_pn_elastic_c1->Divide(2,2);                    //X
  obl_test_pn_elastic_c1->GetPad(1)->SetFrameFillColor(0);//X
  obl_test_pn_elastic_c1->GetPad(2)->SetFrameFillColor(0);//X
  obl_test_pn_elastic_c1->GetPad(3)->SetFrameFillColor(0);//X
  obl_test_pn_elastic_c1->GetPad(4)->SetFrameFillColor(0);//X


  PReaction ela1("0.010","p","n","p n");
  ela1.Print();
  TH1F * theta1 = 
    new TH1F ("theta1","Polar angle at 10 MeV",100,-1,1);
  ela1.Do(theta1,"p1 = [p,1]; p1->Boost([p + n]); _x = p1->CosTheta()");
  ela1.Loop(nev);
  obl_test_pn_elastic_c1->cd(1);       //X
  theta1->SetXTitle("cos (#theta)");//X
  theta1->Draw("");                 //X



  PReaction ela2("0.159","p","n","p n");
  ela2.Print();
  TH1F * theta2 = 
    new TH1F ("theta2","Polar angle at 159 MeV",100,-1,1);
  ela2.Do(theta2,"p1 = [p,1]; p1->Boost([p + n]); _x = p1->CosTheta()");
  ela2.Loop(nev);
  obl_test_pn_elastic_c1->cd(2);       //X
  theta2->SetXTitle("cos (#theta)");//X
  theta2->Draw("");                 //X


    
  PReaction ela3("0.370","p","n","p n");
  ela3.Print();
  TH1F * theta3 = 
    new TH1F ("theta3","Polar angle at 370 MeV",100,-1,1);
  ela3.Do(theta3,"p1 = [p,1]; p1->Boost([p + n]); _x = p1->CosTheta()");
  ela3.Loop(nev);
  obl_test_pn_elastic_c1->cd(3);       //X
  theta3->SetXTitle("cos (#theta)");//X
  theta3->Draw("");                 //X



  PReaction ela4("2","p","n","p n");
  ela4.Print();
  TH1F * theta4 = 
    new TH1F ("theta4","Polar angle at 2 GeV",100,-1,1);
  ela4.Do(theta4,"p1 = [p,1]; p1->Boost([p + n]); _x = p1->CosTheta()");
  ela4.Do(theta4,"_x = 1.");
  ela4.Loop(nev);
  obl_test_pn_elastic_c1->cd(4);       //X
  theta4->SetXTitle("cos (#theta)");//X
  theta4->Draw("");                 //X

  obl_test_pn_elastic_c1->Update();   //X
  obl_test_pn_elastic_c1->Modified(); //X
  obl_test_pn_elastic_c1->Print("obl_test_pn_elastic@Elastic_pn_scattering:_polar_angle.png");           //X

}}
