
{


    TCanvas *obl_test_rho_mass_shape_c1 = new        //X
	TCanvas ("obl_test_rho_mass_shape_c1","rho");//X
    obl_test_rho_mass_shape_c1->SetBorderMode(0);    //X
    obl_test_rho_mass_shape_c1->SetFillColor(0);     //X
    obl_test_rho_mass_shape_c1->SetLogy(1);          //X
    obl_test_rho_mass_shape_c1->SetTicky(1);         //X
    obl_test_rho_mass_shape_c1->SetTickx(1);         //X
    obl_test_rho_mass_shape_c1->SetFrameFillColor(0);//X
    

    PReaction my_reaction("3.5","p","p","p p rho0");
    TH1F * histo1 = new TH1F ("histo1","dilepton mass",100,0.2,1.3);
    histo1->Sumw2();
    histo1->SetMinimum(0.000000001);
    my_reaction.Do(histo1,"_x =  ([rho0])->M()");

    my_reaction.Print();   //The "Print()" statement is optional
    my_reaction.Loop(10000);
    histo1->Draw("e1");
    
    PReaction my_reaction2("3.5","p","p","p p rho0 [e+ e-]");
    TH1F * histo2 = new TH1F ("histo2","dilepton mass",100,0.2,1.3);
    histo2->Sumw2();
    my_reaction2.Do(histo2,"_x =  ([rho0])->M()");

    my_reaction2.Print();   //The "Print()" statement is optional
    my_reaction2.Loop(10000);
    histo2->Draw("samee1");

    histo1->SetXTitle("M [GeV/c^{2}]");//X
    obl_test_rho_mass_shape_c1->Update();   //X
    obl_test_rho_mass_shape_c1->Modified(); //X
    obl_test_rho_mass_shape_c1->Print("obl_test_rho_mass_shape@Rho_mass_shape.png");           //X




}
