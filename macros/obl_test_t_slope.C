//TITLE New t slope with PAnyDistribution
//This macro tests the PAnyDistribution and batch functions

{

    gStyle->SetOptStat(0); //X
    
    PAnyDistribution* decay = new PAnyDistribution("t_slope","A function to add a new t-slope");
    decay->Add("q,     parent");
    decay->Add("p,     daughter");
    decay->Add("rho0,  daughter");

    TH1F * cache  = new TH1F ("cache","Rho0 t cache",400,-4.0,.0);
    
    decay->AddEquation(cache,"beam = _parent->GetBeam(); beam->Boost(_parent) ; t1 = (beam - [rho0])->Mag2(); _x = t1;");

    decay->AddEquation("_f = exp( 7.1*t1 );");

    decay->SetMaxEnhancementFactor(10);

    makeDistributionManager()->Add(decay);

    //Construct the reaction, as usual:
    PReaction my_reaction(3.5,"g","p","p rho0 [pi+ pi-]");

    TH1F * histo1 = new TH1F ("histo1","rho0 t",100,-1,0);
    histo1->Sumw2();

    my_reaction.Do(histo1,"beam2 = ([g+p]->GetBeam()); t12 = (beam2 - [rho0])->Mag2(); _x = t12;");
    my_reaction.Print();   //The "Print()" statement is optional

    //Make a dummy loop to fill the AnyDistribution with some statistics:
    my_reaction.Preheating(100);

    //Final loop
    my_reaction.Loop(500); 

    TCanvas *obl_test_t_slope_c1 = new //X
	TCanvas ("obl_test_t_slope_c1","Multiplicity");//X
    obl_test_t_slope_c1->SetBorderMode(0);    //X
    obl_test_t_slope_c1->SetFillColor(0);     //X
    obl_test_t_slope_c1->SetTicky(1);         //X
    obl_test_t_slope_c1->SetTickx(1);         //X
    obl_test_t_slope_c1->SetFrameFillColor(0);//X
    obl_test_t_slope_c1->SetLogy(); //X
    histo1->SetXTitle("t");//X
    histo1->Draw();  //X Mean should be the selected mean mult.
    obl_test_t_slope_c1->Update();   //X
    obl_test_t_slope_c1->Modified(); //X
    obl_test_t_slope_c1->Print("obl_test_t_slope@T_slope.png");           //X
}
