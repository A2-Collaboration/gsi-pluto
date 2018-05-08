//TITLE Thermal pi0s
//This macro tests the thermal pi0 (in the so-called "multi-mode") 

{
    //X
    gStyle->SetPalette(1,0); //X
    
    Float_t T1    = 0.080;   // temperature in GeV (for pion 2-component spectra)
    Float_t T2    = 0.080;   // temperature in GeV (assume this for thermalized source)
    Float_t frac  = 1.0;     // fraction of pion low-T component (from Jehad's fit to QMD)
    Float_t blast = 0.0;     // radial expansion velocity

    Float_t A2    = 1.0;     // polar distribution (from KaoS pion data)
    Float_t A4    = 0.0;
    Float_t v1    = 0.0;     // side flow
    Float_t v2    = 0.0;     // elliptic flow
    Float_t Eb    = 3.5;     // beam energy in AGeV

    PFireball *source1=new PFireball("pi0",Eb,T1,T2,frac,blast,A2,A4,v1,v2);

    Float_t Mpi0   = 1.;      // pi0 multiplicity... just a test

    Float_t enhancepi=3.;            //enhancement factor for pi
    
    Mpi0 *=enhancepi;

    source1->SetW(1.0/enhancepi);

    PChannel *c1 = source1->makeChannel(20, Mpi0);          // pi0 decay directly
    source1->Print();

    PReaction *r=new PReaction(&c1,"pA",1);

    r->allParticles();        //Must be used, otherwise the PProjector does not consider all particles
    PPlutoBulkDecay *pl = new PPlutoBulkDecay();
    pl->SetRecursiveMode(1);  //Let also the products decay
    pl->SetTauMax(0.001);     //maxTau in ns
    r->AddBulk(pl);

    //Test histos
    TH1F * histo1 = new TH1F ("histo1","pi0 mult",21,-0.5,20.5);

    r->Do("counter = 0;"); 
    r->Do("foreach(pi0); counter = counter +1");
    r->Do(histo1,"_x = counter");

    TH2F * histo2 = new TH2F ("histo2","Rap. vs. Pt",50,-1.5,3.5,50,0,1);
    r->Do(histo2,"foreach(pi0); _x = [pi0]->Rapidity(); _y=[pi0]->Pt(); ");

    TH1F * histo3 = new TH1F ("histo3","ee invariant mass",100,0,0.5);
    r->Do(histo3,"foreach(dilepton); _x = [dilepton]->M(); ");

    TH1F * histo4 = new TH1F ("histo4","g mom ",100,0,1.5);//X
    r->Do("counter = 0;"); //X
    r->Do("foreach(g);  counter = counter +1");//X
    
    r->Print();
    r->loop(20000);

    TCanvas *obl_test_thermal_pi0_c1 = new 
	TCanvas ("obl_test_thermal_pi0_c1","Multiplicity");
    obl_test_thermal_pi0_c1->SetBorderMode(0);    //X
    obl_test_thermal_pi0_c1->SetFillColor(0);     //X
    obl_test_thermal_pi0_c1->SetTicky(1);         //X
    obl_test_thermal_pi0_c1->SetTickx(1);         //X
    obl_test_thermal_pi0_c1->SetFrameFillColor(0);//X
    histo1->SetXTitle("pi0 multiplicity");//X
    histo1->Draw();  //X Mean should be the selected mean mult.

    TCanvas *obl_test_thermal_pi0_c2 = new 
	TCanvas ("obl_test_thermal_pi0_c2","Rap. vs. Pt");
    obl_test_thermal_pi0_c2->SetBorderMode(0);    //X
    obl_test_thermal_pi0_c2->SetFillColor(0);     //X
    obl_test_thermal_pi0_c2->SetTicky(1);         //X
    obl_test_thermal_pi0_c2->SetTickx(1);         //X
    obl_test_thermal_pi0_c2->SetFrameFillColor(0);//X
    histo2->SetXTitle("Rapidity");//X
    histo2->SetYTitle("Pt");//X
    histo2->Draw("colz");//X

    TCanvas *obl_test_thermal_pi0_c3 = new //X
	TCanvas ("obl_test_thermal_pi0_c3","EE Mass");//X
    obl_test_thermal_pi0_c3->SetBorderMode(0);    //X
    obl_test_thermal_pi0_c3->SetFillColor(0);     //X
    obl_test_thermal_pi0_c3->SetLogy(1);          //X
    obl_test_thermal_pi0_c3->SetTicky(1);         //X
    obl_test_thermal_pi0_c3->SetTickx(1);         //X
    obl_test_thermal_pi0_c3->SetFrameFillColor(0);//X
    histo3->SetXTitle("M_{ee}");//X
    histo3->Draw();//X
    obl_test_thermal_pi0_c1->Update();   //X
    obl_test_thermal_pi0_c1->Modified(); //X
    obl_test_thermal_pi0_c1->Print("obl_test_thermal_pi0@Thermal_pi0:_pi0_multiplicity.png");           //X
    obl_test_thermal_pi0_c2->Update();   //X
    obl_test_thermal_pi0_c2->Modified(); //X
    obl_test_thermal_pi0_c2->Print("obl_test_thermal_pi0@Thermal_pi0:_pi0_pt_vs_rapidity.png");           //X
    obl_test_thermal_pi0_c3->Update();   //X
    obl_test_thermal_pi0_c3->Modified(); //X
    obl_test_thermal_pi0_c3->Print("obl_test_thermal_pi0@Thermal_pi0:_pi0_dalitz_decay_dielectron_mass.png");           //X
}
