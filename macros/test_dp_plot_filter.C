
{

    makeDistributionManager()->Unpack("../test/pluto_ee_filter_may07.root");

    // makeDistributionManager()->Startup("_filter_debug=1");

gStyle->SetPalette(8,0);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);


TCanvas *c1 = new TCanvas("bla", "ee invmass",800,800);
c1->SetBorderMode(0);
c1->SetFillColor(0);
c1->SetLeftMargin(0.15);
c1->SetBottomMargin(0.15);
c1->SetTopMargin(0.05);
c1->SetRightMargin(0.1);
c1->SetLogy(1);
c1->SetTicky(1);
c1->SetTickx(1);
 c1->SetFrameFillColor(0);

TH2D * frame = new TH2D("frame","frame",1, 0.0,0.8,1,0.00000000001,0.0001);


frame->GetXaxis()->SetTitleOffset(1.1);
frame->GetYaxis()->SetTitleOffset(1.1);
frame->GetXaxis()->SetTitleSize(.06);
frame->GetYaxis()->SetTitleSize(.06);
frame->GetYaxis()->SetTitleFont(42);
frame->GetXaxis()->SetTitleFont(42);
frame->SetXTitle("m_{ee} [GeV/c^{2}]");

frame->SetYTitle("d#sigma/dm_{ee} [b #upoint GeV ^{-1}c^{2} ]");
 c1->cd();
 frame->Draw();

    Double_t kin_max = 0.95;
    Double_t histo_max = 0.8;
//    Int_t histo_bins = 60;
    Int_t histo_bins = 100;
    
    Int_t num_events = 100000;
             num_events = 1000000;

    makeDistributionManager()->Exec("elementary");
    makeDistributionManager()->Exec("dalitz_mod: krivoruchenko");
    makeDistributionManager()->Exec("dalitz_mod: static_br_thresh=0.100 ; flat_generator");
    //     makeDistributionManager()->Enable("vmd");
    makeDistributionManager()->LinkDB();

    makeDistributionManager()->Exec("brems: kaptari");
    makeDistributionManager()->Exec("brems: sum; weighting");
    makeDistributionManager()->Exec("brems: fsi");
    
    PTCrossWeight * my_cross = 
	new PTCrossWeight("n + p_to_n_p_pi0_pi0/tcross",
			  "Cross section for double pi0 production",-1);
    my_cross->SetCrossSection(0.1 * 0.001 * 3); //0.1mb times 3 (see below) 
    makeDistributionManager()->Add(my_cross);
   PHadesParticleSmearer * smear = new PHadesParticleSmearer();
    smear->SetResolutionFactor(3.);

#if 1
//BREMS
    PReaction *my_reaction2 = new PReaction("2.5","d","p","p n dilepton [e+ e-] (p)",NULL,1,0,0,0);
    my_reaction2->AddBulk(smear);
    //Create my histogram:
    TH1F * pp_sum2 = new TH1F ("pp_sum2","pp DiLepton mass (coherent sum)",histo_bins,0.1,histo_max);
    TH1F * pp_sum = new TH1F ("pp_sum","pp DiLepton mass (coherent sum)",histo_bins,0.,histo_max);
    TH1F * dummy2 = new TH1F ("dummy2","pp DiLepton mass (coherent sum)",10,-1,1);
//    pp_sum2->Sumw2();
    dummy2->Sumw2();
    //Create the container of the histogram list
    PProjector *m2 = new PProjector(); 
    //Dilepton mass
    m2->AddCommand("opang = [e+]->Angle([e-])");
    m2->AddHistogram(pp_sum2,"if opang > (9./180.)*TMath::Pi(); _x=([e+]+[e-])->M()");
    m2->AddHistogram(pp_sum,"if opang > (9./180.)*TMath::Pi(); _x=([e+]+[e-])->M()");
    m2->AddHistogram(dummy2,"_x=cos(([e+]+[e-])->Theta())");
    my_reaction2->AddBulk(m2);
    my_reaction2->Loop(num_events);
    PUtils::correct(pp_sum2); //correct for number of used bins
    PUtils::correct(dummy2);
    
    pp_sum2->SetMarkerStyle(22);
    pp_sum2->SetMarkerColor(4);
    pp_sum2->SetLineWidth(2);
//     pp_sum2->SetLineColor(4);
//     pp_sum2->SetLineStyle(9);
   


#endif

    
    
    PReaction *my_reaction3 = new PReaction("2.5","d","p","n D+ [p dilepton [e+ e-]] (p)",NULL,1,0,0,0);
    my_reaction3->AddBulk(smear);
//    my_reaction3->AddReaction("p n dilepton [e+ e-] p");
    my_reaction3->AddReaction("p D0 [n dilepton [e+ e-]] (p)");

    //Create my histogram:
    TH1F * pp_sum3 = new TH1F ("pp_sum3","pp DiLepton mass (coherent sum)",histo_bins,0.,histo_max);
    TH1F * dummy3 = new TH1F ("dummy3","pp DiLepton mass (coherent sum)",100,0.8,1.6);
    
 
//    pp_sum3->Sumw2();
    dummy3->Sumw2();
 
    //Create the container of the histogram list
    PProjector *m3 = new PProjector(); 
    //Dilepton mass
    m3->AddCommand("opang = [e+]->Angle([e-])");
    m3->AddHistogram(pp_sum3,"if opang > (9./180.)*TMath::Pi(); _x=([e+]+[e-])->M()");
//    m3->AddHistogram(pp_sum,"_x=([e+]+[e-])->M()");
 
    my_reaction3->AddBulk(m3);
    
    my_reaction3->Loop(num_events);
//    my_reaction3->Loop(1);
    my_reaction3->Print();//0.00564527

    PUtils::correct(pp_sum3); //correct for number of used bins

    makeDistributionManager()->Enable("vmd");
     PReaction *my_reaction3a = new PReaction("2.5","d","p","n D+ [p dilepton [e+ e-]] (p)",NULL,1,0,0,0);     
//    my_reaction3->AddReaction("p n dilepton [e+ e-] p");
    my_reaction3a->AddReaction("p D0 [n dilepton [e+ e-]] (p)");
    my_reaction3a->AddBulk(smear);
        //Create my histogram:
    TH1F * pp_sum3a = new TH1F ("pp_sum3a","pp DiLepton mass (coherent sum)",histo_bins,0.,histo_max);
    TH1F * dummy3a = new TH1F ("dummy3a","pp DiLepton mass (coherent sum)",100,0.8,1.6);
    dummy3a->Sumw2();
 
    //Create the container of the histogram list
    PProjector *m3a = new PProjector(); 
    //Dilepton mass
    m3a->AddCommand("opang = [e+]->Angle([e-])");
    m3a->AddHistogram(pp_sum3a,"if opang > (9./180.)*TMath::Pi(); _x=([e+]+[e-])->M()");
    my_reaction3a->AddBulk(m3a);
    
    my_reaction3a->Loop(num_events);



    //pi0 production
    PReaction *my_reaction4 = new PReaction("2.5","d","p",
					    "n D+ [p pi0 [g dilepton [e+ e-]]] (p)",NULL,1,0,0,0);
    my_reaction4->AddReaction("p D0 [n pi0 [g dilepton [e+ e-]]] (p)");
    //double pi0 production
    my_reaction4->AddReaction("p n pi0 pi0 [g dilepton [e+ e-]] (p)");
    //*3 because of p p pi0 pi-
    my_reaction4->AddBulk(smear);
    
    TH1D * pp_sum4 = new TH1D ("pp_sum4","pp DiLepton mass (coherent sum)",histo_bins,0.001,0.14);
    //   pp_sum4->Sumw2();
    
       //Create the container of the histogram list
    PProjector *m4 = new PProjector(); 
    //Dilepton mass
    m4->AddCommand("opang = [e+]->Angle([e-])");
    m4->AddHistogram(pp_sum4,"if opang > (9./180.)*TMath::Pi(); _x=([e+]+[e-])->M()");
    m4->AddHistogram(pp_sum,"if opang > (9./180.)*TMath::Pi(); _x=([e+]+[e-])->M()");

 
    my_reaction4->AddBulk(m4);
  
    my_reaction4->Preheating(10);
    my_reaction4->Loop(num_events);
//    my_reaction4->Loop(1);
    my_reaction4->Print();//0.00564527

    PUtils::correct(pp_sum4); //correct for number of used bins



   //eta production
    PReaction *my_reaction5 = new PReaction("2.5","d","p",
					    "n p eta [g dilepton [e+ e-]] (p)",NULL,1,0,0,0);
    my_reaction5->AddReaction("d eta [g dilepton [e+ e-]] (p)");
    my_reaction5->AddBulk(smear);
    
    TH1D * pp_sum5 = new TH1D ("pp_sum5","pp DiLepton mass (coherent sum)",histo_bins,0.,histo_max);
    //   pp_sum4->Sumw2();
    TH1F * pp_sum5a = new TH1F ("pp_sum5a","pp DiLepton mass (coherent sum)",histo_bins,0.1,histo_max);
    
       //Create the container of the histogram list
    PProjector *m5 = new PProjector(); 
    //Dilepton mass
    m5->AddCommand("opang = [e+]->Angle([e-])");
    m5->AddHistogram(pp_sum5,"if opang > (9./180.)*TMath::Pi(); _x=([e+]+[e-])->M()");
    m5->AddHistogram(pp_sum,"if opang > (9./180.)*TMath::Pi(); _x=([e+]+[e-])->M()");
    m5->AddHistogram(pp_sum5a,"_x=([e+]+[e-])->M()");
 
    my_reaction5->AddBulk(m5);
  
    my_reaction5->Preheating(10);
    my_reaction5->Loop(num_events);
//    my_reaction4->Loop(1);
    my_reaction5->Print();//0.00564527

    PUtils::correct(pp_sum5); //correct for number of used bins
    PUtils::correct(pp_sum5a); //correct for number of used bins
    
    pp_sum3->Add(pp_sum5);
    pp_sum2->Add(pp_sum5a); //diff. axis
    pp_sum2->Draw("Csame");

    pp_sum3->SetMarkerStyle(20);
    pp_sum3->SetMarkerColor(3);
    pp_sum3->SetLineColor(4);
    pp_sum3->SetLineStyle(9);
    pp_sum3->SetLineWidth(2);
    pp_sum3->SetFillColor(10);
    PUtils::correct(pp_sum3a); //correct for number of used bins
    pp_sum3a->Add(pp_sum5);

    pp_sum3a->SetMarkerStyle(20);
    pp_sum3a->Draw("Csame");
    pp_sum3a->SetMarkerColor(3);
    pp_sum3a->SetLineColor(4);
    pp_sum3a->SetLineStyle(7);
    pp_sum3a->SetLineWidth(2);
    pp_sum3a->SetFillColor(17);
    
    pp_sum3->Draw("Csame");

    pp_sum2->Draw("Csame");

    pp_sum4->SetMarkerStyle(27);
    pp_sum4->Draw("Csame");
    pp_sum4->SetMarkerColor(2);
    pp_sum4->SetLineColor(2);
    pp_sum4->SetLineStyle(2);
    pp_sum4->SetLineWidth(2);

    pp_sum5->Draw("Csame");
    pp_sum5->SetMarkerColor(2);
    pp_sum5->SetLineColor(8);
    pp_sum5->SetLineStyle(10);
    pp_sum5->SetLineWidth(2);
    
    
    PUtils::correct(pp_sum);
//    pp_sum->Draw("Csame");
    pp_sum->SetLineWidth(2);

 //OBE:
 TLine *l1 = new TLine(0.3,3.16*1e-04,0.4,3.16*1e-04);
 l1->SetLineWidth(2);
 l1->Draw("same");
 TLatex *x1 = new TLatex(0.45,3.16*1e-04*0.78,"OBE (K&K) + #eta");
 x1->SetTextSize(0.04);
 x1->Draw("same");

 //VMD
 TLine *l2 = new TLine(0.3,1*1e-04,0.4,1*1e-04);
 l2->Draw("same");
 l2->SetLineColor(4);
 l2->SetLineStyle(7);
 l2->SetLineWidth(2);
 TLatex *x2 = new TLatex(0.45,1*1e-04*0.78,"#Delta (2-comp.) + #eta");
 x2->SetTextSize(0.04);
 x2->Draw("same");

 //QED:
 TLine *l3 = new TLine(0.3,3.16e-05,0.4,3.16e-05);
 l3->SetLineWidth(2);
 l3->SetLineStyle(9);
 l3->SetLineColor(4);
 l3->Draw("same");
 TLatex *x3 = new TLatex(0.45,3.16*1e-5*0.78,"#Delta (const.) + #eta");
 x3->SetTextSize(0.04);
 x3->Draw("same");

 //PI
 TLine *l4 = new TLine(0.3,1*1e-05,0.4,1*1e-05);
 l4->Draw("same");
 l4->SetLineColor(2);
 l4->SetLineStyle(2);
 l4->SetLineWidth(2);   
 TLatex *x4 = new TLatex(0.45,1*1e-05*0.78,"#pi^{0}");
 x4->SetTextSize(0.04);
 x4->Draw("same");

  //ETA
 TLine *l4 = new TLine(0.3,3.16*1.0e-06,0.4,3.16*1.0e-06);
 l4->Draw("same");
 l4->SetLineColor(8);
 l4->SetLineStyle(10);
 l4->SetLineWidth(2);   
 TLatex *x4 = new TLatex(0.45,3.16*1e-6*0.78,"#eta");
 x4->SetTextSize(0.04);
 x4->Draw("same");

// Labels
 TLatex *x = new TLatex(0.3,1.8e-03,"dp #rightarrow np e^{+} e^{-} p_{spec}");
 x->Draw("same");

Double_t xAxis[33] = {0/1000., 5/1000., 10/1000., 15/1000., 20/1000., 25/1000., 30/1000., 35/1000., 
		      40/1000., 45/1000., 50/1000., 55/1000., 60/1000., 65/1000., 70/1000., 75/1000., 
		      80/1000., 90/1000., 100/1000., 120/1000., 160/1000., 215/1000., 270/1000., 
		      325/1000., 380/1000., 435/1000., 500/1000., 550/1000., 600/1000., 700/1000., 800/1000., 
		      900/1000., 1000/1000.}; 

 Double_t scale = 1e-6;
 TH1 *hmass_cut33_back_1_sig_norm = new TH1F("hmass_cut33_back_1_sig_norm","",32, xAxis);
   hmass_cut33_back_1_sig_norm->SetBinContent(2,-0.0453748*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(3,6.15357*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(4,22.1315*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(5,31.5663*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(6,33.6503*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(7,31.0997*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(8,26.7148*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(9,23.3826*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(10,21.0868*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(11,18.1938*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(12,15.6738*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(13,14.4015*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(14,11.8571*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(15,9.22774*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(16,7.78852*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(17,6.04399*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(18,2.75297*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(19,1.14051*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(20,0.336649*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(21,0.242977*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(22,0.15181*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(23,0.092777*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(24,0.0712852*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(25,0.0576961*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(26,0.0344135*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(27,0.0132155*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(28,0.000841551*scale);
   hmass_cut33_back_1_sig_norm->SetBinContent(29,0.000131521*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(2,0.0809017*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(3,0.393711*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(4,0.64714*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(5,0.718377*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(6,0.718287*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(7,0.671004*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(8,0.634903*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(9,0.61262*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(10,0.594512*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(11,0.555559*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(12,0.537638*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(13,0.556242*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(14,0.534226*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(15,0.503942*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(16,0.456658*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(17,0.315744*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(18,0.247617*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(19,0.146855*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(20,0.0474574*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(21,0.0222295*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(22,0.0131332*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(23,0.00979851*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(24,0.00680336*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(25,0.00585345*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(26,0.00385863*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(27,0.00262061*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(28,0.000486992*scale);
   hmass_cut33_back_1_sig_norm->SetBinError(29,0.000131521*scale);
   hmass_cut33_back_1_sig_norm->SetEntries(42447);
   hmass_cut33_back_1_sig_norm->SetLineWidth(2);
   hmass_cut33_back_1_sig_norm->SetMarkerStyle(8);
   hmass_cut33_back_1_sig_norm->SetMarkerSize(1.8);
   hmass_cut33_back_1_sig_norm->Draw("same");



c1->Update();
c1->Modified();

c1->GetFrame()->SetBorderMode(0);
c1->Update();
c1->Modified();
c1->cd();
c1->SetSelected(c1);
c1->Print("dp_cocktail_comparison_tmp.eps");

}

