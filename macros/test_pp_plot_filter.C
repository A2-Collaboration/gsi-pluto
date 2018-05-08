
{

    makeDistributionManager()->Unpack("../test/pluto_ee_filter_apr06.root");
    
    //    makeDistributionManager()->Startup("_filter_debug=1");

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

    //TH2D * frame = new TH2D("frame","frame",1, 0.0,0.8,1,0.00000000001,0.01);
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
    Int_t histo_bins = 50;
    
    Int_t num_events = 100000;
    //num_events = 1000000;

    makeDistributionManager()->Exec("elementary");
    makeDistributionManager()->Exec("dalitz_mod: krivoruchenko");
    makeDistributionManager()->Exec("dalitz_mod: static_br_thresh=0.100 ; flat_generator");

    makeDistributionManager()->LinkDB();

    makeDistributionManager()->Exec("brems: kaptari");
    makeDistributionManager()->Exec("brems: sum; weighting; fsi");

    PTCrossWeight * my_cross = 
	new PTCrossWeight("p + p_to_p_p_pi0_pi0/tcross",
			  "Cross section for double pi0 production",-1);
    my_cross->SetCrossSection(0.1 * 0.001 * 2); //0.2mb times 2 (see below) 
    makeDistributionManager()->Add(my_cross);

    PHadesParticleSmearer * smear = new PHadesParticleSmearer();
    smear->SetResolutionFactor(1.);


    //BREMS
    PReaction *my_reaction2 = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);
    my_reaction2->AddBulk(smear);
    
    //Create my histogram:
    TH1F * pp_sum2 = new TH1F ("pp_sum","pp DiLepton mass (coherent sum)",histo_bins,0.1,histo_max);
    TH1F * dummy2 = new TH1F ("dummy2","pp DiLepton mass (coherent sum)",10,-1,1);
    //    pp_sum2->Sumw2();
    dummy2->Sumw2();
    //Create the container of the histogram list
    PProjector *m2 = new PProjector(); 
    //Dilepton mass
    m2->AddCommand("opang = [e+]->Angle([e-])");
    //    m2->AddHistogram(pp_sum2,"if opang > (9./180.)*TMath::Pi(); _x=([e+] + [e-])->M()");
    m2->AddHistogram(pp_sum2,"_x=([e+] + [e-])->M()");
    m2->AddHistogram(dummy2,"_x=cos(([e+] + [e-])->Theta())");
    my_reaction2->AddBulk(m2);
    my_reaction2->Loop(num_events);
    PUtils::correct(pp_sum2); //correct for number of used bins
    PUtils::correct(dummy2);
  
    PBremsstrahlung * brems =  makeDistributionManager()->GetDistribution("p + p_brems_p_p_dilepton");
    brems->SetNeutron(-1);
    brems->SetP2E(1.25);
  
    pp_sum2->SetMarkerStyle(22);
    pp_sum2->SetMarkerColor(4);
    pp_sum2->SetLineWidth(2);
    //    pp_sum2->SetLineColor(4);
    // pp_sum2->SetLineStyle(9);


    //    brems->Draw("same");
    brems->SetLineWidth(1);



    
    makeDistributionManager()->Exec("brems : elastic");

    
    PReaction *my_reaction3 = new PReaction("1.25","p","p","p D+ [p dilepton [e+ e-]]",NULL,1,0,0,0);
    //    my_reaction3->AddReaction("p p dilepton [e+ e-]");
    //    PReaction *my_reaction3 = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);
    //    my_reaction3->AddReaction("p D+ [p dilepton [e+ e-]]");

    //Create my histogram:
    TH1F * pp_sum3 = new TH1F ("pp_sum3","pp DiLepton mass (coherent sum)",histo_bins,0.,histo_max);
    TH1F * dummy3 = new TH1F ("dummy3","pp DiLepton mass (coherent sum)",100,0.8,1.6);
    my_reaction3->AddBulk(smear);
 
    //    pp_sum3->Sumw2();
    dummy3->Sumw2();
 
    //Create the container of the histogram list
    PProjector *m3 = new PProjector(); 
    //Dilepton mass
    m3->AddCommand("opang = [e+]->Angle([e-])");
    //    m3->AddHistogram(pp_sum3,"if opang > (9./180.)*TMath::Pi(); _x=([e+] + [e-])->M()");
    m3->AddHistogram(pp_sum3,"_x=([e+] + [e-])->M()");
 
    my_reaction3->AddBulk(m3);
    
    my_reaction3->Loop(num_events);
    //    my_reaction3->Loop(1);
    my_reaction3->Print();//0.00564527


    makeDistributionManager()->Enable("vmd");
    PReaction *my_reaction3a = new PReaction("1.25","p","p","p D+ [p dilepton [e+ e-]]",NULL,1,0,0,0);
    //    my_reaction3->AddReaction("p p dilepton [e+ e-]");
    //    PReaction *my_reaction3 = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);
    //    my_reaction3->AddReaction("p D+ [p dilepton [e+ e-]]");

    //Create my histogram:
    TH1F * pp_sum3a = new TH1F ("pp_sum3a","pp DiLepton mass (coherent sum)",histo_bins,0.,histo_max);
    TH1F * dummy3a = new TH1F ("dummy3a","pp DiLepton mass (coherent sum)",100,0.8,1.6);
    dummy3a->Sumw2();
    my_reaction3a->AddBulk(smear);

    //Create the container of the histogram list
    PProjector *m3a = new PProjector(); 
    //Dilepton mass
    m3a->AddCommand("opang = [e+]->Angle([e-])");
    //    m3a->AddHistogram(pp_sum3a,"if opang > (9./180.)*TMath::Pi(); _x=([e+] + [e-])->M()");
    m3a->AddHistogram(pp_sum3a,"_x=([e+] + [e-])->M()");
    my_reaction3a->AddBulk(m3a);
    
    my_reaction3a->Loop(num_events);


    PUtils::correct(pp_sum3); //correct for number of used bins

    pp_sum3->SetMarkerStyle(20);
    pp_sum3->SetMarkerColor(3);
    pp_sum3->SetLineColor(4);
    pp_sum3->SetLineStyle(9);
    pp_sum3->SetLineWidth(2);
    pp_sum3->SetFillColor(10);
    PUtils::correct(pp_sum3a); //correct for number of used bins

    pp_sum3a->SetMarkerStyle(20);
    pp_sum3a->Draw("Csame");
    pp_sum3a->SetMarkerColor(3);
    pp_sum3a->SetLineColor(4);
    pp_sum3a->SetLineStyle(7);
    pp_sum3a->SetLineWidth(2);
    pp_sum3a->SetFillColor(17);
    
    pp_sum3->Draw("Csame");

    pp_sum2->Draw("Csame");

    //pi0 production
    PReaction *my_reaction4 = new PReaction("1.25","p","p","p D+ [p pi0 [g dilepton [e+ e-]]]","test_filter_pi0",1,0,0,0);
    //double pi0 production
    //my_reaction4->AddReaction("p p pi0 pi0 [g dilepton [e+ e-]]");
    my_reaction4->AddBulk(smear);

    //TH1D * pp_sum4 = new TH1D ("pp_sum4","pp DiLepton mass (coherent sum)",histo_bins,0.001,0.14);
    TH1D * pp_sum4 = new TH1D ("pp_sum4","pp DiLepton mass (coherent sum)",20,0.001,0.14);
   //  pp_sum4->Sumw2();
    
    //Create the container of the histogram list
    PProjector *m4 = new PProjector(); 
    //Dilepton mass
    m4->AddCommand("opang = [e+]->Angle([e-])");
    //m4->AddHistogram(pp_sum4,"if opang > (9./180.)*TMath::Pi(); _x=([e+] + [e-])->M()");
    m4->AddHistogram(pp_sum4," _x=([e+] + [e-])->M()");

 
    my_reaction4->AddBulk(m4);
  
    my_reaction4->Preheating(10);
    my_reaction4->Loop(num_events);
    //my_reaction4->Loop(10);
    my_reaction4->Print();//0.00564527

    PUtils::correct(pp_sum4); //correct for number of used bins

    //    pp_sum4->SetMarkerStyle(27);
    pp_sum4->Draw("Csame");
    pp_sum4->SetMarkerColor(2);
    pp_sum4->SetLineColor(2);
    //    pp_sum4->SetLineStyle(10);
    pp_sum4->SetLineStyle(2);
    pp_sum4->SetLineWidth(2);


    //OBE:
    TLine *l1 = new TLine(0.3,3.16*1.0e-06,0.4,3.16*1.0e-06);
    l1->SetLineWidth(2);
    l1->Draw("same");
    TLatex *x1 = new TLatex(0.45,3.16*1e-06*0.78,"OBE (K&K)");
    x1->SetTextSize(0.04);
    x1->Draw("same");

    //VMD
    TLine *l2 = new TLine(0.3,1.*1e-06,0.4,1.*1e-06);
    l2->Draw("same");
    l2->SetLineColor(4);
    l2->SetLineStyle(7);
    l2->SetLineWidth(2);
    TLatex *x2 = new TLatex(0.45,1.*1e-06*0.78,"#Delta (2-comp.)");
    x2->SetTextSize(0.04);
    x2->Draw("same");

    //QED:
    TLine *l3 = new TLine(0.3,3.16*1.0e-07,0.4,3.16*1.0e-07);
    l3->SetLineWidth(2);
    l3->SetLineStyle(9);
    l3->SetLineColor(4);
    l3->Draw("same");
    TLatex *x3 = new TLatex(0.45,3.16*1e-7*0.78,"#Delta (const.)");
    x3->SetTextSize(0.04);
    x3->Draw("same");

    //PI
    TLine *l4 = new TLine(0.3,1.*1e-07,0.4,1.*1e-07);
    l4->Draw("same");
    l4->SetLineColor(2);
    l4->SetLineStyle(2);
    l4->SetLineWidth(2);   
    TLatex *x4 = new TLatex(0.45,1.*1e-07*0.78,"#pi^{0}");
    x4->SetTextSize(0.04);
    x4->Draw("same");


    // Labels
    TLatex *x = new TLatex(0.3,1.8e-03,"pp #rightarrow pp e^{+} e^{-} ");
    //    x->Draw("same");

    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetTitleOffset(1.1);
    frame->GetXaxis()->SetTitleSize(.06);
    frame->GetYaxis()->SetTitleSize(.06);
    frame->GetYaxis()->SetTitleFont(42);
    frame->GetXaxis()->SetTitleFont(42);
    frame->SetXTitle("m_{ee} [GeV/c^{2}]");



       Double_t xAxis[33] = {0/1000., 5/1000., 10/1000., 15/1000., 20/1000., 25/1000., 30/1000., 
			     35/1000., 40/1000., 45/1000., 50/1000., 55/1000., 60/1000., 65/1000., 
			     70/1000., 75/1000., 80/1000., 90/1000., 100/1000., 120/1000., 160/1000., 
			     215/1000., 270/1000., 325/1000., 380/1000., 435/1000., 500/1000., 550/1000., 
			     600/1000., 700/1000., 800/1000., 900/1000., 1000/1000.}; 
       Double_t scale = 1e-6;
   TH1 *hmass_cut20_back_1_sig_norm = new TH1F("hmass_cut20_back_1_sig_norm","",32, xAxis);
   hmass_cut20_back_1_sig_norm->SetBinContent(2,0.0197184*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(3,3.08897*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(4,12.323*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(5,16.1652*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(6,17.6273*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(7,16.1717*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(8,14.7628*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(9,12.7399*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(10,12.0236*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(11,10.1641*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(12,8.89723*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(13,8.35003*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(14,6.9851*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(15,5.70519*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(16,4.91632*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(17,3.81329*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(18,2.13384*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(19,0.877159*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(20,0.0695648*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(21,0.036028*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(22,0.0249973*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(23,0.00915055*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(24,0.00387053*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(25,0.00144041*scale);
   hmass_cut20_back_1_sig_norm->SetBinContent(26,0.000655378*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(2,0.0536937*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(3,0.239983*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(4,0.428645*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(5,0.457265*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(6,0.461813*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(7,0.399751*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(8,0.372479*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(9,0.353626*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(10,0.343067*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(11,0.331567*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(12,0.318666*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(13,0.350142*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(14,0.293553*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(15,0.277189*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(16,0.279299*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(17,0.180449*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(18,0.139262*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(19,0.068018*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(20,0.0201587*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(21,0.00838737*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(22,0.00440742*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(23,0.00219634*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(24,0.00116418*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(25,0.000641552*scale);
   hmass_cut20_back_1_sig_norm->SetBinError(26,0.000331567*scale);
   hmass_cut20_back_1_sig_norm->SetEntries(25);
   hmass_cut20_back_1_sig_norm->SetMarkerStyle(20);
   hmass_cut20_back_1_sig_norm->SetMarkerSize(1.8);
   hmass_cut20_back_1_sig_norm->Draw("same");


  TGraphAsymmErrors *grae = new TGraphAsymmErrors(32);
   grae->SetName("Graph");
   grae->SetTitle("Graph");
   grae->SetFillColor(17);
   grae->SetLineColor(2);
   grae->SetLineWidth(2);
   grae->SetMarkerSize(0.1);
   grae->SetPoint(0,2.5,0);
   grae->SetPointError(0,0,0,0,0);
   grae->SetPoint(1,7.5,0.0197185*scale);
   grae->SetPointError(1,0,0,0.00440918,0.00440918);
   grae->SetPoint(2,12.5,3.08898*scale);
   grae->SetPointError(2,0,0,0.690717,0.690717);
   grae->SetPoint(3,17.5,12.3231*scale);
   grae->SetPointError(3,0,0,2.75552,2.75552);
   grae->SetPoint(4,22.5,16.1652*scale);
   grae->SetPointError(4,0,0,3.61464,3.61464);
   grae->SetPoint(5,27.5,17.6273*scale);
   grae->SetPointError(5,0,0,3.94158,3.94158);
   grae->SetPoint(6,32.5,16.1717*scale);
   grae->SetPointError(6,0,0,3.6161,3.6161);
   grae->SetPoint(7,37.5,14.7628*scale);
   grae->SetPointError(7,0,0,3.30107,3.30107);
   grae->SetPoint(8,42.5,12.7399*scale);
   grae->SetPointError(8,0,0,2.84873,2.84873);
   grae->SetPoint(9,47.5,12.0236*scale);
   grae->SetPointError(9,0,0,2.68855,2.68855);
   grae->SetPoint(10,52.5,10.1641*scale);
   grae->SetPointError(10,0,0,2.27276,2.27276);
   grae->SetPoint(11,57.5,8.89723*scale);
   grae->SetPointError(11,0,0,1.98948,1.98948);
   grae->SetPoint(12,62.5,8.35002*scale);
   grae->SetPointError(12,0,0,1.86712,1.86712);
   grae->SetPoint(13,67.5,6.9851*scale);
   grae->SetPointError(13,0,0,1.56192,1.56192);
   grae->SetPoint(14,72.5,5.7052*scale);
   grae->SetPointError(14,0,0,1.27572,1.27572);
   grae->SetPoint(15,77.5,4.91633*scale);
   grae->SetPointError(15,0,0,1.09932,1.09932);
   grae->SetPoint(16,85,3.81329*scale);
   grae->SetPointError(16,0,0,0.852678,0.852678);
   grae->SetPoint(17,95,2.13384*scale);
   grae->SetPointError(17,0,0,0.477141,0.477141);
   grae->SetPoint(18,110,0.87716*scale);
   grae->SetPointError(18,0,0,0.196139,0.196139);
   grae->SetPoint(19,140,0.0695647*scale);
   grae->SetPointError(19,0,0,0.0155551,0.0155551);
   grae->SetPoint(20,187.5,0.036028*scale);
   grae->SetPointError(20,0,0,0.0080561,0.0080561);
   grae->SetPoint(21,242.5,0.0249974*scale);
   grae->SetPointError(21,0,0,0.00558959,0.00558959);
   grae->SetPoint(22,297.5,0.00915055*scale);
   grae->SetPointError(22,0,0,0.00204613,0.00204613);
   grae->SetPoint(23,352.5,0.00387052*scale);
   grae->SetPointError(23,0,0,0.000865475,0.000865475);
   grae->SetPoint(24,407.5,0.00144041*scale);
   grae->SetPointError(24,0,0,0.000322086,0.000322086);
   grae->SetPoint(25,467.5,0.000655379*scale);
   grae->SetPointError(25,0,0,0.000146547,0.000146547);
   grae->SetPoint(26,525,0);
   grae->SetPointError(26,0,0,0,0);
   grae->SetPoint(27,575,0);
   grae->SetPointError(27,0,0,0,0);
   grae->SetPoint(28,650,0);
   grae->SetPointError(28,0,0,0,0);
   grae->SetPoint(29,750,0);
   grae->SetPointError(29,0,0,0,0);
   grae->SetPoint(30,850,0);
   grae->SetPointError(30,0,0,0,0);
   grae->SetPoint(31,950,0);
   grae->SetPointError(31,0,0,0,0);
   grae->Draw("[]");

    c1->Update();
    c1->Modified();

    TFile *filter_old = new TFile("test_pp_plott_filter_old_hist.root");
    PUtils::correct(invmasdilac9deg);
     PUtils::correct(invmasdilac);
     // invmasdilac9deg->Draw("same");
    //invmasdilac->Draw("same");
    
c1->GetFrame()->SetBorderMode(0);
c1->Update();
c1->Modified();
c1->cd();
c1->SetSelected(c1);
c1->Print("pp_cocktail_comparison_tmp.eps");



}

