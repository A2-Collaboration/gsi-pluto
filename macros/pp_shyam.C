
{

    gSystem->CompileMacro("./HadesAkzeptanz.C");
    TFile matrix1=TFile("/d/jspc03/cernuser/HParamsTat/matricesEffSingle.PairCode.all.cuts.ep.may07_v2.root");
    TH3F * ep_acc = (TH3F*) gROOT->FindObject("acce3DPosi");
    TH3F * ep_eff = (TH3F*) gROOT->FindObject("effi3DPosiAllCut");
    
    matrix2=TFile("/d/jspc03/cernuser/HParamsTat/matricesEffSingle.PairCode.all.cuts.em.may07_v2.root");
    TH3F * em_acc = (TH3F*) gROOT->FindObject("acce3DEle");
    TH3F * em_eff = (TH3F*) gROOT->FindObject("effi3DEleAllCut");

    HadesAcceptance * ak = new HadesAcceptance();
    //    ak->SetEpAcc(ep_acc);
    // ak->SetEmAcc(em_acc);

    PBatch x;
    // x.AddCommand("my9deg = (9./180.)*TMath::Pi();");
    x.AddCommand("my9deg = 0.;");
    x.Execute();

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
    
    TH2D * frame = new TH2D("frame","frame",1, 0.0,0.8,1,0.00000000001,0.0001);
    
    
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetTitleOffset(1.1);
    frame->GetXaxis()->SetTitleSize(.06);
    frame->GetYaxis()->SetTitleSize(.06);
    frame->GetYaxis()->SetTitleFont(42);
    frame->GetXaxis()->SetTitleFont(42);
    frame->SetXTitle("m_{ee} [GeV/c^{2}]");
    
    frame->SetYTitle("d#sigma/dM [b #upoint GeV ^{-1}c^{2} ]");
    c1->cd();
    frame->Draw();
    
    Double_t kin_max = 0.95;
    Double_t histo_max = 0.8;
    Int_t histo_bins = 50;
    
    Int_t num_events = 10000;
    //num_events = 50000;
     num_events = 500000;

    makeDistributionManager()->Exec("brems: shyam");
    makeDistributionManager()->Exec("brems: sum; weighting");
    makeDistributionManager()->Exec("dalitz_mod: krivoruchenko");
    makeDistributionManager()->Exec("dalitz_mod: static_br_thresh=0.100 ; flat_generator");    
    makeDistributionManager()->Enable("vmd");

    //TOTAL SUM
    PReaction *my_reaction2 = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);
    my_reaction2->AddBulk(ak);
    //Create my histogram:
    TH1F * pp_sum2 = new TH1F ("pp_sum","pp DiLepton mass (coherent sum)",histo_bins,0.01,histo_max);
    my_reaction2->Do(pp_sum2,"if (_acc_valid); if (([e+]->Angle([e-])) > my9deg); _x=[dilepton]->M()");
    my_reaction2->Print();
    my_reaction2->Loop(num_events);
    pp_sum2->SetLineWidth(2);
    PUtils::correct(pp_sum2); //correct for number of used bins
    pp_sum2->Draw("Csame");

    //compared with KK
    makeDistributionManager()->Exec("brems: kaptari");
    makeDistributionManager()->Exec("brems: sum; weighting");

    PReaction *my_reaction3 = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);
    //Create my histogram:
    TH1F * pp_sum3 = new TH1F ("pp_sum3","pp DiLepton mass (coherent sum)",histo_bins,0.1,histo_max);
    my_reaction3->AddBulk(ak);
    // my_reaction3->Do("_acc_valid->Print()");
    my_reaction3->Do(pp_sum3,"if (_acc_valid); if (([e+]->Angle([e-])) > my9deg); _x=[dilepton]->M()");
    my_reaction3->Print();
    my_reaction3->Loop(num_events);
    PUtils::correct(pp_sum3); //correct for number of used bins
    pp_sum3->Draw("Csame");
    pp_sum3->SetLineStyle(2);
    pp_sum3->SetLineWidth(2);
    
    //ELASTIC
    makeDistributionManager()->Exec("brems: shyam; elastic; weighting");
    PReaction *my_reaction2_elastic = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);
    my_reaction2_elastic->AddBulk(ak);
    my_reaction2_elastic->Print();
    //Create my histogram:
    TH1F * pp_sum2_elastic = new TH1F ("pp_sum2_elastic","pp DiLepton mass (coherent sum)",histo_bins,0.01,histo_max);
    my_reaction2_elastic->Do(pp_sum2_elastic,"if (_acc_valid); if (([e+]->Angle([e-])) > my9deg); _x=[dilepton]->M()");
    my_reaction2_elastic->Loop(num_events);
    pp_sum2_elastic->SetLineColor(2);
    pp_sum2_elastic->SetLineWidth(2);
    PUtils::correct(pp_sum2_elastic); //correct for number of used bins
    pp_sum2_elastic->Draw("Csame");

    //compared with KK
    makeDistributionManager()->Exec("brems: kaptari; elastic; weighting");

    PReaction *my_reaction3_elastic = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);
    my_reaction3_elastic->AddBulk(ak);
    my_reaction3_elastic->Print();
    //Create my histogram:
    TH1F * pp_sum3_elastic = new TH1F ("pp_sum3_elastic","pp DiLepton mass (coherent sum)",histo_bins,0.1,histo_max);
    my_reaction3_elastic->Do(pp_sum3_elastic,"if (_acc_valid); if (([e+]->Angle([e-])) > my9deg); _x=[dilepton]->M()");
    my_reaction3_elastic->Loop(num_events);
    PUtils::correct(pp_sum3_elastic); //correct for number of used bins
    pp_sum3_elastic->Draw("Csame");
    pp_sum3_elastic->SetLineColor(2);
    pp_sum3_elastic->SetLineStyle(2);
    pp_sum3_elastic->SetLineWidth(2);

    //DELTA
    makeDistributionManager()->Exec("brems: shyam; delta; weighting");
    PReaction *my_reaction2_delta = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);
    my_reaction2_delta->AddBulk(ak);
    my_reaction2_delta->Print();
    //Create my histogram:
    TH1F * pp_sum2_delta = new TH1F ("pp_sum2_delta","pp DiLepton mass (coherent sum)",histo_bins,0.01,histo_max);
    my_reaction2_delta->Do(pp_sum2_delta,"if (_acc_valid); if (([e+]->Angle([e-])) > my9deg); _x=[dilepton]->M()");
    my_reaction2_delta->Loop(num_events);
    pp_sum2_delta->SetLineColor(8);
    pp_sum2_delta->SetLineWidth(2);
    PUtils::correct(pp_sum2_delta); //correct for number of used bins
    pp_sum2_delta->Draw("Csame");

    //compared with KK
    makeDistributionManager()->Exec("brems: kaptari; delta; weighting");

    PReaction *my_reaction3_delta = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);
    my_reaction3_delta->AddBulk(ak);
    my_reaction3_delta->Print();
    //Create my histogram:
    TH1F * pp_sum3_delta = new TH1F ("pp_sum3_delta","pp DiLepton mass (coherent sum)",histo_bins,0.1,histo_max);
    my_reaction3_delta->Do(pp_sum3_delta,"if (_acc_valid); if (([e+]->Angle([e-])) > my9deg); _x=[dilepton]->M()");
    my_reaction3_delta->Loop(num_events);
    PUtils::correct(pp_sum3_delta); //correct for number of used bins
    pp_sum3_delta->Draw("Csame");
    pp_sum3_delta->SetLineColor(8);
    pp_sum3_delta->SetLineStyle(2);
    pp_sum3_delta->SetLineWidth(2);

    //Finally PLUTO delta
    PReaction *my_reaction3_pluto = new PReaction("1.25","p","p","p D+ [p dilepton [e+ e-]]");
    my_reaction3_pluto->AddBulk(ak);
    TH1F * pp_sum3_pluto = new TH1F ("pp_sum3_pluto","pp DiLepton mass (coherent sum)",histo_bins,0.01,histo_max);
    my_reaction3_pluto->Do(pp_sum3_pluto,"if (_acc_valid); if (([e+]->Angle([e-])) > my9deg); _x=[dilepton]->M()");
    my_reaction3_pluto->Print();
    my_reaction3_pluto->Loop(num_events);
    pp_sum3_pluto->SetLineColor(8);
    pp_sum3_pluto->SetLineStyle(3);
    pp_sum3_pluto->SetLineWidth(2);
    PUtils::correct(pp_sum3_pluto); //correct for number of used bins
    pp_sum3_pluto->Draw("Csame");


    //N1520
    makeDistributionManager()->Exec("brems: shyam; n1520; weighting");
    PReaction *my_reaction2_n1520 = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);
    my_reaction2_n1520->AddBulk(ak);
    my_reaction2_n1520->Print();
    //Create my histogram:
    TH1F * pp_sum2_n1520 = new TH1F ("pp_sum2_n1520","pp DiLepton mass (coherent sum)",histo_bins,0.01,histo_max);
    my_reaction2_n1520->Do(pp_sum2_n1520,"if (_acc_valid); if (([e+]->Angle([e-])) > my9deg); _x=[dilepton]->M()");
    my_reaction2_n1520->Loop(num_events);
    pp_sum2_n1520->SetLineColor(4);
    pp_sum2_n1520->SetLineWidth(2);
    PUtils::correct(pp_sum2_n1520); //correct for number of used bins
    pp_sum2_n1520->Draw("Csame");




    c1->Update();
    c1->Modified();
    
    c1->GetFrame()->SetBorderMode(0);
    c1->Update();
    c1->Modified();
    c1->cd();
    c1->SetSelected(c1);
    c1->Print("pp_shyam_tmp.eps");
    
}

