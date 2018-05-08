
{

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


TH2D * frame = new TH2D("frame","frame",1, 0.0,0.6,1,0.000000000001,0.0001);


frame->GetXaxis()->SetTitleOffset(1.1);
frame->GetYaxis()->SetTitleOffset(1.1);
frame->GetXaxis()->SetTitleSize(.06);
frame->GetYaxis()->SetTitleSize(.06);
frame->GetYaxis()->SetTitleFont(42);
frame->GetXaxis()->SetTitleFont(42);
frame->SetXTitle("M [GeV/c^{2}]");

frame->SetYTitle("d#sigma/dM [b]");
 frame->Draw();

    Double_t kin_max = 0.95;
    Double_t histo_max = 0.6;
//    Int_t histo_bins = 60;
    Int_t histo_bins = 100;
    
    Int_t num_events = 50000;
    num_events = 10000;
#if 1
    makeDistributionManager()->Exec("elementary");
    makeDistributionManager()->LinkDB();
    //Now the pp case
    PReaction *my_reaction = new PReaction("1.25","p","p","p D+ [p dilepton [e+ e-]]",NULL,1,0,0,0);
    my_reaction->Print();
    //Create my histogram:
    TH1F * pp_sum = new TH1F ("pp_sum","pp DiLepton mass (coherent sum)",histo_bins,0.,histo_max);
    TH1F * dummy = new TH1F ("dummy","pp DiLepton mass (coherent sum)",10,-1,1);
    pp_sum->Sumw2();
    dummy->Sumw2();
    //Create the container of the histogram list
    PProjector *m1 = new PProjector(); 
    //Dilepton mass
    m1->AddHistogram(pp_sum,"_x=[dilepton]->M()");
    m1->AddHistogram(dummy,"_x=cos([dilepton]->Theta())");
    my_reaction->AddBulk(m1);
    my_reaction->Loop(10000);
    PUtils::correct(pp_sum); //correct for number of used bins
//    PUtils::correct(dummy);

    pp_sum->Draw("same");


    makeDistributionManager()->Exec("brems_kaptari : kin_max=1. ; weighting");
    
//     TF1 *flat = new TF1("flat","1*(x>0.1) + 0.0001*(x<=0.1)",0,0.9);
//     pp_gen = 
// 	new PInclusiveModel("flatx@p + p_brems_p_p_dilepton/generator","Dilepton generator",-1);
    
//     //The distribution template:
//     pp_gen->Add("q,parent");
//     pp_gen->Add("p,daughter");
//     pp_gen->Add("p,daughter");
//     pp_gen->Add("dilepton,daughter,primary");
//     pp_gen->SetSampleFunction(flat);
//     //Enable distribution as a generator
//     pp_gen->EnableGenerator();    
//     makeDistributionManager()->Add(pp_gen);
    

//    makeDistributionManager()->GetDistribution("brems_kaptari")->Exec("fsi");
    makeDistributionManager()->Exec("brems_kaptari : delta");

    PReaction *my_reaction2 = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);

    //Create my histogram:
    TH1F * pp_sum2 = new TH1F ("pp_sum","pp DiLepton mass (coherent sum)",histo_bins,0.1,histo_max);
    TH1F * dummy2 = new TH1F ("dummy2","pp DiLepton mass (coherent sum)",10,-1,1);
    pp_sum2->Sumw2();
    dummy2->Sumw2();
    //Create the container of the histogram list
    PProjector *m2 = new PProjector(); 
    //Dilepton mass
    m2->AddHistogram(pp_sum2,"_x=[dilepton]->M()");
    m2->AddHistogram(dummy2,"_x=cos([dilepton]->Theta())");
    my_reaction2->AddBulk(m2);
    my_reaction2->Loop(num_events);
    PUtils::correct(pp_sum2); //correct for number of used bins
    PUtils::correct(dummy2);
  
    PBremsstrahlung * brems =  makeDistributionManager()->GetDistribution("p + p_brems_p_p_dilepton");
    brems->SetNeutron(-1);
    brems->SetP2E(1.25);
  
    pp_sum2->SetMarkerStyle(22);
    pp_sum2->Draw("same");

    brems->Draw("same");

    //now the classical for weighting
#endif
    TF1 *flat = new TF1("flat","1",0,1);
    
    PInclusiveModel *dilepton_generator = new PInclusiveModel("D+_flat@D+_to_p_dilepton/generator",
							      "Dilepton generator",-1);    
    dilepton_generator->Set("dilepton,primary");
    dilepton_generator->SetSampleFunction(flat);
    dilepton_generator->EnableGenerator();
    
    makeDistributionManager()->Add(dilepton_generator);

    makeDistributionManager()->GetDistribution("D+_dalitz")->EnableWeighting();
    //makeDistributionManager()->GetDistribution("D+_dalitz")->SetExpectedWeightMean(-1);
    //Use Delta Dalitz model for weighting
    
    PReaction *my_reaction3 = new PReaction("1.25","p","p","p D+ [p dilepton [e+ e-]]",NULL,1,0,0,0);
    my_reaction3->Print();

    //Create my histogram:
    TH1F * pp_sum3 = new TH1F ("pp_sum3","pp DiLepton mass (coherent sum)",histo_bins,0.,histo_max);
 
    pp_sum3->Sumw2();
 
    //Create the container of the histogram list
    PProjector *m3 = new PProjector(); 
    //Dilepton mass
    m3->AddHistogram(pp_sum3,"_x=[dilepton]->M()");
 
    my_reaction3->AddBulk(m3);
    my_reaction3->Loop(num_events);
    PUtils::correct(pp_sum3); //correct for number of used bins

    pp_sum3->SetMarkerStyle(21);
    pp_sum3->Draw("same");



    //Now the kriv
    gSystem->CompileMacro( "../macros/PDeltaDalitzKrivoruchenko.C");

    PDeltaDalitzKrivoruchenko * newmodel2 = new PDeltaDalitzKrivoruchenko("kriv@D+_to_p_dilepton",
									  "dgdm from Krivoruchenko",-1);
    newmodel2->EnableWeighting();
    newmodel2->SetExpectedWeightMean(-1);
	
    TF1 *flat3 = new TF1("flat3","1",0,1);
    
    PInclusiveModel *dilepton_generator = new PInclusiveModel("flat@D+_to_p_dilepton/generator",
							      "Dilepton generator",-1);    
    dilepton_generator->Set("dilepton,primary");
    dilepton_generator->SetSampleFunction(flat);
    dilepton_generator->EnableGenerator();
    
    makeDistributionManager()->Add(dilepton_generator);
    makeDistributionManager()->Add(newmodel2);
    
    
    PReaction *my_reaction4 = new PReaction("1.25","p","p","p D+ [p dilepton [e+ e-]]",NULL,1,0,0,0);
    my_reaction4->Print();

    //Create my histogram:
    TH1F * pp_sum4 = new TH1F ("pp_sum4","pp DiLepton mass (coherent sum)",histo_bins,0.,histo_max);
 
    pp_sum4->Sumw2();
 
    //Create the container of the histogram list
    PProjector *m4 = new PProjector(); 
    //Dilepton mass
    m4->AddHistogram(pp_sum4,"_x=[dilepton]->M()");
 
    my_reaction4->AddBulk(m4);
    my_reaction4->Loop(num_events);
    //cout<< pp_sum4->Integral() << endl;
    
    PUtils::correct(pp_sum4); //correct for number of used bins

    pp_sum4->SetMarkerStyle(20);
    pp_sum4->Draw("same");
    pp_sum4->SetMarkerColor(3);

}

