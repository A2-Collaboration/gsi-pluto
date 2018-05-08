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

TH2D * frame = new TH2D("frame","frame",1, 0.0,0.6,1,0.0000000000001,0.0001);


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

    Double_t histo_max = 0.6;
//    Int_t histo_bins = 60;
    Int_t histo_bins = 100;
    
    Int_t num_events = 50000;
    num_events = 100000;


makeDistributionManager()->Exec("brems_kaptari : sum; weighting");

PReaction *my_reaction2 = new PReaction("1.25","p","n","p n dilepton [e+ e-]",NULL,1,0,0,0);

//Create my histogram:
TH1F * pn_sum2 = new TH1F ("pn_sum","pn DiLepton mass (coherent sum)",histo_bins,0.1,histo_max);

//pn_sum2->Sumw2();

//Create the container of the histogram list
PProjector *m2 = new PProjector(); 
//Dilepton mass
m2->AddHistogram(pn_sum2,"_x=[dilepton]->M()");

my_reaction2->AddBulk(m2);
my_reaction2->Loop(num_events);
PUtils::correct(pn_sum2); //correct for number of used bins

//now the incoherent sum:

//here we have to add 2 independent reactions into the SAME histogram
makeDistributionManager()->Exec("brems_kaptari : delta");
TH1F * pn_inc2 = new TH1F ("pn_inc","pn DiLepton mass (incoherent sum)",histo_bins,0.1,histo_max);
pn_inc2->Sumw2();
PReaction *my_reaction3 = new PReaction("1.25","p","n","p n dilepton [e+ e-]",NULL,1,0,0,0);

//Create the container of the histogram list
PProjector *m3 = new PProjector(); 
//Dilepton mass
m3->AddHistogram(pn_inc2,"_x=[dilepton]->M()");

my_reaction3->AddBulk(m3);
my_reaction3->Loop(num_events);



makeDistributionManager()->Exec("brems_kaptari : elastic");
PReaction *my_reaction3a = new PReaction("1.25","p","n","p n dilepton [e+ e-]",NULL,1,0,0,0);
TH1F * pn_el2 = new TH1F ("pn_el","pn DiLepton mass (elastic)",histo_bins,0.1,histo_max);
pn_el2->Sumw2();

//Create the container of the histogram list
PProjector *m3a = new PProjector(); 
//Dilepton mass
m3a->AddHistogram(pn_inc2,"_x=[dilepton]->M()");
m3a->AddHistogram(pn_el2,"_x=[dilepton]->M()");

my_reaction3a->AddBulk(m3a);
my_reaction3a->Loop(num_events);


PUtils::correct(pn_inc2); //correct for number of used bins
PUtils::correct(pn_el2);

pn_inc2->Draw("same");
pn_sum2->Draw("Csame");
pn_el2->Draw("same");


c1->Update();
c1->Modified();

c1->GetFrame()->SetBorderMode(0);
c1->Update();
c1->Modified();
c1->cd();
c1->SetSelected(c1);
c1->Print("kk_incoherent_pn.eps");


}
