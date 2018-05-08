{
TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pp_delta_dalitz_sam.root");
if (!f) f = new TFile("pp_delta_dalitz_sam.root");
TTree *Reaction = (TTree*)gDirectory->Get("data");
TClonesArray *evt=new TClonesArray("PParticle",7);
Reaction->SetBranchAddress("Particles",&evt);
PParticle *par[7];

Int_t nentries = Reaction->GetEntries();
TH1F *mass=new TH1F("mass","Dilepton mass",100,0,1.);
mass->Sumw2();

Double_t sum=0;
Double_t sum1=0,sum2=0,sum3=0,;

PParticle *dil,*d;

for (Int_t i=0; i<nentries;i++) {
    
    Reaction->GetEntry(i);
    
    for (int j=0;j<7;j++) {
	par[j]=(PParticle*)evt.At(j);
	if (par[j]->is("dilepton")) dil    =(PParticle*)evt.At(j);
	if (par[j]->is("D+")) d    =(PParticle*)evt.At(j);
    }
    sum1+=dil->W();
    mass->Fill(dil->M(),(Stat_t) dil->W());
    
}
mass->Draw("e1");

TFile *f2 = (TFile*)gROOT->GetListOfFiles()->FindObject("pp_delta_dalitz_wei.root");
if (!f2) f2 = new TFile("pp_delta_dalitz_wei.root");
Reaction = (TTree*)gDirectory->Get("data");

Reaction->SetBranchAddress("Particles",&evt);

nentries = Reaction->GetEntries();
TH1F *mass2=new TH1F("mass2","Dilepton mass",100,0,1.);
mass2->Sumw2();

for (Int_t i=0; i<nentries;i++) {
    
    Reaction->GetEntry(i);
    
    for (int j=0;j<7;j++) {
	par[j]=(PParticle*)evt.At(j);
	if (par[j]->is("dilepton")) dil    =(PParticle*)evt.At(j);
	if (par[j]->is("D+")) d    =(PParticle*)evt.At(j);
    }
//     sum+=dil->W();
//     cout << dil->W() << endl;
    sum2+=dil->W();
    mass2->Fill(dil->M(),(Stat_t) dil->W());
    
}

//mass2->Scale(180000);
//mass2->Divide(mass);
//mass2->Draw("e1");
mass2->Draw("samee1");
mass2->SetLineColor(2);

TFile *f22 = (TFile*)gROOT->GetListOfFiles()->FindObject("pp_delta_dalitz_wei2.root");
if (!f22) f22 = new TFile("pp_delta_dalitz_wei2.root");
Reaction = (TTree*)gDirectory->Get("data");

Reaction->SetBranchAddress("Particles",&evt);

nentries = Reaction->GetEntries();
TH1F *mass22=new TH1F("mass22","Dilepton mass",100,0,1.);
mass22->Sumw2();
//nentries=100;
for (Int_t i=0; i<nentries;i++) {
    
    Reaction->GetEntry(i);
    
    for (int j=0;j<7;j++) {
	par[j]=(PParticle*)evt.At(j);
	if (par[j]->is("dilepton")) dil    =(PParticle*)evt.At(j);
	if (par[j]->is("D+")) d    =(PParticle*)evt.At(j);
    }
    mass22->Fill(dil->M(),(Stat_t) dil->W());
//    cout << dil->W() << endl;

    sum3+=dil->W();
//    mass22->Fill(dil->M());
   
}

//mass22->Scale(0.3);
//mass22->Divide(mass2);
//mass22->Draw("e1");
mass22->Draw("samee1");
mass22->SetLineColor(3);

cout << sum1 << ":" << sum2 << ":" << sum3 << endl;

}
