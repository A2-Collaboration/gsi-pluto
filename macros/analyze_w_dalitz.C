{
// Reset ROOT and connect tree file
//TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("w_dalitz.root");
// set up tree and branch addresses
//TTree *Reaction = (TTree*)gDirectory->Get("data"); // get 1st tree
TChain *Reaction = new TChain("data");
Reaction->Add("w_dalitz.root");
TClonesArray *evt=new TClonesArray("PParticle",8);        // must create TClonesArray
TArrayI *Status=new TArrayI(3);                           // overall filter status
Reaction->SetBranchAddress("Filters",Status->GetArray()); // recover the filter status
Reaction->SetBranchAddress("Particles",&evt);             // recover PParticles

// intermediate variables
const double r2d=180./3.14159265358979323846; 
TVector3 v1,v2;
PParticle *par[8];                         // pointer to array of PParticles

// set up histograms
TH1F *hf1=new TH1F("hf1","a) #pi^{-} beam p_{z} profile",200,1.18,1.28);
TH1F *hf2=new TH1F("hf2","b) #omega mass",200,0.72,0.84);
TH1F *hf3=new TH1F("hf3","c) dilepton mass   ",200,0.,0.8);
TH1F *hf4=new TH1F("hf4","d) polar angle all e^{+} and e^{-} (deg)",200,0,180);
TH1F *hf5=new TH1F("hf5","e) polar angle accepted e^{+} and e^{-} (deg)",200,0,180);
TH1F *hf6=new TH1F("hf6","f) accepted e^{+}e^{-} opening angle (deg)",200,0,180);

// display output
c1 = new TCanvas("c1","pi- + p -> n + w -> n + pi0 + gamma* -> n + pi0 + e- + e+",
                 200,10,900,900);
c1->Divide(3,2);
c1_1->SetTitle("c1_1");
c1_2->SetTitle("c1_2");
c1_3->SetLogy();
c1_3->SetTitle("c1_3");
c1_4->SetTitle("c1_4");
c1_5->SetTitle("c1_5");
c1_6->SetLogy();
c1_6->SetTitle("c1_6");
c1_1->Draw();
c1_2->Draw();
c1_3->Draw();
c1_4->Draw();
c1_5->Draw();
c1_6->Draw();

// event loop
Long64_t nentries = Reaction->GetEntries();      // get number of entries
for (Int_t i=0; i<nentries;i++) {                // enter event loop
  Reaction->GetEntry(i);                         // current event
//  par[0]=(PParticle*)evt.At(0);                  // retrieve pi- beam particle
  par[3]=(PParticle*)evt.At(2);                  // retrieve omega meson
  par[5]=(PParticle*)evt.At(4);                  // retrieve virtual photon
  par[6]=(PParticle*)evt.At(5);                  // retrieve electron
  par[7]=(PParticle*)evt.At(6);                  // retrieve positron
//  hf1->Fill(par[0]->Pz());                       // beam Pz momentum profile
  hf2->Fill(par[3]->M());                        // omega-meson mass
  hf3->Fill(par[5]->M());                        // dilepton mass
  hf4->Fill(r2d*(par[6]->Theta()));              // electron polar angle (deg)
  hf4->Fill(r2d*(par[7]->Theta()));              // positron polar angle (deg)
  if ( (*Status)[0]==0 ) {                       // accepted dileptons only
    hf5->Fill(r2d*(par[6]->Theta()));            // electron polar angle (deg)
    hf5->Fill(r2d*(par[7]->Theta()));            // positron polar angle (deg)
    v1=par[6]->Vect();                           // e- momentum vector
    v2=par[7]->Vect();                           // e+ momentum vector
    hf6->Fill(r2d*(acos(v1.Dot(v2)/(v1.Mag()*v2.Mag())))); // e+e- angle (deg)
    }
}
c1_1->cd();
hf1->Draw();
c1_2->cd();
hf2->Draw();
c1_3->cd();
hf3->Draw();
c1_4->cd();
hf4->Draw();
c1_5->cd();
hf5->Draw();
c1_6->cd();
hf6->Draw();
}


