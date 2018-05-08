{


//First define our histograms
TH2F * histo = new TH2F ("histo","Dalitz",30,2.,3.7,30,2.,3.7);
histo->Sumw2();
TH1F * histo1 = new TH1F ("histo1","p/pi0 missing mass",100,0.1,2.0);
TH1F * histo2 = new TH1F ("histo2","cos theta of pp",20,-1.,1.);
TH1F * histo3 = new TH1F ("histo3","cos theta of D+ decay",20,-1.,1.);

makeDistributionManager();

PNNFSI *bla = new PNNFSI("p + p_to_p_p_eta/nnfsi","bla",-1);
bla->EnableWeighting();
bla->SetExpectedWeightMean(-1);
makeDistributionManager()->Add(bla);

//Define the reaction as usual
PReaction my_reaction(3.13,"p","p","p p eta", NULL,1,0,0,0);

//Create the container of the histogram list
PProjector *m1 = new PProjector(); 

//Combine the masses of p,1 and pi0, p,2, pi0 to the Dalitz plot
m1->AddHistogram(histo,"_x = ([p,1] + [eta])->M2() ; _y = ([p,2] + [eta])->M2() ");

my_reaction.AddBulk(m1);

my_reaction.Loop(50000);
//my_reaction.Loop(10);

gStyle->SetPalette(1,0);
histo->Draw("colz");

}
