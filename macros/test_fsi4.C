{


makeDistributionManager()->GetDistribution("brems")->Exec("fsi");

PReaction *my_reaction;
my_reaction = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);

//Create my histogram:
TH1F * pn_sum = new TH1F ("pn_sum","pn DiLepton mass (coherent sum)",30,0.,0.5);
pn_sum->Sumw2();
//Create the container of the histogram list
PProjector *m1 = new PProjector(); 
//Dilepton mass
m1->AddHistogram(pn_sum,"_x=[dilepton]->M()");
my_reaction->AddBulk(m1);


my_reaction->Loop(100000);

pn_sum->Draw();

//now with wei

makeDistributionManager()->GetDistribution("brems")->Exec("kin_max=0.9 ; weighting");


PReaction *my_reaction2;
my_reaction2 = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);

//Create my histogram:
TH1F * pn_sum2 = new TH1F ("pn_sum2","pn DiLepton mass (coherent sum)",30,0.,0.5);
pn_sum2->Sumw2();
//Create the container of the histogram list
PProjector *m2 = new PProjector(); 
//Dilepton mass
m2->AddHistogram(pn_sum2,"_x=[dilepton]->M()");
my_reaction2->AddBulk(m2);


my_reaction2->Loop(100000);

pn_sum2->Draw("e1");
pn_sum->Draw("same");

pn_sum2->Divide(pn_sum);
pn_sum2->Draw("e1");

}
