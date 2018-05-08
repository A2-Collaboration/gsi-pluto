{


//makeDistributionManager()->GetDistribution("brems")->Exec("kin_max=0.9 ; weighting");
makeDistributionManager()->GetDistribution("brems")->Exec("fsi");

PReaction *my_reaction;
my_reaction = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);

//Create my histogram:
TH2F * pn_sum = new TH2F ("pn_sum","pn DiLepton mass (coherent sum)",30,0.9,1.8,30,0.9,1.8);
pn_sum->Sumw2();
//Create the container of the histogram list
PProjector *m1 = new PProjector(); 
//Dilepton mass
m1->AddHistogram(pn_sum,"_x=([dilepton] + [p,1])->M() ; _y=([dilepton] + [p,2])->M()");
my_reaction->AddBulk(m1);


my_reaction->Loop(100000);

pn_sum->Draw("box");

}
