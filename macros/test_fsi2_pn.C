{


makeDistributionManager()->Exec("brems_kaptari : weighting");
//makeDistributionManager()->Exec("brems_kaptari : fsi");
PReaction *my_reaction;
my_reaction = new PReaction("1.25","p","n","p n dilepton [e+ e-]",NULL,1,0,0,0);

//Create my histogram:
TH1F * pn_sum = new TH1F ("pn_sum","pn DiLepton mass (coherent sum)",30,0.,0.6);
TH1F * pn_m = new TH1F ("pn_m","pn DiLepton mass (coherent sum)",30,1.8,2.5);
TH2F * pn_2d = new TH2F ("pn_2d","pn DiLepton mass (coherent sum)",10,0.,0.6,10,1.8,2.5);
pn_sum->Sumw2();
pn_2d->Sumw2();
pn_m->Sumw2();
//Create the container of the histogram list
PProjector *m1 = new PProjector(); 
//Dilepton mass
m1->AddHistogram(pn_sum,"_x=[dilepton]->M()");
m1->AddHistogram(pn_2d,"_x=[dilepton]->M() ; _y=([p,1] + [n,1])->M()");
m1->AddHistogram(pn_m,"_x=([p,1] + [n,1])->M()");
my_reaction->AddBulk(m1);


my_reaction->Loop(100000);

pn_sum->Draw();

//now with FSI

makeDistributionManager()->Exec("brems_kaptari : fsi");

PReaction *my_reaction2;
my_reaction2 = new PReaction("1.25","p","n","p n dilepton [e+ e-]",NULL,1,0,0,0);

//Create my histogram:
TH1F * pn_sum2 = new TH1F ("pn_sum2","pn DiLepton mass (coherent sum)",30,0.,0.6);
TH2F * pn_2d2 = new TH2F ("pn_2d2","pn DiLepton mass (coherent sum)",10,0.,0.6,10,1.8,2.5);
TH1F * pn_m2 = new TH1F ("pn_m2","pn DiLepton mass (coherent sum)",30,1.8,2.5);
pn_sum2->Sumw2();
pn_2d2->Sumw2();
pn_m2->Sumw2();
//Create the container of the histogram list
PProjector *m2 = new PProjector(); 
//Dilepton mass
m2->AddHistogram(pn_sum2,"_x=[dilepton]->M()");
//m2->AddHistogram(pn_2d2,"_x=[dilepton]->M() ; _y=([p,1] + [p,2])->M(); _y->Print()");
m2->AddHistogram(pn_2d2,"_x=[dilepton]->M() ; _y=([p,1] + [n,1])->M()");
m2->AddHistogram(pn_m2,"_x=([p,1] + [n,1])->M()");
my_reaction2->AddBulk(m2);


my_reaction2->Loop(100000);

pn_sum2->Draw("e1");
pn_sum->Draw("same");

pn_sum2->Divide(pn_sum);
pn_sum2->Draw("e1");

pn_2d2->Divide(pn_2d);
//pn_2d2->Draw("box");
pn_m2->Divide(pn_m);

}
