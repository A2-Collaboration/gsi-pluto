{


makeDistributionManager()->Exec("brems_kaptari : weighting");
//makeDistributionManager()->Exec("brems_kaptari : fsi");
PReaction *my_reaction;
my_reaction = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);

//Create my histogram:
TH1F * pp_sum = new TH1F ("pp_sum","pp DiLepton mass (coherent sum)",30,0.,0.6);
TH2F * pp_2d = new TH2F ("pp_2d","pp DiLepton mass (coherent sum)",10,0.,0.6,10,1.8,2.5);
TH1F * pp_m = new TH1F ("pp_m","pp DiLepton mass (coherent sum)",30,1.8,2.5);
pp_sum->Sumw2();
pp_2d->Sumw2();
pp_m->Sumw2();
//Create the container of the histogram list
PProjector *m1 = new PProjector(); 
//Dilepton mass
m1->AddHistogram(pp_sum,"_x=[dilepton]->M()");
m1->AddHistogram(pp_2d,"_x=[dilepton]->M() ; _y=([p,1] + [p,2])->M()");
m1->AddHistogram(pp_m,"_x=([p,1] + [p,2])->M()");
my_reaction->AddBulk(m1);


my_reaction->Loop(100000);

pp_sum->Draw();

//now with FSI

makeDistributionManager()->Exec("brems_kaptari : fsi");

PReaction *my_reaction2;
my_reaction2 = new PReaction("1.25","p","p","p p dilepton [e+ e-]",NULL,1,0,0,0);

//Create my histogram:
TH1F * pp_sum2 = new TH1F ("pp_sum2","pp DiLepton mass (coherent sum)",30,0.,0.6);
TH2F * pp_2d2 = new TH2F ("pp_2d2","pp DiLepton mass (coherent sum)",10,0.,0.6,10,1.8,2.5);
TH1F * pp_m2 = new TH1F ("pp_m2","pp DiLepton mass (coherent sum)",30,1.8,2.5);
pp_sum2->Sumw2();
pp_2d2->Sumw2();
pp_m2->Sumw2();
//Create the container of the histogram list
PProjector *m2 = new PProjector(); 
//Dilepton mass
m2->AddHistogram(pp_sum2,"_x=[dilepton]->M()");
//m2->AddHistogram(pp_2d2,"_x=[dilepton]->M() ; _y=([p,1] + [p,2])->M(); _y->Print()");
m2->AddHistogram(pp_2d2,"_x=[dilepton]->M() ; _y=([p,1] + [p,2])->M()");
m2->AddHistogram(pp_m2,"_x=([p,1] + [p,2])->M()");
my_reaction2->AddBulk(m2);


my_reaction2->Loop(100000);

pp_sum2->Draw("e1");
pp_sum->Draw("same");

pp_sum2->Divide(pp_sum);
pp_sum2->Draw("e1");

pp_2d2->Divide(pp_2d);
//pp_2d2->Draw("box");

pp_m2->Divide(pp_m);

}
