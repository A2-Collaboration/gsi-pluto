Double_t angle(Double_t *x, Double_t *par) {
   return (4 + 3.1*0.5*(3*x[0]*x[0]-1) + 2.0*0.125*(35*x[0]*x[0]*x[0]*x[0] - 30* x[0]*x[0] +3)) / 9.;
}


omega_3pi(){

gStyle->SetPalette(8,0);
gSystem->CompileMacro( "../macros/PHadesAcc.C");


TF1 *f1 = new TF1("myfunc",angle,-1,1,0);
//f1->Draw();
 PAngularDistribution * pp_omega_prod_angle = 
     new PAngularDistribution("pp_omega_prod_angle",
			      "Omega polar angles (DISTO, 2.85GeV)");
 pp_omega_prod_angle->Add("w,  daughter, primary");
 pp_omega_prod_angle->Add("p,    daughter");
 pp_omega_prod_angle->Add("p,    daughter");
 pp_omega_prod_angle->Add("q,    parent,   reference");
 pp_omega_prod_angle->SetRotate(kFALSE);
 pp_omega_prod_angle->SetAngleFunction(f1);
 makeDistributionManager()->Add(pp_omega_prod_angle);

TH2F * histo2 = new TH2F ("histo2","Dalitz",20,0.05,.5,20,0.05,.5);
histo2->Sumw2();
TH2F * histo1 = new TH2F ("histo1","Dalitz /acc",20,0.05,.5,20,0.05,.5);
histo1->Sumw2();

 TH1F * histo3 = new TH1F ("histo3","ang",20,-1,1);
histo3->Sumw2();

PReaction my_reaction(4.338,"p","p","p p w [pi+ pi- pi0]");

my_reaction.AddBulk(new PHadesAcc());

PProjector *m1 = new PProjector(); 
m1->AddHistogram(histo2,"_x=([pi0]+[pi+])->M2();_y=([pi0]+[pi-])->M2()");
m1->AddHistogram(histo1,"if(_hadacc)");
 m1->AddHistogram(histo3,"_w=([pi0]+[pi+]+[pi0]);_w->Boost([p + p]);_x=cos(_w->Theta())");
my_reaction.AddBulk(m1);



my_reaction.Print();
my_reaction.Loop(500000);


histo1->Draw("colz");
 histo2->Draw("boxsame");
 
 new TCanvas();
 histo3->Draw();

}
