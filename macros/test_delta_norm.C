{

Double_t beam_energy = 1.25;
//Double_t beam_energy = 2.2;
Char_t beam_energy_c[10];
sprintf(beam_energy_c,"%f",beam_energy);

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

makeDistributionManager()->Exec("dalitz_mod: krivoruchenko");
makeDistributionManager()->Exec("dalitz_mod: static_br_thresh=0.100 ; flat_generator");

makeDistributionManager()->Exec("elementary");

//new FF
gSystem->CompileMacro( "../macros/PDeltaDalitzFF.C");

PDeltaDalitzFF * newmodel = new PDeltaDalitzFF("iachello@D+_to_p_dilepton/formfactor",
					       "Iachello ff for D+ -> p e+e-",-1);
newmodel->SetQED(1);
//newmodel->SetCC(2.7,0,0);
newmodel->SetCC(3.0,0,0);

makeDistributionManager()->Add(newmodel);
makeDistributionManager()->LinkDB();

//static BR from Kriv.
PChannelModel* kriv = ((PChannelModel*) makeDistributionManager()->GetDistribution("D+_krivoruchenko"));
kriv->Draw();
//cout <<  kriv->Integral(0.,0.5)/ 0.120<<endl;

cout << "Full width BR:" << endl;
cout <<  kriv->GetWidth(1.232) / 0.120 <<endl;


//return;

PReaction *my_reaction;
my_reaction = new PReaction(beam_energy_c,"p","p",
			    "p D+ [dilepton [e+ e-] p]",
			    NULL,1,0,0,0);

//Create my histogram:
TH1F * pn_sum = new TH1F ("pn_sum","pn DiLepton mass (coherent sum)",30,0.,0.8);
TH1F * pn_sum2 = new TH1F ("pn_sum2","pn DiLepton mass (coherent sum)",30,0.,0.8);
pn_sum->Sumw2();
pn_sum2->Sumw2();
//Create the container of the histogram list
PProjector *m1 = new PProjector(); 
//Dilepton mass
m1->AddHistogram(pn_sum,"_x=[dilepton,1]->M()");

my_reaction->AddBulk(m1);


my_reaction->Loop(100000);


newmodel->SetQED(0);


PReaction *my_reaction2;
my_reaction2 = new PReaction(beam_energy_c,"p","p",
			    "p D+ [dilepton [e+ e-] p]",
			    NULL,1,0,0,0);

//Create the container of the histogram list
PProjector *m2 = new PProjector(); 
//Dilepton mass

m2->AddHistogram(pn_sum2,"_x=[dilepton,1]->M()");
my_reaction2->AddBulk(m2);


my_reaction2->Loop(100000);


//Get the total cross section for our beam energy
PParticle p("p");
PParticle p2("p",beam_energy);
PParticle p3= p+p2;
PChannelModel* t = ((PChannelModel*) makeDistributionManager()->GetDistribution("p + p_to_D+_p/tcross"));
cout << "QED:" << endl;
cout << pn_sum->Integral() / t->Eval(p3.M()) << endl;

cout << "VMD:" << endl;
cout << pn_sum2->Integral() / t->Eval(p3.M()) << endl;

PUtils::correct(pn_sum); //correct for number of used bins
PUtils::correct(pn_sum2); //correct for number of used bins




pn_sum->Draw("");
//pn_sum->SetMaximum(0.01);
pn_sum2->Draw("same");
pn_sum2->SetLineColor(2);
}
