{
//Add the unknown decay
makeStaticData()->AddDecay(-1,"eta -> dilepton + dilepton", "eta", "dilepton,dilepton",0.0000001);
//--->The last number is the BR, replace it

gSystem->CompileMacro( "./PEtaDoubleDalitz.C");
gSystem->CompileMacro( "./PEtaDoubleDalitzFF.C");

//***** Add a model for di-lepton sampling ******************************
PEtaDoubleDalitz *dilepton_generator =
new PEtaDoubleDalitz("blabla@eta_to_dilepton_dilepton",
		     "Dilepton generator for eta -> dilepton + dilepton",-1);

//Activate the model
makeDistributionManager()->Add(dilepton_generator);
makeDistributionManager()->Print();
//***********************************************************************

PEtaDoubleDalitzFF * newmodel 
= new PEtaDoubleDalitzFF("double_ff@eta_to_dilepton_dilepton/formfactor",
			 "Eta DD FF",-1);


makeDistributionManager()->Add(newmodel);
makeDistributionManager()->Print();
    



PReaction my_reaction(2.2,"p","p","p p eta [dilepton [e+ e-] dilepton [e+ e-]]","etae+e-e+e-");

//This is for on-line debugging
TH1F * histo1 = new TH1F ("histo1","ee invariant mass",100,0.,0.4);
TH2F * histo2 = new TH2F ("histo2","ee1 vs. ee2",20,0.,0.4,20,0.,0.4);

my_reaction.Do(histo1,"_x = ([e+,1] + [e-,1])->M() ");
my_reaction.Do(histo1,"_x = ([e+,2] + [e-,2])->M() ");

my_reaction.Do(histo2,"_x = ([e+,1] + [e-,1])->M() ; _y =  ([e+,2] + [e-,2])->M()");


my_reaction.Print();
my_reaction.Loop(100000);



//Now the same but without form factor

makeDistributionManager()->Disable("double_ff");

PReaction my_reaction2(2.2,"p","p","p p eta [dilepton [e+ e-] dilepton [e+ e-]]","etae+e-e+e-");

//This is for on-line debugging
TH1F * histo1a = new TH1F ("histo1","ee invariant mass",100,0.,0.4);

my_reaction2.Do(histo1a,"_x = ([e+,1] + [e-,1])->M() ");
my_reaction2.Do(histo1a,"_x = ([e+,2] + [e-,2])->M() ");


my_reaction2.Print();
my_reaction2.Loop(100000);





c1 = new TCanvas();

histo1->Draw("");
histo1a->Draw("same");

c2 = new TCanvas();

histo2->Draw("colz");

}
