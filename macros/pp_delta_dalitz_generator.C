// p + p -> p + Delta+ -> p + p + gamma* -> p + p + e- + e+
// using the generator

{
gROOT->Reset();

Int_t version=0; //d-dalitz generator
//version=1; //flat generator
 //version=2; //x2 generator
 
 if ((version==1) || (version==2)) {
 
     makeDistributionManager()->GetDistribution("D+_dalitz")->EnableWeighting();
     makeDistributionManager()->GetDistribution("eta_dalitz")->EnableWeighting();
     //makeDistributionManager()->GetDistribution("D+_dalitz")->SetExpectedWeightMean(-1);
     //makeDistributionManager()->GetDistribution("D+_dalitz")->DisableSampling();

     TF1 *flat;
//A flat distribution can be done via the "inclusive model"
     if (version==1) 
	 flat=new TF1("flat","1",0,1);
     else
	 //flat=new TF1("flat","(x+1)*(x+1)",0,1);
     //flat=new TF1("flat","0.5+x",0,1);
	 flat=new TF1("flat","0.2+x*x",0,1);



     //PInclusiveModel *dilepton_generator = new PInclusiveModel("flat@D+_to_p_dilepton","Dilepton generator",-1);
     PInclusiveModel *dilepton_generator = new PInclusiveModel("flat@eta_to_g_dilepton","Dilepton generator",-1);
     //dilepton_generator->Add("D+,parent");
     //dilepton_generator->Add("p,daughter");
     dilepton_generator->Add("eta,parent");
     dilepton_generator->Add("g,daughter");
     dilepton_generator->Add("dilepton,daughter,primary");
     dilepton_generator->SetSampleFunction(flat);
     dilepton_generator->EnableGenerator();
     
     makeDistributionManager()->Add(dilepton_generator);
     //makeDistributionManager()->Print("decay_models");
     
 }

//3.13,"p","p","p p eta [dilepton [e+ e-] g]", "eta_dalitz",1,0,0,0

 PReaction *r;

 if (version==0) 
     //r=new PReaction(3.13,"p","p","p D+ [dilepton [e+ e-] p]", "pp_delta_dalitz_sam",1,0,0,0);
     r=new PReaction(3.13,"p","p","p p eta [dilepton [e+ e-] g]", "pp_delta_dalitz_sam",1,0,0,0);
 if (version==1) 
     //r=new PReaction(3.13,"p","p","p D+ [dilepton [e+ e-] p]", "pp_delta_dalitz_wei",1,0,0,0);
     r=new PReaction(3.13,"p","p","p p eta [dilepton [e+ e-] g]", "pp_delta_dalitz_wei",1,0,0,0);
 if (version==2) 
     //r=new PReaction(3.13,"p","p","p D+ [dilepton [e+ e-] p]", "pp_delta_dalitz_wei2",1,0,0,0);
     r=new PReaction(3.13,"p","p","p p eta [dilepton [e+ e-] g]", "pp_delta_dalitz_wei2",1,0,0,0);

r->Print();
 r->Preheating(100);
r->loop(10000);

 r->PrintReport();

}


//10000: 4.28412541397422963e-05
//30000: 1.42804152004180791e-05
//300000: 1.42804150282578059e-06
