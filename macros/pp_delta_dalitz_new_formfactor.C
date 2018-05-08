
#include "../doc/make_delta_dalitz_fig.C"

pp_delta_dalitz_new_formfactor() {
    gROOT->Reset();
    gStyle->SetOptLogy(1);

    Int_t use_weighting=kTRUE;
    //use_weighting=kFALSE;

    //
    // We use the weighting option with a flat di-lepton generator
    //
    
    if (use_weighting) {
	
	makeDistributionManager()->GetDistribution("D+_dalitz")->EnableWeighting();
	//Use Delta Dalitz model for weighting
	
	TF1 *flat = new TF1("flat","1",0,1);
    
	PInclusiveModel *dilepton_generator = new PInclusiveModel("flat@D+_to_p_dilepton/generator",
								  "Dilepton generator",-1);    
	dilepton_generator->Set("dilepton,primary");
	dilepton_generator->SetSampleFunction(flat);
	dilepton_generator->EnableGenerator();
	
	makeDistributionManager()->Add(dilepton_generator);

    }
    
    //
    // Compile and add the new form factor model
    //
    
    gSystem->CompileMacro( "./PDeltaDalitzFF.C");

    PDeltaDalitzFF * newmodel = new PDeltaDalitzFF("iachello@D+_to_p_dilepton/formfactor",
						   "Iachello ff for D+ -> p e+e-",-1);

    
    makeDistributionManager()->Add(newmodel);
    makeDistributionManager()->Disable("iachello");
    makeDistributionManager()->Print();
    

    //
    // One can draw it!
    // 

    TCanvas *c1b = new TCanvas("c1b","Facteur de forme time like de Iachello",400,300);
    newmodel->SetRange(0.,2.);
    newmodel->Draw();


    //
    // Go to the event loop
    // File output is supressed, the projector is used instead
    // to save the analysis macro
    // 

    PReaction *r = new PReaction(3.13,"p","p","p D+ [dilepton [e+ e-] p]");
    
    //Create my histogram:
    TH1F * histo1 = new TH1F ("histo1","DiLepton mass",100,0.,2.0);
    histo1->Sumw2();
    //Create the container of the histogram list
    PProjector *m1 = new PProjector(); 
    //Dilepton mass
    m1->AddHistogram(histo1,"_x=[dilepton]->M()");


    
    r->AddBulk(m1);


    r->Print();
    r->Preheating(500);
    r->loop(5000);
    
    r->PrintReport();
     
    TCanvas *c4 = new TCanvas("c4"," Dalitz decay ",400,300);

    histo1->Draw("e1");
    histo1->SetMinimum(0.00000000001);
    cout << histo1->Integral(0,1000) << endl;

    make_delta_dalitz_fig();
    
    
}

