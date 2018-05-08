
#include "../doc/make_delta_dalitz_fig.C"
#include "../doc/make_delta_partial_contribution_fig.C"

pp_delta_dalitz_krivoruchenko() {
    gROOT->Reset();
    gStyle->SetOptLogy(1);

    

    Int_t use_weighting=kTRUE;
    //use_weighting=kFALSE;

    //
    // We use the weighting option with a flat di-lepton generator
    //
    
    makeDistributionManager()->Exec("elementary");
    
    //
    // Compile and add the new form factor model
    //
    
    gSystem->CompileMacro( "./PDeltaDalitzFF.C");

    PDeltaDalitzFF * newmodel = new PDeltaDalitzFF("iachello@D+_to_p_dilepton/formfactor",
						   "Iachello ff for D+ -> p e+e-",-1);

    
    gSystem->CompileMacro( "./PDeltaDalitzKrivoruchenko.C");

    PDeltaDalitzKrivoruchenko * newmodel2 = new PDeltaDalitzKrivoruchenko("kriv@D+_to_p_dilepton",
									  "dgdm from Krivoruchenko",-1);

    if (use_weighting) {
	
	newmodel2->EnableWeighting();
	newmodel2->SetExpectedWeightMean(-1);
	makeDistributionManager()->GetDistribution("D+_dalitz")->EnableWeighting();
	//Use Delta Dalitz model for weighting (also the old one if new one is disabled)
	
	TF1 *flat = new TF1("flat","1",0,1);
    
	PInclusiveModel *dilepton_generator = new PInclusiveModel("flat@D+_to_p_dilepton/generator",
								  "Dilepton generator",-1);    
	dilepton_generator->Set("dilepton,primary");
	dilepton_generator->SetSampleFunction(flat);
	dilepton_generator->EnableGenerator();
	
	makeDistributionManager()->Add(dilepton_generator);

    }

    makeDistributionManager()->Add(newmodel);
    makeDistributionManager()->Add(newmodel2);
    makeDistributionManager()->Disable("iachello");
    makeDistributionManager()->Print();
    

    //
    // One can draw it!
    // 

//     TCanvas *c1b = new TCanvas("c1b","Facteur de forme time like de Iachello",400,300);
//     newmodel->SetRange(0.,2.);
//     newmodel->Draw();


    //
    // Go to the event loop
    // File output is supressed, the projector is used instead
    // to save the analysis macro
    // 

    PReaction *r = new PReaction("1.25","p","p","p D+ [dilepton [e+ e-] p]");
    
    //Create my histogram:
    TH1F * histo1 = new TH1F ("histo1","DiLepton mass",100,0.,0.6);
    histo1->Sumw2();
    //Create the container of the histogram list
    PProjector *m1 = new PProjector(); 
    //Dilepton mass
    m1->AddHistogram(histo1,"_x=[dilepton]->M()");


    
    r->AddBulk(m1);


    r->Print();

#if 0
    r->Preheating(500);
    r->loop(10000);
    
    r->PrintReport();
     
    TCanvas *c4 = new TCanvas("c4"," Dalitz decay ",400,300);

    PUtils::correct(histo1);

    histo1->Draw("e1");
//    histo1->SetMinimum(0.0000001);


#endif

//The following lines are not needed, they just compare the old vs. new spectral functions
    

#if 1
    make_delta_dalitz_fig();

     newmodel2->SetNpx(200);
     newmodel2->SetLineColor(3);

     newmodel2->SetParameter(1,-1);
     newmodel2->SetMinimum(0.000000000001);
     
     PDeltaDalitzKrivoruchenko * newmodel2x= newmodel2->Clone();
     newmodel2x->SetLineStyle(10);
     newmodel2x->SetUseQED(1);
     
     cout << newmodel2x->Integral(0.,1.) << endl;
//newmodel2x->Draw("same");
     
     PDeltaDalitzKrivoruchenko * newmodel23= newmodel2->Clone();
     newmodel23->SetLineStyle(7);
     newmodel23->draw_parent_mass=1.8;
     
     newmodel23->Draw("same");
     
     PDeltaDalitzKrivoruchenko * newmodel22= newmodel2->Clone();
     newmodel22->SetLineStyle(9);
     newmodel22->draw_parent_mass=1.5;
     newmodel22->Draw("same");
     
     
     
     
     newmodel2->Draw("same");
#endif
     //   make_delta_partial_contribution_fig();


    
}

