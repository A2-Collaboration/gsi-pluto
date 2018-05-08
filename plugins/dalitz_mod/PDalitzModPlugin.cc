/////////////////////////////////////////////////////////////////////
//Plugin to modify the Dalitz decays
//
//This plugin is part of the distribution manager
//It can be activated with 
//makeDistributionManager()->Exec("dalitz_mod: command");
//
//where the plugin supports the following commands:
//static_br_thresh=value   Threshold for disabling static br in GeV
//flat_generator           Enables a flat generator for all registered 
//                         Dalitz decays. Instead of the static br we
//                         use the dG/dm directly, if the parent res
//                         width is > static_br_thresh
//krivoruchenko            Use Delta Dalitz from Krivoruchenko [L1]
//
//                         Author:  I. Froehlich
//                         Written: 17.9.2008
//
//References:
//[L1] M.I. Krivoruchenko, A. Faessler (Tubingen U.), Phys.Rev.D65:017502,2002, nucl-th/0104045
//[L2] M. Zetenyi and Gy. Wolf, Heavy Ion Phys.17:27-39,2003, arXiv:nucl-th/0202047v1
//
//////////////////////////////////////////////////////////////////////

#include "PDalitzModPlugin.h"

#include "PDeltaDalitzFF.h"
#include "PResonanceDalitz.h"

PDalitzModPlugin::PDalitzModPlugin(const Char_t *id, const Char_t *de):
    PDistributionCollection(id, de) {
    
    resonances_done = kFALSE;
    RequiresPlugin("pdg:extend_resonances");
}

Bool_t PDalitzModPlugin::Activate(void) {

    //    PDistributionManagerUtil * pdmutil = makeDistributionManagerUtil();
    static_br_thresh = 1111.;  //something very large -> disabled
    PDistributionManagerUtil *pdmutil = makeDistributionManagerUtil();
    kriv1 = kriv2 = NULL;
    pdmutil->SetGroup("decay_models");
    return kTRUE;
}

PDalitzModPlugin::~PDalitzModPlugin() {
}

Bool_t PDalitzModPlugin::ExecCommand(const char *command, Double_t value) {


    //Resonances:
    if (!resonances_done) {
	PDistributionManagerUtil *pdmutil = makeDistributionManagerUtil();

	//model from Ref. [L2]

	pdmutil->AddSubGroup("dalitz_zetenyi_wolf", "Dalitz decay from Zetenyi/Wolf, nucl-th/0202047v1", "decay_models");
	pdmutil->SetGroup("dalitz_zetenyi_wolf");

	makeStaticData()->AddDecay(-1, "N*(1520)+ -> p + dilepton", "ND13+", "p, dilepton", 4.0e-5);
	
	PResonanceDalitz *newmodel = 
	    new PResonanceDalitz("ND13+_dalitz@ND13+_to_p_dilepton",
				 "dgdm from Zetenyi/Wolf", -1);
	newmodel->setGm(0.793);
	pdmutil->Add(newmodel);

	makeStaticData()->AddDecay(-1, "N*(1440)+ -> p + dilepton", "N*(1440)+", "p, dilepton", 1.1e-5);
	newmodel = new PResonanceDalitz("NP11+_dalitz@NP11+_to_p_dilepton",
					"dgdm from Zetenyi/Wolf", -1);
	newmodel->setGm(0.139);
	pdmutil->Add(newmodel);


	makeStaticData()->AddDecay(-1,"N*(1535)+ -> p + dilepton", "N*(1535)+", "p, dilepton", 1.6e-5);
	newmodel = new PResonanceDalitz("NS11+_dalitz@NS11+_to_p_dilepton",
					"dgdm from Zetenyi/Wolf", -1);
	newmodel->setGm(0.634);
	pdmutil->Add(newmodel);


	makeStaticData()->AddDecay(-1, "N*(1520)0 -> n + dilepton", "ND130", "n, dilepton", 4.0e-5);
	newmodel = new PResonanceDalitz("ND130_dalitz@ND130_to_n_dilepton",
					"dgdm from Zetenyi/Wolf", -1);
	newmodel->setGm(0.719);
	pdmutil->Add(newmodel);

	makeStaticData()->AddDecay(-1, "N*(1440)0 -> n + dilepton", "N*(1440)0", "n, dilepton", 1.1e-5);
	newmodel = new PResonanceDalitz("NP110_dalitz@NP110_to_n_dilepton",
					"dgdm from Zetenyi/Wolf", -1);
	newmodel->setGm(0.098); 
	pdmutil->Add(newmodel);


	makeStaticData()->AddDecay(-1, "N*(1535)0 -> n + dilepton", "N*(1535)0", "n, dilepton", 1.6e-5);
	newmodel = new PResonanceDalitz("NS110_dalitz@NS110_to_n_dilepton",
					"dgdm from Zetenyi/Wolf", -1);
	newmodel->setGm(0.549);
	pdmutil->Add(newmodel);

	


	// setGm(0 , 1.98  , 1.980, 1, 3);              // D1232
	// setGm(1 , 0.098 , 0.139, 1, 1);              // N1440
	// setGm(2 , 0.719 , 0.793,-1, 3);              // N1520
	// setGm(3 , 0.549 , 0.634,-1, 1);              // N1535
	// setGm(4 , 0.285 , 0.315,-1, 1);              // N1650
	// setGm(5 , 1.52  , 0.678,-1, 5);              // N1675
	// setGm(6 , 0.971 , 2.74 , 1, 5);              // N1680
	// setGm(7 , 0.193 , 0.126,-1, 3);              // N1700
	// setGm(8 , 0.019 , 0.030, 1, 1);              // N1710
	// setGm(9 , 0.386 , 0.193, 1, 3);              // N1720
	// setGm(10, 0.162 , 0.162, 1, 3);              // D1600
	// setGm(11, 0.162 , 0.162,-1, 1);              // D1620
	// setGm(12, 0.549 , 0.549,-1, 3);              // D1700
	// setGm(13, 0.713 , 0.713, 1, 5);              // D1905
	// setGm(14, 0.066 , 0.066, 1, 1);              // D1910
	// setGm(15, 0.479 , 0.479,-1, 5);              // D1930 
	
	resonances_done = kTRUE;
    }


    //only by command:
    if (strcmp (command,"flat_generator") == 0) {

	PDistributionManagerUtil *pdmutil = makeDistributionManagerUtil();
	pdmutil->LinkDB();

	pdmutil->AddSubGroup("generators", "Generator models", "root");
	pdmutil->SetGroup("generators");

	TF1 *flat = new TF1("flat", "1", 0, 1);

	Int_t header = makeDataBase()->GetEntry("std_set");
	Int_t particlekey = -1;
	//      Int_t generator_key=makeStaticData()->MakeDirectoryEntry("modeldef","generator");

	//loop over particles
	while (makeDataBase()->MakeListIterator(header, "snpart", "slink", &particlekey)) {
	
	    Int_t decaykey = -1;
	    Int_t pid = makeStaticData()->GetParticleIDByKey(particlekey);

	    //loop over decays
	    while ((makeStaticData()->GetParticleNChannelsByKey(particlekey)>0) && 
		   (makeDataBase()->MakeListIterator(particlekey, "pnmodes", "link", &decaykey))) {
		Int_t tid[11];
		tid[0] = 10; 
		makeStaticData()->GetDecayModeByKey(decaykey, tid); // retrieve current mode info

		if (PData::IsDalitz(pid, tid[1], tid[2])) {
	    
		    TString *id = new TString(makeStaticData()->GetParticleName(pid));
		    id->Append("_generator_");
		    for (int p=1; p<=tid[0]; p++) {
			id->Append(makeStaticData()->GetParticleName(tid[p]));
			if (p != tid[0]) id->Append("_");
		    }
		    id->Append("@");
		    id->Append(makeStaticData()->GetParticleName(pid));
		    id->Append("#generator#");
		    for (int p=1; p<=tid[0]; p++) {
			id->Append(makeStaticData()->GetParticleName(tid[p]));
			if (p != tid[0]) id->Append("#");
		    }
		    id->Append("&generator");

		    PInclusiveModel *dilepton_generator = 
			new PInclusiveModel((char*)id->Data(), "Dilepton generator", -2);    
		    dilepton_generator->Set("dilepton,primary");
		    dilepton_generator->SetSampleFunction(flat);
		    dilepton_generator->EnableGenerator();
	    
		    pdmutil->Add(dilepton_generator);

		    PChannelModel *pmodel = makeDynamicData()->GetDecayModelByKey(decaykey);
		    if (pmodel) {
			pmodel->EnableWeighting();
			pmodel->SetVersionFlag(VERSION_GENERATOR_MC);
			//disable re-scaling if parent width is larger
			//then static_br_thresh
			if (makeStaticData()->GetParticleTotalWidth(pid) > static_br_thresh) {
			    pmodel->SetExpectedWeightMean(-1);
			    Info("ExecCommand","Model <%s> uses dGamma/dM for the branching ratio",
				 pmodel->GetIdentifier());
			} 
		    } else {
			Warning("ExecCommand", "Primary model not found");
		    }
		} //isDalitz
	    }
	}
      
	return kTRUE;

    } else if (strcmp (command,"static_br_thresh") == 0) {

	static_br_thresh = value;
	return kTRUE;

    } else if (strcmp (command,"krivoruchenko") == 0) {

	if (kriv1) return kTRUE;
	PDistributionManagerUtil *pdmutil = makeDistributionManagerUtil();
	pdmutil->LinkDB();
	pdmutil->SetGroup("decay_models");
	kriv1 = new PDeltaDalitzKrivoruchenko("D+_krivoruchenko@D+_to_p_dilepton",
					      "dgdm from Krivoruchenko", -1);

	kriv2 = new PDeltaDalitzKrivoruchenko("D0_krivoruchenko@D0_to_n_dilepton",
					      "dgdm from Krivoruchenko", -1);
	pdmutil->Add(kriv1);
	pdmutil->Add(kriv2);

	//Add VMD-FF in the vmd group
	pdmutil->AddSubGroup("vmd", "VMD form factors", "root");
	pdmutil->SetGroup("vmd");

	PDeltaDalitzFF *vmd_newmodel = new PDeltaDalitzFF("D+_iachello_ff@D+_to_p_dilepton/formfactor",
							  "Iachello ff for D+ -> p e+e-", -1);
	pdmutil->Add(vmd_newmodel);
	vmd_newmodel = new PDeltaDalitzFF("D0_iachello_ff@D0_to_n_dilepton/formfactor",
					  "Iachello ff for D0 -> n e+e-", -1);
	pdmutil->Add(vmd_newmodel);
	

	//Add QED-FF in the qed group
	pdmutil->AddSubGroup("qed", "QED form factors", "root");
	pdmutil->SetGroup("qed");

	PDeltaDalitzFF *qed_newmodel = new PDeltaDalitzFF("D+_qed_ff@D+_to_p_dilepton/formfactor",
							  "QED ff for D+ -> p e+e-", -1);
	qed_newmodel->SetQED(1);
	qed_newmodel->SetCC(3.0, 0, 0);
	pdmutil->Add(qed_newmodel);
	qed_newmodel = new PDeltaDalitzFF("D0_qed_ff@D0_to_n_dilepton/formfactor",
					  "QED ff for D0 -> n e+e-", -1);
	qed_newmodel->SetQED(1);
	qed_newmodel->SetCC(3.0,0,0);
	pdmutil->Add(qed_newmodel);
	return kTRUE;
    } 
	
    return kFALSE;
}



ClassImp(PDalitzModPlugin)



    
