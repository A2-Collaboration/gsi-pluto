{

    //Open the root file with the matrices:
    TFile *filters_ep = new TFile("../../matricesEffSingle.PairCode.all.cuts.ep.may07_v2.root");
    //Get the histos we need:
    TH3F *efficiency_eplus = effi3DPosiAllCut->Clone("efficiency_eplus");
    TH3F *acceptance_eplus = acce3DPosi->Clone("acceptance_eplus");


    TFile *filters_ep = new TFile("../../matricesEffSingle.PairCode.all.cuts.em.may07_v2.root");
    TH3F *efficiency_eminus = effi3DEleAllCut->Clone("efficiency_eminus");
    TH3F *acceptance_eminus = acce3DEle->Clone("acceptance_eminus");




    TFile *f = new TFile("pluto_ee_filter.root","RECREATE");

    PCommandList *p = new PCommandList("_main", "echo **** HADES ee acceptance/efficiency filter for May07, v2");
    
    p->AddCommand("echo ********************************************************************************");
    p->AddCommand("echo Output:     calculates the variables eff, bool_eff, acc, bool_acc");
    p->AddCommand("echo ********************************************************************************");
    p->AddCommand("echo Usage:      This filter works in 2 modi:");
    p->AddCommand("echo Usage:      1.) the event filter mode:");
    p->AddCommand("echo Usage:          e.g. #eefilter=bool_eff");
    p->AddCommand("echo Usage:      2.) as a particle remover:");
    p->AddCommand("echo Usage:          makeDistributionManager()->Startup(\"_filter_remove_particles=1\");");
    p->AddCommand("echo ********************************************************************************");
    p->AddCommand("echo Debug mode:     makeDistributionManager()->Startup(\"_filter_debug=1\");");
    p->AddCommand("echo ********************************************************************************");
    p->AddCommand("_filter_debug=0");
    p->AddCommand("_filter_remove_particles=0");

    p->Write();

    
    PCommandList *s = 
	new PCommandList("_startup", 
			 "echo ********************************************************************************");
    
    s->AddCommand("echo This event loop uses the HADES ee acceptance/efficiency filter for May07, v2 ");
    s->AddCommand("echo ********************************************************************************");
    s->Write();


    //Labels should be unique. Avoid common names
    PCommandList *l = new PCommandList("_loop");
    l->AddCommand("if (_filter_debug); echo *************** new event");
    l->AddCommand("if (_filter_remove_particles); goto _hades_may07_filter_remove_particles;");
    l->AddCommand("bool_eff = 1");
    l->AddCommand("bool_acc = 1");

    //event remover:
    l->AddCommand("_hades_may07_filter_label_eplus:");
    l->AddCommand("_x = [e+]->Phi(); if (_x > TMath::Pi()); _x = _x - (2.*TMath::Pi())");
    l->AddCommand("_y = [e+]->Theta(); ");
    l->AddCommand("_z = (([e+]->P()) * 1000.); if ([e+]->P() > 1.9) ; _z = 1900");
    l->AddCommand("if (_filter_debug); foreach e+; echo e+: $_x, $_y, $_z, acc: $acc_ep");
    l->AddCommand(efficiency_eplus,"eff_ep = Eval(); rnd = [e+]->Sample(); if (rnd>eff_ep); bool_eff = 0");
    l->AddCommand(acceptance_eplus,"acc_ep = Eval(); rnd = [e+]->Sample(); if (rnd>acc_ep); bool_acc = 0");
    l->AddCommand("if (_y < 0.1) || (_y > 1.5)  ; bool_acc = 0 ");
    l->AddCommand("if (_filter_debug); foreach e+; echo e+: $_x, $_y, $_z eff: $eff_ep, acc: $acc_ep");
    l->AddCommand("formore e+; echo This filter supports only one eplus");

    l->AddCommand("_hades_may07_filter_label_eminus:");
    l->AddCommand("_x = [e-]->Phi(); if (_x > TMath::Pi()); _x = _x - (2*TMath::Pi())");
    l->AddCommand("_y = [e-]->Theta(); ");
    l->AddCommand("_z = [e-]->P() * 1000.; if ([e-]->P() > 1.8) ; _z= 1800");
    l->AddCommand(efficiency_eminus,"eff_em = Eval(); rnd = [e-]->Sample(); if (rnd>eff_em); bool_eff = 0");
    l->AddCommand(acceptance_eminus,"acc_em = Eval(); rnd = [e-]->Sample(); if (rnd>acc_em); bool_acc = 0");
    l->AddCommand("if (_y < 0.1) || (_y > 1.5)  ; bool_acc = 0 ");
    l->AddCommand("if (_filter_debug); foreach e-; echo e-: $_x, $_y, $_z eff: $eff_em, acc: $acc_em");
    l->AddCommand("formore e-; echo This filter supports only one eminus");
    l->AddCommand("goto _filter_end");
    
    //particle remover:
    l->AddCommand("_hades_may07_filter_remove_particles:");
    l->AddCommand("_x = [e-]->Phi(); if (_x > TMath::Pi()); _x = _x - (2*TMath::Pi())");
    l->AddCommand("_y = [e-]->Theta(); ");
    l->AddCommand("_z = [e-]->P() * 1000.; if ([e-]->P() > 1.8) ; _z = 1800");
    l->AddCommand(acceptance_eminus,"acc_em = Eval(); rnd = [e-]->Sample(); if (rnd>acc_em); [e-]->SetInActive(); if (_filter_debug); echo removed e- (acc)");    
    //l->AddCommand(efficiency_eminus,"eff_em = Eval(); rnd = [e-]->Sample(); if (rnd>eff_em); [e-]->SetInActive(); if (_filter_debug); echo removed e- (eff)");
    l->AddCommand("if (_filter_debug); echo e-: $_x, $_y, $_z eff: $eff_em, acc: $acc_em");
    l->AddCommand("formore e-; goto _hades_may07_filter_remove_particles");

    l->AddCommand("_hades_may07_filter_remove_particles_eplus:");
    l->AddCommand("_x = [e+]->Phi(); if (_x > TMath::Pi()); _x = _x - (2*TMath::Pi())");
    l->AddCommand("_y = [e+]->Theta(); ");
    l->AddCommand("_z = [e+]->P() * 1000.; if ([e+]->P() > 1.9) ; _z = 1900");
    //l->AddCommand("_z = 0.; ");
    l->AddCommand(acceptance_eplus,"acc_ep = Eval(); rnd = [e+]->Sample(); if (rnd>acc_ep); [e+]->SetInActive(); if (_filter_debug); echo removed e+ (acc)");
    //l->AddCommand(efficiency_eplus,"eff_ep = Eval(); rnd = [e+]->Sample(); if (rnd>eff_ep); [e+]->SetInActive(); if (_filter_debug); echo removed e+ (eff)");
    l->AddCommand("if (_filter_debug); echo e+: $_x, $_y, $_z eff: $eff_ep, acc: $acc_ep");
    l->AddCommand("formore e+; goto _hades_may07_filter_remove_particles_eplus");

    //for pairs -> only single e+e-
    l->AddCommand("opang = [e+]->Angle([e-])");
    //l->AddCommand("if opang < (9./180.)*TMath::Pi(); [e+]->SetInActive(); [e-]->SetInActive();  if (_filter_debug); echo removed pair");


    l->AddCommand("_filter_end:");

    l->Write();
    efficiency_eplus->Write();
    acceptance_eplus->Write();
    efficiency_eminus->Write();
    acceptance_eminus->Write();
    f->Write();
 

}
