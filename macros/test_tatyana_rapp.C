{
    
    //Open the root file with the matrices:
    TFile *matrices = new TFile("pbpb20gev.LMR.root");
    TH2D  *h2_lmr_vac_copy = h2_lmr_vac;
    TFile *matrices2 = new TFile("mom_dist_rho0_auau20gev_new.root");
    TH2D  *pt_vs_pz;
    matrices2->GetObject("sum pt vs pz (rho0)", pt_vs_pz);
    //pt_vs_pz = new TH2D("pt_vs_pz", "pt pz", 4, -5, 5, 4, 0, 5);pt_vs_pz->SetBinContent(1, 1, 1);pt_vs_pz->SetBinContent(2, 2, 1);pt_vs_pz->SetBinContent(3, 3, 1);pt_vs_pz->SetBinContent(4, 4, 1);

    PReaction my_reaction("out");

    TH2D *filled  = new TH2D("filled",  "m pt",  90, 0, 1.59999, 101, 0, 5.02499);
    TH2D *filled2 = new TH2D("filled2", "pt pz", 50, -5, 5, 50, 0, 5);
    
    PProjector *input = new PProjector();
    input->AddHistogram(h2_lmr_vac_copy, "m=0; pt=0; GetRandom(m, pt)");
    input->AddHistogram(pt_vs_pz, "pz = GetRandomX(pt)");
    //construct the dilepton:
    input->Do("dilepton = P3M(pt, 0., pz, m);");
    input->Do("push(dilepton)");
    my_reaction.AddPrologueBulk(input); //The "prolog" is done before the decay

    my_reaction.SetDecayAll(1.);

    my_reaction.Do(filled, "_x = m; _y = pt");
    my_reaction.Do(filled2, "_x = pz; _y = pt");
    //my_reaction.Do("echo $m, $pt, $pz");


    
    


    my_reaction.Loop(1000);

}
