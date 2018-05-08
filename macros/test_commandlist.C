{

    TH1F *histo = new TH1F("histo", "ee mass", 20, 0., 0.6);
    histo->Fill(0.3);

    TFile *f = new TFile("list.root","RECREATE");

    PCommandList *p = new PCommandList("_main", "echo **** Demo acceptance file");

    p->AddCommand("echo Output: .....");
    p->AddCommand("echo Usage: .....");

    p->Write();

    
    PCommandList *s = new PCommandList("_startup", "echo do not forget, that a filter is present...");
    s->Write();


    //Labels should be unique. Avoid common names
    PCommandList *l = new PCommandList("_loop_start");
    l->AddCommand("#protonfilter = 1");
    l->AddCommand("_demo_filter_label:");
    l->AddCommand("if [p]->Theta() < (18*TMath::Pi() / 180) ; #protonfilter = 0; formore p ; goto _demo_filter_label");


    l->AddCommand("#efilter = 1");
    l->AddCommand(histo,"_x = [dilepton]->M(); val = Eval(); if (val<0.5); #efilter = 0");

    l->Write();
    histo->Write();

    f->Write();
 

}
