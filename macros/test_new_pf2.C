//TITLE Macro to test the PF2 extension
{

    Double_t x[10]= { 0,   0,   0.8, 0.8, 0.3,  0.7, 0.13, 0.523, 0.4234, 0.8 };
    Double_t y[10]= { 0,   0.8, 0,   0.8, 0.32, 0.2, 0.53, 0.132, 0.78,   0.523};
    Double_t z[10]= { 0.2, 0.4, 0.2, 0.7, 0.5,  0.6, 0.7,  0.3,   0.8,    1};

    TGraph2D * gra = new TGraph2D(10,x,y,z);
    gra->SetNpx(300);
    gra->SetNpy(300);
    //gra->Draw("surf4");

    Double_t xd[11]= { 0, 0.1, 0.2, 0.3, 0.4, 0.5,  0.55, 0.6, 0.65, 0.7, 0.8 };
    TH1D * hist1 = new TH1D("hist1", "Hist 1", 10, xd);
    hist1->SetBinContent(1,0.1);
    hist1->SetBinContent(2,0.2);
    hist1->SetBinContent(3,0.3);
    hist1->SetBinContent(4,0.4);
    hist1->SetBinContent(5,0.3);
    hist1->SetBinContent(6,0.2);
    hist1->SetBinContent(7,0.25);
    hist1->SetBinContent(8,0.3);
    hist1->SetBinContent(9,0.4);
    hist1->SetBinContent(10,0.5);

    TH2D * hist2 = new TH2D("hist2", "Hist 2", 3, 0, 0.8, 3, 0, 0.8);
    hist2->SetBinContent(1,1,0.1);
    hist2->SetBinContent(2,1,0.2);
    hist2->SetBinContent(3,1,0.3);
    hist2->SetBinContent(1,2,0.2);
    hist2->SetBinContent(2,2,0.3);
    hist2->SetBinContent(3,2,0.2);
    hist2->SetBinContent(1,3,0.25);
    hist2->SetBinContent(2,3,0.35);
    hist2->SetBinContent(3,3,0.3);
    

    

    PF2 * only_tgraph = new PF2("Only TGraph2D", 0, 0.8, 0, 0.8);
    only_tgraph->Add(gra, "_f = Eval(_x,_y)");
    only_tgraph->SetNpx(100);
    only_tgraph->SetNpy(100);
    only_tgraph->Draw("surf4");

    
    PF2 * tgraph_and_hist1 = new PF2("TGraph2D mixed with a TH1", 0, 0.8, 0, 0.8);
    tgraph_and_hist1->Add(gra, "_f = Eval(_x,_y)");
    tgraph_and_hist1->Add(hist1, "if (_x < 0.3) _f = Eval(_y)");
    tgraph_and_hist1->SetNpx(100);
    tgraph_and_hist1->SetNpy(100);
    tgraph_and_hist1->Draw("surf4");

    PF2 * tgraph_and_hist1a = new PF2("TGraph2D mixed with a TH1 and some dummy calculations", 0, 0.8, 0, 0.8);
    tgraph_and_hist1a->Add(gra, "_f = Eval(_x,_y)");
    tgraph_and_hist1a->Add("my_x = _x - _y; my_y = (_x + _y) * 0.5; ");
    tgraph_and_hist1a->Add(hist1, "if (my_x > 0.1 && my_x < 0.35) _f = Eval(my_y)");
    tgraph_and_hist1a->SetNpx(300);
    tgraph_and_hist1a->SetNpy(300);
    tgraph_and_hist1a->Draw("surf4");

    
    PF2 * tgraph_and_hist2 = new PF2("TGraph2D merged with a TH2", 0, 0.8, 0, 0.8);
    tgraph_and_hist2->Add(gra, "_f = Eval(_x,_y)");
    tgraph_and_hist2->Add(hist2, "_f = _f * Eval(_x, _y)");
    tgraph_and_hist2->SetNpx(300);
    tgraph_and_hist2->SetNpy(300);
    tgraph_and_hist2->Draw("surf4");
    

}
