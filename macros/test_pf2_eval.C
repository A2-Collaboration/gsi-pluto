{

    TFile *f = new TFile("histo.root");
    TH2F *hf1old= new TH2F(*hf1);

    makeDynamicData()->SetBatchHistogram("testhisto",hf1old);
 
    PF2EvalBatch *p = new PF2EvalBatch("p",0,0.3,0,0.3);

    //p->AddEquation("_f = _x + _y");
    //p->AddHistogram(hf1old,"_f = Eval(_x,_y)");
    p->AddEquation("_f = testhisto->Eval(_x,_y,0,4)");

    p->Draw("colz");


}
