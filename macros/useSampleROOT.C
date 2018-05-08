//TITLE Use the dummy sampling model

double f_my_function(double *x, double * par)
{
 
    return pow(sin(x[0]*10+x[1]*x[1]*10)*x[0]*x[1],2);
}


useSampleROOT(){


 TH2F *hf=new TH2F("hf","X vs. Y",150,0.,1.,150,0.,1.);

TF2 *f2=new TF2("f2",(void*) f_my_function,0.,1.,0.,1.);
 f2->SetNpx(200);
 f2->SetNpy(200);

TStopwatch timer;
timer.Start();

 Double_t x,y;

 for (int i=0;i<3000000;i++) {
     
     f2->GetRandom2(x,y);
     hf->Fill(x,y);
     //cout << bla.GetArrayValue(0) << ":" << bla.GetArrayValue(1) << endl; 
 }



timer.Stop();
printf("CPU time=%f\n",timer.CpuTime());


 
 hf->Draw("surf4");

}
