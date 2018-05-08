//TITLE Use the dummy sampling model

{

//Load the model into the root session
//Here, only compiled code works, the ".L" does not work
gSystem->CompileMacro( "./PSampleModel.C");

//Construct the model from the class
//Use the parser
PSampleModel *newmodel = new PSampleModel("flat@dilepton_x", "My private model",-1);

//make it known to the Pluto world:
//makeDistributionManager()->Add(newmodel);
//makeDistributionManager()->Print("user");

 TH2F *hf=new TH2F("hf","X vs. Y",100,0.,1.,100,0.,1.);

PAdaptiveMeshN bla(0x3, 2, newmodel, 1.);
 bla.SetRange(0,0.,1.);
 bla.SetRange(1,0.,1.);
 bla.SetThreshold(1.1,0.01);
 bla.SetMCPoints(100000);
 bla.Divide(5,2);
 //bla.Divide(2,0);



TStopwatch timer;
timer.Start();

//bla.PrintMesh();exit(1);

// for (int i=0;i<3000000;i++) {

 cout << "starting loop" << endl;

 for (int i=0;i<10000;i++) {
     cout << i << endl;
     bla.GetRandom();
     hf->Fill(bla.GetArrayValue(0),bla.GetArrayValue(1));
     //cout << bla.GetArrayValue(0) << ":" << bla.GetArrayValue(1) << endl; 
 }

// bla.PrintMesh();

timer.Stop();
printf("CPU time=%f\n",timer.CpuTime());


 
 hf->Draw("surf4");

}
