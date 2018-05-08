double f_my_function(double *x, double * par)
{
    if (x[0] > 0.99) return 1;
    return 0;

  // (1+x*x)/2
    return (1+x[0]*x[0])/2;
    // It is important that the test fuction does not exceed 1 (but comes close to it)
}


useAngularDistribution_test_kaon() {

    Double_t x[3], y[3];
		 


y[0]=0.3;y[1]=0.6;y[2]=1.;

x[0]=-1.;x[1]=0;x[2]=1.;

angles3 = new TGraph(3,x,y);

TF1 *angles2=new TF1("angles2",f_my_function,-1,1,1);



PAngularDistribution *ang = new PAngularDistribution("my_angle","My angular dist");
ang->Add("q,parent,reference");
ang->Add("n,daughter");
ang->Add("K+,daughter,primary");
ang->Add("n,grandparent,base_reference");
//ang->SetAngleFunction(angles3,kFALSE,kFALSE);
ang->SetAngleFunction(angles2);
ang->Print();

makeDistributionManager()->Add(ang);
TH1F * histo1 = new TH1F ("histo1","cos theta of pp",20,-1.,1.);

//PReaction my_reaction (3,"K+","n","n K+ ");
 PParticle q =  PParticle("K+",20,0,1) + PParticle("n",0,30,0);
PReaction my_reaction (&q,"n K+ ");

//first copy the particles:
 my_reaction.Do("n=[n]; q = [K+ + n]; K=[K+]; parent=n+K");
 // my_reaction.Do("parent->Print(); K ->Print()");
 //Reconstruct the K+ beam
 my_reaction.Do("K_beam = P3M(20.000000,0.000000,1.000000,20.047014)");
 //boost into parent rest frame
 my_reaction.Do("K_beam->Boost(parent); q->Boost(parent); K->Boost(parent); n->Boost(parent);");
 //rotate such that K_beam point to z-axis
 //my_reaction.Do("K->Rot(K_beam); n->Rot(K_beam);");
 // my_reaction.Do(histo1,"_x=((-1)*(cos(K->Theta())))");
 my_reaction.Do("angle = K->Angle(K_beam);");
my_reaction.Do(histo1,"_x= cos(angle)");


makeDistributionManager()->Print();
my_reaction.Print();
my_reaction.Loop(10000);
histo1->Draw();
 histo1->SetLineColor(2);

TH1F * histo2 = new TH1F ("histo2","cos theta of pp",20,-1.,1.);
PReaction my_reaction2 (3,"K+","n","n K+ ");
my_reaction2.Do(histo2,"p1=[K+,1]; q = [K+ + n]; p1->Boost(q); _x=cos(p1->Theta())");



//my_reaction2.Loop(20000);
histo2->Draw("same");
}
