
{
    gStyle->SetPalette(1,0);

    Float_t Eb   = 10;
    //  Float_t T1    = 0.150;   // temperature in GeV
    Float_t T1    = 0.255;   // temperature in GeV --------> Goran
    Float_t T2    = 0.;      // temperature in GeV
    Float_t blast = 0.;      // radial expansion velocity
    Int_t const num_of_react=40000;
    
    //PUtils::SetSeed(0);	 ---> not longer needed, setted to systime by Pluto, see "Info in <PUtilsREngine::PUtilsREngine>: Random seed set to...."
    
    PFireball *source_JPsi=new PFireball("J/Psi",Eb,T1,T2,1.,blast,0.,0.,0.,0.5);
    source_JPsi->setSigma(0.8);
    //source_JPsi->setSigma(0.23); //  --------> Goran
    //setting flow (v1,v2)
    source_JPsi->setFlow(0.2,0.3);
    //source_JPsi->setFlow(0.,0.);

    cout<<source_JPsi->getV2()<<endl;
    
    source_JPsi->Print();
    
    PParticle *JPsi = new PParticle("J/Psi");
    //  JPsi->SetM(3.686); //Psi'
    PParticle *mumJPsi = new PParticle("mu-");
    PParticle *mupJPsi = new PParticle("mu+");
    
    PParticle* s_JPsi[]    = {source_JPsi,JPsi};
    PChannel* c_sJPsi   = new PChannel(s_JPsi,1,1);
    
    PParticle *s_JPsidimu[]  ={JPsi,mumJPsi,mupJPsi};
    PChannel  *c_JPsidimu  = new PChannel(s_JPsidimu,2,1);
    PChannel  *cc_JPsi[]   = {c_sJPsi,c_JPsidimu};
    
    PReaction *r_JPsi = new PReaction(cc_JPsi,"Jpsi_1k_10GeV",sizeof(cc_JPsi)/sizeof(cc_JPsi[0]),0,0,0,0);
 
    //The following lines provide some inline histogram production
    r_JPsi->allParticles();  //required to look for unstables
    TH2F * histo2 = new TH2F ("histo2","Rap. vs. Pt",50,-1.5,3.5,50,0,2.);
    TH1F * histo1 = new TH1F ("histo1","Phi",100,-3.5,3.5);
    r_JPsi->Do(histo2,"foreach(J/Psi); _x = [J/Psi]->Rapidity(); _y=[J/Psi]->Pt(); ");
    r_JPsi->Do(histo1,"foreach(J/Psi); _x = [J/Psi]->Phi();");
   
    r_JPsi->Print();
    r_JPsi->setHGeant(0);
    r_JPsi->loop(num_of_react);
    //histo2->Draw("colz");
    histo1->Draw("");


}


