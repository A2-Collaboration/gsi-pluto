
void smear(PParticle * p, double gamma) {
  // smearing of PParticle 4-vector due to the HADES resolution

const long double pi=3.14159265358979323846; 

  double ptot=(p->Vect()).Mag();
  if (ptot<=0.) return;
  if (gamma<=0.) return;
  double th=p->Theta(), ph=p->Phi(), sth =sin(th), mass=p->M();
  if (fabs(sth)<1.e-6) sth = (sth>0.) ? 1.e-6 : -1.e-6;
  double smeared_ptot = fabs( PUtils::sampleGaus(ptot, gamma ) );
  double smeared_theta = PUtils::sampleGaus(th, 0.0023 );
  double smeared_phi = ( PUtils::sampleGaus(sth*ph, 0.0023 ) )/sth;
  double theta_check=smeared_theta-2.*pi*int(smeared_theta/(2.*pi));
  if (theta_check>pi||theta_check<0.) smeared_phi -= pi;
  p->SetRho(smeared_ptot);                    // keeping theta, phi constant
  p->SetTheta(smeared_theta);                 // keeping |p|, phi constant
  p->SetPhi(smeared_phi);                     // keeping |p|, theta constant
  double newp=p->P();
  p->SetE(sqrt(newp*newp+mass*mass));         // reset energy
}


double GetResolution(PParticle *pPart, int iSetup)
{
int pid = pPart->ID();
double p, mass, res;

    {
        //p = (pPart->Vect()).Mag()/1000.;
        p = (pPart->Vect()).Mag();

        //if(pid == ID_Em || pid == ID_Ep) printf("p %7.4f\n", p);

        mass = pPart->M();
    
        switch(iSetup)
        {
            case 2:
                res = 2.0 + 13.63 * p;
                break;

            case 3:
                res = 1.0 + 3.6 * p;
		break;

           case 4:
	       res = 0.5 + 1.36 * p;  
                break;

            default:
                return 0.1;
        }

//	cout<<"iSetup= "<<iSetup<< "res= "<<res<<endl;

        if(pid ==2 || pid==3) return 0.01*p*res;

        double possq = res*res - 0.6*0.6;         // position resolution
        double betasq = 1.+(mass/p)*(mass/p);
        double scatsq = 0.6*0.6*betasq;           // scale multiple scatt
        res = sqrt(possq + scatsq);

        return 0.01*p*res;
    }

    return 0.0;
}


int selectAc(PParticle* p)
{ // Lepton acceptance & efficiency function

  static Int_t count=0;

  static TH3F *accEplus;
  static TH3F *accEminus;
  

  TFile *iFile = NULL;

   if(count==0) 
   {
     // printf("count: %d \n",count);                                                      
     //iFile =  new TFile("matrices/Acc3D_MediumField_HighRes.root");
     //iFile =  new TFile("matrices/acceptances.root"); 
     //iFile =  new TFile("matrices/Acc3D_FullField_HighRes.root");

     //accEplus = (TH3F*)(gROOT->FindObject("EPlusAcc"));   // read CC acc matrices
     //accEminus = (TH3F*)(gROOT->FindObject("EMinusAcc"));
     //accPiplus = (TH3F*)(gROOT->FindObject("PiPlusAcc"));   // read acc matrices
     //accPiminus = (TH3F*)(gROOT->FindObject("PiMinusAcc"));
     //accProton = (TH3F*)(gROOT->FindObject("ProtonAcc"));

//    pp1.25 GeV specific

     cout<<"read matrices"<<endl;
    iFile =  new TFile("../../matricesEffSingle.PairCode.all.cuts.em.apr06.root");
    accEminus = (TH3F*)(gROOT->FindObject("acce3DEle"));  
    iFile =  new TFile("../../matricesEffSingle.PairCode.all.cuts.ep.apr06.root");
    accEplus =  (TH3F*)(gROOT->FindObject("acce3DPosi"));   // read acc e+/e- matrices             
    cout<<"done"<<endl;
//
   }
   count++;

  if ((int)p==-1)
  {
    count++;
    //printf("= ok <<<\n\n");
    return -1;
  }

  Int_t id = p->ID();     // get particle ID

  if (p->IsActive())         //assume to be same for all particles
  {
      Float_t theta,phi,mom;
      theta=0.0;phi=0.0;mom=0.0;

      if (id==14) {
/* convention used by Witek CC/pp */

    //printf("= active <<<\n\n");            
    theta = 57.29578*(p->Theta());  // polar angle in degree
    phi = 90.-57.29578*(p->Phi());  // azimuthal angle in degree
    if (phi<0.) phi += 360.;                // TLorentzVector uses ATan2(x,y) and returns -pi < Phi() < +pi
    mom = 1000.*(p->P());           // momentum in MeV/c
    
    phi -=  ((Int_t)(phi/60.))*60.;         // move to sector system

    if (theta>90 || theta<5) return 0;

      }
/* convention used by Tetyana & Jochen  pp @ 1.25 GeV*/

      if (id==2 || id==3) {

    theta = p->Theta();      // polar angle in rad
    phi = p->Phi();         // azimuthal angle in rad  
    mom = 1000.*(p->P());  // momentum in MeV/c
  
    if(theta>1.5 || theta<0.1) return 0;
				   }
Float_t testAcc=0.;

    if (id==2)
    {  // check posive charge
      //printf("= active 1 <<<\n\n");
           if (mom>1900.) mom=1900.;
           //cout<<"positron"<<"mom"<<mom<<"phi"<<phi<<"theta"<<theta<<endl;
//     	   testAcc = accEplus->GetBinContent(accEplus->FindBin(mom,theta,phi));  // convention used by Witek
           testAcc = accEplus->GetBinContent(accEplus->FindBin(phi,theta,mom));  // convention used by Tetyana
           //cout<< "positron: acceptance"<< testAcc << endl;

      //      testAcc=1.;
    }

  
    else if (id==3 )
    {  // check negative charge
      //printf("= active 2 <<<\n\n");
	if (mom>1800.) mom=1800.;
//	testAcc = accEminus->GetBinContent(accEminus->FindBin(mom,theta,phi)); // convention used by Witek
//        cout<<" electron"<<"phi"<<phi <<"theta"<< theta<<"mom" << mom << endl;
	testAcc = accEminus->GetBinContent(accEminus->FindBin(phi,theta,mom)); // convention used by Tetyana              //       cout<<"acceptance"<<testAcc<<endl;
       //testAcc=1.;
    }


    if (PUtils::sampleFlat()<testAcc)
    {
   mom /= 1000.; 
      return 1;
    }
  }
  return 0;
}


Float_t getEfficiencyFactor (float phi1, float theta1, float mom1, float chrg1, float phi2, float theta2, float mom2, float chrg2)
{ 
  static Int_t count1=0;
  static TH3F *p3DeffEleLow;
  static TH3F *p3DeffPosiLow;
  Float_t pi=3.14159;

  if (count1==0) { 

   // Init efficiency matrices pp1.25 GeV Tetyana & Jochen

  TFile *pLowRangeEffFile = new TFile("../../matricesEffSingle.PairCode.all.cuts.em.apr06.root","READ");
   if (pLowRangeEffFile)  
   { 
       pLowRangeEffFile->cd(); 
       p3DeffEleLow = (TH3F*) pLowRangeEffFile->Get("effi3DEleAllCut");  
   } 
   else Error("Init","pointer to eff matrix file is NULL"); 
 
   TFile *pHighRangeEffFile = new TFile("../../matricesEffSingle.PairCode.all.cuts.ep.apr06.root","READ");
   if (pHighRangeEffFile)  
   { 
       pHighRangeEffFile->cd(); 
       p3DeffPosiLow = (TH3F*) pHighRangeEffFile->Get("effi3DPosiAllCut");  
   } 
   else Error("Init","pointer to eff matrix file is NULL"); 
      }

 
  count1++;

  mom1=mom1*1000;
  mom2=mom2*1000;
 
  if(mom1>1990) mom1=1990;
  if(mom2>1990) mom2=1990;
 
    Float_t fEff1 = 1.; 
    Float_t fEff2 = 1.;  
    Float_t r2d = TMath::RadToDeg(); 

    Float_t phi1_tmp = phi1; 
    Float_t phi2_tmp = phi2; 
  
    if (chrg1==2) // positron 
    { 
 
	    fEff1 = p3DeffPosiLow-> 
	    GetBinContent(p3DeffPosiLow->FindBin(phi1_tmp,theta1,mom1)); 
//	    *eff_p=fEff1;
// if (mom1>400) cout << " phi1_tmp " << phi1_tmp  <<" theta "<< theta1*59.1 <<"  mom1  "<<mom1<< " eff1  " << fEff1 << endl; 
   } 
    else if (chrg1==3) // electron 
    { 
 
	    fEff1 = p3DeffEleLow-> 
		GetBinContent(p3DeffEleLow->FindBin(phi1_tmp,theta1,mom1)); 
    } 
  
 
    // second leg of the pair 
    if (chrg2==2) // positron 
    { 
	 
	    fEff2 = p3DeffPosiLow-> 
		GetBinContent(p3DeffPosiLow->FindBin(phi2_tmp,theta2,mom2)); 
//	    *eff_p=fEff2;
 
    } 
    else if (chrg2==3) // electron 
    { 
 
	    fEff2 = p3DeffEleLow->GetBinContent(p3DeffEleLow->FindBin(phi2_tmp,theta2,mom2)); 
 
    } 
    // if (mom2>400) cout << " phi2_tmp " << phi2_tmp  <<" theta "<< theta2*59.1 <<"  mom2  "<<mom2<< " eff2  " << fEff2 << endl; 
 
    if (fEff1>0.00 && fEff2>0.00) return fEff1*fEff2; 
    else { // * efficiency hole * // ?? 
	if (mom2<2000 &&mom1<2000){ 
//	    cout<<" phi2_tmp "<< r2d*phi2_tmp<<"  theta2  "<< r2d*theta2<<"   mom2  "<<mom2<<endl; 
//	    cout<<" phi1_tmp "<< r2d*phi1_tmp<<"  theta1  "<< r2d*theta1<<"   mom1  "<<mom1<<endl; 
	} 
	return 0.; 
    } 
} 

Float_t getEfficiencyFactorP (float phi1, float theta1, float mom1)
{ 
  static Int_t count1=0;
  static TH3F *p3Deff;
  Float_t pi=3.14159;

  if (count1==0) { 

   // Init efficiency matrices pp1.25 GeV Tetyana & Jochen & Witek convention

  TFile *pFile = new TFile("matrices/proton_eff_Apr06_realgeom.root","READ"); 
  if (pFile) 
   { 
       pFile->cd(); 
       p3Deff = (TH3F*) pFile->Get("ProtonEff");  
   } 
   else Error("Init","pointer to eff matrix file is NULL"); 
  
 
  count1++;

  }
 
  mom1=mom1*1000;
 
  if(mom1>1990) mom1=1990;
 
    Float_t fEff1 = 1.; 
    Float_t r2d = TMath::RadToDeg(); 


    Float_t phi1_tmp = phi1; 
   
 
	    fEff1 = p3Deff-> 
	    GetBinContent(p3Deff->FindBin(phi1_tmp,theta1,mom1)); 
//	    *eff_p=fEff1;
// cout << " phi1_tmp " << phi1_tmp*r2d  <<" theta "<< theta1*r2d <<"  mom1  "<<mom1<< " eff1  " << fEff1 << endl; 
  

    if (fEff1>0.00) return fEff1; 
    else { // * efficiency hole * // ?? 
	if (mom1<2000){ 
//	    cout<<" phi2_tmp "<< r2d*phi2_tmp<<"  theta2  "<< r2d*theta2<<"   mom2  "<<mom2<<endl; 
//	    cout<<" phi1_tmp "<< r2d*phi1_tmp<<"  theta1  "<< r2d*theta1<<"   mom1  "<<mom1<<endl; 
	} 
	return 0.; 
    } 
} 
  

//void analyse(const Char_t *inFile,  Float_t Ebeam , const Char_t *outFile, Int_t nev=0) {
void anapp_dev() {

    makeDistributionManager()->Unpack("../test/pluto_ee_filter_apr06.root");
    makeDistributionManager()->Startup("_filter_remove_particles=1");
    
    Char_t *inFile = "test_filter_pi0_to_p_D+.p_pi0.g_dilepton.e+_e-";
    Float_t Ebeam=1.25;
    const Char_t *outFile="test_pp_plott_filter_old";
    Int_t nev=0;
    

// input file

 TChain data("data");
    
   TString outfile(outFile);
   outfile+="_hist.root";
   cout<<outfile.Data()<<endl;
   TString infile(inFile);
   infile+=".root";
   cout<<infile.Data()<<endl;
   TFile file(infile, "READ");
   if (!file.IsZombie()) {
     data.Add(infile);
   }
 
TClonesArray *evt=new TClonesArray("PParticle",10); // create TClonesArray

 data.SetBranchAddress("Particles",&evt);           // and connect it to chain3

Float_t r2d =180.0/3.14159;
const Float_t bPar=0.0;
const Float_t PHISEC=30.0/r2d;


 
  TH1F *invmasdil =new TH1F("invmasdil","Invarian mass of  e+ and e- ",100,0.,1.);
  TH1F *inv_del =new TH1F("inv_del","Invarian mass of delta ",100,1.,2.);
  TH2F *invmasdildelta =new TH2F("invmasdildelta","Invariant mass of dil/delta ",100,0.,1.,100,1.,2.);
  TH2F *p1tp2t =new TH2F("p1tp2t","p1 vs p2",90,0.,90.0,90,0.,90.);
  TH1F *opangle =new TH1F("opangle","op. angle ",100,0.,360.);
  TH1F *misspi0 =new TH1F("misspi0","pi0 missing mass",600,-0.1,1.0); 
  TH1F *invmasdilac =new TH1F("invmasdilac","Invarian mass of  e+ and e- accepted",100,0.,1.);
  TH1F *invmasdilac9deg =new TH1F("invmasdilac9deg","Invarian mass of  e+ and e- accepted, 9deg opang",100,0.,1.);
  TH1F *invmasdilac9deg_mev =new TH1F("invmasdilac9deg_mev","Invariant mass of  e+ and e- accepted, 9deg opang mev",100,0.,1000.);
  TH1F *invmaspi0ac9deg =new TH1F("invmaspi0ac9deg","Invarian mass of e+,e- pi0,  accepted, 9deg opang",100,0.,1.);
  TH1F *invmasdeltac9deg =new TH1F("invmasdeltac9deg","Invarian mass of  e+, e- delta accepted, 9deg opang",100,0.,1.);
 TH1F *invmasetac9deg =new TH1F("invmasetac9deg","Invarian mass of  e+, e- eta accepted, 9deg opang",100,0.,1.);
 TH1F *invmasomegac9deg =new TH1F("invmasomegac9deg","Invarian mass of  e+, e- omega accepted, 9deg opang",100,0.,1.); 
 TH1F *invmasomegdalac9deg =new TH1F("invmasomegdalac9deg","Invarian mass of  e+, e- omega dalitz accepted, 9deg opang",100,0.,1.); 
 TH1F *invmasrhoac9deg =new TH1F("invmasrhoac9deg","Invarian mass of  e+, e- rho accepted, 9deg opang",100,0.,1.);
  TH1F *invmasdilacef =new TH1F("invmasdilacef","Invariant mass of  e+ and e- accepted, rec. eff",100,0.,1.);
  TH1F *ptac9deg =new TH1F("ptac9deg","Pt of  e+ and e- accepted, 9deg opang",50,0.,1.);
  TH1F *rapac9deg =new TH1F("rapac9deg","Rap of  e+ and e- accepted, 9deg opang",75,-1.,2.);
  TH1F *ptac9deg_hm =new TH1F("ptac9deg_hm","Pt of  e+ and e- accepted, 9deg opang",50,0.,1.);
  TH1F *rapac9deg_hm =new TH1F("rapac9deg_hm","Rap of  e+ and e- accepted, 9deg opang",75,-1.,2.);
  TH1F *invmasdilac9degef =new TH1F("invmasdilac9degef","Invarian mass of  e+ and e- accepted, 9deg opang, rec. eff",100,0.,1.);
 TH1F *invmasdilac9deg_ex =new TH1F("invmasdilac9deg_ex","Invarian mass of  e+ and e- accepted, 9deg opang, rec. ex",100,0.,1.);
  TH1F *invmasom =new TH1F("invmasom","Invarian mass of  omega ",100,0.,1.);
  TH1F *misspp2 =new TH1F("misspp2","pp missing mass",100,0.,1.); 
  TH1F *misspp2ac =new TH1F("misspp2ac","pp missing mass accepted",100,0.,1.); 
  TH2F *dalitzpom =new TH2F("dalitzpom","Dalitz plot proton-omega ",100,2.,5.0,100,2.,5.0); 
  TH2F *dalitzpomac =new TH2F("dalitzpomac","Dalitz plot proton-omega accepted",100,2.,5.0,100,2.,5.0);  
  TH2F *dalitzpomac1p =new TH2F("dalitzpomac1p","Dalitz plot proton-omega 1 proton accepted",100,2.,5.0,100,2.,5.0);
  TH2F *p1momtheta =new TH2F("p1momtheta","proton1  momentum vs theta",100,0.,4.0,45,0.,90.);
  TH2F *p2momtheta =new TH2F("p2momtheta","proton2  momentum vs theta",100,0.,4.0,45,0.,90.);
  TH2F *emomtheta =new TH2F("emomtheta","electron  momentum vs theta",100,0.,1000,45,0.,90.);
  TH2F *pmomtheta =new TH2F("pmomtheta","positron  momentum vs theta",100,0.,1000,45,0.,90.);
  TH2F *emomthetaef =new TH2F("emomthetaef","electron  momentum vs theta, rec.eff",50,0.,1.0,45,0.,90.);
  TH2F *pmomthetaef =new TH2F("pmomthetaef","positron  momentum vs theta, rec.eff",50,0.,1.0,45,0.,90.);
  TH1F *delta_mass1 =new TH1F("delta_mass1","pe+e- inv. mass-proton accepted",100,0.8,1.6);
  TH1F *delta_mass1c =new TH1F("delta_mass1c","pe+e- inv. mass-proton accepted correct",100,.8,1.6);
  TH1F *delta_costheta1= new TH1F("delta_costheta1","delta cos theta-proton accepted",250,-1,1);
  TH1F *delta_costheta2= new TH1F("delta_costheta2","delta cos theta-proton missing",250,-1,1);
  TH1F *omega_costheta= new TH1F("omega_costheta","omega cos theta",100,-1,1);
  TH2F *ev4_st_costheta= new TH2F("ev4_st_costheta","e- cos theta in gamma* frame",20,-1.,1.,50,0.,1.0);
  TH1F *gamma_st_costheta= new TH1F("gamma_st_costheta","gamma* cos theta vs inv mass",100,-1.,1.);
  TH1F *ev4_st_costheta_ex= new TH1F("ev4_st_costheta_ex","e- cos theta in gamma* frame",20,-1,1);
  TH1F *gamma_st_costheta_ex= new TH1F("gamma_st_costheta_ex","gamma* cos theta vs inv mass",100,-1.,1.);
  TH1F *cmcosTheta_del= new TH1F("cmcosTheta_del","delta cos theta",500,-1,1);
  TH1F *pid=new TH1F("pid","particle PID",100,0,100); 


   Float_t xbins[] =
{0.000,0.005,0.010,0.015,0.020,0.025,0.030,0.035,0.040,0.045,0.050,0.055,0.070,0.085,0.100,0.120,0.140,0.180,0.230,0.280,0.350,0.420,0.500};
   Int_t nbins = sizeof(xbins)/sizeof(Float_t);

   TH1F *sig_all_var = new TH1F("sig_all_var","sig_all_var",nbins-1,xbins);

PParticle *beam   = new PParticle(14,Ebeam);      // beam  14-proton
PParticle *target = new PParticle(14);          // target
PParticle *cm      = new PParticle(*beam + *target);


 PParticle *par[30],*pPart;
  
TLorentzVector p1v4,p2v4,p3v4,piv4,ev4,pv4,dilep,omega,p1om,p2om,elpi,ppi,misspp,cm4,delv4,gamma_st,ev4_st,pv4_st,gamma_st_pi0;	  

Int_t nentries = (Int_t)(data.GetEntries()); // get number of entries
if (nev<=0) nev = nentries;

printf("+++ nentries : %d +++ nev : %d \n", nentries, nev);

TStopwatch timer; // timer
timer.Start();

for (Int_t iev=0; iev<nev; iev++) {          // enter event loop
   cout << iev << endl;
    Float_t p1p,p1t,p2t,p2p,et,pt,ep,pp,omega_cos,gamma_st_cos,ev4_st_cos,pv4_st_cos;

     data.GetEntry(iev);

  Int_t npart = (Int_t)(evt->GetEntries());

  Int_t ip=0;
  Int_t pd,IDal=0;
  Int_t part[7];
  PParticle *fpart[100];
  double weight=0.0;
  Float_t angle;
  for (Int_t i=0;i<7;i++) {part[i]=0;}

  for (Int_t j=0;j<npart;j++) 
   {     
     pPart =(PParticle*)(evt->At(j));
     fpart[j]=pPart;
     pd = pPart->ID();
     weight=pPart->W();
     pid->Fill(pd,weight);
     if (pd==14) 
      {
	  par[ip]=(PParticle*)(evt->At(j));     //protons
          part[ip]=1;
          ip=ip+1;
      }
     else if(pd==7 || pd==1) {
	 par[4]=(PParticle*)(evt->At(j));  // neutral particle:pi0 or gamma for Dalitz decays: pi0(gamma)e+e-
         part[4]=1;
       IDal=1;
     }
     else if(pd==2) {part[2]=1; par[2]=(PParticle*)(evt->At(j));} //e+
     else if(pd==3) {part[3]=1; par[3]=(PParticle*)(evt->At(j));} //e-
     else if(pd==8) {part[5]=1; par[5]=(PParticle*)(evt->At(j));} //pi+
     else if(pd==9) {part[6]=1; par[6]=(PParticle*)(evt->At(j));} //pi-
   } 

/*
     par[0]=(PParticle*)(evt->At(0)); // get particles from Tree proton
     par[1]=(PParticle*)(evt->At(1));     //proton
     par[2]=(PParticle*)(evt->At(2));     //pi0 or gamma
     par[3]=(PParticle*)(evt->At(3));     //e+
     par[4]=(PParticle*)(evt->At(4));     //e-
*/
     for (Int_t i=0; i<6; i++){
     if(part[i]==1) par[i]->SetActive();
     }                       

   Int_t test=0; 

/* Acceptance of HADES */

#if 1
   for (Int_t i=0;i<6; i++){ 
       if(part[i]==1){
       test=selectAc(par[i]);
       if(test!=1) par[i]->SetInActive();
       }
   }
#endif

//cout << "npart:" << npart << endl;
//makeDistributionManager()->GetLoopFilter()->Modify(fpart,NULL,&npart,npart);

/* Resolution */

    for (Int_t i=0;i<6; i++){ 
       if(part[i]==1){    
       double smear_factor= GetResolution(par[i],3);
       //smear_factor= 0;
       smear(par[i],smear_factor);   // finite resolution for 4 MDC only
       }
    }
/* hadrons example 2 proton analysis in HADES acceptance*/ 

       if(part[0]==1) p1v4=par[0]->Vect4();                //proton 1
       if(part[1]==1) p2v4=par[1]->Vect4();                //proton 2
       cm4=cm->Vect4();
//       cout<<p1v4.P()<<"next"<<p2v4.P()<<endl;
      if(part[0]==1 &&  part[1]==1){
        delv4=cm4-p1v4;            //second proton from resonance
       if(par[0]->IsActive() && par[1]->IsActive()){          //HADES acceptance 
//	{
//       p3v4=par[5]->Vect4();                
       angle = r2d*(p1v4.Phi()-p2v4.Phi());  // cooplanarity angle cut
       opangle->Fill(angle,weight);
       p1tp2t->Fill(p1v4.Theta()*r2d,p2v4.Theta()*r2d,weight);
//       delv4=p3v4+p1v4;      //p pi0
//       delv4=p1v4;             // elastic
       misspi0->Fill((cm4-p1v4-p2v4).Mag2(),weight); 
       TLorentzVector cm_del=TLorentzVector(delv4);
       cm_del.Boost(-cm4.BoostVector());
       cmcosTheta_del->Fill(cm_del.Pz()/cm_del.P(),weight); 
       inv_del->Fill(delv4.Mag(),weight);
       }
       }
/* Dalitz decay reconstruction */
        gamma_st_cos=0;
        ev4_st_cos=0; 
       if(part[2]==1 && part[3]==1) {
       if(IDal) piv4=par[4]->Vect4();                //pi0 or gamma from Dalitz
       pv4=par[2]->Vect4();                 //e+
       ev4=par[3]->Vect4();                 //e-
       dilep=ev4+pv4;
       invmasdil->Fill(dilep.Mag(),weight);
       if(IDal) omega=ev4+pv4+piv4;            //pi0 or omega or eta 4 vector
       else { omega=ev4+pv4;
       invmasdildelta->Fill(omega.Mag(),delv4.Mag(),weight); // Delta dalitz decay
       }

// only if omega vector is delta rest frame ! 

     if (PUtils::sampleFlat()<0.5) omega = delv4; 
     else omega=p2v4+ev4+pv4;

/* cm angular distributions */
       TLorentzVector omega_cm=TLorentzVector(omega);
       omega_cm.Boost(-cm4.BoostVector());
       omega_cos=omega_cm.Pz()/omega_cm.P(); 
       omega_costheta->Fill(omega_cos,weight);
       
       gamma_st=ev4+pv4;
       gamma_st.Boost(-cm4.BoostVector());   // gamma* 4vector in CM
       gamma_st_cos=gamma_st.Pz()/gamma_st.P();
       gamma_st_costheta->Fill(gamma_st_cos,weight);
       
       gamma_st_pi0=ev4+pv4;                                 // e- 4vector in gamma* rest frame
       gamma_st_pi0.Boost(-omega.BoostVector());              //pi0 rest frame
       ev4_st = ev4;
       pv4_st = pv4;

       ev4_st.Boost(-omega.BoostVector());
       ev4_st.Boost(-gamma_st_pi0.BoostVector());
       ev4_st_cos = cos(ev4_st.Angle(gamma_st_pi0.Vect()));
       ev4_st_costheta->Fill(ev4_st_cos,dilep.Mag(),weight);

       pv4_st.Boost(-omega.BoostVector());
       pv4_st.Boost(-gamma_st_pi0.BoostVector());
       pv4_st_cos = cos(pv4_st.Angle(gamma_st_pi0.Vect()));

       misspp= cm4-p1v4-p2v4;
       p1om=omega+p1v4;
       p2om=omega+p2v4;

       p1p=p1v4.P();
       p2p=p2v4.P();
       ep=ev4.P();
       pp=pv4.P();
       et = ev4.Theta()*r2d;
       pt = pv4.Theta()*r2d;
//      cout<<"ID p1,p2,pi,el,pos       
       p1t = p1v4.Theta()*r2d;
       p2t = p2v4.Theta()*r2d;

       if(IDal) invmasom->Fill(omega.Mag(),weight);
       if(IDal) dalitzpom->Fill(p1om.Mag2(),p2om.Mag2(),weight);
       misspp2->Fill(misspp.Mag2(),weight);
       Float_t angle = r2d*ev4.Angle(pv4.Vect());  // opening angle cut      


// check if particles are in the HADES acceptance

       if(par[2]->IsActive() &&  par[3]->IsActive() ) 
       {                                                         // e+ e- Acceptance !
	   Int_t procesId1=par[2]->GetParentId();
           Int_t procesId2=par[3]->GetParentId();
//           cout<<"yes,yes"<<procesId1<<procesId2<<endl;
             // reconstruction efficiency electrons
  float eff_s=getEfficiencyFactor(par[2]->Phi(),par[2]->Theta(),par[2]->P(),par[2]->ID(),par[3]->Phi(),par[3]->Theta(),par[3]->P(), par[3]->ID()); 
//               eff_s=1.0;                                      // no efficiency
               if(eff_s>0.0) {                                   // efficiency >0.05
               invmasdilac->Fill(dilep.Mag(),weight);
               emomtheta->Fill(ep*1000.,et,weight);
               pmomtheta->Fill(pp*1000.,pt,weight);
               if(angle>9) {
                 invmasdilac9deg->Fill(dilep.Mag(),weight);
                 invmasdilac9deg_mev->Fill(dilep.Mag()*1000,weight);
                 ptac9deg->Fill(dilep.Pt(),weight);
                 rapac9deg->Fill(dilep.Rapidity(),weight);
                 if(dilep.Mag()>0.14) { ptac9deg_hm->Fill(dilep.Pt(),weight);
                    rapac9deg_hm->Fill(dilep.Rapidity(),weight);}
                 if (procesId1==7051 && procesId2==7051) invmaspi0ac9deg->Fill(dilep.Mag(),weight);
                 if (procesId1==36051 && procesId2==36051) invmasdeltac9deg->Fill(dilep.Mag(),weight);
                 if (procesId1==17051 && procesId2==17051) invmasetac9deg->Fill(dilep.Mag(),weight);
                 if (procesId1==52 && procesId2==52) invmasomegac9deg->Fill(dilep.Mag(),weight);
                 if (procesId1==41 && procesId2==41) invmasrhoac9deg->Fill(dilep.Mag(),weight);
                 if (procesId1==52051 && procesId2==52051) invmasomegdalac9deg->Fill(dilep.Mag(),weight);
               	       }
// rec efficiency included in weights;
// dielectrons
 	        eff_s=eff_s*weight;
                invmasdilacef->Fill(dilep.Mag(),eff_s);
                emomthetaef->Fill(ep,et,eff_s);
                pmomthetaef->Fill(pp,pt,eff_s);
                if(angle>9) invmasdilac9degef->Fill(dilep.Mag(),eff_s);
		if((par[0]->IsActive() || par[1]->IsActive()) && IDal==0) dalitzpomac1p->Fill(p1om.Mag2(),p2om.Mag2(),eff_s);
// only for exclusive pe+e- reconstruction	

		if((par[0]->IsActive()  ||  par[1]->IsActive()) && angle>9)  // pp acceptance 
	      {
// proton
                  float eff_p=1.0;
                  if(par[0]->IsActive()) eff_p=getEfficiencyFactorP(par[0]->Phi(),par[0]->Theta(),par[0]->P());
		  //     eff_s=eff_s*eff_p;
                  if(par[1]->IsActive()) eff_p=getEfficiencyFactorP(par[1]->Phi(),par[1]->Theta(),par[1]->P());
                  eff_s=eff_s*eff_p;
                  TLorentzVector delta_v41,delta_v42,delta_cm1,delta_cm2,delta_corr1; 
                  if (dilep.Mag()>0.14) {  // select mass region && fiducial volume
                      delta_v41=cm4-p1v4;delta_corr1=cm4-p1v4;delta_v42=cm4-p2v4;
		   delta_mass1->Fill(delta_v41.Mag(),eff_s);
                   delta_mass1c->Fill(delta_corr1.Mag(),eff_s);
                   delta_mass1->Fill(delta_v42.Mag(),eff_s);
                   delta_cm1=delta_v41;
                   delta_cm2=delta_v42;
                   delta_cm1.Boost(-cm4.BoostVector());
                   delta_cm2.Boost(-cm4.BoostVector());
                   delta_costheta1->Fill(delta_cm1.Pz()/delta_cm1.P(),eff_s);
                   delta_costheta2->Fill(delta_cm2.Pz()/delta_cm2.P(),eff_s);

                     if (ep>0.1 && pp>0.1 && et>26 && et<82 && pt>26 && pt<82) {    
                     gamma_st_costheta_ex->Fill(gamma_st_cos,eff_s);
                     ev4_st_costheta_ex->Fill(ev4_st_cos,eff_s);
                     ev4_st_costheta_ex->Fill(pv4_st_cos,eff_s);
                     }
                  }
//
              dalitzpomac->Fill(p1om.Mag2(),p2om.Mag2(),eff_s);
              misspp2ac->Fill(misspp.Mag2(),eff_s);
              p1momtheta->Fill(p1p,p1t,eff_s);
              p2momtheta->Fill(p2p,p2t,eff_s);
              invmasdilac9deg_ex->Fill(dilep.Mag(),eff_s);
              sig_all_var->Fill(dilep.Mag(),eff_s);
	      //            cout<<"exlusive ev4_st_cos="<<ev4_st_cos<<endl;
	      } //pp acceptance 
	   }     // efficiency >0.05	   
       }         //e+,e- acceptance 
      }          // e+e- in the event
}	         //event loop  

//  TFile myfile("ana_test.root","RECREATE");

  TFile myfile(outfile,"RECREATE");
  pid->Write();
  invmasdil->Write();
  inv_del->Write();
  invmasdildelta->Write();
  p1tp2t->Write(); 
  opangle->Write();
  misspi0->Write();
  invmasdilac->Write();
  invmasdilac9deg->Write();
  invmasdilac9deg_mev->Write();
  invmaspi0ac9deg->Write();
  invmasdeltac9deg->Write();
  invmasomegac9deg->Write();
  invmasomegdalac9deg->Write();
  invmasrhoac9deg->Write();
  invmasetac9deg->Write();
  ptac9deg->Write();
  rapac9deg->Write();
  ptac9deg_hm->Write();
  rapac9deg_hm->Write();
  invmasdilac9degef->Write();
  invmasdilacef->Write();
  invmasdilac9deg_ex->Write();
  sig_all_var->Write();
  invmasom->Write();
  dalitzpom->Write();
  dalitzpomac->Write();
  dalitzpomac1p->Write();
  p1momtheta->Write();
  p2momtheta->Write();
  emomtheta->Write();
  pmomtheta->Write();
  emomthetaef->Write();
  pmomthetaef->Write();
  misspp2->Write();
  misspp2ac->Write();
  omega_costheta->Write();
  gamma_st_costheta->Write();
  ev4_st_costheta->Write(); 
  gamma_st_costheta_ex->Write();
  ev4_st_costheta_ex->Write(); 
  cmcosTheta_del->Write(); 
  delta_mass1->Write();
  delta_mass1c->Write();
  delta_costheta1->Write();
  delta_costheta2->Write();
  myfile.Close();
}
