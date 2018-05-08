#include <cstdlib>
#include <string>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <iomanip> 
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <math.h>
#include <fstream>
#include "iostream"
#include "math.h"
#include "stdio.h" 
#include "stdlib.h"
#include "time.h"  
#include "string"  
#include "unistd.h"
#include "string.h"
#include "fstream"
//#include "Math.h"
#include "Math/Math.h"
#include "TROOT.h"  
#include "TTree.h"
#include <TFile.h>  
#include "TChain.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TArc.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"

//#include "header.h"

using namespace std;

int piony()

{
	    
      long double T11[51],T12[51],T13[51],T14[51],T16[51];
      long double T111[51],T112[51],T113[51],T114[51],T116[51];
      long double T122[51],T123[51],T124[51],T126[51];
      long double T133[51],T134[51],T136[51];
      long double T144[51],T146[51];
      long double T166[51],along[51];

      long double T21[51],T22[51],T23[51],T24[51],T26[51];      
      long double T211[51],T212[51],T213[51],T214[51],T216[51];
      long double T222[51],T223[51],T224[51],T226[51];
      long double T233[51],T234[51],T236[51];
      long double T244[51],T246[51];
      long double T266[51];      

      long double T31[51],T32[51],T33[51],T34[51],T36[51];
      long double T311[51],T312[51],T313[51],T314[51],T316[51];
      long double T322[51],T323[51],T324[51],T326[51];
      long double T333[51],T334[51],T336[51];
      long double T344[51],T346[51];
      long double T366[51];

      long double T41[51],T42[51],T43[51],T44[51],T46[51]; 
      long double T411[51],T412[51],T413[51],T414[51],T416[51];
      long double T422[51],T423[51],T424[51],T426[51];
      long double T433[51],T434[51],T436[51];
      long double T444[51],T446[51];
      long double T466[51];
           
      Int_t ipass[51];
      Int_t kpass[51],rpass[51];
     char tm[11];
      char tc[5];
      char tl[31];
      char txt14[15];
      char txt22[23];
      char txt2[3];
      char txt5[6];
      char txt6[7];
      
  long double rapport = 0;
  long double offset = 0;
  long double rx = 0;
  Int_t itotpaass = 0;
  Int_t itotq9 = 0;
  Int_t idet1 = 0;
  Int_t idet2 = 0;
  long double dist1 = 0;
  long double dist2 = 0;
  long double ouvrescaleh = 0;
  long double ouvrescalev = 0;
  long double r1 = 0;
  Int_t ii1 = 0;
  Int_t ii2 = 0;
  long double r2 = 0;
  long double r3 = 0;
  Int_t ii3 = 0;
  Int_t ii4 = 0;
  long double r4 = 0;
  long double rh = 0;
  long double rv = 0;
  long double rres = 0;
  long double rdm = 0;
  long double shiftx = 0;
  long double shifty = 0;
  long double dx = 0;
  long double dxp = 0;
  long double dy = 0;
  long double dyp = 0;
  long double dpspmax = 0;
  long double dpsppas = 0;
  Int_t ntir = 0;
  long double ppion = 0;
  long double esi = 0;
  long double resoldet = 0;
  long double betapion = 0;
  long double diffmul = 0;
  long double resoldmx = 0;
  long double resoldmy = 0;
  Int_t idineg = 0;
  long double xmax = 0;
  long double ymax = 0;
  long double xpmax = 0;
  long double ypmax = 0;
  long double dpmax = 0;
  //Int_t k = 0;
  long double relem = 0;
  long double z1 = 0;
  long double z2 = 0;
  long double z3 = 0;
  long double z4 = 0;
  long double z5 = 0;
  long double z6 = 0;
  long double z7 = 0;
  long double z8 = 0;
  long double z9 = 0;
  long double z10 = 0;
  long double z11 = 0;
  long double z12 = 0;
  long double z13 = 0;
  long double z14 = 0;
  long double z15 = 0;
  Int_t itir = 1;
  long double yy = 0;
  long double xi = 0;
  long double xpi = 0;
  long double yi = 0;
  long double ypi = 0;
  long double dp = 0;
  long double dd = 0;
  long double ndp = 0;
  long double uns = 0;
  Int_t kii = 0;
  long double seuil = 0;
  Int_t iok = 0;
  long double xf = 0;
  long double yf = 0;
  long double xpf = 0;
  long double ypf = 0;
  long double xf1 = 0;
  long double yf1 = 0;
  long double xpf1 = 0;
  long double ypf1 = 0;
  long double xf2 = 0;
  long double xpf2 = 0;
  long double yf2 = 0;
  long double ypf2 = 0;
  long double xdet1 = 0;
  long double ydet1 = 0;
  long double xdet2 = 0;
  long double ydet2 = 0;
  Int_t kk = 0;
  long double deltax = 0;
  long double deltay = 0;
  long double xhad = 0;
  long double yhad = 0;
  long double radius = 0;
  Int_t ijk = 0;
  Int_t itotpass = 0;
  long double xdet1sav = 0;
  long double xca = 0;
  Int_t nca = 0;
  long double xdet2sav = 0;
  long double ydet1sav = 0;
  long double ydet2sav = 0;
  long double yca = 0;
  long double alfa = 0;
  long double beta = 0;
  long double gamma = 0;
  long double discri = 0;
  long double galfa = 0;
  long double gbeta = 0;
  long double ggamma = 0;
  long double gdelta = 0;
  long double aa1 = 0;
  long double aa2 = 0;
  long double aa3 = 0;
  long double aa4 = 0;
  long double sol = 0;
  Int_t iter = 1;
  long double yyy1 = 0;
  long double yyp1 = 0;
  long double yyy2 = 0;
  long double yyp2 = 0;
  long double yyy3 = 0;
  long double yyp3 = 0;
  long double yyy4 = 0;
  long double yyp4 = 0;
  long double yyy5 = 0;
  long double yyp5 = 0;
  long double teta = 0;
  long double bb = 0;
  long double cc = 0;
  long double aap = 0;
  long double bbp = 0;
  long double ccp = 0;
  long double deter = 0;
  long double yci = 0;
  long double phi = 0;
  long double poids = 0;
  long double xfth = 0;
  long double yfth = 0;
  long double xfcal = 0;
  long double yfcal = 0;
  long double rr = 0;
  Int_t ipi = 0;
  Int_t ips2 = 0;
  Int_t idif = 0;
  long double xabs = 0;
  long double yabs = 0;
  Int_t ki = 0;
  long double reste = 0;
  long double ri = 0;
  Int_t kij = 0;
  long double alossfraction = 0;
  
 // string txt2 = 0;
  long double aa = 0;
  long double sig11 = 0;
//  string txt5 = 0;
  long double sig22 = 0;
  long double r12 = 0;
  long double sig33 = 0;
  long double r13 = 0;
  long double r23 = 0;
  long double sig44 = 0;
  long double r14 = 0;
  long double r24 = 0;
  long double r34 = 0;
  long double sig55 = 0;
  long double r15 = 0;
  long double r25 = 0;
  long double r35 = 0;
  long double r45 = 0;
  long double sig66 = 0;
  long double r16 = 0;
  long double r26 = 0;
  long double r36 = 0;
  long double r46 = 0;
  long double r56 = 0;
 // string txt14 = 0;
  Int_t ii = 0;
  Int_t ind1 = 0;
  Int_t ind23 = 0;
  double rreste = 0;
  Int_t retape = 0;
      
     Int_t  ouvh[51] = {0,6.,6.,6.,6.,7.,7.,7.,7.,7.,6.,
      6.,6.,6.,6.,6.,5.,5.,6.,6.,6., 
      6.,6.,6.,6.,6.,6.,6.,6.,6.,6.,
      6.,6.,0.6,5.,5.,5.,5.,5.,5.,5.,
      5.,5.,5.,5.,5.,5.,5.,5.,5.,5.  };
     Int_t ouvv[51] =  {0,6.,6.,6.,6.,3.5,3.5,3.5,3.5,3.5,3.5,
      3.5,6.,6.,6.,6.,6.,4.,4.,6.,6.,
      6.,3.5,3.5,6.,6.,6.,6.,6.,6.,6.,
      6.,6.,0.6,5.,5.,5.,5.,5.,5.,5.,
      5.,5.,5.,5.,5.,5.,5.,5.,5.,5.};   
     Int_t kind[11]= {0, 2, 4, 6, 10, 13, 17, 23, 26, 29, 33 };
    // string str;
   // std::cin.fixed;
   // cin.fixed;
   // cout.fixed;
    std::ifstream plik1;
    plik1.open("pibeam_set6.data");
   // getline(plik1,str);
   // std::cout<<"dupa"<<str<<endl;
    //ifstream plik1("pibeam_set3.data");
    std::ofstream plik2;
    plik2.open("pibeam_set3.out");
    
  //  cout.precision(10);
  //  cin.precision(10);
    
    //ios::precision(); 
   // ifstream.precision(10);
  // ofstream.precision (std::numeric_limits<double>::digits10 + 1);
    Int_t  nelem=34;
    long double renorgauss=sqrt(3./5.);
    Int_t x0si=22.;
    long double  rhosi=2.33;
    long double amasspi=0.13958;
    //Int_t dupa;
    //char asdf[3];
    
        for(Int_t iii=0; iii<51; iii++) 
    {
    ipass[iii] = 0;
    kpass[iii] = 0;
    rpass[iii] = 0.;
  }
  
  //--TTree--- branch variables----------------------
  Int_t sector =0;
  double rT166 = 0,rT266 = 0,rT366 = 0,rT466 = 0,raa = 0,ralong = 0;
  
   TFile *output = new TFile("piony.root", "RECREATE");
     
     TTree *tree = new TTree("Pion","Pionbeam");
     tree->Branch("sectors", &sector);
     tree->Branch("T166", &rT166);
     tree->Branch("T266", &rT266);
     tree->Branch("T366", &rT366);
     tree->Branch("T466", &rT466);
     tree->Branch("aa", &raa);
     tree->Branch("along", &ralong);
     
     sector = 1;
    string str;
    for(Int_t i=1; i<nelem; i++)
    {
		//i= i + 1;
	str.clear();	
	//getline(plik1,str);
		//stringstream plik1(str);
		//std::cout<<ouvh[i]<<endl;
    //C
    //C      plik1.read (5,*) tm,typ,txt6,al,txt5
    //C  100 format (A10)
    //plik1.read(plik1), along[i], txt2, aa, sig11, txt5;
   
    
   plik1>>along[i]>>txt2>>aa>>sig11>>txt5;
   plik1>>aa>>sig22>>txt5>>r12;
   plik1>>aa>>sig33>>txt5>>r13>>r23;
   plik1>>aa>>sig44>>txt5>>r14>>r24>>r34;
   plik1>>aa>>sig55>>txt5>>r15>>r25>>r35>>r45;
   plik1>>aa>>sig66>>txt5>>r16>>r26>>r36>>r46>>r56;

  

 getline(plik1,str,'\n'); 
 getline(plik1,str,'\n');

    
    plik1>>T11[i]>>T12[i]>>T13[i]>>T14[i]>>aa>>T16[i];

    cout << T11[i] << endl;

    plik1>>T21[i]>>T22[i]>>T23[i]>>T24[i]>>aa>>T26[i];
    plik1>>T31[i]>>T32[i]>>T33[i]>>T34[i]>>aa>>T36[i]>>std::fixed;
   // std::cout<<T31[i]<<endl;
    plik1>>T41[i]>>T42[i]>>T43[i]>>T44[i]>>aa>>T46[i];
    plik1>>aa>>aa>>aa>>aa>>aa>>aa;
    plik1>>aa>>aa>>aa>>aa>>aa>>aa;
    getline(plik1,str);
    

    //C
    //C     lecture second ordre T1ij
    //C
   getline(plik1,str);
   plik1>>ind1>>ind23>>T111[i];
   plik1>>ind1>>ind23>>T112[i]>>ind1>>ind23>>T122[i];
   plik1>>ind1>>ind23>>T113[i]>>ind1>>ind23>>T123[i]>>ind1>>ind23>>T133[i];
   plik1>>ind1>>ind23>>T114[i]>>ind1>>ind23>>T124[i]>>ind1>>ind23>>T134[i]>>ind1>>ind23>>T144[i];
   plik1>>ind1>>ind23>>aa>>ind1>>ind23>>aa>>ind1>>ind23>>aa>>ind1>>ind23>>aa>>ind1>>ind23>>aa;
   plik1>>ind1>>ind23>>T116[i]>>ind1>>ind23>>T126[i]>>ind1>>ind23>>T136[i]>>ind1>>ind23>>T146[i]>>ind1>>ind23>>aa>>ind1>>ind23>>T166[i];
   getline(plik1,str);
   plik2<<".."<<i<<".."<<T166[i]<<endl;
   
    //C
    //C     lecture second ordre T2ij
    //C
    plik1>>ii>>ii>>T211[i];
    plik1>>ii>>ii>>T212[i]>>ii>>ii>>T222[i];
    plik1>>ii>>ii>>T213[i]>>ii>>ii>>T223[i]>>ii>>ii>>T233[i];
    plik1>>ii>>ii>>T214[i]>>ii>>ii>>T224[i]>>ii>>ii>>T234[i]>>ii>>ii>>T244[i];
    plik1>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa;
    plik1>>ii>>ii>>T216[i]>>ii>>ii>>T226[i]>>ii>>ii>>T236[i]>>ii>>ii>>T246[i]>>ii>>ii>>aa>>ii>>ii>>T266[i]; 
    getline(plik1,str);
    plik2<<".."<<i<<".."<<T266[i]<<endl;
    rT266 = T266[i];
    
    //C
    //C     lecture second ordre T3ij
    //C
   plik1>>ind1>>ind23>>T311[i];
   plik1>>ind1>>ind23>>T312[i]>>ind1>>ind23>>T322[i];
   plik1>>ind1>>ind23>>T313[i]>>ind1>>ind23>>T323[i]>>ind1>>ind23>>T333[i];
   plik1>>ind1>>ind23>>T314[i]>>ind1>>ind23>>T324[i]>>ind1>>ind23>>T334[i]>>ind1>>ind23>>T344[i];
   plik1>>ind1>>ind23>>aa>>ind1>>ind23>>aa>>ind1>>ind23>>aa>>ind1>>ind23>>aa>>ind1>>ind23>>aa;
   plik1>>ind1>>ind23>>T316[i]>>ind1>>ind23>>T326[i]>>ind1>>ind23>>T336[i]>>ind1>>ind23>>T346[i]>>ind1>>ind23>>aa>>ind1>>ind23>>T366[i];
   getline(plik1,str);
   plik2<<".."<<i<<".."<<T366[i]<<endl;
   rT366 = T366[i];
   
    //C
    //C     lecture second ordre T4ij
    //C
   
   plik1>>ii>>ii>>T411[i];
   plik1>>ii>>ii>>T412[i]>>ii>>ii>>T422[i];
   plik1>>ii>>ii>>T413[i]>>ii>>ii>>T423[i]>>ii>>ii>>T433[i];
   plik1>>ii>>ii>>T414[i]>>ii>>ii>>T424[i]>>ii>>ii>>T434[i]>>ii>>ii>>T444[i];
   plik1>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa;
   plik1>>ii>>ii>>T416[i]>>ii>>ii>>T426[i]>>ii>>ii>>T436[i]>>ii>>ii>>T446[i]>>ii>>ii>>aa>>ii>>ii>>T466[i];
   getline(plik1,str);
   plik2<<".."<<i<<".."<<T466[i]<<endl;
   rT466 = T466[i];
    //C
    //C     lecture second ordre T5ij
    //C
  
  plik1>>ii>>ii>>aa;
  plik1>>ii>>ii>>aa>>ii>>ii>>aa;
  plik1>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa;
  plik1>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa;
  plik1>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa;
  plik1>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa>>ii>>ii>>aa;
  getline(plik1,str);
  plik2<<".."<<i<<".."<<aa<<endl;
  raa = aa;

 // cout<<"nowa petla"<<endl;
 // getline(plik1,str);
 // std::cout<<"pusto"<<asdf<<"???????????????"<<endl;
  
    //C
    //C     fin de lecture des Tijk
    //C
    rapport = 0;
    //C  10 cm :offset pour digitaliser la position
    offset = 10;
    //C     sur les detecteurs
    if (T126[i] != 0) {
      rapport = -100 * T12[i] / T126[i];
    }
    if (rapport <  - 15) {
      rapport = -15;
    }
    if (rapport >   15) {
      rapport = 15;
    }
   // rx = fem::ffloat(i) + 0.01f;
   rx = float(i) + 0.01;
    //C     call hfill (40,rx,0.,rapport+20.)
    //C      write (6,*) i,rapport
   plik2<<"..along ="<<along[i]<<"..T126 ="<<T126[i]<<endl;
   ralong = along[i];
   getline(plik1,str,'\n');
  // getline(plik1,str);
  rT166 = T166[i];
   sector = i + 1;
 //if (sector > 0) 
 tree->Fill();
  }
  
  plik1.close();
  
  itotpaass = 0;
  itotq9 = 0;
  //C plan intermediaire
  idet1 = 17;
  //C      idet2=20  ! entree Q5
  //C      idet2=22    !  entree dipole 2
  //C      idet2=24      ! entre Q7
  //C  point intermediaire entre Q7 et Q8
  idet2 = 26;
  //C      idet2=27    ! entree Q8
  //C      idet2=28     ! sortie Q8
  //C      idet2=29  ! entree Q9
  //C      idet2=30   !  sortie Q9
  //C      idet2=31   !  1 metre avant cible HADES
  //C   calcul des T1jk a une distance dist2 du point standard pour detecteur 2
  //C distance a mettre en metre
  dist1 = 0;
  //C      write (6,*) t11[idet2],t12[idet2],t13[idet2],t14[idet2],t16[idet2]
 //std::setprecision(20);
  // is it really necessary? dist1=0.0 !
  /*
      T11[idet1]=T11[idet1]+0.1*dist1*T21[idet1];
      T12[idet1]=T12[idet1]+0.1*dist1*T22[idet1];
      T13[idet1]=T13[idet1]+0.1*dist1*T23[idet1];
      T14[idet1]=T14[idet1]+0.1*dist1*T24[idet1];
      T16[idet1]=T16[idet1]+0.1*dist1*T26[idet1];
      T111[idet1]=T111[idet1]+0.1*dist1*T211[idet1];
      T112[idet1]=T112[idet1]+0.1*dist1*T212[idet1];
      T113[idet1]=T113[idet1]+0.1*dist1*T213[idet1];
      T114[idet1]=T114[idet1]+0.1*dist1*T214[idet1];
      T116[idet1]=T116[idet1]+0.1*dist1*T216[idet1];
      T122[idet1]=T122[idet1]+0.1*dist1*T222[idet1];
      T123[idet1]=T123[idet1]+0.1*dist1*T223[idet1];
      T124[idet1]=T124[idet1]+0.1*dist1*T224[idet1];
      T126[idet1]=T126[idet1]+0.1*dist1*T226[idet1];
      T133[idet1]=T133[idet1]+0.1*dist1*T233[idet1];
      T134[idet1]=T134[idet1]+0.1*dist1*T234[idet1];
      T136[idet1]=T136[idet1]+0.1*dist1*T236[idet1];
      T144[idet1]=T144[idet1]+0.1*dist1*T244[idet1];
      T146[idet1]=T146[idet1]+0.1*dist1*T246[idet1];
      T166[idet1]=T166[idet1]+0.1*dist1*T266[idet1];
      T31[idet1]=T31[idet1]+0.1*dist1*T41[idet1];
      T32[idet1]=T32[idet1]+0.1*dist1*T42[idet1];
      T33[idet1]=T33[idet1]+0.1*dist1*T43[idet1];
      T34[idet1]=T34[idet1]+0.1*dist1*T44[idet1];
      T36[idet1]=T36[idet1]+0.1*dist1*T46[idet1];
      T311[idet1]=T311[idet1]+0.1*dist1*T411[idet1];
      T312[idet1]=T312[idet1]+0.1*dist1*T412[idet1];
      T313[idet1]=T313[idet1]+0.1*dist1*T413[idet1];
      T314[idet1]=T314[idet1]+0.1*dist1*T414[idet1];
      T316[idet1]=T316[idet1]+0.1*dist1*T416[idet1];
      T322[idet1]=T322[idet1]+0.1*dist1*T422[idet1];
      T323[idet1]=T323[idet1]+0.1*dist1*T423[idet1];
      T324[idet1]=T324[idet1]+0.1*dist1*T424[idet1];
      T326[idet1]=T326[idet1]+0.1*dist1*T426[idet1];
      T333[idet1]=T333[idet1]+0.1*dist1*T433[idet1];
      T334[idet1]=T334[idet1]+0.1*dist1*T434[idet1];
      T336[idet1]=T336[idet1]+0.1*dist1*T436[idet1];
      T344[idet1]=T344[idet1]+0.1*dist1*T444[idet1];
      T346[idet1]=T346[idet1]+0.1*dist1*T446[idet1];
      T366[idet1]=T366[idet1]+0.1*dist1*T466[idet1];
      */
      //----------------------------------------------
      // dist2=-1.6;  // ! distance a mettre en metre 
      // T11[idet2]=T11[idet2]+0.1*dist2*T21[idet2];
      // T12[idet2]=T12[idet2]+0.1*dist2*T22[idet2];
      // T13[idet2]=T13[idet2]+0.1*dist2*T23[idet2];
      // T14[idet2]=T14[idet2]+0.1*dist2*T24[idet2];
      // T16[idet2]=T16[idet2]+0.1*dist2*T26[idet2];
      // T111[idet2]=T111[idet2]+0.1*dist2*T211[idet2];
      // T112[idet2]=T112[idet2]+0.1*dist2*T212[idet2];
      // T113[idet2]=T113[idet2]+0.1*dist2*T213[idet2];
      // T114[idet2]=T114[idet2]+0.1*dist2*T214[idet2];
      // T116[idet2]=T116[idet2]+0.1*dist2*T216[idet2];
      // T122[idet2]=T122[idet2]+0.1*dist2*T222[idet2];
      // T123[idet2]=T123[idet2]+0.1*dist2*T223[idet2];
      // T124[idet2]=T124[idet2]+0.1*dist2*T224[idet2];
      // T126[idet2]=T126[idet2]+0.1*dist2*T226[idet2];
      // T133[idet2]=T133[idet2]+0.1*dist2*T233[idet2];
      // T134[idet2]=T134[idet2]+0.1*dist2*T234[idet2];
      // T136[idet2]=T136[idet2]+0.1*dist2*T236[idet2];
      // T144[idet2]=T144[idet2]+0.1*dist2*T244[idet2];
      // T146[idet2]=T146[idet2]+0.1*dist2*T246[idet2];
      // T166[idet2]=T166[idet2]+0.1*dist2*T266[idet2];
      // T31[idet2]=T31[idet2]+0.1*dist2*T41[idet2];
      // T32[idet2]=T32[idet2]+0.1*dist2*T42[idet2];
      // T33[idet2]=T33[idet2]+0.1*dist2*T43[idet2];
      // T34[idet2]=T34[idet2]+0.1*dist2*T44[idet2];
      // T36[idet2]=T36[idet2]+0.1*dist2*T46[idet2];
      // T311[idet2]=T311[idet2]+0.1*dist2*T411[idet2];
      // T312[idet2]=T312[idet2]+0.1*dist2*T412[idet2];
      // T313[idet2]=T313[idet2]+0.1*dist2*T413[idet2];
      // T314[idet2]=T314[idet2]+0.1*dist2*T414[idet2];
      // T316[idet2]=T316[idet2]+0.1*dist2*T416[idet2];
      // T322[idet2]=T322[idet2]+0.1*dist2*T422[idet2];
      // T323[idet2]=T323[idet2]+0.1*dist2*T423[idet2];
      // T324[idet2]=T324[idet2]+0.1*dist2*T424[idet2];
      // T326[idet2]=T326[idet2]+0.1*dist2*T426[idet2];
      // T333[idet2]=T333[idet2]+0.1*dist2*T433[idet2];
      // T334[idet2]=T334[idet2]+0.1*dist2*T434[idet2];
      // T336[idet2]=T336[idet2]+0.1*dist2*T436[idet2];
      // T344[idet2]=T344[idet2]+0.1*dist2*T444[idet2];
      // T346[idet2]=T346[idet2]+0.1*dist2*T446[idet2];
      // T366[idet2]=T366[idet2]+0.1*dist2*T466[idet2];
 
  
  plik2<<".."<<T11[idet2 - 1]<<".."<<T12[idet2 - 1]<<".."<<T13[idet2]<<".."<<T16[idet2 - 1]<<endl;
  plik2<<".."<<T11[idet2]<<".."<<T12[idet2]<<".."<<T13[idet2]<<".."<<T14[idet2]<<".."<<T16[idet2]<<endl;;
  plik2<<".."<<T31[idet2 - 1]<<".."<<T32[idet2 - 1]<<".."<<T33[idet2]<<".."<<T36[idet2 - 1]<<endl;
  plik2<<".."<<T31[idet2]<<".."<<T32[idet2]<<".."<<T33[idet2]<<".."<<T34[idet2]<<".."<<T36[idet2]<<endl;
  

  
   TFile *output1 = new TFile("piony1.root", "RECREATE");
   TTree *treex = new TTree("Pion1","Pionbeam1");
   treex->Branch("reste", &rreste);
   treex->Branch("etape", &retape);
   TH1F *random = new TH1F("random","random",3000,-1.5,1.5);
  //C
  //C zone utile detecteurs det1 et det2
  ouvh[idet1] = 5.0;
  //C zone utile detecteurs det1 et det2
  ouvv[idet1] = 5.0;
  //C zone utile detecteurs det1 et det2
  ouvh[idet2] = 5.0;
  //C zone utile detecteurs det1 et det2
  ouvv[idet2] = 5.0;
  //C
  ouvh[33] = 4.6;
  ouvv[33] = 4.6;
  //C
  ouvrescaleh = 1.;
  ouvrescalev = 1.;
  //float ide2,ii1,ii2,idet2,

  
  r1 = float(idet2) + 0.001;
  r1 = r1 / 10.;
  ii1 = r1;
  r1 = float(ii1) + 0.001;
  ii2 = idet1 - 10 * ii1;
  r2 = float(ii2) + 0.001;
  r3 = float(idet2) + 0.001;
  r3 = r3 / 10;
  ii3 = r3;
  r3 = float(ii3) + 0.001;
  ii4 = idet2 - 10 * ii3;
  r4 = float(ii4) + 0.001;
  
  //std::cout<<r4<<endl;
  
  //C
  //C     r1 a r4 quand le numero d'element ou sont mis les detecteurs
  //C     rh,rv,rres,rdm donne le facteur de reduction sur l'effet de
  //C     la taille h de l'objet
  //C     la taille v de l'objet
  //C     la resolution des detecteurs
  //C     la diffusion multiple
  //C     une valeur =0 indique que c'est pris en compte totalement
  //C     une valeur non nulle indique qu'on divise l'effet
  //C     correspondant par 10**rh, 10**rv, 10**rres ou 10**rdm
  //C
  //C on multiplie la taille hor de l'objet par 2
  rh = -0.30103;
  rh = 0.;
  rv = 0.;
  //C groupement des pistes par 4
  rres = -0.60206;
  //C  groupement des pistes par 2
  rres = -0.30103;
  //C non groupement des pistes
  rres = 0.;
  //C      rres=0.30103   ! multiplication du nombre de pistes par 2
  //C multiplication du nombre de pistes par 2**3
  rres = 0.30103 * 3.;
  if (rres > 5.) {
    rres = 5.;
  }
  rdm = 0.;

  //C  beam offset (in cm) at prod target in x
  shiftx = 0.0;
  //C  beam offset (in cm) at prod target in y
  shifty = 0.0;
  dx = 0.05 / pow(10., rh);
  //C      change  change change
  dxp = 10.;
  dy = 0.05 / pow(10., rv);
  dyp = 50.;
  //C doit etre une valeur entiere en %
  dpspmax = 6.;
  dpsppas = 1.0;


  //ntir = 100000;
  ntir = 1;




  //C  p en gev/c
  ppion = 1.3;
  //C      ppion=0.7     !  p en gev/c
  //C      call hfill (41,9.001,0.,10.*ppion)
  //C epaisseur 300 microns Si -> g.cm-2
  esi = 300. / 10000. * rhosi;
  //C  =800 microns carre
  resoldet = 800. / pow(10, rres) / 10000.;
  //C
  //C    *********  fin des donnees definissant le probleme
  //C
  betapion = ppion / sqrt(pow(ppion,2) + pow(amasspi,2));
  diffmul = 13.6 / ppion / betapion * sqrt(esi / x0si);
  //Cdm raporte a det2
  resoldmx = T12[idet2] * diffmul / T22[idet1];
  resoldmx = resoldmx / pow(10., rdm);
  //Cdm raporte a det2
  resoldmy = T34[idet2] * diffmul / T44[idet1];
  resoldmy = resoldmy / pow(10., rdm);
  //write(plik2,"(' diff mult=',f5.2,' mrad     resol equiv det2=',f5.2,' mm')"),
 // write (6,777) diffmul,10.*resoldet
  plik2<<"diff mul = "<<diffmul<<"resol="<<resoldet*10<<endl;
 // diffmul, 10. * resoldet;
  idineg = 0;
  //C
  //C     ********  debut boucle sur les tirages
  //C
  xmax = 0.05;
  ymax = 0.05;
  xpmax = 10.;
  ypmax = 50.;
  //C utilise uniquement pour le calcul des contributions ma
  dpmax = 5.;
  //FEM_DO_SAFE(k, 1, nelem) 
  
  for(Int_t k = 1; k < nelem; k++){
   // relem = fem::ffloat[k] + 0.1f;
   relem = float(k) + 0.1;

   
    z1 = 10. * fabs(T31[k] * xmax);
    z2 = 10. * fabs(T32[k] * xpmax);
    z3 = 10. * fabs(T33[k] * ymax);
    z4 = 10. * fabs(T34[k] * ypmax);
    z5 = 10. * fabs(T36[k] * dpmax);
    z6 = 10. * fabs(T322[k] * xpmax * xpmax);
    z7 = 10. * fabs(T323[k] * xpmax * ymax);
    z8 = 10. * fabs(T324[k] * xpmax * ypmax);
    z9 = 10. * fabs(T326[k] * xpmax * dpmax);
    z10 = 10. * fabs(T333[k] * ymax * ymax);
    z11 = 10. * fabs(T334[k] * ymax * ypmax);
    z12 = 10. * fabs(T336[k] * ymax * dpmax);
    z13 = 10. * fabs(T344[k] * ypmax * ypmax);
    z14 = 10. * fabs(T346[k] * ypmax * dpmax);
    z15 = 10. * fabs(T366[k] * dpmax * dpmax);



  }
  
  double rxf = 0,ryf = 0,rdp = 0,rxpf = 0,rypf = 0, relement = 0;
    
   TFile *output_new = new TFile("piony_new.root", "RECREATE");
   TTree *treec = new TTree("Pion5","Pionbeam5");
   treec->Branch("xi", &rxf);
   treec->Branch("yi", &ryf);
   treec->Branch("dp", &rdp);
   treec->Branch("xpi", &rxpf);
   treec->Branch("ypi", &rypf);
   treec->Branch("elemnt", &relement);
   
    
     for(Int_t itir = 0; itir < ntir; itir++)
     {
    //C
    //C    tirage gaussien dans l'extension objet en x et y
    //C
    //C    tirage uniforme en theta et phi
    //C
    //C    dp par pas de 1%  de +dpspmax a -dpspmax
    //C
    // yy = (2. * float(rand())/RAND_MAX - 1.);
    //  //cout<<yy<<endl;
    //  random->Fill(yy);
    // yy = yy + (2. * float(rand())/RAND_MAX - 1.);
    // yy = yy + (2. * float(rand())/RAND_MAX - 1.);
    // yy = yy + (2. * float(rand())/RAND_MAX - 1.);
    // yy = yy + (2. * float(rand())/RAND_MAX - 1.);
    // yy = yy * renorgauss;
    // xi = dx * yy + shiftx;
    // xpi = dxp * (2. * float(rand())/RAND_MAX - 1.);
    // yy = 2. * float(rand())/RAND_MAX - 1.;
    // yy = yy + (2. * float(rand())/RAND_MAX - 1.);
    // yy = yy + (2. * float(rand())/RAND_MAX - 1.);
    // yy = yy + (2. * float(rand())/RAND_MAX - 1.);
    // yy = yy + (2. * float(rand())/RAND_MAX - 1.);
    // yy = yy * renorgauss;
    // yi = dy * yy + shifty;
    // ypi = dyp * (2. * float(rand())/RAND_MAX - 1.);

	 xi = 0;
	 yi = 0;
	 xpi = 0;
	 ypi = 0.05;
	 

   
    //C
    //C   on genere de facon equiprobable des delta(p)/p  de dpspmax
    //C
    //C    entre  -dpspmax  et  +dpspmax  (bornes comprises)
    //C
    //C          par pas de dpsppas
    //C
    // dp = dpspmax * 2. * float(rand())/RAND_MAX - 1.;
    // dd = float(rand())/RAND_MAX;
    // ndp = 2 * int((dpspmax + 0.01) / dpsppas) + 1;
    // dp = float((ndp - 1) / 2) * dpsppas;
    // uns = 1. / float(ndp);
  //  dp = 1;

	 dp = -0.08;
  //  dd = 1;
    
    // //FEM_DO_SAFE(kii, 1, ndp) 
    // for(Int_t kii = 1; kii < ndp; kii++){
    //   seuil = uns * float(ndp - kii);
    //   if (dd < (seuil)) {
    //     dp = dp - dpsppas;
    //   }
    // }

    //C      xi=1.
    //C      xpi=0.yf=T31(k)*xi+T32(k)*xpi+T33(k)*yi+T34(k)*ypi+T36(k)*dp
    //C      yi=1.
    //C      ypi=20.
    //C      dp=-6.
    //C
    //C     ********  debut boucle sur les elements magnetiques
    //C
    iok = 1;
    //FEM_DO_SAFE(k, 1, nelem) 
    for(Int_t k = 1; k < nelem; k++){ //xf is OK!!!!!!!!!!!
	xf =0;
	//	yf =0;
	//std::setprecision(20);
	//std::cin.scientific
      z1=T11[k]*xi;
      z2=T12[k]*xpi;
      z3=T13[k]*yi;
      z4=T14[k]*ypi;
      z5=T16[k]*dp;
      //      cout << xf << endl;
      
      //cout << z1 << ":" << z2 << ":" << z3 << ":" << z4 << ":" << z5 << endl;

     // cout<<z5<<endl;
      xf=z1+z2+z3+z4+z5;
      yf=T31[k]*xi+T32[k]*xpi+T33[k]*yi+T34[k]*ypi+T36[k]*dp;  //yf=T31(k)*xi+T32(k)*xpi+T33(k)*yi+T34(k)*ypi+T36(k)*dp
    //  std::cout<<"T31 = .. "<<T31[k]<<"..T32 = .. "<<T32[k]<<"..T33 = .. "<<T33[k]<<"..T34 = .. "<<T34[k]<<"..T36 = .. "<<T36[k]<<std::fixed<<endl;
      xpf=T21[k]*xi+T22[k]*xpi+T23[k]*yi+T24[k]*ypi+T26[k]*dp;
      ypf=T41[k]*xi+T42[k]*xpi+T43[k]*yi+T44[k]*ypi+T46[k]*dp;
      xf1=xf;
      yf1=yf;
      xpf1=xpf;
      ypf1=ypf;
       
      //cout << xf << endl;

      xf=xf1+T111[k]*xi*xi+T112[k]*xi*xpi+T113[k]*xi*yi;
      xf=xf+T114[k]*xi*ypi+T116[k]*xi*dp;
      xf=xf+T122[k]*xpi*xpi+T123[k]*xpi*yi+T124[k]*xpi*ypi;
      xf=xf+T126[k]*xpi*dp;
      xf=xf+T133[k]*yi*yi+T134[k]*yi*ypi+T136[k]*yi*dp;
      xf=xf+T144[k]*ypi*ypi+T146[k]*ypi*dp;
      xf=xf+T166[k]*dp*dp;

      //cout << xf << endl;

     // xf=xf1+T126[k]*xpi*dp+T166[k]*dp*dp;
    //  xf=xf1+T126[k]*xpi*dp;

      xpf=xpf1+T211[k]*xi*xi+T212[k]*xi*xpi+T213[k]*xi*yi;
      xpf=xpf+T214[k]*xi*ypi+T216[k]*xi*dp;
      xpf=xpf+T222[k]*xpi*xpi+T223[k]*xpi*yi+T224[k]*xpi*ypi;
      xpf=xpf+T226[k]*xpi*dp;
      xpf=xpf+T233[k]*yi*yi+T234[k]*yi*ypi+T236[k]*yi*dp;
      xpf=xpf+T244[k]*ypi*ypi+T246[k]*ypi*dp;
      xpf=xpf+T266[k]*dp*dp;

      yf=yf1+T311[k]*xi*xi+T312[k]*xi*xpi+T313[k]*xi*yi;
      yf=yf+T314[k]*xi*ypi+T316[k]*xi*dp;
      yf=yf+T322[k]*xpi*xpi+T323[k]*xpi*yi+T324[k]*xpi*ypi;
      yf=yf+T326[k]*xpi*dp;
      yf=yf+T333[k]*yi*yi+T334[k]*yi*ypi+T336[k]*yi*dp;
      yf=yf+T344[k]*ypi*ypi+T346[k]*ypi*dp;
      yf=yf+T366[k]*dp*dp;

      ypf=ypf1+T411[k]*xi*xi+T412[k]*xi*xpi+T413[k]*xi*yi;
      ypf=ypf+T414[k]*xi*ypi+T416[k]*xi*dp;
      ypf=ypf+T422[k]*xpi*xpi+T423[k]*xpi*yi+T424[k]*xpi*ypi;
      ypf=ypf+T426[k]*xpi*dp;
      ypf=ypf+T433[k]*yi*yi+T434[k]*yi*ypi+T436[k]*yi*dp;
      ypf=ypf+T444[k]*ypi*ypi+T446[k]*ypi*dp;
      ypf=ypf+T466[k]*dp*dp;
      //C
      //C      yf=yf1+T346[k]*ypi*dp+T366[k]*dp*dp
      //C      yf=yf1+T346[k]*ypi*dp
     // rxf = 
      xf2 = xf;
      xpf2 = xpf;
      yf2 = yf;
      ypf2 = ypf;

       //xf is ok
    //   std::cout<<xf<<endl;
      if (fabs(xf) > (ouvrescaleh * ouvh[k])) { iok = 0;}
      if (fabs(yf) > (ouvrescalev * ouvv[k])) { iok = 0;}

      if (iok > 0) {
        ipass[k]++;
      }
      deltax = 0;
      deltay = 0;
      if (xf1 != 0) {
        deltax = fabs((xf2 - xf1) / xf1 * 100);
      }
      if (yf1 != 0) {
        deltay = fabs((yf2 - yf1) / yf1 * 100);
      }
      xhad = xf2;
      yhad = yf2;
      
      if ((k == 30) && (iok == 1)) {
        itotq9++;
      
        xf = xf1 + T111[33] * xi * xi + T112[33] * xi * xpi + T113[33] * xi * yi;
        xf = xf + T114[33] * xi * ypi + T116[33] * xi * dp;
        xf = xf + T122[33] * xpi * xpi + T123[33] * xpi * yi + T124[33] * xpi * ypi;
        xf = xf + T126[33] * xpi * dp;
        xf = xf + T133[33] * yi * yi + T134[33] * yi * ypi + T136[33] * yi * dp;
        xf = xf +T144[33] * ypi * ypi + T146[33] * ypi * dp;
        xf = xf +T166[33] * dp * dp;
        yf = yf1 + T311[33] * xi * xi + T312[33] * xi * xpi + T313[33] * xi * yi;
        yf = yf + T314[33] * xi * ypi + T316[33] * xi * dp;
        yf = yf + T322[33] * xpi * xpi + T323[33] * xpi * yi + T324[33] * xpi * ypi;
        yf = yf + T326[33] * xpi * dp;
        yf = yf + T333[33] * yi * yi + T334[33] * yi * ypi + T336[33] * yi * dp;
        yf = yf + T344[33] * ypi * ypi + T346[33] * ypi * dp;
        yf = yf + T366[33] * dp * dp;
        radius = sqrt(pow(xf,2) + pow(yf,2));
        //C      call hfill (42,radius,0.,1./282450.)
        //C      call hfill (42,radius,0.,1./240539.)
      }
      rxf = xhad;
      ryf =yhad;
      rdp = dp;
      rxpf = xpf2;
      rypf = ypf2;
      relement = k + 1;
      treec ->Fill();

      cout << xf2 << ":" << yf2 << ":" << xpf2 << ":" << ypf2 << endl;
      
      
    }
    //C
    //C     ********  fin boucle sur les elements magnetiques
    //C
    ijk = (dp + 10.01);

    if (iok > 0) {
      kpass[ijk]++;
    }
   //  if (iok == 1) {
   //    itotpass++;
    
   //    //C        call hfill (36,ydet1,ydet2,1.)
   //    //C
   //    //C       discretisation du detecteur 1
   //    xdet1sav = xdet1;
   //    xca = (xdet1 + 0.5 * resoldet + offset) / resoldet;
   //    nca = xca;
   //    xdet1 = float(nca) * resoldet - offset;
   //    //C
   //    //C       discretisation du detecteur 2
   //    xca = (xdet2 + 0.5 * resoldet + offset) / resoldet;
   //    nca = xca;
   //    xdet2 = float(nca) * resoldet - offset;
   //    xdet2sav = xdet2;
   //    //C
   //    //C       diffusion multiple
   //    yy = (2. * float(rand())/RAND_MAX - 1.);
   //    yy = yy + (2. * float(rand())/RAND_MAX - 1.);
   //    yy = yy + (2. * float(rand())/RAND_MAX - 1.);
   //    yy = yy + (2. * float(rand())/RAND_MAX - 1.);
   //    yy = yy + (2. * float(rand())/RAND_MAX - 1.);
   //    yy = yy * renorgauss;
   //    //C diff mult gaussien
   //    xdet2 = xdet2 + resoldmx * yy;
   //    //C
   //    //C       discretisation du detecteur 1 en y
   //    //C
   //    ydet1sav = ydet1;
   //    ydet2sav = ydet2;

   //    alfa = -(T16[idet1] * T126[idet2] - T16[idet2] * T126[idet1]);
   //    beta = T16[idet2] * T12[idet1] - xdet2 * T126[idet1];
   //    beta = beta - T16[idet1] * T12[idet2] + xdet1 * T126[idet2];
   //    gamma = xdet1 * T12[idet2] - xdet2 * T12[idet1];
   //    discri = pow(beta,2) - 4. * alfa * gamma;
   //    galfa = T126[idet1] * T166[idet2] - T126[idet2] * T166[idet1];
   //    gbeta = T12[idet1] * T166[idet2] - T12[idet2] * T166[idet1];
   //    gbeta = gbeta - T16[idet1] * T126[idet2] + T16[idet2] * T126[idet1];
   //    ggamma = T12[idet1] * T16[idet2] - T12[idet2] * T16[idet1];
   //    ggamma = ggamma + xdet1 * T126[idet2] - xdet2 * T126[idet1];
   //    gdelta = xdet1 * T12[idet2] - xdet2 * T12[idet1];
   //    aa1 = 1.;
   //    aa2 = gbeta / galfa;
   //    aa3 = ggamma / galfa;
   //    aa4 = gdelta / galfa;
   //    sol = 3.;
   //    iter = 0;
   //    //C        write (6,755) iter,sol,dp
   //    yyy1 = galfa * pow(sol,3) + gbeta * pow(sol,2) +
   //      ggamma * sol + gdelta;
   //    yyp1 = 3. * galfa * pow(sol,2) + 2. * gbeta * sol + ggamma;
   //    sol = sol - yyy1 / yyp1;
   //    iter = iter + 1;
   //    //C       write (6,755) iter,sol,dp
   //    yyy2 = galfa * pow(sol,3) + gbeta * pow(sol,2) +
   //      ggamma * sol + gdelta;
   //    yyp2 = 3. * galfa * pow(sol,2) + 2. * gbeta * sol + ggamma;
   //    sol = sol - yyy2 / yyp2;
   //    iter = iter + 1;
   //    //C        write (6,755) iter,sol,dp
   //    yyy3 = galfa * pow(sol,3) + gbeta * pow(sol,2) +
   //      ggamma * sol + gdelta;
   //    yyp3 = 3. * galfa * pow(sol,2) + 2. * gbeta * sol + ggamma;
   //    sol = sol - yyy3 / yyp3;
   //    iter = iter + 1;
   //    //C        write (6,755) iter,sol,dp
   //    yyy4 = galfa * pow(sol,3) + gbeta * pow(sol,2) +
   //      ggamma * sol + gdelta;
   //    yyp4 = 3. * galfa * pow(sol,2) + 2. * gbeta * sol + ggamma;
   //    sol = sol - yyy4 / yyp4;
   //    iter = iter + 1;
   //    //C        write (6,755) iter,sol,dp
   //    yyy5 = galfa * pow(sol,3) + gbeta * pow(sol,2) +
   //      ggamma * sol + gdelta;
   //    yyp5 = 3. * galfa * pow(sol,2) + 2. * gbeta * sol + ggamma;
   //    sol = sol - yyy5 / yyp5;
   //    iter = iter + 1;
   //    //C        write (6,755) iter,sol,dp
   //    teta = (xdet1 - T16[idet1] * sol - T166[idet1] * sol * sol);
   //    teta = teta / (T12[idet1] + T126[idet1] * sol);
     
   //  //  goto statment_999;
   //    aa = T33[idet1];
   //    bb = T34[idet1] + T346[idet1] * sol + T324[idet1] * teta;
   //    cc = ydet1 - T36[idet1] * sol;
   //    cc = cc - T326[idet1] * teta * sol - T366[idet1] * sol * sol;
   //    aap = T33[idet2];
   //    bbp = T34[idet2] + T346[idet2] * sol + T324[idet2] * teta;
   //    ccp = ydet2 - T36[idet2] * sol;
   //    ccp = ccp - T326[idet2] * teta * sol - T366[idet2] * sol * sol;
   //    deter = aa * bbp - aap * bb;
   //    yci = (cc * bbp - bb * ccp) / deter;
   //    phi = (aa * ccp - cc * aap) / deter;
   // //   statment_999:

   //    poids = 100. / float(ntir) * float(ndp);
   //    poids = 100. / 102962.;

   //    xf = xf1 + T111[33] * xi * xi + T112[33] * xi * xpi + T113[33] * xi * yi;
   //    xf = xf + T114[33] * xi * ypi + T116[33] * xi * dp;
   //    xf = xf + T122[33] * xpi * xpi + T123[33] * xpi * yi + T124[33] * xpi * ypi;
   //    xf = xf + T126[33] * xpi * dp;
   //    xf = xf + T133[33] * yi * yi + T134[33] * yi * ypi + T136[33] * yi * dp;
   //    xf = xf + T144[33] * ypi * ypi + T146[33] * ypi * dp;
   //    xfth = xf + T166[33] * dp * dp;
   //    yf = yf1 + T311[33] * xi * xi + T312[33] * xi * xpi + T313[33] * xi * yi;
   //    yf = yf + T314[33] * xi * ypi + T316[33] * xi * dp;
   //    yf = yf + T322[33] * xpi * xpi + T323[33] * xpi * yi + T324[33] * xpi * ypi;
   //    yf = yf + T326[33] * xpi * dp;
   //    yf = yf + T333[33] * yi * yi + T334[33] * yi * ypi + T336[33] * yi * dp;
   //    yf = yf + T344[33] * ypi * ypi + T346[33] * ypi * dp;
   //    yfth = yf + T366[33] * dp * dp;
   //    //C
   //    xf = xf1 + T111[33] * 0. * 0. + T112[33] * 0. * teta + T113[33] * 0. * yci;
   //    xf = xf + T114[33] * 0. * phi + T116[33] * 0. * sol;
   //    xf = xf + T122[33] * teta * teta + T123[33] * teta * yci + T124[33] * teta * phi;
   //    xf = xf + T126[33] * teta * sol;
   //    xf = xf + T133[33] * yci * yci + T134[33] * yci * phi + T136[33] * yci * sol;
   //    xf = xf + T144[33] * phi * phi + T146[33] * phi * sol;
   //    xfcal = xf + T166[33] * sol * sol;
   //    yf = yf1 + T311[33] * 0. * 0. + T312[33] * 0. * teta + T313[33] * 0. * yci;
   //    yf = yf + T314[33] * 0. * phi + T316[33] * 0. * sol;
   //    yf = yf + T322[33] * teta * teta + T323[33] * teta * yci + T324[33] * teta * phi;
   //    yf = yf + T326[33] * teta * sol;
   //    yf = yf + T333[33] * yci * yci + T334[33] * yci * phi + T336[33] * yci * sol;
   //    yf = yf + T344[33] * phi * phi + T346[33] * phi * sol;
   //    yfcal = yf + T366[33] * sol * sol;
   //    rr = sqrt(pow((xfth - xfcal),2) + pow((yfth - yfcal),2)); // matlab.hpp
   //    ipi = dp + 6.001;
   //    ips2 = ipi / 1;
   //    idif = ipi - ips2 * 1;
   //    //C        if (idif.eq.0) call hfill (230+ips2,10.*yfth,0.,1.)
   //    xabs = 10. * xfth;
   //    yabs = 10. * (xfcal - xfth) + dp * 10.;

   //    xabs = 10. * yfth;
   //    yabs = 10. * (yfcal - yfth) + dp * 5.;

   //  }
    
  }
  random->Write();
  double rrpass = 0;
  Int_t dpsp = 0;
  
   TFile *output2 = new TFile("piony2.root", "RECREATE");
   TTree *treez = new TTree("Pion2","Pionbeam2");
   treez->Branch("alossfraction", &rreste);
   treez->Branch("dpsp", &retape);
 // output1->Write();
  /////////////////////////////////////////////////////////

  
  for(Int_t ki=1; ki<nelem; ki++){
                                reste = 100. * float(ipass[ki]) / float(ntir);
                                plik2<<"etape.."<<ki<<"..reste.."<<reste<<endl;
                                ri = float(ki) + 0.01;
				rreste = reste;
				retape = ki+1;
				if (retape > 0) treex->Fill();
				

  }

  for(Int_t kij = 1; kij<16; kij++){
                                   rpass[kij] = float(kpass[kij]) * 15. / float(ntir) * 100.;
                                   plik2<<"dpsp.."<<kij - 7<<"..transsmision.."<<rpass[kij]<<endl;
				   rrpass = rpass[kij];
				   dpsp = (kij - 7)+1;
				   if (dpsp > 0) treez->Fill();
  }
                                   plik2<<"idineg...."<<idineg<<endl;
                                   alossfraction = float(itotq9 - itotpass) / float(itotq9);
                                   plik2<<"fraction hors cible.."<<alossfraction * 100.;    
     
   
     plik2.close();
     output->Write();
     output1->Write();
     output2->Write();
     output_new->Write();
}
