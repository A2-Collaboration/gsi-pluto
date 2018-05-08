////////////////////////////////////////////////////////
// Read in urqmd f15 file (collision history) and put the particles
// on the PLUTO bulk stack. Either dileptons (D0,D+,w,rho0,phi,eta,eta'
// and pi0) can be decayed via shining (translated for original UrQMD
// fortran routines, authors: Katharina Schmidt, Sascha Vogel,
// Christian Sturm 2007-2009) into e+/e- or by PLUTO bulk decay.
// Optional stable non dileptons can be written out. PHUrReader adds
// a TClonesArray "Addon" to the PLUTO Tree that contains original
// UrQMD infos about creation and absorbtion times, densities etc.
//
// USAGE: see urqmd_f15_input.C
////////////////////////////////////////////////////////

#include "PHUrReader.h"
#include "TSystem.h"
#include "PStaticData.h"
#include "PDynamicData.h"

#include "PDistributionManager.h"

ClassImp(PHUrReader)


TClonesArray *PHUrReader::fAdd   = 0;
TClonesArray *PHUrReader::fEvent = 0;

PHUrReader::PHUrReader(TString inputfile) {

    Input(inputfile);
    collisionsIn .reserve(20000) ;  // reserve enough space to avoid copying of objects during dynamic growing
    collisionsOut.reserve(20000);   // reserve enough space to avoid copying of objects during dynamic growing
    collisions   .reserve(2000);    // reserve enough space to avoid copying of objects during dynamic growing
    printEvtHeader = kTRUE;
    printColHeader = kFALSE;
    printParticle  = kFALSE;
    fEvtCt         = 0;
    fOutputAscii   = NULL;
    outputNonDiLeptons       = kFALSE;
    outputLeptons            = kTRUE;
    outputFreezeoutDiLeptons = kFALSE;
    fNumMax  = 10000;
    /*
     Taken from Unigen:
     UrQMD PIDs are in fact composite - a particle is fully defined by the
     type specifier (ityp), the charge (ichg) and in case of baryons, the
     third component isospin (iso3; ignored by us at present). For
     simplicity, our conversion tables collapse these into a single number
     as follows:
      - propagate the sign of ityp (particle-antiparticle distinction for
        baryons, strangeness-antistrangeness distinction for mesons) to that
        of the stored value;
      - shift the ichg range from -2..2 (UrQMD does not support nuclear
        fragments other than protons and neutrons so all particles it
        produces fall in this range) to 0..4 to make sure it doesn't
        interfere with the above;
      - multiply shifted charge by +/-1000 and add it to type. The latter
        is guaranteed to be smaller than 1000 (baryon types are one- or
        two-digit, meson types three-digit) so that way no ambiguities
        occur.
  Int_t id;
  if (ityp >= 0)
    id = 1000 * (ichg + 2) + ityp;
  else
    id = -1000 * (ichg + 2) + ityp;
    */
    //
    pdg_param = 0;
    pid_param = 0;

    mUrqmdToPdg[1017   ] = 1114;
    mUrqmdToPdg[1018   ] =31114;
    mUrqmdToPdg[1019   ] =1112 ;
    mUrqmdToPdg[1020   ] =11114;
    mUrqmdToPdg[1021   ] =11112;
    mUrqmdToPdg[1022   ] =1116 ;
    mUrqmdToPdg[1023   ] =21112;
    mUrqmdToPdg[1024   ] =21114;
    mUrqmdToPdg[1025   ] =11116;
    mUrqmdToPdg[1026   ] =1118 ;
    mUrqmdToPdg[1040   ] =3112 ;
    mUrqmdToPdg[1041   ] =3114 ;
    mUrqmdToPdg[1042   ] =13112;
    mUrqmdToPdg[1043   ] =13114;
    mUrqmdToPdg[1044   ] =23112;
    mUrqmdToPdg[1045   ] =3116 ;
    mUrqmdToPdg[1046   ] =13116;
    mUrqmdToPdg[1047   ] =23114;
    mUrqmdToPdg[1048   ] =3118 ;
    mUrqmdToPdg[1049   ] =3312 ;
    mUrqmdToPdg[1050   ] =3314 ;
    mUrqmdToPdg[1051   ] =23314;
    mUrqmdToPdg[1052   ] =13314;
    mUrqmdToPdg[1053   ] =33314;
    mUrqmdToPdg[1054   ] =13316;
    mUrqmdToPdg[1055   ] =3334 ;
    mUrqmdToPdg[1101   ] =-211 ;
    mUrqmdToPdg[1104   ] =-213 ;
    mUrqmdToPdg[1111   ] =-9000211;
    mUrqmdToPdg[1114   ] =-20213 ;
    mUrqmdToPdg[1118   ] =-215   ;
    mUrqmdToPdg[1122   ] =-10213 ;
    mUrqmdToPdg[1126   ] =-100213;
    mUrqmdToPdg[1130   ] =-30213 ;
    mUrqmdToPdg[2001   ] =2112   ;
    mUrqmdToPdg[2002   ] =12112  ;
    mUrqmdToPdg[2003   ] =1214   ;
    mUrqmdToPdg[2004   ] =22112  ;
    mUrqmdToPdg[2005   ] =32112  ;
    mUrqmdToPdg[2006   ] =2116   ;
    mUrqmdToPdg[2007   ] =12116  ;
    mUrqmdToPdg[2008   ] =21214  ;
    mUrqmdToPdg[2009   ] =42112  ;
    mUrqmdToPdg[2010   ] =31214  ;
    mUrqmdToPdg[2011   ] =41214  ;
    mUrqmdToPdg[2012   ] =12118  ;
    mUrqmdToPdg[2013   ] =52114  ;
    mUrqmdToPdg[2016   ] =100012110;
    mUrqmdToPdg[2017   ] =2114    ;
    mUrqmdToPdg[2018   ] =32114   ;
    mUrqmdToPdg[2019   ] =1212    ;
    mUrqmdToPdg[2020   ] =12114   ;
    mUrqmdToPdg[2021   ] =11212   ;
    mUrqmdToPdg[2022   ] =1216    ;
    mUrqmdToPdg[2023   ] =21212   ;
    mUrqmdToPdg[2024   ] =22114   ;
    mUrqmdToPdg[2025   ] =11216   ;
    mUrqmdToPdg[2026   ] =2118    ;
    mUrqmdToPdg[2027   ] =3122    ;
    mUrqmdToPdg[2028   ] =13122   ;
    mUrqmdToPdg[2029   ] =3124    ;
    mUrqmdToPdg[2030   ] =23122   ;
    mUrqmdToPdg[2031   ] =33122   ;
    mUrqmdToPdg[2032   ] =13124   ;
    mUrqmdToPdg[2033   ] =43122   ;
    mUrqmdToPdg[2034   ] =53122   ;
    mUrqmdToPdg[2035   ] =3126    ;
    mUrqmdToPdg[2036   ] =13126   ;
    mUrqmdToPdg[2037   ] =23124   ;
    mUrqmdToPdg[2038   ] =3128    ;
    mUrqmdToPdg[2039   ] =23126   ;
    mUrqmdToPdg[2040   ] =3212    ;
    mUrqmdToPdg[2041   ] =3214    ;
    mUrqmdToPdg[2042   ] =13212   ;
    mUrqmdToPdg[2043   ] =13214   ;
    mUrqmdToPdg[2044   ] =23212   ;
    mUrqmdToPdg[2045   ] =3216    ;
    mUrqmdToPdg[2046   ] =13216   ;
    mUrqmdToPdg[2047   ] =23214   ;
    mUrqmdToPdg[2048   ] =3218    ;
    mUrqmdToPdg[2049   ] =3322    ;
    mUrqmdToPdg[2050   ] =3324    ;
    mUrqmdToPdg[2051   ] =23324   ;
    mUrqmdToPdg[2052   ] =13324   ;
    mUrqmdToPdg[2053   ] =33324   ;
    mUrqmdToPdg[2054   ] =13326   ;
    mUrqmdToPdg[2100   ] =22      ;
    mUrqmdToPdg[2101   ] =111     ;
    mUrqmdToPdg[2102   ] =221     ;
    mUrqmdToPdg[2103   ] =223     ;
    mUrqmdToPdg[2104   ] =113     ;
    mUrqmdToPdg[2105   ] =9000221 ;
    mUrqmdToPdg[2106   ] =311     ;
    mUrqmdToPdg[2107   ] =331     ;
    mUrqmdToPdg[2108   ] =313     ;
    mUrqmdToPdg[2109   ] =333     ;
    mUrqmdToPdg[2110   ] =333     ;
    mUrqmdToPdg[2111   ] =9000111 ;
    mUrqmdToPdg[2112   ] =10221   ;
    mUrqmdToPdg[2113   ] =20313   ;
    mUrqmdToPdg[2114   ] =20113   ;
    mUrqmdToPdg[2115   ] =20223   ;
    mUrqmdToPdg[2116   ] =20333   ;
    mUrqmdToPdg[2117   ] =315     ;
    mUrqmdToPdg[2118   ] =115     ;
    mUrqmdToPdg[2119   ] =225     ;
    mUrqmdToPdg[2120   ] =335     ;
    mUrqmdToPdg[2121   ] =10313   ;
    mUrqmdToPdg[2122   ] =10113   ;
    mUrqmdToPdg[2123   ] =10223   ;
    mUrqmdToPdg[2124   ] =10333   ;
    mUrqmdToPdg[2125   ] =100313  ;
    mUrqmdToPdg[2126   ] =100113  ;
    mUrqmdToPdg[2127   ] =100223  ;
    mUrqmdToPdg[2128   ] =100333  ;
    mUrqmdToPdg[2129   ] =30313   ;
    mUrqmdToPdg[2130   ] =30113   ;
    mUrqmdToPdg[2131   ] =30223   ;
    mUrqmdToPdg[2132   ] =337     ;
    mUrqmdToPdg[3001   ] =2212    ;
    mUrqmdToPdg[3002   ] =12212   ;
    mUrqmdToPdg[3003   ] =2124    ;
    mUrqmdToPdg[3004   ] =22212   ;
    mUrqmdToPdg[3005   ] =32212   ;
    mUrqmdToPdg[3006   ] =2216    ;
    mUrqmdToPdg[3007   ] =12216   ;
    mUrqmdToPdg[3008   ] =22124   ;
    mUrqmdToPdg[3009   ] =42212   ;
    mUrqmdToPdg[3010   ] =32124   ;
    mUrqmdToPdg[3011   ] =42124   ;
    mUrqmdToPdg[3012   ] =12218   ;
    mUrqmdToPdg[3013   ] =52214   ;
    mUrqmdToPdg[3016   ] =100012210;
    mUrqmdToPdg[3017   ] =2214     ;
    mUrqmdToPdg[3018   ] =32214    ;
    mUrqmdToPdg[3019   ] =2122     ;
    mUrqmdToPdg[3020   ] =12214    ;
    mUrqmdToPdg[3021   ] =12122    ;
    mUrqmdToPdg[3022   ] =2126     ;
    mUrqmdToPdg[3023   ] =22122    ;
    mUrqmdToPdg[3024   ] =22214    ;
    mUrqmdToPdg[3025   ] =12126    ;
    mUrqmdToPdg[3026   ] =2218     ;
    mUrqmdToPdg[3040   ] =3222     ;
    mUrqmdToPdg[3041   ] =3224     ;
    mUrqmdToPdg[3042   ] =13222    ;
    mUrqmdToPdg[3043   ] =13224    ;
    mUrqmdToPdg[3044   ] =23222    ;
    mUrqmdToPdg[3045   ] =3226     ;
    mUrqmdToPdg[3046   ] =13226    ;
    mUrqmdToPdg[3047   ] =23224    ;
    mUrqmdToPdg[3048   ] =3228     ;
    mUrqmdToPdg[3101   ] =211      ;
    mUrqmdToPdg[3104   ] =213      ;
    mUrqmdToPdg[3106   ] =321      ;
    mUrqmdToPdg[3108   ] =323      ;
    mUrqmdToPdg[3110   ] =10321    ;
    mUrqmdToPdg[3111   ] =9000211  ;
    mUrqmdToPdg[3113   ] =20323    ;
    mUrqmdToPdg[3114   ] =20213    ;
    mUrqmdToPdg[3117   ] =325      ;
    mUrqmdToPdg[3118   ] =215      ;
    mUrqmdToPdg[3121   ] =10323    ;
    mUrqmdToPdg[3122   ] =10213    ;
    mUrqmdToPdg[3125   ] =100323   ;
    mUrqmdToPdg[3126   ] =100213   ;
    mUrqmdToPdg[3129   ] =30323    ;
    mUrqmdToPdg[3130   ] =30213    ;
    mUrqmdToPdg[4017   ] =2224     ;
    mUrqmdToPdg[4018   ] =32224    ;
    mUrqmdToPdg[4019   ] =2222     ;
    mUrqmdToPdg[4020   ] =12224    ;
    mUrqmdToPdg[4021   ] =12222    ;
    mUrqmdToPdg[4022   ] =2226     ;
    mUrqmdToPdg[4023   ] =22222    ;
    mUrqmdToPdg[4024   ] =22224    ;
    mUrqmdToPdg[4025   ] =12226    ;
    mUrqmdToPdg[4026   ] =2228     ;
    mUrqmdToPdg[-3055  ] = -3334   ;
    mUrqmdToPdg[-3054  ] = -13316  ;
    mUrqmdToPdg[-3053  ] = -33314  ;
    mUrqmdToPdg[-3052  ] = -13314  ;
    mUrqmdToPdg[-3051  ] = -23314  ;
    mUrqmdToPdg[-3050  ] = -3314   ;
    mUrqmdToPdg[-3049  ] = -3312   ;
    mUrqmdToPdg[-3048  ] = -3118   ;
    mUrqmdToPdg[-3047  ] = -23114  ;
    mUrqmdToPdg[-3046  ] = -13116  ;
    mUrqmdToPdg[-3045  ] = -3116   ;
    mUrqmdToPdg[-3044  ] = -23112  ;
    mUrqmdToPdg[-3043  ] = -13114  ;
    mUrqmdToPdg[-3042  ] = -13112  ;
    mUrqmdToPdg[-3041  ] = -3114   ;
    mUrqmdToPdg[-3040  ] = -3112   ;
    mUrqmdToPdg[-3026  ] = -1118   ;
    mUrqmdToPdg[-3025  ] = -11116  ;
    mUrqmdToPdg[-3024  ] = -21114  ;
    mUrqmdToPdg[-3023  ] = -21112  ;
    mUrqmdToPdg[-3022  ] = -1116   ;
    mUrqmdToPdg[-3021  ] = -11112  ;
    mUrqmdToPdg[-3020  ] = -11114  ;
    mUrqmdToPdg[-3019  ] = -1112   ;
    mUrqmdToPdg[-3018  ] = -31114  ;
    mUrqmdToPdg[-3017  ] = -1114   ;
    mUrqmdToPdg[-2129  ] = -30313  ;
    mUrqmdToPdg[-2125  ] = -100313 ;
    mUrqmdToPdg[-2121  ] = -10313  ;
    mUrqmdToPdg[-2117  ] = -315    ;
    mUrqmdToPdg[-2113  ] = -20313  ;
    mUrqmdToPdg[-2110  ] = -10311  ;
    mUrqmdToPdg[-2108  ] = -313    ;
    mUrqmdToPdg[-2106  ] = -311    ;
    mUrqmdToPdg[-2055  ] = -3334   ;
    mUrqmdToPdg[-2054  ] = -13326  ;
    mUrqmdToPdg[-2053  ] = -33324  ;
    mUrqmdToPdg[-2052  ] = -13324  ;
    mUrqmdToPdg[-2051  ] = -23324  ;
    mUrqmdToPdg[-2050  ] = -3324   ;
    mUrqmdToPdg[-2049  ] = -3322   ;
    mUrqmdToPdg[-2048  ] = -3218   ;
    mUrqmdToPdg[-2047  ] = -23214  ;
    mUrqmdToPdg[-2046  ] = -13216  ;
    mUrqmdToPdg[-2045  ] = -3216   ;
    mUrqmdToPdg[-2044  ] = -23212  ;
    mUrqmdToPdg[-2043  ] = -13214  ;
    mUrqmdToPdg[-2042  ] = -13212  ;
    mUrqmdToPdg[-2041  ] = -3214   ;
    mUrqmdToPdg[-2040  ] = -3212   ;
    mUrqmdToPdg[-2039  ] = -23126  ;
    mUrqmdToPdg[-2038  ] = -3128   ;
    mUrqmdToPdg[-2037  ] = -23124  ;
    mUrqmdToPdg[-2036  ] = -13126  ;
    mUrqmdToPdg[-2035  ] = -3126   ;
    mUrqmdToPdg[-2034  ] = -53122  ;
    mUrqmdToPdg[-2033  ] = -43122  ;
    mUrqmdToPdg[-2032  ] = -13124  ;
    mUrqmdToPdg[-2031  ] = -33122  ;
    mUrqmdToPdg[-2030  ] = -23122  ;
    mUrqmdToPdg[-2029  ] = -3124   ;
    mUrqmdToPdg[-2028  ] = -13122  ;
    mUrqmdToPdg[-2027  ] = -3122   ;
    mUrqmdToPdg[-2026  ] = -2118   ;
    mUrqmdToPdg[-2025  ] = -11216  ;
    mUrqmdToPdg[-2024  ] = -22114  ;
    mUrqmdToPdg[-2023  ] = -21212  ;
    mUrqmdToPdg[-2022  ] = -1216   ;
    mUrqmdToPdg[-2021  ] = -11212  ;
    mUrqmdToPdg[-2020  ] = -12114  ;
    mUrqmdToPdg[-2019  ] = -1212   ;
    mUrqmdToPdg[-2018  ] = -32114  ;
    mUrqmdToPdg[-2017  ] = -2114   ;
    mUrqmdToPdg[-2016  ] = -100012110;
    mUrqmdToPdg[-2013  ] = -52114    ;
    mUrqmdToPdg[-2012  ] = -12118    ;
    mUrqmdToPdg[-2011  ] = -41214    ;
    mUrqmdToPdg[-2010  ] = -31214    ;
    mUrqmdToPdg[-2009  ] = -42112    ;
    mUrqmdToPdg[-2008  ] = -21214    ;
    mUrqmdToPdg[-2007  ] = -12116    ;
    mUrqmdToPdg[-2006  ] = -2116     ;
    mUrqmdToPdg[-2005  ] = -32112    ;
    mUrqmdToPdg[-2004  ] = -22112    ;
    mUrqmdToPdg[-2003  ] = -1214     ;
    mUrqmdToPdg[-2002  ] = -12112    ;
    mUrqmdToPdg[-2001  ] = -2112     ;
    mUrqmdToPdg[-1129  ] = -30323    ;
    mUrqmdToPdg[-1125  ] = -100323   ;
    mUrqmdToPdg[-1121  ] = -10323    ;
    mUrqmdToPdg[-1117  ] = -325      ;
    mUrqmdToPdg[-1113  ] = -20323    ;
    mUrqmdToPdg[-1110  ] = -10321    ;
    mUrqmdToPdg[-1108  ] = -323      ;
    mUrqmdToPdg[-1106  ] = -321      ;
    mUrqmdToPdg[-1048  ] = -3228     ;
    mUrqmdToPdg[-1047  ] = -23224    ;
    mUrqmdToPdg[-1046  ] = -13226    ;
    mUrqmdToPdg[-1045  ] = -3226     ;
    mUrqmdToPdg[-1044  ] = -23222    ;
    mUrqmdToPdg[-1043  ] = -13224    ;
    mUrqmdToPdg[-1042  ] = -13222    ;
    mUrqmdToPdg[-1041  ] = -3224     ;
    mUrqmdToPdg[-1040  ] = -3222     ;
    mUrqmdToPdg[-1026  ] = -2218     ;
    mUrqmdToPdg[-1025  ] = -12126    ;
    mUrqmdToPdg[-1024  ] = -22214    ;
    mUrqmdToPdg[-1023  ] = -22122    ;
    mUrqmdToPdg[-1022  ] = -2126     ;
    mUrqmdToPdg[-1021  ] = -12122    ;
    mUrqmdToPdg[-1020  ] = -12214    ;
    mUrqmdToPdg[-1019  ] = -2122     ;
    mUrqmdToPdg[-1018  ] = -32214    ;
    mUrqmdToPdg[-1017  ] = -2214     ;
    mUrqmdToPdg[-1016  ] = -100012210;
    mUrqmdToPdg[-1013  ] = -52214    ;
    mUrqmdToPdg[-1012  ] = -12218    ;
    mUrqmdToPdg[-1011  ] = -42124    ;
    mUrqmdToPdg[-1010  ] = -32124    ;
    mUrqmdToPdg[-1009  ] = -42212    ;
    mUrqmdToPdg[-1008  ] = -22124    ;
    mUrqmdToPdg[-1007  ] = -12216    ;
    mUrqmdToPdg[-1006  ] = -2216     ;
    mUrqmdToPdg[-1005  ] = -32212    ;
    mUrqmdToPdg[-1004  ] = -22212    ;
    mUrqmdToPdg[-1003  ] = -2124     ;
    mUrqmdToPdg[-1002  ] = -12212    ;
    mUrqmdToPdg[-1001  ] = -2212     ;
    mUrqmdToPdg[-26   ] =-2228       ;
    mUrqmdToPdg[-25   ] =-12226      ;
    mUrqmdToPdg[-24   ] =-22224      ;
    mUrqmdToPdg[-23   ] =-22222      ;
    mUrqmdToPdg[-22   ] =-2226       ;
    mUrqmdToPdg[-21   ] =-12222      ;
    mUrqmdToPdg[-20   ] =-12224      ;
    mUrqmdToPdg[-19   ] =-2222       ;
    mUrqmdToPdg[-18   ] =-32224      ;
    mUrqmdToPdg[-17   ] =-2224       ;


    dilep.mUrqmdToPdg = &mUrqmdToPdg;

    mUrqmdProcess[1]  = "NN->ND";
    mUrqmdProcess[2]  = "NN->NN*";
    mUrqmdProcess[3]  = "NN->ND*";
    mUrqmdProcess[4]  = "NN->DD";
    mUrqmdProcess[5]  = "NN->DN*";
    mUrqmdProcess[6]  = "NN->DD*";
    mUrqmdProcess[7]  = "NN->N*N*,N*D*,D*D*";
    mUrqmdProcess[8]  = "ND->DD";
    mUrqmdProcess[10] = "MB->B'";
    mUrqmdProcess[11] = "MM->M'";
    mUrqmdProcess[13] =  "BB (but not pp,pn) elastic scattering";
    mUrqmdProcess[14] = "inelastic scattering (no string excitation)";
    mUrqmdProcess[15] = "BB->2 strings";
    mUrqmdProcess[17] = "pn-elastic";
    mUrqmdProcess[19] = "pp-elastic";
    mUrqmdProcess[20] = "decay";
    mUrqmdProcess[22] = "BBar elastic";
    mUrqmdProcess[23] = "BBar annihilation->1 string";
    mUrqmdProcess[24] = "BBar diffractive->2 strings";
    mUrqmdProcess[26] = "MB elastic scattering";
    mUrqmdProcess[27] = "MB,MM->1 string";
    mUrqmdProcess[28] = "MB,MM->2 string";
    mUrqmdProcess[30] = "ND->NN";
    mUrqmdProcess[31] = "DD->DN";
    mUrqmdProcess[32] = "DD->NN";
    mUrqmdProcess[35] = "ND inelastic";
    mUrqmdProcess[36] = "Danielewicz forward delay (MB->B')";
    mUrqmdProcess[37] = "Danielewicz forward delay (MM->M')";
    mUrqmdProcess[38] = "MM elastic scattering";
    mUrqmdProcess[39] = "BBar inelastic scattering (no annihilation)";
}

PHUrReader::~PHUrReader() {
    if(fInputAscii.is_open()) fInputAscii.close();
}

const Char_t *PHUrReader::GetUrQMDProcess(Int_t id) {
    map<Int_t,TString>::iterator it = mUrqmdProcess.find(id);
    if(it!=mUrqmdProcess.end()) return it->second.Data();
    else  Error("GetUrQMDProcess()", "Process %i unknown!", id);
    return 0;
}

void PHUrReader::ClearVector(vector<vector<PHUrParticle> >& v) {

    for(UInt_t i = 0; i < v.size(); i++){
	v[i].clear();
    }
    v.clear();
}

void PHUrReader::Input(TString filename) {

    if(filename.CompareTo("") != 0) {

	if(gSystem->AccessPathName(filename.Data()) != 0){
	    Error("Input()", "File %s not does not exist!",filename.Data());
	} else {
	    if(fInputAscii.is_open()) fInputAscii.close();
	    fInputName = filename;
	    fInputAscii.open(filename.Data());
	}
    } else Error("Input()", "File name not specified !");
}

void PHUrReader::Output(TString filename) {
    if(filename.CompareTo("") != 0) {
	fOutputName = filename;
    } else Error("Output()", "File name not specified !");
}

void PHUrReader::SetMapID(Int_t Id_in, Int_t Id_out) {
    // remap UrQMD id to another id
    if(mapIDs.find(Id_in) == mapIDs.end()){
	mapIDs[Id_in] = Id_out;
    } else {
	Error("SetMapID()", "ID=%i already mapped! Will be ignored ....", Id_in);
    }
}

Bool_t PHUrReader::ReadEvent() {
    // read function for UrQMD f15 files

    if(!fInputAscii.is_open() || !fInputAscii.good()) return kFALSE;

    if(fEvtCt == 0){
	// only in first event. In next
	// events the eventheader will be read
        // by the previous call.
	evtheader.Clear(NULL);
	if(!evtheader.Read(fInputAscii)){
	    Error("ReadEvent()", "Could not read event header!");
	    return kFALSE;
	}
    } else {
        evtheader = evtheaderCP;
    }

    if(evtheader.id == -1 || evtheader.id == 0) {     // detected  eventheader

	if(printEvtHeader) { evtheader.Print(); }

	ClearVector(collisionsIn);  // clear collisions of previous event
	ClearVector(collisionsOut);
	collisions.clear();

	while ( !fInputAscii.eof() ) {
            colheader.Clear(NULL);
	    if(!colheader.Read(fInputAscii)) return kFALSE;   // read next collision header

            if(colheader.n_in == -1 || colheader.n_in == 0){  // detected  eventheader
		evtheaderCP.Copy(colheader);                  // move previous evtheader to current eventheader
                fEvtCt++;
		return kTRUE;
	    }

	    collisions.push_back(colheader);

	    if(printColHeader) colheader.Print();

	    //-------------------------------------------------------------------------------
	    // read in-going particles
	    vector< PHUrParticle > InParticles;

	    for(Int_t i = 0; i < colheader.n_in; i++) {
		PHUrParticle part;
		if(!part.Read(fInputAscii)) { 
		    Error("ReadEvent()", "Could not read in-going particle!"); 
		    return kFALSE; 
		}

                Int_t pdg = 1000;
                Int_t id  = part.pdg;
		if(mUrqmdToPdg.find(id) != mUrqmdToPdg.end()) pdg = mUrqmdToPdg[id];
		part.pdg = pdg;

		if(printParticle) part.Print();
		InParticles.push_back(part);
	    }
	    collisionsIn.push_back(InParticles);
	    //-------------------------------------------------------------------------------

	    //-------------------------------------------------------------------------------
	    // read out-going particles
	    vector< PHUrParticle > OutParticles;

	    for(Int_t i = 0; i < colheader.n_out; i++) {
		PHUrParticle part;
		if(!part.Read(fInputAscii)) { 
		    Error("ReadEvent()", "Could not read out-going particle!"); 
		    return kFALSE; 
		}

		Int_t pdg = 1000;
                Int_t id  = part.pdg;
		if(mUrqmdToPdg.find(id) != mUrqmdToPdg.end()) pdg = mUrqmdToPdg[id];
		part.pdg = pdg;

		if(printParticle)part.Print();
		OutParticles.push_back(part);
	    }
	    collisionsOut.push_back(OutParticles);
	    //-------------------------------------------------------------------------------

	}
        fEvtCt++; // last event only
        return kTRUE;

    } else {
	Error("ReadEvent()", "This is not an event header!");
	return kFALSE;
    }

}

Bool_t PHUrReader::Modify(PParticle **mstack, int *, int *num, int stacksize) {

    if(!ReadEvent()) return kFALSE;
    

    //----------------------------------------------------------------------------------
    // for the first call : create event structure
    if(fEvtCt == 1) {

	makeDistributionManager()->Exec("pdg:init");
	pdg_param = makeDataBase()->GetParamInt("pdg");
	if (pdg_param < 0) return kFALSE;
	pid_param = makeDataBase()->GetParamInt("pid");

	TTree *T = GetTree();

	 if(T) {
             dilep.outputLeptons = outputLeptons;
             dilep.reader        = this;
	     if(outputLeptons) dilep.Output(fInputName,fOutputName);

	     if(stacksize == 0){
		 Error("Modify()", "FATAL ERROR : stacksize == 0!, Check your statcksize setiings for overrun!");
                 exit(1);
	     }

	     CreateParticleArray(T,stacksize );
	     CreateAddonArray   (T,stacksize );

	 } else {
	     Error("Modify()", "retrieved NULL pointer for storage tree!");
             return kFALSE;
	 }
    }
    //----------------------------------------------------------------------------------

    //###########################################################################################
    // loop over urqmd event
    //c.. other
    Int_t flag       = -1 ;
    Int_t flag_noLep = -1 ;
    //c.. start analysis
    //c.. loop over all interactions ww

    UInt_t numberCol = collisions.size();
    fNumPart      = 0;
    fNumPartAddon = 0;
    fEvent->Clear();
    fAdd  ->Clear();

    //----------------------------------------------------------------------------------
    // map the occurence of one particle by index (for expample elastic scattering does not
    // change the index of the particle, but momentum)


    vector<PHUrParticle*> particles;
    map   <Int_t,PHUrParticle*> mLastParticle;

    vector<Int_t> particlesCt;
    particles  .resize(20000,0);
    particlesCt.resize(20000,0);

    Int_t nparticle = 0;
    for(UInt_t icol = 0; icol < numberCol; icol++) {  // collision
	for(UInt_t iin = 0; iin < collisionsIn[icol].size(); iin++) {   // in-going particles
	    PHUrParticle *particle = &collisionsIn [icol][iin];
	    if(particles[particle->ind] == 0) {
		particles[particle->ind] = particle;
		nparticle++;
		particle->first = 0;
                mLastParticle[particle->ind] = particle;
	    } else {
                mLastParticle[particle->ind] = particle;
	    }
            particlesCt[particle->ind]++;
	    particle->instance = particlesCt[particle->ind];

	    particles[particle->ind]->listInCollIndex.push_back(icol);
	}

	for(UInt_t iout = 0; iout <  collisionsOut[icol].size(); iout++) { // out-going particles
	    PHUrParticle *particle = &collisionsOut [icol][iout];
	    if(particles[particle->ind] == 0) {
		particles[particle->ind] = particle;
                nparticle++;
		particlesCt[particle->ind]++;
		particle->first = 1;
		mLastParticle[particle->ind] = particle;
	    } else {
		mLastParticle[particle->ind] = particle;
	    }
	    particle->instance = particlesCt[particle->ind];
	    particles[particle->ind]->listOutCollIndex.push_back(icol);
	}
    }

    for(map< Int_t, PHUrParticle*>::iterator iter = mLastParticle.begin(); iter != mLastParticle.end(); ++iter ) {
	if((*iter).second ) (*iter).second->last=1;
    }

    particles  .erase(particles  .begin()+nparticle,  particles.end());
    particlesCt.erase(particlesCt.begin()+nparticle,particlesCt.end());

    //----------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------
    if(outputNonDiLeptons == 2) {
        Int_t ct = 0;
	for(UInt_t icol = 0; icol < numberCol; icol++) {  // collision
	    for(UInt_t iin = 0; iin < collisionsIn[icol].size(); iin++) {   // in-going particles
		PHUrParticle &particle = collisionsIn [icol][iin];
		PHUrCollisionHeader &collision = collisions[icol];

		if(particle.id == dilep.omega ||
		   particle.id == dilep.delta ||   // omega || delta || rho || phi
		   particle.id == dilep.rho   ||
		   particle.id == dilep.phi   ||
		   particle.id == dilep.pion     ||
		   particle.id == dilep.etaprime ||
		   particle.id == dilep.eta
		   ) {  // dilepton interest
		    continue;
		}
		
		if(ct == fNumMax-50) {
                    continue;
		}

		Int_t index = 0;
		PParticle *pMoth    = CreateParticle(index);
		PHUrAddon *pMothAdd = CreateAddon(index);

		pMoth->Reset(particle.pdg, particle.Px(), particle.Py(), particle.Pz(), particle.E(), 1.);

		pMothAdd->t_cre       = particle.t;
		pMothAdd->t_abs       = particle.t;
		pMothAdd->dens_cre    = collision.baryon_density;
		pMothAdd->dens_abs    = collision.baryon_density;
		pMothAdd->fr_cre      = particle.fr;
		pMothAdd->fr_abs      = particle.fr;
		pMothAdd->stable      = 0;
		pMothAdd->n_in        = collision.n_in;
		pMothAdd->n_out       = collision.n_out;
		pMothAdd->process_id  = collision.process_id;
		pMothAdd->n_collision = collision.n_collision;
		pMothAdd->instance    = particle.instance;
		pMothAdd->first       = particle.first;
		pMothAdd->last        = particle.last;

	    }

	    for(UInt_t iout = 0; iout <  collisionsOut[icol].size(); iout++) { // out-going particles
	    
		PHUrParticle &particle = collisionsOut [icol][iout];
		PHUrCollisionHeader &collision = collisions[icol];

		if(particle.id == dilep.omega ||
		   particle.id == dilep.delta ||   // omega || delta || rho || phi
		   particle.id == dilep.rho   ||
		   particle.id == dilep.phi   ||
		   particle.id == dilep.pion     ||
		   particle.id == dilep.etaprime ||
		   particle.id == dilep.eta
		   ) {  // dilepton interest
		    continue;
		}

		if(ct == fNumMax-50) {
                    continue;
		}

		Int_t index = 0;
		PParticle *pMoth    = CreateParticle(index);
		PHUrAddon *pMothAdd = CreateAddon(index);

		pMoth->Reset(particle.pdg, particle.Px(), particle.Py(), particle.Pz(), particle.E(), 1.);

		pMothAdd->t_cre       = particle.t;
		pMothAdd->t_abs       = particle .t;
		pMothAdd->dens_cre    = collision.baryon_density;
		pMothAdd->dens_abs    = collision.baryon_density;
		pMothAdd->fr_cre      = particle .fr;
		pMothAdd->fr_abs      = particle .fr;
		pMothAdd->stable      = 0;
		pMothAdd->n_in        = collision.n_in;
		pMothAdd->n_out       = collision.n_out;
		pMothAdd->process_id  = collision.process_id;
		pMothAdd->n_collision = collision.n_collision;
		pMothAdd->instance    = particle.instance;
		pMothAdd->first       = particle.first;
		pMothAdd->last        = particle.last;

	    }
	}
    }
    //----------------------------------------------------------------------------------



    for(UInt_t icol = 0; icol < numberCol; icol++) {  // collision
	for(UInt_t iout = 0; iout <  collisionsOut[icol].size(); iout++) { // out-going particles

            //----------------------------------------------------------------------------------
	    // try to find the particle as incoming particle in the
            // following collisions
	    for(UInt_t icol2 = icol+1; icol2 < numberCol; icol2++) {   // following collisions
		for(UInt_t iin = 0; iin < collisionsIn[icol2].size(); iin++) {   // in-going particles
	
		    PHUrParticle &particleIn  = collisionsIn [icol2][iin];
		    PHUrParticle &particleOut = collisionsOut[icol ][iout];
		    
		    if (particleIn.id == particleOut.id && particleIn.Pz() == particleOut.Pz()) {
			// same particle in following collision
			
			// found the same particle
			Bool_t dilepton = kFALSE;
			
			if(particleIn.id == dilep.omega ||
			   particleIn.id == dilep.delta ||   // omega || delta || rho || phi
			   particleIn.id == dilep.rho   ||
			   particleIn.id == dilep.phi) {  // dilepton interest
			    
			    dilep.evtheader    = &evtheader;
			    dilep.particleIn   = &particleIn;
			    dilep.particleOut  = &particleOut;
			    dilep.collisionIn  = &collisions[icol2];
			    dilep.collisionOut = &collisions[icol ];
			    
			    if(!outputFreezeoutDiLeptons) { 
				dilep.Dilep(); 
			    } else    
				flag = 0;             // not final state
			    
			    dilepton = kTRUE;
			    
			} // if dilepton interest
			
			if(particleIn.id == dilep.pion     ||
			   particleIn.id == dilep.etaprime ||
			   particleIn.id == dilep.eta) {  // pion || etaprime || eta
			    flag = 0;
			    dilepton = kTRUE;
			}
			
			if(!dilepton) flag_noLep = 0;			
			
		    } // if same particle in following collision
		} // for iin
	    } // for col2
            //----------------------------------------------------------------------------------
	    
	    PHUrParticle &particleOut = collisionsOut[icol ][iout];
	    
            //----------------------------------------------------------------------------------
	    // in case of freeze out dileptons we have to
            // take care
	    if(particleOut.id == dilep.omega ||
	       particleOut.id == dilep.delta ||   // omega || delta || rho || phi
	       particleOut.id == dilep.rho   ||
	       particleOut.id == dilep.phi) {  // dilepton interest
		
		if(outputFreezeoutDiLeptons && flag != 0)  {
		    
		    dilep.evtheader    = &evtheader;
		    dilep.particleIn   = &particleOut;
		    dilep.particleOut  = &particleOut;
		    dilep.collisionIn  = &collisions[icol];
		    dilep.collisionOut = &collisions[icol];
		    
		    dilep.Dilep();
		}
	    }
            //----------------------------------------------------------------------------------


	    if(particleOut.id == dilep.pion     ||
	       particleOut.id == dilep.etaprime ||
	       particleOut.id == dilep.eta) {  // pion || etaprime || eta

		if(flag != 0) { //  particle has not undergone any following collsion (last state)

		    dilep.evtheader    = &evtheader;
		    dilep.particleIn   = &particleOut;
		    dilep.particleOut  = &particleOut;
		    dilep.collisionIn  = &collisions[icol];
		    dilep.collisionOut = &collisions[icol];

		    dilep.Dilep();
		} // if flag != 0
	    } else {

		if( !(particleOut.id == dilep.omega ||
		      particleOut.id == dilep.delta ||   // omega || delta || rho || phi
		      particleOut.id == dilep.rho   ||
		      particleOut.id == dilep.phi )
		    
		    ) {
		    
		    if (flag_noLep != 0) {  // final out non lepton

			if (outputNonDiLeptons == 1) {

			    PHUrCollisionHeader &collisionOut = collisions[icol];
                            Int_t index = 0;
			    PParticle *pMoth    = CreateParticle(index);
			    PHUrAddon *pMothAdd = CreateAddon(index);

			    pMoth->Reset(particleOut.pdg, particleOut.Px(), particleOut.Py(),
					 particleOut.Pz(), particleOut.E(), 1.);

			    pMothAdd->t_cre       = particleOut.t;
			    pMothAdd->t_abs       = particleOut .t;
			    pMothAdd->dens_cre    = collisionOut.baryon_density;
			    pMothAdd->dens_abs    = collisionOut.baryon_density;
			    pMothAdd->fr_cre      = particleOut .fr;
			    pMothAdd->fr_abs      = particleOut .fr;
			    pMothAdd->stable      = 2;
			    pMothAdd->n_in        = collisionOut.n_in;
                            pMothAdd->n_out       = collisionOut.n_out;
                            pMothAdd->process_id  = collisionOut.process_id;
                            pMothAdd->n_collision = collisionOut.n_collision;
			    pMothAdd->instance    = particleOut.instance;
			    pMothAdd->first       = particleOut.first;
                            pMothAdd->last        = particleOut.last;

			}
		    }
		}
	    }

	    flag       = 1;
            flag_noLep = 1;
	}  // for iout
    } // for icol

    // end loop over urqmd event
    //###########################################################################################


    Int_t nParticle = fEvent->GetEntries();
    Int_t startInd  = *num;
    Int_t cur       = startInd;
    Int_t *i_result = 0;
    for (Int_t i = 0; i < nParticle; i++) {

	if (i == stacksize) {
	    Warning("Modify", "Stack size too small, increase '_system_particle_stacksize'");
	    return kTRUE;
	}

	PParticle *part = (PParticle*) fEvent->At(i);
	Int_t id = part->ID();   // already urqmd -> pdg code
	
	if (!makeDataBase()->GetParamInt (pdg_param, id, pid_param, &i_result)) {

	    if(find(funknownIDs.begin(),funknownIDs.end(),id) == funknownIDs.end()){
		funknownIDs.push_back(id);
		Warning("Modify()",
			"Particle with PDG code %i is not defined in Pluto, further particles will get the 'dummy' pid", id);
	    }

	} else {

            Int_t id = *i_result;
	    if(mapIDs.find(*i_result) != mapIDs.end()){ id = mapIDs[*i_result]; }
	    part->SetID(id);

	}

	*(mstack[cur]) = *part;

	cur++;
    }

    *num += cur-startInd;

    return kTRUE;
}


void PHUrReader::CreateAddonArray(TTree *T, Int_t stacksize) {
   // create Addon TClonesArray and add it to the PLUTO tree

    fAdd = new TClonesArray("PHUrAddon", stacksize);
    T->Branch("Addon", &fAdd, 5000, 99);
    fNumMax = stacksize;
}

PHUrAddon *PHUrReader::CreateAddon(Int_t &index) {
    // create an Addon object at index
    if (fNumPartAddon == fNumMax){
	Warning("CreateAddon()", "Reached max stacksize, no addon created!");
        return 0;
    }

    new((*fAdd)[fNumPartAddon++]) PHUrAddon();
    index = fNumPartAddon - 1;
    return (PHUrAddon*)((*fAdd)[fNumPartAddon - 1]);
}

void PHUrReader::CreateParticleArray(TTree *, Int_t stacksize) {
   // create UrQMD PParticle local TClonesArray
    fEvent  = new TClonesArray( "PParticle", stacksize);
    //T->Branch("Event", &fEvent,8000,99);
    fNumMax = stacksize;
}

PParticle *PHUrReader::CreateParticle(Int_t &index) {
    // create PParticle  object in local event at index

    if (fNumPart == fNumMax){
	Warning("CreateAddon()", "Reached max stacksize, no addon created!");
        return 0;
    }

    new((*fEvent)[fNumPart++]) PParticle();
    index = fNumPart - 1;
    return (PParticle*)((*fEvent)[fNumPart - 1]);
}

Int_t PHUrReader::UserAnalysis(PParticle**,Int_t npart) {

    cout << "userAna nparticle = " << npart << endl;

    return 1;
}

