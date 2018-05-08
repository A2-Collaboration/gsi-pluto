/////////////////////////////////////////////////////////////////////
//
//Plugin to add the pdg code to the Pluto data base
//Required for the unigen input source
//
//////////////////////////////////////////////////////////////////////

#include "PPDGPlugin.h"


PPDGPlugin::PPDGPlugin(const Char_t *id, const Char_t *de):
    PDistributionCollection(id, de) {
    is_initialized       = 0;
    is_extend_resonances = 0;
}

Bool_t PPDGPlugin::Activate(void) {
    return kTRUE;
}

PPDGPlugin::~PPDGPlugin() {
}


Bool_t PPDGPlugin::ExecCommand(const char *command, Double_t) {

    if (strcmp (command,"extend_resonances") == 0) {
	if (!is_extend_resonances) {
	    is_extend_resonances = 1;
	    //============================= N*(1680)0 ==================================
	    makeStaticData()->AddParticle(80,"NF150", 1.68);   //PID, name, mass
	    makeStaticData()->AddAlias("NF150","N*(1680)0");
	    makeStaticData()->SetParticleTotalWidth("NF150", 0.13);
	    makeStaticData()->SetParticleBaryon("NF150", 1);
	    makeStaticData()->SetParticleSpin("NF150", 5);
	    makeStaticData()->SetParticleCharge("NF150", 0);
	    makeStaticData()->SetParticleIsospin("NF150", 1);
	    makeStaticData()->SetParticleParity("NF150", 1);
	    //-------------------- N,pi ---------------
	    makeStaticData()->AddDecay("NF150 --> p + pi-", "NF150",
				       "p,pi-", 0.65/3. );   ////PDG - {0.65-0.7}
	    makeStaticData()->AddDecay("NF150 --> n + pi+", "NF150",
				       "n,pi0", 0.65*2./3. );
	    //------------------ D1232,pi -------------
	    makeStaticData()->AddDecay("NF150 --> Delta+ + pi-", "NF150",
				       "D+,pi-", 0.097482/2. );	//PDG - {0.05-0.15}
	    makeStaticData()->AddDecay("NF150 --> Delta0 + pi0", "NF150",
				       "D0,pi0", 0.097482/3. );
	    makeStaticData()->AddDecay("NF150 --> Delta- + pi+", "NF150",
				       "D-,pi+", 0.097482/6. );
	    //------------------ N,rho -------------
	    makeStaticData()->AddDecay("NF150 --> p + rho-", "NF150",
				       "p,rho-", 0.1/3. );		//PDG - {0.03-0.15}		
	    makeStaticData()->AddDecay("NF150 --> n + rho0", "NF150",
				       "n,rho0", 0.1*2./3. );
	    //------------------ p,pi,pi -------------
	    makeStaticData()->AddDecay("NF150 --> n + pi+ + pi-", "NF150",
				       "p,pi+,pi-", 0.1475*2./3. );		//PDG - {0.05-0.20}
	    makeStaticData()->AddDecay("NF150 --> n + pi0 + pi0", "NF150",
				       "p,pi0,pi0", 0.1475/3. );
	    //------------------ p,gamma -------------
	    makeStaticData()->AddDecay("NF150 --> n + g", "NF150",
				       "p,g", 0.0025 );	//PDG - {0.0021-0.0032}
	    //0.0025/137 = 1.8e-5 left for Dalitz Decay

	    //============================= N*(1680)+ ==================================
	    makeStaticData()->AddParticle(81,"NF15+", 1.68);   //PID, name, mass
	    makeStaticData()->AddAlias("NF15+","N*(1680)+");
	    makeStaticData()->SetParticleTotalWidth("NF15+", 0.13);
	    makeStaticData()->SetParticleBaryon("NF15+", 1);
	    makeStaticData()->SetParticleSpin("NF15+", 5);
	    makeStaticData()->SetParticleCharge("NF15+", 1);
	    makeStaticData()->SetParticleIsospin("NF15+", 1);
	    makeStaticData()->SetParticleParity("NF15+", 1);
	    //-------------------- N,pi ---------------
	    makeStaticData()->AddDecay("NF15+ --> p + pi0", "NF15+",
				       "p,pi0", 0.65/3. );   ////PDG - {0.65-0.7}
	    makeStaticData()->AddDecay("NF15+ --> n + pi+", "NF15+",
				       "n,pi+", 0.65*2./3. );
	    //------------------ D1232,pi -------------
	    makeStaticData()->AddDecay("NF15+ --> Delta++ + pi-", "NF15+",
				       "D++,pi-", 0.097482/2. );	//PDG - {0.05-0.15}
	    makeStaticData()->AddDecay("NF15+ --> Delta+ + pi0", "NF15+",
				       "D+,pi0", 0.097482/3. );
	    makeStaticData()->AddDecay("NF15+ --> Delta0 + pi+", "NF15+",
				       "D0,pi+", 0.097482/6. );
	    //------------------ N,rho -------------
	    makeStaticData()->AddDecay("NF15+ --> p + rho0", "NF15+",
				       "p,rho0", 0.1/3. );		//PDG - {0.03-0.15}		
	    makeStaticData()->AddDecay("NF15+ --> n + rho+", "NF15+",
				       "n,rho+", 0.1*2./3. );
	    //------------------ p,pi,pi -------------
	    makeStaticData()->AddDecay("NF15+ --> p + pi+ + pi-", "NF15+",
				       "p,pi+,pi-", 0.1475*2./3. );		//PDG - {0.05-0.20}
	    makeStaticData()->AddDecay("NF15+ --> p + pi0 + pi0", "NF15+",
				       "p,pi0,pi0", 0.1475/3. );
	    //------------------ p,gamma -------------
	    makeStaticData()->AddDecay("NF15+ --> p + g", "NF15+",
				       "p,g", 0.0025 );	//PDG - {0.0021-0.0032}
	    //0.0025/137 = 1.8e-5 left for Dalitz Decay

	    //=========================================================================
	    //=========================================================================
	    //=========================================================================

	    //============================= D(1700)- ==================================
	    makeStaticData()->AddParticle(82,"DD33-", 1.7);   //PID, name, mass
	    makeStaticData()->AddAlias("DD33-","Delta(1700)+");
	    makeStaticData()->SetParticleTotalWidth("DD33-", 0.3);
	    makeStaticData()->SetParticleBaryon("DD33-", 1);
	    makeStaticData()->SetParticleSpin("DD33-", 3);
	    makeStaticData()->SetParticleCharge("DD33-", -1);
	    makeStaticData()->SetParticleIsospin("DD33-", 3);
	    makeStaticData()->SetParticleParity("DD33-", -1);
	    //-------------------- N,pi ---------------
	    makeStaticData()->AddDecay("DD33- --> n + pi-", "DD33-",
				       "n,pi-", 0.15 );		//PDG - {0.1-0.2}
	    //------------------ D1232,pi -------------
	    makeStaticData()->AddDecay("DD33- --> Delta0 + pi-", "DD33-",
				       "D0,pi-", 4.495854e-01*2./5. );//PDG - {0.3-0.6}
	    makeStaticData()->AddDecay("DD33- --> Delta- + pi0", "DD33-",
				       "D-,pi0", 4.495854e-01*3./5. );
	    //------------------ N,rho -------------
	    makeStaticData()->AddDecay("DD33- --> n + rho-", "DD33-",
				       "n,rho-", 0.45 );	//PDG - {0.3-0.55}

	    //============================= D(1700)0 ==================================
	    makeStaticData()->AddParticle(83,"DD330", 1.7);   //PID, name, mass
	    makeStaticData()->AddAlias("DD330","Delta(1700)+");
	    makeStaticData()->SetParticleTotalWidth("DD330", 0.3);
	    makeStaticData()->SetParticleBaryon("DD330", 1);
	    makeStaticData()->SetParticleSpin("DD330", 3);
	    makeStaticData()->SetParticleCharge("DD330", 0);
	    makeStaticData()->SetParticleIsospin("DD330", 3);
	    makeStaticData()->SetParticleParity("DD330", -1);
	    //-------------------- N,pi ---------------
	    makeStaticData()->AddDecay("DD330 --> p + pi-", "DD330",
				       "p,pi-", 0.15/3. );		//PDG - {0.1-0.2}
	    makeStaticData()->AddDecay("DD330 --> n + pi0", "DD330",
				       "n,pi0", 0.15*2./3. );
	    //------------------ D1232,pi -------------
	    makeStaticData()->AddDecay("DD330 --> Delta+ + pi-", "DD330",
				       "D+,pi-", 4.495854e-01*8./15. );//PDG - {0.3-0.6}
	    makeStaticData()->AddDecay("DD330 --> Delta0 + pi0", "DD330",
				       "D0,pi0", 4.495854e-01/15. );
	    makeStaticData()->AddDecay("DD330 --> Delta- + pi+", "DD330",
				       "D-,pi+", 4.495854e-01*2./5. );
	    //"D0,pi+", 0.5897*8./15. );
	    //------------------ N,rho -------------
	    makeStaticData()->AddDecay("DD330 --> n + rho0", "DD330",
				       "n,rho0", 0.45*2./3. );	//PDG - {0.3-0.55}
	    makeStaticData()->AddDecay("DD330 --> p + rho-", "DD330",
				       "p,rho-", 0.45/3. );
	    //------------------ p,gamma -------------
	    makeStaticData()->AddDecay("DD330 --> n + g", "DD330",
				       "n,g", 0.002 );


	    //============================= D(1700)+ ==================================
	    makeStaticData()->AddParticle(84,"DD33+", 1.7);   //PID, name, mass
	    makeStaticData()->AddAlias("DD33+","Delta(1700)+");
	    makeStaticData()->SetParticleTotalWidth("DD33+", 0.3);
	    makeStaticData()->SetParticleBaryon("DD33+", 1);
	    makeStaticData()->SetParticleSpin("DD33+", 3);
	    makeStaticData()->SetParticleCharge("DD33+", 1);
	    makeStaticData()->SetParticleIsospin("DD33+", 3);
	    makeStaticData()->SetParticleParity("DD33+", -1);
	    //-------------------- N,pi ---------------
	    makeStaticData()->AddDecay("DD33+ --> p + pi0", "DD33+",
				       "p,pi0", 0.15*2./3. );		//PDG - {0.1-0.2}
	    makeStaticData()->AddDecay("DD33+ --> n + pi+", "DD33+",
				       "n,pi+", 0.15/3. );
	    //------------------ D1232,pi -------------
	    makeStaticData()->AddDecay("DD33+ --> Delta++ + pi-", "DD33+",
				       "D++,pi-", 4.495854e-01*2./5. );//PDG - {0.3-0.6}
	    makeStaticData()->AddDecay("DD33+ --> Delta+ + pi0", "DD33+",
				       "D+,pi0", 4.495854e-01/15. );
	    makeStaticData()->AddDecay("DD33+ --> Delta0 + pi+", "DD33+",
				       "D0,pi+", 4.495854e-01*8./15. );
	    //------------------ N,rho -------------
	    makeStaticData()->AddDecay("DD33+ --> p + rho0", "DD33+",
				       "p,rho0", 0.45*2./3. );	//PDG - {0.3-0.55}
	    makeStaticData()->AddDecay("DD33+ --> n + rho+", "DD33+",
				       "n,rho+", 0.45/3. );
	    //------------------ p,gamma -------------
	    makeStaticData()->AddDecay("DD33+ --> p + g", "DD33+",
				       "p,g", 0.002 );
	    //0.002/137 = 1.46e-5 left for Dalitz Decay

	    //============================= D(1700)++ ==================================
	    makeStaticData()->AddParticle(85,"DD33++", 1.7);   //PID, name, mass
	    makeStaticData()->AddAlias("DD33++","Delta(1700)++");
	    makeStaticData()->SetParticleTotalWidth("DD33++", 0.3);
	    makeStaticData()->SetParticleBaryon("DD33++", 1);
	    makeStaticData()->SetParticleSpin("DD33++", 3);
	    makeStaticData()->SetParticleCharge("DD33++", 1);
	    makeStaticData()->SetParticleIsospin("DD33++", 3);
	    makeStaticData()->SetParticleParity("DD33++", -1);
	    //-------------------- N,pi ---------------
	    makeStaticData()->AddDecay("DD33++ --> p + pi+", "DD33++",
				       "p,pi+", 0.15 );		//PDG - {0.1-0.2}
	    //------------------ D1232,pi -------------
	    makeStaticData()->AddDecay("DD33++ --> Delta++ + pi0", "DD33++",
				       "D++,pi0", 4.495854e-01*3./5. );//PDG - {0.3-0.6}
	    makeStaticData()->AddDecay("DD33++ --> Delta+ + pi+", "DD33++",
				       "D+,pi+", 4.495854e-01*2./5. );
	    //------------------ N,rho -------------
	    makeStaticData()->AddDecay("DD33++ --> p + rho+", "DD33++",
				       "p,rho+", 0.45 );	//PDG - {0.3-0.55}
	    //0.002/137 = 1.46e-5 left for Dalitz Decay

	    //=========================================================================
	    //=========================================================================
	    //=========================================================================

	    //============================= D(1905)- ==================================
	    makeStaticData()->AddParticle(86,"DF35-", 1.89);   //PID, name, mass
	    makeStaticData()->AddAlias("DF35-","Delta(1905)+");
	    makeStaticData()->SetParticleTotalWidth("DF35-", 0.33);
	    makeStaticData()->SetParticleBaryon("DF35-", 1);
	    makeStaticData()->SetParticleSpin("DF35-", 5);
	    makeStaticData()->SetParticleCharge("DF35-", -1);
	    makeStaticData()->SetParticleIsospin("DF35-", 3);
	    makeStaticData()->SetParticleParity("DF35-", 1);
	    //-------------------- N,pi ---------------
	    makeStaticData()->AddDecay("DF35- --> n + pi-", "DF35-",
				       "n,pi-", 0.12 );		//PDG - {0.09-0.15}
	    //------------------ D1232,pi -------------
	    makeStaticData()->AddDecay("DF35- --> Delta0 + pi-", "DF35-",
				       "D0,pi-", 1.479854e-01*2./5. );//PDG - {0.3-0.6}
	    makeStaticData()->AddDecay("DF35- --> Delta- + pi0", "DF35-",
				       "D-,pi0", 1.479854e-01*3./5. );
	    //------------------ N,rho -------------
	    makeStaticData()->AddDecay("DF35- --> n + rho-", "DF35-",
				       "n,rho-", 0.73 );	//PDG - {0.3-0.55}

	    //============================= D(1905)0 ==================================
	    makeStaticData()->AddParticle(87,"DF350", 1.89);   //PID, name, mass
	    makeStaticData()->AddAlias("DF350","Delta(1905)+");
	    makeStaticData()->SetParticleTotalWidth("DF350", 0.33);
	    makeStaticData()->SetParticleBaryon("DF350", 1);
	    makeStaticData()->SetParticleSpin("DF350", 5);
	    makeStaticData()->SetParticleCharge("DF350", 0);
	    makeStaticData()->SetParticleIsospin("DF350", 3);
	    makeStaticData()->SetParticleParity("DF350", 1);
	    //-------------------- N,pi ---------------
	    makeStaticData()->AddDecay("DF350 --> p + pi-", "DF350",
				       "p,pi-", 0.12/3. );		//PDG - {0.09-0.15}
	    makeStaticData()->AddDecay("DF350 --> n + pi0", "DF350",
				       "n,pi0", 0.12*2./3. );		
	    //------------------ D1232,pi -------------
	    makeStaticData()->AddDecay("DF350 --> Delta+ + pi-", "DF350",
				       "D+,pi-", 1.479854e-01*8./15. );//PDG - {0.3-0.6}
	    makeStaticData()->AddDecay("DF350 --> Delta0 + pi0", "DF350",
				       "D0,pi0", 1.479854e-01/15. );
	    makeStaticData()->AddDecay("DF350 --> Delta- + pi+", "DF350",
				       "D-,pi+", 1.479854e-01*2./5. );
	    //"D0,pi+", 0.5897*8./15. );
	    //------------------ N,rho -------------
	    makeStaticData()->AddDecay("DF350 --> n + rho0", "DF350",
				       "n,rho0", 0.73*2./3. );	//PDG - {0.3-0.55}
	    makeStaticData()->AddDecay("DF350 --> p + rho-", "DF350",
				       "p,rho-", 0.73/3. );
	    //------------------ p,gamma -------------
	    makeStaticData()->AddDecay("DF350 --> n + g", "DF350",
				       "n,g", 0.002 );


	    //============================= D(1905)+ ==================================
	    makeStaticData()->AddParticle(88,"DF35+", 1.89);   //PID, name, mass
	    makeStaticData()->AddAlias("DF35+","Delta(1905)+");
	    makeStaticData()->SetParticleTotalWidth("DF35+", 0.33);
	    makeStaticData()->SetParticleBaryon("DF35+", 1);
	    makeStaticData()->SetParticleSpin("DF35+", 5);
	    makeStaticData()->SetParticleCharge("DF35+", 1);
	    makeStaticData()->SetParticleIsospin("DF35+", 3);
	    makeStaticData()->SetParticleParity("DF35+", 1);
	    //-------------------- N,pi ---------------
	    makeStaticData()->AddDecay("DF35+ --> p + pi0", "DF35+",
				       "p,pi0", 0.12*2./3. );		//PDG - {0.09-0.15}
	    makeStaticData()->AddDecay("DF35+ --> n + pi+", "DF35+",
				       "n,pi+", 0.12/3. );		
	    //------------------ D1232,pi -------------
	    makeStaticData()->AddDecay("DF35+ --> Delta++ + pi-", "DF35+",
				       "D++,pi-", 1.479854e-01*2./5. );//PDG - {<25%}
	    makeStaticData()->AddDecay("DF35+ --> Delta+ + pi0", "DF35+",
				       "D+,pi0", 1.479854e-01/15. );
	    makeStaticData()->AddDecay("DF35+ --> Delta0 + pi+", "DF35+",
				       "D0,pi+", 1.479854e-01*8./15. );
	    //------------------ N,rho -------------
	    makeStaticData()->AddDecay("DF35+ --> p + rho0", "DF35+",
				       "p,rho0", 0.73*2./3. );	//PDG - {>60%}
	    makeStaticData()->AddDecay("DF35+ --> n + rho+", "DF35+",
				       "n,rho+", 0.73/3. );
	    //------------------ p,gamma -------------
	    makeStaticData()->AddDecay("DF35+ --> p + g", "DF35+",
				       "p,g", 0.002 );
	    //0.002/137 = 1.46e-5 left for Dalitz Decay

	    //============================= D(1905)++ ==================================
	    makeStaticData()->AddParticle(89,"DF35++", 1.89);   //PID, name, mass
	    makeStaticData()->AddAlias("DF35++","Delta(1905)++");
	    makeStaticData()->SetParticleTotalWidth("DF35++", 0.33);
	    makeStaticData()->SetParticleBaryon("DF35++", 1);
	    makeStaticData()->SetParticleSpin("DF35++", 5);
	    makeStaticData()->SetParticleCharge("DF35++", 1);
	    makeStaticData()->SetParticleIsospin("DF35++", 3);
	    makeStaticData()->SetParticleParity("DF35++", 1);
	    //-------------------- N,pi ---------------
	    makeStaticData()->AddDecay("DF35++ --> p + pi+", "DF35++",
				       "p,pi+", 0.12 );		//PDG - {0.09-0.15}
	    //------------------ D1232,pi -------------
	    makeStaticData()->AddDecay("DF35++ --> Delta++ + pi0", "DF35++",
				       "D++,pi0", 1.479854e-01*3./5. );//PDG - {0.3-0.6}
	    makeStaticData()->AddDecay("DF35++ --> Delta+ + pi+", "DF35++",
				       "D+,pi+", 1.479854e-01*2./5. );
	    //------------------ N,rho -------------
	    makeStaticData()->AddDecay("DF35++ --> p + rho+", "DF35++",
				       "p,rho+", 0.73 );	//PDG - {0.3-0.55}
	    //0.002/137 = 1.46e-5 left for Dalitz Decay
	    return kTRUE;
	}
    }

    if (strcmp (command,"init") == 0) {
	if (!is_initialized) {
	    is_initialized = 1;

	    makeDataBase()->MakeParamInt("pdg", "PDG code");

	    //this arry was copied from UniGen (pluto_pdg.dat)

	    Int_t pdg[70] = {0,22,
			     -11,
			     11,
			     12,
			     -13,
			     13,
			     111,
			     211,
			     -211,
			     130,
			     321,
			     -321,
			     2112,
			     2212,
			     -2212,
			     310,
			     221,
			     3122,
			     3222,
			     3212,
			     3112,
			     3322,
			     3312,
			     3334,
			     -2112,
			     -3122,
			     -3112,
			     -3212,
			     -3222,
			     -3322,
			     -3312,
			     -3334 ,
			     0,
			     2114,
			     2224,
			     2214,
			     1114,
			     0,
			     0,
			     0,
			     113,
			     213,
			     -213,
			     0,
			     0,
			     0,
			     0,
			     0,
			     0,
			     0,
			     0,
			     223,
			     331,
			     0,
			     333,
			     0,
			     0,
			     0,
			     0,
			     0,
			     0,
			     0,
			     0,
			     0,
			     0,
			     0,
			     443,
			     100443,
			     0};
	    
	    for (int i=0; i<71; i++) { //loop over known particles
		
		int pkey = makeStaticData()->GetParticleKey(i);
		
		if (pkey > 0) {
		    if (!makeDataBase()->SetParamInt (pkey, "pdg", new Int_t(pdg[i])))
			return kFALSE;
// 		    else
// 			cout << pkey << ":" << pdg[i] << endl;
		}
		
	    }
	}
	return kTRUE;
    }

    
    return kFALSE;
}



ClassImp(PPDGPlugin)



    
