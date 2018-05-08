/////////////////////////////////////////////////////////////////////
// 
// Dalitz decay for N* resonances
// 
//Reference:
//Dilepton decays of baryon resonances.
//M. Zetenyi, Gy. Wolf, (Budapest, RMKI) . Feb 2002. 16pp.
//Published in Heavy Ion Phys.17:27-39,2003.
//e-Print: nucl-th/0202047
//
//
//                                  Author:  A. Rustamov
/////////////////////////////////////////////////////////////////////

#include "PResonanceDalitz.h"

PDistribution* PResonanceDalitz::Clone(const char*) const {
    //clone the object
    return new PResonanceDalitz((const PResonanceDalitz &)* this);
};

PResonanceDalitz::PResonanceDalitz(const Char_t *id, const Char_t *de, Int_t key) : PDalitzDecay(id, de, key) {
    //Constructor
    mass_x = mass_p;
    spin   = makeStaticData()->GetParticleSpin(parent_id);
    par    = makeStaticData()->GetParticleParity(parent_id);
    charge = makeStaticData()->GetParticleCharge(parent_id);
};

double PResonanceDalitz::GetMatrixT(int &spin, int &parity, int&, const double& ecm, const double& M) {
    float matr;
    if (spin == 1) {
	if (parity == 1) {
	    matr = TMath::Power((ecm2 - mn2),2)*(ecmPmn2 - M2);
	    return matr/2./mn4;
	} else {
	    matr = TMath::Power((ecm2 - mn2),2)*(ecmMmn2 - M2);
	    return matr/2./mn4;
	}
    } else if (spin == 3) {
	if (parity == 1) {
	    matr = (ecmMmn2 - M*M )*(3.*ecm4 + 6.*ecm3*mn + 4.*ecm2*mn2 +
				     2.*ecm*mn3 + mn4 - 2.*ecm*mn*M2 - 2.*mn2*M2 + M4);
	    return matr/12./ecm2/mn2;
	} else {
	    matr = (ecmPmn2 - M*M )*(3.*ecm4 - 6.*ecm3*mn + 4.*ecm2*mn2 -
				     2.*ecm*mn3 + mn4 + 2.*ecm*mn*M2 - 2.*mn2*M2 + M4);
	    return matr/12./ecm2/mn2;
	}
    } else if (spin == 5) {
	if (parity == 1) {
	    matr = (ecmMmn2 - M*M )*(ecmPmn2 - M*M )*(ecmPmn2 - M*M )*
		(2.*ecm4 - 4.*ecm3*mn + 3.*ecm2*mn2 - 2.*ecm*mn3 + mn4 +
		 2.*ecm*mn*M2 - 2.*mn2*M2 + M4);
	    return matr/480./ecm4/mn4;
        } else {
	    matr = (ecmMmn2 - M*M )*(ecmMmn2 - M*M )*(ecmPmn2 - M*M )*
		(2.*ecm4 + 4.*ecm3*mn + 3.*ecm2*mn2 + 2.*ecm*mn3 + mn4 -
		 2.*ecm*mn*M2 - 2.*mn2*M2 + M4);
	    return matr/480./ecm4/mn4;
        }
    } else {
	cout << "enter resonable spin " << endl;
	return 0;
    }
}

double PResonanceDalitz::GetMatrixL(int &spin, int &parity, int&, const double&, const double& M) {
    float matr;
    if (spin == 1) {
	if (parity == 1) {
	    matr = M2*ecmPmn2*(ecmPmn2 - M2);
	    return matr/2./mn4;
	} else {
	    matr = M2*ecmMmn2*(ecmMmn2 - M2);
	    return matr/2./mn4;
	}
    } else if (spin == 3) {
	if (parity == 1) {
	    matr = M2*(ecmMmn2 - M*M );
	  return matr/3/mn2;
	} else {
	    matr = M2*(ecmPmn2 - M*M );
	    return matr/3./mn2;
	}
    } else if (spin == 5) {
	if (parity == 1) {
	    matr = M2*(ecmMmn2 - M*M )*(ecmPmn2 - M*M)*(ecmPmn2 - M*M);
	    return matr/120/ecm2/mn4;
	} else {
	    matr =M2*(ecmMmn2 - M*M )*(ecmMmn2 - M*M )*(ecmPmn2 - M*M);
	    return matr/120/ecm2/mn4;
	}
    } else {
	cout<<"enter resonable spin "<<endl;
	return 0;
    }
}


double PResonanceDalitz::getLambda(double a, double b, double c) {
    return a*a + b*b + c*c - 2*(a*b + b*c + a*c);
}


double PResonanceDalitz::dGdM(const int&, const double& m, const double& ecm) {
  
    mn      = 0.938273;
    mn2     = mn*mn;
    mn3     = mn*mn*mn;
    mn4     = mn*mn*mn*mn;
    ecm2    = TMath::Power(ecm,2);
    ecm4    = TMath::Power(ecm,4);
    ecm3    = TMath::Power(ecm,3);
    ecmPmn2 = TMath::Power((ecm + mn),2);
    ecmMmn2 = TMath::Power((ecm - mn),2);
    M2      = m*m;
    M4      = m*m*m*m;   
    double lam  = getLambda(ecm*ecm,mn*mn,m*m);
    double alfa = 0.007297;
    double pi   = 3.14159265; 
    
    double dgdm = alfa*alfa*g_Em*g_Em*sqrt(lam)*
	(2.*GetMatrixT(spin,par,charge,ecm,m) +
	 GetMatrixL(spin,par,charge,ecm,m))/6./pi/ecm/ecm/ecm/m;
    
    if (dgdm > 0) return dgdm;
    return 0;
}



ClassImp(PResonanceDalitz)

