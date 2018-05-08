#include "../src/PChannelModel.h"
#include "TRandom2.h"

//Class definition

class PEtaDoubleDalitz : public PChannelModel  {
  
public:
  
  PEtaDoubleDalitz(Char_t *id, Char_t *de, Int_t key);
  PDistribution* Clone(const char*delme=NULL) const;
  
  Bool_t Init(void);
  
  using PChannelModel::SampleMass;
  Bool_t SampleMass(void);
  
private:
  
  Double_t Gen2lepton1(Double_t m);
  TRandom2 *gRand;
  PParticle *dil1,*dil2,*parent;
  
  ClassDef(PEtaDoubleDalitz,0)  //Just a dummy model for mass  sampling, NO PHYSICS!
};

PDistribution* PEtaDoubleDalitz::Clone(const char*) const {
  //clone the object
  return new PEtaDoubleDalitz((const PEtaDoubleDalitz &)* this);
};

PEtaDoubleDalitz::PEtaDoubleDalitz(Char_t *id, Char_t *de, Int_t key) : PChannelModel(id, de,key) {
  //Constructor
  gRand = new TRandom2();
  dil1=dil2=parent=NULL;
} ;

Double_t PEtaDoubleDalitz::Gen2lepton1(Double_t m)
{
    Double_t MEl=PData::LMass("e+");
    Double_t MEl2=2*MEl;
    Double_t m2=m*m;
    Double_t ARAND[2];
    Double_t WTD;
    Double_t AM;
    do{
        gRand->RndmArray(2,ARAND);
        AM=MEl2*TMath::Exp(TMath::Log(m/MEl2)*ARAND[0]);
        Double_t AM2=AM*AM;
        Double_t RE=MEl2*MEl2/AM2;
        WTD =TMath::Sqrt(1-RE)*(1+RE/2)*TMath::Power(1-AM2/m2,3);
    }while(ARAND[1]>WTD);// ???????????? while(ARAND[1]>WTD)
    return AM;
    cout << AM << endl ;
    cout << WTD << endl ;
}


Bool_t PEtaDoubleDalitz::Init(void) {
  
  parent = GetParticle("parent");
  if (!parent) {
    Warning("Init","Parent not found");
    return kFALSE;
  }

  dil1 = GetParticle("dilepton");
  dil2 = GetParticle("dilepton");
  
  if (!dil1 || !dil2) {
    Warning("Init","Dileptons not found");
    return kFALSE;
  }
  
  return kTRUE;
}

Bool_t PEtaDoubleDalitz::SampleMass(void) {

  Double_t MEta = parent->M();
  
  
  Double_t mVV[2];
  do {
    mVV[0]=Gen2lepton1(MEta);//IM squared of the 1-st pair e+e-
    mVV[1]=Gen2lepton1(MEta);//IM squared of the 2-nd pair e+e-
  } while(mVV[0]+mVV[1]>MEta);

  
  dil1->SetM(mVV[0]);
  dil2->SetM(mVV[1]);
  
  //do whatever you want here!!!!!!!!!!!!!
  //dil1->SetM(0.121323);
  //dil2->SetM(0.243424);
  
  return kTRUE;
};

ClassImp(PEtaDoubleDalitz)
