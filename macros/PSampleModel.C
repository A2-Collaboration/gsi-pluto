//TITLE A dummy model for multi-dimensional sampling

#include "../src/PChannelModel.h"

//Class definition

class PSampleModel : public PChannelModel  {
  
 public:

    using PDistribution::GetWeight;
    PSampleModel(Char_t *id, Char_t *de, Int_t key);
    PDistribution* Clone(const char*delme=NULL) const;
    Double_t GetWeight(Double_t *mass, Int_t *didx=NULL);

 private:

    ClassDef(PSampleModel,0)  //Just a dummy model for mass sampling, NO PHYSICS!
};

PDistribution* PSampleModel::Clone(const char*) const {
    //clone the object
    return new PSampleModel((const PSampleModel &)* this);
};

PSampleModel::PSampleModel(Char_t *id, Char_t *de, Int_t key) : PChannelModel(id, de,key) {
    //Constructor
} ;

Double_t PSampleModel::GetWeight(Double_t *mass, Int_t *didx) {
    //cout << "GetWeight" << mass[0] << ":" << mass[1]<<endl;
    //return pow(sin(mass[0]*10+mass[1]*mass[1]*10)*mass[0]*mass[1],2);
    // return 1/((mass[0]-0.5)+0.00000001 + (mass[1]*mass[1]-.5)); --> sehr langsam
    //return 1/((mass[0]*mass[0]*mass[0]*mass[0])+0.00000001 + (mass[1]*mass[1]*mass[1]*mass[1])); //--> sauschnell

    //if ((fabs(mass[0]-0.5)<0.001) && (fabs(mass[1]-0.5)<0.001)) return 1.;
    if ((fabs(mass[0]-0.1)<0.001) && (fabs(mass[1]-0.1)<0.001)) return 1.;
    


    //if ((mass[0]<0.5) || (mass[1]<0.5)) return 1.;
    return 0.;
					    

}


ClassImp(PSampleModel)


