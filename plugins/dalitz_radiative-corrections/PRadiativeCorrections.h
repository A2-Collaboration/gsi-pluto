#ifndef _PRADIATIVECORRECTIONS_H_
#define _PRADIATIVECORRECTIONS_H_

#include <algorithm>

#include "PChannelModel.h"

#include "TGraph2D.h"
#include "TH2D.h"

class PRadiativeCorrections : public PChannelModel {

 public:
    PRadiativeCorrections();
    PRadiativeCorrections(const Char_t *id, const Char_t *de, Int_t key);

    Bool_t Init();

    using PDistribution::GetWeight;
    using PChannelModel::GetWeight;

    virtual Double_t GetWeight() = 0;

    Bool_t IsValid();

 protected:
    PParticle *parent, *meson, *lp, *lm;

    PChannelModel *rad_corrections;

    virtual Double_t ApproximateValue(const double, const double);
    virtual void SetLimits() = 0;
    bool limits_set;

    Double_t weight_max;
    // vectors containing data points to determine the limits of the correction table
    std::vector<double>* y_vals = nullptr;
    std::vector<double>* x_tuples = nullptr;
    // pointer to the x and y values as well as the corrections
    // used as a fallback to approximate the correction if the delauny interpolation fails
    std::vector<double>* x_vec = nullptr;
    std::vector<double>* y_vec = nullptr;
    std::vector<double>* z_vec = nullptr;

  private:
    double interpolate_linear(const double value, const double x_low, const double x_up, const double y_low, const double y_up);
    double interpolate_bilinear(const double x, const double y, const double x_low, const double x_up, const double y_low, const double y_up,
                                const double val_x_low_y_low, const double val_x_up_y_low, const double val_x_low_y_up, const double val_x_up_y_up);


    ClassDef(PRadiativeCorrections,0)
};

#endif
