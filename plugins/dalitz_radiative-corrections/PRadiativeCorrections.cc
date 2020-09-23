#include "PRadiativeCorrections.h"

#include <iostream>

using namespace std;


ClassImp(PRadiativeCorrections)

PRadiativeCorrections::PRadiativeCorrections() {
}

PRadiativeCorrections::PRadiativeCorrections(const Char_t *id, const Char_t *de, Int_t key) :
    PChannelModel(id, de, key) {

    if (is_channel<0)
        Warning("PRadiativeCorrections","This model should be bound to CHANNELS only");

    parent = nullptr;
    meson = nullptr;
    lp = nullptr;
    lm = nullptr;

    limits_set = false;
    weight_max = 1.;
}

Bool_t PRadiativeCorrections::Init() {
    //Init function called once for each PChannel

    //looking for parent. This is mandatory
    parent = GetParticle("parent");
    if (!parent) {
        Warning("Init","Parent not found");
        return kFALSE;
    }
    bool muon;
    if (parent->Is("dilepton"))
        muon = false;
    else if (parent->Is("dimuon"))
        muon = true;
    else {
        Warning("Init","Channel is not a dilepton/dimuon decay");
        return kFALSE;
    }

    lp = muon ? GetParticle("mu+") : GetParticle("e+");
    if (!lp) {
        Warning("Init","Positron / mu+ not found");
        return kFALSE;
    }
    lm = muon ? GetParticle("mu-") : GetParticle("e-");
    if (!lm) {
        Warning("Init","Electron / mu- not found");
        return kFALSE;
    }

    rad_corrections = GetSecondaryModel("corrections");

    if (!rad_corrections) {
        Warning("Init","Model for radiative corrections not found");
        return kFALSE;
    }

    return kTRUE;
}

Double_t PRadiativeCorrections::ApproximateValue(const double x, const double y)
{
    double y_lower, y_upper;
    auto it = upper_bound(y_vals->begin(), y_vals->end(), y);
    y_lower = *(it-1);
    y_upper = *it;
    // if the iterator is at end() of the vector this means the y value is larger than the largest y value available
    bool y_max = it == y_vals->end() ? true : false;
    // check for the special case that we have an exact y value which is covered in the correction table
    // note: lower_bound and upper_bound might return unintuitive results; if the returned iterators are the same
    // this means the current value which is used sits between values covered in the container / correction table
    bool y_exact = it != lower_bound(y_vals->begin(), y_vals->end(), y);
    // determined lower and upper y values, use those to get the corresponding x value ranges
    // first for the lower y value
    auto it_y_lower_begin = lower_bound(y_vec->begin(), y_vec->end(), y_lower);
    int pos_y_lower_begin = it_y_lower_begin - y_vec->begin();
    auto it_y_lower_end = upper_bound(it_y_lower_begin, y_vec->end(), y_lower);
    int pos_y_lower_end = it_y_lower_end - y_vec->begin();
    // get the lower and upper x value within this range; take into account that end is outside the range of interest and that upper_bound
    // returns a pointer to the position above the value of interest which means if x is below the range it returns begin
    auto it_x_upper_y_low = upper_bound(x_vec->begin()+pos_y_lower_begin, x_vec->begin()+pos_y_lower_end, x);
    double x_lower_y_low = it_x_upper_y_low == (x_vec->begin()+pos_y_lower_begin) ? *it_x_upper_y_low : *(it_x_upper_y_low-1);
    double x_upper_y_low = it_x_upper_y_low == (x_vec->begin()+pos_y_lower_end) ? x_lower_y_low : *it_x_upper_y_low;
    // if the iterator returned by the upper_bound call above is the beginning of the search range means the x value is smaller the the smallest one
    bool x_min = x < x_lower_y_low ? true : false;
    // in case the iterator points to the end of the range
    bool x_max = x > x_upper_y_low ? true : false;
    int pos_x_low_y_low = it_x_upper_y_low - x_vec->begin() - 1;
    if (x_min)  // in case x is smaller than the smallest value in the range, the position is at begin and subtracting one leads to a wrong position
        pos_x_low_y_low++;
    int pos_x_up_y_low = it_x_upper_y_low - x_vec->begin();
    if (x_max)  // if the x is larger than the highest value, the position is at end which is one step above the range of interest
        pos_x_up_y_low--;
    bool xmin_low = x_min;
    bool xmax_low = x_max;
    // check if an exact x value has been calculated which is covered by the correction table
    bool x_exact = it_x_upper_y_low != lower_bound(x_vec->begin()+pos_y_lower_begin, x_vec->begin()+pos_y_lower_end, x);
    // if the y value is the max value, we're at the edge of the correction table
    // just do a linear interpolation between the two found points basically along the correction table edge
    // the same holds true if the calculated y value is an exact value covered in the table, do an interpolation along y in this case as well
    // note: below the min y value is not possible since y is a kinematic variable >= 0 and 0 is inlcuded in the corrections
    if (y_max || y_exact) {
        // there might be the special case that we're outside the range for both x and y, in this case use the closest edge as the approx. value
        if (x_min)
            return z_vec->at(pos_x_up_y_low);
        else if (x_max)
            return z_vec->at(pos_x_low_y_low);
        // there might be the very rare case that both x and y are exactly covered in the correction table, return the exact value in this case
        if (y_exact && x_exact)
            return z_vec->at(pos_x_up_y_low);
        // if the edge cases above do not apply, perform the aforementioned linear interpolation
        return interpolate_linear(x, x_lower_y_low, x_upper_y_low, z_vec->at(pos_x_low_y_low), z_vec->at(pos_x_up_y_low));
    }
    // now do the same for the upper y value
    // y is within the correction table, continue search for x value for the upper y value
    auto it_y_upper_begin = lower_bound(y_vec->begin(), y_vec->end(), y_upper);
    int pos_y_upper_begin = it_y_upper_begin - y_vec->begin();
    auto it_y_upper_end = upper_bound(it_y_upper_begin, y_vec->end(), y_upper);
    int pos_y_upper_end = it_y_upper_end - y_vec->begin();
    // get the lower and upper x value within this range
    auto it_x_upper_y_up = upper_bound(x_vec->begin()+pos_y_upper_begin, x_vec->begin()+pos_y_upper_end, x);
    double x_lower_y_up = it_x_upper_y_up == (x_vec->begin()+pos_y_upper_begin) ? *it_x_upper_y_up : *(it_x_upper_y_up-1);
    double x_upper_y_up = it_x_upper_y_up == (x_vec->begin()+pos_y_upper_end) ? x_lower_y_up : *it_x_upper_y_up;
    bool xmin_up = x < x_lower_y_up ? true : false;
    bool xmax_up = x_upper_y_up < x ? true : false;
    int pos_x_low_y_up = it_x_upper_y_up - x_vec->begin() - 1;
    if (xmin_up)  // in case x is smaller than the smallest value in the range, the position is at begin and subtracting one leads to a wrong position
        pos_x_low_y_up++;
    int pos_x_up_y_up = it_x_upper_y_up - x_vec->begin();
    if (xmax_up)  // if the x is larger than the highest value, the position is at end which is one step above the range of interest
        pos_x_up_y_up--;
    // check if the x value is an exact value which is covered --> no bilinear interpolation needed, just a linear one is sufficient
    if (x_exact ||
            (it_x_upper_y_up != lower_bound(x_vec->begin()+pos_y_upper_begin, x_vec->begin()+pos_y_upper_end, x)))
        return interpolate_linear(y, y_lower, y_upper, z_vec->at(pos_x_low_y_low), z_vec->at(pos_x_low_y_up));
    // check if xmin or xmax has been triggered, only interpolate y along the edge of the correction table
    x_min |= xmin_up;
    x_max |= xmax_up;
    if (x_min || x_max)
        return interpolate_linear(y, y_lower, y_upper, z_vec->at(pos_x_low_y_low), z_vec->at(pos_x_low_y_up));
    // all special cases (outside correction table, at the edge, exact values) where only a linear interpolation is needed are checked now
    // at this point we have only the case remaining that both the x and y values are covered and sitting between precalculated corrections
    // do a bilinear interpolation between those four surrounding correction values
    return interpolate_bilinear(x, y, x_lower_y_up, x_upper_y_up, y_lower, y_upper,
                                z_vec->at(pos_x_low_y_low), z_vec->at(pos_x_up_y_low),
                                z_vec->at(pos_x_low_y_up), z_vec->at(pos_x_up_y_up));
}

double PRadiativeCorrections::interpolate_linear(const double x, const double x0, const double x1,
                                                 const double y0, const double y1)
{
    // y = y_{0}\left(1-{\frac {x-x_{0}}{x_{1}-x_{0}}}\right) + y_{1}\left({\frac {x-x_{0}}{x_{1}-x_{0}}}\right)

    double w = (x-x0)/(x1-x0);
    return y0*(1-w) + y1*w;
}

double PRadiativeCorrections::interpolate_bilinear(const double x, const double y,
                                                   const double x0, const double x1,
                                                   const double y0, const double y1,
                                                   const double val_x0_y0, const double val_x1_y0,
                                                   const double val_x0_y1, const double val_x1_y1)
{
    // f(x,y) = \frac{1}{(x_{2}-x_{1})(y_{2}-y_{1})} \left(f(x_{1},y_{1})(x_{2}-x)(y_{2}-y) + f(x_{2},y_{1})(x-x_{1})(y_{2}-y) + f(x_{1},y_{2})(x_{2}-x)(y-y_{1}) + f(x_{2},y_{2})(x-x_{1})(y-y_{1})\right)

    return 1./((x1-x0)*(y1-y0))
            * (val_x0_y0*(x1-x)*(y1-y) + val_x1_y0*(x-x0)*(y1-y)
               + val_x0_y1*(x1-x)*(y-y0) + val_x1_y1*(x-x0)*(y-y0));
}

Bool_t PRadiativeCorrections::IsValid() {
    //Use rejection mode...

    if (GetVersionFlag() & VERSION_WEIGHTING) return kTRUE;
    //...but not if weighting enabled.

    meson = parent->GetParent();
    double q2 = parent->M2();
    double im2 = meson->M2();
    double x = q2/im2;
    double y = meson->Vect4().Dot(lp->Vect4()-lm->Vect4());
    y = 2*abs(y)/meson->M2()/(1-x);
    double nu2 = 4*lm->M2()/im2;
    double beta = sqrt(1-nu2/x);
    if ((x < nu2) || (x > 1.))
        return kFALSE;
    if ((y < 0.) || (y > beta))
        return kFALSE;

    // valid x and y values, do rejection now
    double weight = GetWeight();

    if ((weight/weight_max) > PUtils::sampleFlat())
        return kTRUE;

    return kFALSE;
}
