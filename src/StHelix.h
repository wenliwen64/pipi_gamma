#ifndef ST_HELIX_H
#define ST_HELIX_H
#include "Rtypes.h"
#include "StThreeVector.h"

/*
#ifdef __ROOT__
#include "Rtypes.h"
#endif
#ifndef __CINT__
#include <iostream>
#include <math.h>
#include "StThreeVector.h"
#ifdef GNU_GCC
#    include <stddef.h>
#endif
#if defined (__SUNPRO_CC) && __SUNPRO_CC < 0x500
#    include <stdcomp.h>
#endif
#ifndef ST_NO_EXCEPTIONS
#    include <stdexcept>
#    if !defined(ST_NO_NAMESPACES)
using std::out_of_range;
#    endif
#endif
#endif // __CINT__
*/

class StHelix {
public:
    StHelix();
    StHelix(double c, double d, double phase, const StThreeVector<double>& o, int h);
    virtual ~StHelix();

    void setParameters(double c, double dip, double phase, const StThreeVector<double>& o, int h);
    double       dipAngle()   const;

protected:
    void setCurvature(double);  /// performs also various checks   
    void setPhase(double);              
    void setDipAngle(double);

protected:
    bool                   mSingularity;        // true for straight line case (B=0)
    StThreeVector<double>  mOrigin;
    double                 mDipAngle;
    double                 mCurvature;
    double                 mPhase;
    int                    mH;                  // -sign(q*B);

    double                 mCosDipAngle;
    double                 mSinDipAngle;
    double                 mCosPhase;
    double                 mSinPhase;

//#ifdef __ROOT__
  ClassDef(StHelix,1)
//#endif /* __ROOT__ */
};



inline double StHelix::dipAngle() const {return mDipAngle;}

inline void StHelix::setDipAngle(double val)
{
    mDipAngle    = val;
    mCosDipAngle = cos(mDipAngle);
    mSinDipAngle = sin(mDipAngle);
}












#endif
