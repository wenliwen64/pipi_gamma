/***************************************************************************
 *
 * $Id: StLorentzVectorF.hh,v 1.6 2005/07/06 18:49:56 fisyak Exp $
 * $Log: StLorentzVectorF.hh,v $
 * Revision 1.6  2005/07/06 18:49:56  fisyak
 * Replace StHelixD, StLorentzVectorD,StLorentzVectorF,StMatrixD,StMatrixF,StPhysicalHelixD,StThreeVectorD,StThreeVectorF by templated version
 *

****************************************************************************/
#ifndef ST_LORENTZ_VECTOR_F_H
#define ST_LORENTZ_VECTOR_F_H
#include "StThreeVectorF.h"
#include "StThreeVectorD.h"
#include "StLorentzVector.h"
typedef StLorentzVector<float> StLorentzVectorF;
#endif
