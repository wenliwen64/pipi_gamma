/***************************************************************************
 *
 * $Id: StLorentzVectorD.hh,v 1.6 2005/07/06 18:49:56 fisyak Exp $
 * $Log: StLorentzVectorD.hh,v $
 * Revision 1.6  2005/07/06 18:49:56  fisyak
 * Replace StHelixD, StLorentzVectorD,StLorentzVectorF,StMatrixD,StMatrixF,StPhysicalHelixD,StThreeVectorD,StThreeVectorF by templated version
 *

****************************************************************************/
#ifndef ST_LORENTZ_VECTOR_D_H
#define ST_LORENTZ_VECTOR_D_H
#include "StThreeVectorF.h"
#include "StThreeVectorD.h"
#include "StLorentzVector.h"
typedef StLorentzVector<double> StLorentzVectorD;
#endif
