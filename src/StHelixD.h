/***************************************************************************
 *
 * $Id: StHelixD.hh,v 1.8 2005/07/06 18:49:56 fisyak Exp $
 * $Log: StHelixD.hh,v $
 * Revision 1.8  2005/07/06 18:49:56  fisyak
 * Replace StHelixD, StLorentzVectorD,StLorentzVectorF,StMatrixD,StMatrixF,StPhysicalHelixD,StThreeVectorD,StThreeVectorF by templated version
 *

****************************************************************************/
#ifndef ST_HELIX_D_H
#define ST_HELIX_D_H
#include "StThreeVectorD.h"
#include "StHelix.h"
#include <utility>
typedef StHelix StHelixD;
typedef pair<double,double> pairD;
#endif
