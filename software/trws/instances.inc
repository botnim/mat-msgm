#ifndef __INSTANCES_H_
#define __INSTANCES_H_

#include "instances.h"
#include "MRFEnergy.h"

template class MRFEnergy<TypeBinary>;
template class MRFEnergy<TypeBinaryFast>;
template class MRFEnergy<TypePotts>;
template class MRFEnergy<TypeGeneral>;
template class MRFEnergy<typeGeneralCompact>;
template class MRFEnergy<TypeTruncatedLinear>;
template class MRFEnergy<TypeTruncatedQuadratic>;
template class MRFEnergy<TypeTruncatedLinear2D>;
template class MRFEnergy<TypeTruncatedQuadratic2D>;


#endif // __INSTANCES_H_
