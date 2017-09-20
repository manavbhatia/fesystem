// $Id: InterpolationHex8.h,v 1.7 2006-10-29 05:08:27 manav Exp $

#ifndef __fesystem_interpolation_hex8_h__
#define __fesystem_interpolation_hex8_h__

// C++ includes


// FESystem includes
#include "Interpolation/FEInterpolationElem.h"

// libMesh includes

/**
*	a 1-D conduction element
 */


#ifndef FE_INTERPOLATION_HEX8_ENUM_ID
#define FE_INTERPOLATION_HEX8_ENUM_ID 14
#else
#error
#endif


#ifndef FE_INTERPOLATION_HEX8_ENUM_NAME
#define FE_INTERPOLATION_HEX8_ENUM_NAME "FE_INTERPOLATION_HEX8"
#else
#error
#endif


DeclareFESystemElemTypeEnumName(FE_INTERPOLATION_HEX8,
                                FE_INTERPOLATION_HEX8_ENUM_ID,
                                FE_INTERPOLATION_HEX8_ENUM_NAME,
                                HEX8)


class InterpolationHex8: public FEInterpolationElem
{
public:
	InterpolationHex8(Discipline::AnalysisDisciplineBase& discipline);

	~InterpolationHex8();
	
private:
	
};


#endif //__fesystem_interpolation_hex8_h__
