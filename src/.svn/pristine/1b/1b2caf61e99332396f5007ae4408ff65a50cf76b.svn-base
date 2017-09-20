// $Id: BendingElemPostProcessQty.h,v 1.2 2006-09-05 20:41:50 manav Exp $

#ifndef __bending_elem_post_process_qty_h__
#define __bending_elem_post_process_qty_h__ 


// C++ includes



// FESystem includes
#include "ElemPostProcessQty.h"


// libMesh includes


class BendingElemPostProcessQty : public ElemPostProcessQty
{
public:
	BendingElemPostProcessQty();
	
	~BendingElemPostProcessQty();

	
	/// returns the location of the quantity, +ve indicates above the mid layer
	unsigned int getLayer() const;
	
protected:
	
	/// the layer of the element, +1 indicates upper and -1 indicates lower.
	unsigned int layer;
};


#endif // __bending_elem_post_process_qty_h__