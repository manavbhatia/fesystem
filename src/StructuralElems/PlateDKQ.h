// $Id: PlateDKQ.h,v 1.1 2006-09-05 21:11:44 manav Exp $

#ifndef __plate_dkq_h__
#define __plate_dkq_h__

// C++ includes

// FESystem includes
#include "structural_elem.h"



/**
 *	a class definig the spring structural element
 */

class PlateDKQ: public StructuralElem
{
 public:
  PlateDKQ(const unsigned int, const Elem* , AnalysisDriver* );
	
  ~PlateDKQ();
	

  /// this function returns a vector of the post process quanties for this element. The quantities will 
  /// consist of the strain, stress and their sensitivities. 
  /// @param vector of load case IDs for which the post process quantities need to be calculated
  /// @param vector of DVs for which the quantities need to be calculated
  virtual std::auto_ptr<ElemPostProcessQty> 
    getElementPostProcessQty(std::vector<unsigned int> , std::vector<DesignData::DesignParameter*>);
  	   
	   
 protected:
		   
  virtual bool calculate_M(DenseMatrix<double>*, 
			   bool sensitivity_calculation = false);
	   
  virtual bool calculate_K(DenseMatrix<double>*, 
			   bool sensitivity_calculation = false);
	   
	   
  virtual bool calculate_F_T(DenseVector<double>*, 
			     bool sensitivity_calculation = false);

  void calculateShearStiffnessFactor(DenseMatrix<double>* matrix);


  /// this is an  method to calculate the strain operator matrices
  /// @param the vector where all the strain operators will be returned
  /// @param vector of points, defined in the element local coordinate system where the operators will be 
  /// calculated
  //virtual void calculateStrainOperator(std::vector<DenseMatrix<double> >* , const std::vector<Point>* );

};


#endif // __plate_dkq_h__
