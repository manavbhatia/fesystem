// $Id: thermal_elem.h,v 1.23.4.2 2007-05-14 16:45:52 manav Exp $

#ifndef __fesystem_thermal_elem_h__
#define __fesystem_thermal_elem_h__

// C++ includes
#include <memory>

// FESystem includes
#include "FESystem/FESystemElem.h"
#include "Utilities/NameEnumHandler.h"
#include "Properties/ElemDataCard.h"

// libMesh includes
#include "numerics/dense_matrix.h"
#include "numerics/dense_vector.h"


// Forward declerations
class Load;
namespace Discipline
{
  class AnalysisDisciplineBase;
}


#ifndef THERMAL_C_MATRIX_ENUM_ID
#define THERMAL_C_MATRIX_ENUM_ID 1
#else
#error
#endif

#ifndef THERMAL_C_MATRIX_ENUM_NAME
#define THERMAL_C_MATRIX_ENUM_NAME "THERMAL_C_MATRIX"
#else
#error
#endif


#ifndef THERMAL_K_C_MATRIX_ENUM_ID
#define THERMAL_K_C_MATRIX_ENUM_ID 2
#else
#error
#endif

#ifndef THERMAL_K_C_MATRIX_ENUM_NAME
#define THERMAL_K_C_MATRIX_ENUM_NAME "THERMAL_K_C_MATRIX"
#else
#error
#endif

#ifndef THERMAL_K_H_MATRIX_ENUM_ID
#define THERMAL_K_H_MATRIX_ENUM_ID 3
#else
#error
#endif

#ifndef THERMAL_K_H_MATRIX_ENUM_NAME
#define THERMAL_K_H_MATRIX_ENUM_NAME "THERMAL_K_H_MATRIX"
#else
#error
#endif


#ifndef THERMAL_F_H_VECTOR_ENUM_ID
#define THERMAL_F_H_VECTOR_ENUM_ID 4
#else
#error
#endif

#ifndef THERMAL_F_H_VECTOR_ENUM_NAME
#define THERMAL_F_H_VECTOR_ENUM_NAME "THERMAL_F_H_VECTOR"
#else
#error
#endif

#ifndef THERMAL_F_VOL_VECTOR_ENUM_ID
#define THERMAL_F_VOL_VECTOR_ENUM_ID 5
#else
#error
#endif

#ifndef THERMAL_F_VOL_VECTOR_ENUM_NAME
#define THERMAL_F_VOL_VECTOR_ENUM_NAME "THERMAL_F_VOL_VECTOR"
#else
#error
#endif

#ifndef THERMAL_F_Q_SURF_VECTOR_ENUM_ID
#define THERMAL_F_Q_SURF_VECTOR_ENUM_ID 6
#else
#error
#endif

#ifndef THERMAL_F_Q_SURF_VECTOR_ENUM_NAME
#define THERMAL_F_Q_SURF_VECTOR_ENUM_NAME "THERMAL_F_Q_SURF_VECTOR"
#else
#error
#endif

#ifndef THERMAL_F_EMITTED_RAD_VECTOR_ENUM_ID
#define THERMAL_F_EMITTED_RAD_VECTOR_ENUM_ID 7
#else
#error
#endif

#ifndef THERMAL_F_EMITTED_RAD_VECTOR_ENUM_NAME
#define THERMAL_F_EMITTED_RAD_VECTOR_ENUM_NAME "THERMAL_F_EMITTED_RAD_VECTOR"
#else
#error
#endif


#ifndef THERMAL_EMITTED_RAD_JAC_MATRIX_ENUM_ID
#define THERMAL_EMITTED_RAD_JAC_MATRIX_ENUM_ID 8
#else
#error
#endif

#ifndef THERMAL_EMITTED_RAD_JAC_MATRIX_ENUM_NAME
#define THERMAL_EMITTED_RAD_JAC_MATRIX_ENUM_NAME "THERMAL_EMITTED_RAD_JAC_MATRIX"
#else
#error
#endif



#ifndef THERMAL_F_SIGMA_FACTOR_ENUM_ID
#define THERMAL_F_SIGMA_FACTOR_ENUM_ID 9
#else
#error
#endif

#ifndef THERMAL_F_SIGMA_FACTOR_ENUM_NAME
#define THERMAL_F_SIGMA_FACTOR_ENUM_NAME "THERMAL_F_SIGMA_FACTOR"
#else
#error
#endif


#ifndef THERMAL_K_C_JAC_MATRIX_ENUM_ID
#define THERMAL_K_C_JAC_MATRIX_ENUM_ID 10
#else
#error
#endif

#ifndef THERMAL_K_C_JAC_MATRIX_ENUM_NAME
#define THERMAL_K_C_JAC_MATRIX_ENUM_NAME "THERMAL_K_C_JAC_MATRIX"
#else
#error
#endif


#ifndef THERMAL_C_JAC_MATRIX_ENUM_ID
#define THERMAL_C_JAC_MATRIX_ENUM_ID 11
#else
#error
#endif

#ifndef THERMAL_C_JAC_MATRIX_ENUM_NAME
#define THERMAL_C_JAC_MATRIX_ENUM_NAME "THERMAL_C_JAC_MATRIX"
#else
#error
#endif



namespace FESystemElem
{
  
  DeclareEnumClass(ThermalElemQtyEnum);
  
  DeclareEnumName(THERMAL_C_MATRIX, FESystemElem::ThermalElemQtyEnum,
                  THERMAL_C_MATRIX_ENUM_ID, THERMAL_C_MATRIX_ENUM_NAME);
  
  DeclareEnumName(THERMAL_K_C_MATRIX, FESystemElem::ThermalElemQtyEnum,
                  THERMAL_K_C_MATRIX_ENUM_ID, THERMAL_K_C_MATRIX_ENUM_NAME);
  
  DeclareEnumName(THERMAL_K_H_MATRIX, FESystemElem::ThermalElemQtyEnum,
                  THERMAL_K_H_MATRIX_ENUM_ID, THERMAL_K_H_MATRIX_ENUM_NAME);
  
  DeclareEnumName(THERMAL_F_H_VECTOR, FESystemElem::ThermalElemQtyEnum,
                  THERMAL_F_H_VECTOR_ENUM_ID, THERMAL_F_H_VECTOR_ENUM_NAME);
  
  DeclareEnumName(THERMAL_F_VOL_VECTOR, FESystemElem::ThermalElemQtyEnum,
                  THERMAL_F_VOL_VECTOR_ENUM_ID, THERMAL_F_VOL_VECTOR_ENUM_NAME);
  
  DeclareEnumName(THERMAL_F_Q_SURF_VECTOR, FESystemElem::ThermalElemQtyEnum,
                  THERMAL_F_Q_SURF_VECTOR_ENUM_ID, THERMAL_F_Q_SURF_VECTOR_ENUM_NAME);
  
  DeclareEnumName(THERMAL_F_EMITTED_RAD_VECTOR, FESystemElem::ThermalElemQtyEnum,
                  THERMAL_F_EMITTED_RAD_VECTOR_ENUM_ID, THERMAL_F_EMITTED_RAD_VECTOR_ENUM_NAME);
  
  DeclareEnumName(THERMAL_EMITTED_RAD_JAC_MATRIX, FESystemElem::ThermalElemQtyEnum,
                  THERMAL_EMITTED_RAD_JAC_MATRIX_ENUM_ID, THERMAL_EMITTED_RAD_JAC_MATRIX_ENUM_NAME);
  
  DeclareEnumName(THERMAL_F_SIGMA_FACTOR, FESystemElem::ThermalElemQtyEnum,
                  THERMAL_F_SIGMA_FACTOR_ENUM_ID, THERMAL_F_SIGMA_FACTOR_ENUM_NAME);
  
  DeclareEnumName(THERMAL_K_C_JAC_MATRIX, FESystemElem::ThermalElemQtyEnum,
                  THERMAL_K_C_JAC_MATRIX_ENUM_ID, THERMAL_K_C_JAC_MATRIX_ENUM_NAME);

  DeclareEnumName(THERMAL_C_JAC_MATRIX, FESystemElem::ThermalElemQtyEnum,
                  THERMAL_C_JAC_MATRIX_ENUM_ID, THERMAL_C_JAC_MATRIX_ENUM_NAME);

  /**
    *	a base class defining the heat conduction elements
   */
  
  class ThermalElem: public FESystemElem::FESystemElemBase
    {
public:
      /// constructor
      ThermalElem(const unsigned int dim, 
                  const unsigned int elem_enum_ID,
                  Discipline::AnalysisDisciplineBase& discipline);
      
      /// @returns the transient system order for the elem equations
      virtual unsigned int getTransientSystemOrder() const;

      /// destructor
      virtual ~ThermalElem();
      
      
      /// returns the assembled vector quantity for the element
      virtual void getElementAssembledQty(const unsigned int name_enum_ID,
                                          DenseVector<double>* qty); 
      
      /// returns the assembled matrix quantity for the element
      virtual void getElementAssembledQty(const unsigned int name_enum_ID,
                                          DenseMatrix<double>* qty);
      
      /// returns the assembled matrix quantity sensitivity for the element
      virtual void  getElementAssembledQtySensitivity(const unsigned int name_enum_ID,
                                                      DenseVector<double>* qty);
      
      /// returns the assembled matrix quantity sensitivity for the element
      virtual void  getElementAssembledQtySensitivity(const unsigned int name_enum_ID,
                                                      DenseMatrix<double>* qty);
      
      
      /// this function returns the element post process quantity
      /// @param vector of load case IDs for which the quantitis should be calculated
      /// @param vector of DVs for which the quantities should be calculated
      virtual std::auto_ptr<ElemPostProcessQty>
        getElementPostProcessQty(std::vector<unsigned int> load_cases,
                                 std::vector<DesignData::DesignParameter*> DV_vector);
      
      
      
protected:
        
        /// inserts the finite element types to be used for this element in the input 
        /// parameters
        virtual void getFETypes(std::vector<FEType>& fetypes);
      
      
      /// inserts the quadrature rules to be used for the different elements
      virtual void getQuadratureRules
        (std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures);
      
      
      /// function to direct the call to appropriate calculation routine based on the 
      /// quantity seeked. The arguements are:
      /// @param  pointer of quantity to return the data in
      /// @param quantity name
      /// @param bool specifying if sensitivity has been asked
      void calculateAssembledQty(DenseMatrix<double>* ,
                                 const unsigned int qty_enum_ID,
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation);
      
      
      /// function to direct the call to appropriate calculation routine based on the 
      /// quantity seeked. The arguements are:
      /// @param  pointer of quantity to return the data in
      /// @param quantity name
      /// @param bool specifying if sensitivity has been asked
      void calculateAssembledQty(DenseVector<double>* ,
                                 const unsigned int qty_enum_ID,
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation);
      
      
      /// this method returns the quantity specified in the arguement list. This has 
      /// been templatized so that the template can be used in derived element types.
      /// This element class will provide its element quantity and calculation routines.
      /// the sensitivity calculated from this is only shape sensitivity, and not property 
      /// sensitivity
      /// @param name  name of the quantity to be returned
      /// @param qty data structure in which the quantity will be removed
      /// @param side_num number of the side if this is a side quantity
      template <class QtyType >
        void getThermalElemQty(const unsigned int name,
                               QtyType* qty,
                               const unsigned int design_point_enum_ID,
                               const unsigned int domain,
                               const bool sensitivity);
      
      
      
      
      /// this method is an interface to the specific functions to calculate element
      /// quantities. It calls the method that performs the calculations and returns the
      /// quantity in the first arguement. 
      void calculateThermalElemQty(DenseVector<double>* qty, 
                                   const unsigned int qty_enum_ID, 
                                   const unsigned int design_point_enum_ID,
                                   const unsigned int domain_enum_ID);
      
      
      /// this method is an interface to the specific functions to calculate element
      /// quantities. It calls the method that performs the calculations and returns the
      /// quantity in the first arguement. 
      void calculateThermalElemQty(DenseMatrix<double>* qty, 
                                   const unsigned int qty_enum_ID, 
                                   const unsigned int design_point_enum_ID,
                                   const unsigned int domain_enum_ID);
      
      
      /**
        *	function to calculate the element capacitance matrix. The first arguement is the matrix
       *	that will be resized to the appropriate size, and filled with the matrix elements. If 
       *	sensitivity of the matrix wrt to a property is to be calculated, then second arguement 
       *	with the name of the property should be provided. 
       */
      virtual void calculate_C(DenseMatrix<double>* matrix, 
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation);
      
      virtual void calculate_C_Jac(DenseMatrix<double>* matrix, 
                                   const unsigned int design_point_enum_ID);
      /**
        *	function to calculate the element conductance matrix. The first arguement is the matrix
       *	that will be resized to the appropriate size, and filled with the matrix elements. If 
       *	sensitivity of the matrix wrt to a property is to be calculated, then second arguement 
       *	with the name of the property should be provided. 
       */
      virtual void calculate_K_c(DenseMatrix<double>* matrix, 
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation);
      
      /**
        *	function to calculate the conductance of element conductance matrix. 
       * The first arguement is the matrix
       *	that will be resized to the appropriate size, and filled with the matrix elements. 
       * This quantity exists only when the material properties are temperature dependent. 
       * Otherwise, it is zero
       */
      virtual void calculate_K_c_Jac(DenseMatrix<double>* matrix,
                                     const unsigned int design_point_enum_ID);
      
      
      /**
        *	function to calculate the element capacitance matrix contribution due to 
       *	surface convection load. The first arguement is the matrix
       *	that will be resized to the appropriate size, and filled with the matrix elements. 
       *  The second arguement should be the vector of surface convection loads for the element.
       *	If sensitivity of the matrix wrt to a property is to be calculated, then third arguement 
       *	with the name of the property should be provided. 
       */
      virtual void calculate_K_h(DenseMatrix<double>* matrix, 
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation);
      
      /**
        *	function to calculate the element surface  heat flux load vector.
       *	The first arguement is the vector
       *	that will be resized to the appropriate size, and filled with the vector elements. 
       *  The second arguement should be the vector of surface flux loads for the element.
       *	If sensitivity of the vector wrt to a property is to be calculated, then third arguement 
       *	with the name of the property should be provided. 
       */
      virtual void calculate_F_qsurf(DenseVector<double>* vector, 
                                     const unsigned int design_point_enum_ID,
                                     bool sensitivity_calculation);
      
      /**
        *	function to calculate the element volume heat load vector.
       *	 The first arguement is the vector
       *	that will be resized to the appropriate size, and filled with the vector elements. 
       *  The second arguement should be the vector of volume heat loads for the element.
       *	If sensitivity of the vector wrt to a property is to be calculated, then third arguement 
       *	with the name of the property should be provided. 
       */
      virtual void calculate_F_Qvol(DenseVector<double>* vector, 
                                    const unsigned int design_point_enum_ID,
                                    bool sensitivity_calculation);
      
      
      /**
        *	function to calculate the element surface convection load vector.
       *  The first arguement is the vector
       *	that will be resized to the appropriate size, and filled with the vector elements. 
       *  The second arguement should be the vector of surface convection loads for the element.
       *	If sensitivity of the vector wrt to a property is to be calculated, then third arguement 
       *	with the name of the property should be provided. 
       */
      virtual void calculate_F_h(DenseVector<double>* vector, 
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation);
      
      
      /**
        *	function to calculate the element emitted radiation load vector.
       *	The first arguement is the vector
       *	that will be resized to the appropriate size, and filled with the vector elements. 
       *  The third arguement should be the vector of emitted radiation load vector for the element.
       *	If sensitivity of the vector wrt to a property is to be calculated, then third arguement 
       *	with the name of the property should be provided. 
       */
      virtual void calculate_F_sigma(DenseVector<double>* vector, 
                                     const unsigned int design_point_enum_ID,
                                     bool sensitivity_calculation);
      
      /// calculated and returns the jacobian contribution from the radiation loads
      virtual void calculate_F_sigma_Jac(DenseMatrix<double>* matrix,
                                         const unsigned int design_point_enum_ID);
      
      
      
      
      /// this calculates the factor for emitted radiation factor
      void calculateEmittedRadiationLoadFactor(DenseVector<double>* load_factor, 
                                               const unsigned int design_point_enum_ID,
                                               const unsigned int domain);
      
      
      /// initializes the property card at element local temperature, 
      /// which is taken as the average of all nodal temperatures. This is needed for temperature 
      /// dependent properties. This is an approximate way of using parameter dependent 
      /// properties, since ideally they should be calculated at each quadrature point. 
      /// That is, however, not done due to memory and CPU constraints.
      virtual void initPropertyCard();
      
      /// this clears the temperature initialization only for the property cards.
      virtual void clearPropertyCardInitialization();
    };
}




inline
unsigned int FESystemElem::ThermalElem::getTransientSystemOrder() const
{
  return 1;
}


inline
void 
FESystemElem::ThermalElem::getElementAssembledQty(const unsigned int name_enum_ID,
                                                  DenseVector<double>* qty)
{
  Assert(qty != NULL, ExcInternalError());
  qty->zero();
  
  // get the quantity at the base design point
  this->calculateAssembledQty(qty, name_enum_ID, FESystemElem::BASE_ELEM::num(), false);
}


inline
void
FESystemElem::ThermalElem::getElementAssembledQty(const unsigned int name_enum_ID,
                                                  DenseMatrix<double>* qty)
{
  Assert(qty != NULL, ExcInternalError());
  qty->zero();

  // get the quantity
  this->calculateAssembledQty(qty, name_enum_ID, FESystemElem::BASE_ELEM::num(), false);
}



inline
void
FESystemElem::ThermalElem::getElementAssembledQtySensitivity(const unsigned int name_enum_ID,
                                                             DenseVector<double>* qty)
{
  assert (qty != NULL);
  qty->zero();
  
  this->calculateAssembledQty(qty, name_enum_ID, FESystemElem::BASE_ELEM::num(),
                              true);
}



inline
void
FESystemElem::ThermalElem::getElementAssembledQtySensitivity(const unsigned int name_enum_ID,
                                                             DenseMatrix<double>* qty)
{
  assert (qty != NULL);
  qty->zero();
  
  this->calculateAssembledQty(qty, name_enum_ID, FESystemElem::BASE_ELEM::num(),
                              true);
}





template <class QtyType >
void 
FESystemElem::ThermalElem::getThermalElemQty(const unsigned int name,
                                             QtyType* qty,
                                             const unsigned int design_point_enum_ID,
                                             const unsigned int domain,
                                             const bool sensitivity)
{
  // the quantities are not stored in the database.
  assert (qty != NULL);
  qty->zero();
  
  if (!sensitivity)
    this->calculateThermalElemQty(qty, name, design_point_enum_ID, domain);
  else
    {
    // this is only for shape sensitivity, since the local quantities specific to this 
    // element are assumed not to be property dependent
    Assert (this->sensitivity_parameter == DesignData::SHAPE_PARAMETER::num(), 
            ExcInternalError());

    static QtyType  base_qty;
    
    // next, obtain the quantity for the perturbed element
    this->calculateThermalElemQty(qty, name,
                                  FESystemElem::BASE_PLUS_DELTA_ELEM::num(), domain);
    
    // also, obtain the quantity for the base element. The base quantity also needs to
    // be evaluated at the same load case, dof_value for sensitivity calculation
    base_qty.zero();
    this->getThermalElemQty(name, &base_qty,  FESystemElem::BASE_ELEM::num(), domain, false);
    
    qty->add(-1.0, base_qty);
    qty->scale(1.0/this->perturbation);
    }
}
  
  
#endif // __fesystem_thermal_elem_h__
