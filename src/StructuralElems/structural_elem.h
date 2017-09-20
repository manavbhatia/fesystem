// $Id: structural_elem.h,v 1.28.6.1 2007-05-14 16:45:45 manav Exp $

#ifndef __fesystem_structural_elem_h__
#define __fesystem_structural_elem_h__

// C++ includes
#include <memory>


// FESystem includes
#include "FESystem/FESystemElem.h"
#include "Utilities/NameEnumHandler.h"
#include "Properties/ElemDataCard.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "Mesh/FEMeshData.h"

// libMesh includes
#include "numerics/dense_matrix.h"
#include "numerics/dense_vector.h"


// Forward declerations
class Load;
namespace Discipline
{
  class AnalysisDisciplineBase;
}




#ifndef STRUCTURAL_M_MATRIX_ENUM_ID
#define STRUCTURAL_M_MATRIX_ENUM_ID 1
#else
#error
#endif

#ifndef STRUCTURAL_M_MATRIX_ENUM_NAME
#define STRUCTURAL_M_MATRIX_ENUM_NAME "STRUCTURAL_M_MATRIX"
#else
#error
#endif

#ifndef STRUCTURAL_C_MATRIX_ENUM_ID
#define STRUCTURAL_C_MATRIX_ENUM_ID 2
#else
#error
#endif

#ifndef STRUCTURAL_C_MATRIX_ENUM_NAME
#define STRUCTURAL_C_MATRIX_ENUM_NAME "STRUCTURAL_C_MATRIX"
#else
#error
#endif

#ifndef STRUCTURAL_K_MATRIX_ENUM_ID
#define STRUCTURAL_K_MATRIX_ENUM_ID 3
#else
#error
#endif

#ifndef STRUCTURAL_K_MATRIX_ENUM_NAME
#define STRUCTURAL_K_MATRIX_ENUM_NAME "STRUCTURAL_K_MATRIX"
#else
#error
#endif

#ifndef STRUCTURAL_K_G_MATRIX_ENUM_ID
#define STRUCTURAL_K_G_MATRIX_ENUM_ID 7
#else
#error
#endif

#ifndef STRUCTURAL_K_G_MATRIX_ENUM_NAME
#define STRUCTURAL_K_G_MATRIX_ENUM_NAME "STRUCTURAL_K_G_MATRIX"
#else
#error
#endif


#ifndef STRUCTURAL_F_T_VECTOR_ENUM_ID
#define STRUCTURAL_F_T_VECTOR_ENUM_ID 4
#else
#error
#endif

#ifndef STRUCTURAL_F_T_VECTOR_ENUM_NAME
#define STRUCTURAL_F_T_VECTOR_ENUM_NAME "STRUCTURAL_F_T_VECTOR"
#else
#error
#endif

#ifndef STRUCTURAL_F_PRESSURE_VECTOR_ENUM_ID
#define STRUCTURAL_F_PRESSURE_VECTOR_ENUM_ID 5
#else
#error
#endif

#ifndef STRUCTURAL_F_PRESSURE_VECTOR_ENUM_NAME
#define STRUCTURAL_F_PRESSURE_VECTOR_ENUM_NAME "STRUCTURAL_F_PRESSURE_VECTOR"
#else
#error
#endif


#ifndef STRUCTURAL_STRAIN_OPERATOR_ENUM_ID
#define STRUCTURAL_STRAIN_OPERATOR_ENUM_ID 6
#else
#error
#endif

#ifndef STRUCTURAL_STRAIN_OPERATOR_ENUM_NAME
#define STRUCTURAL_STRAIN_OPERATOR_ENUM_NAME "STRUCTURAL_STRAIN_OPERATOR"
#else
#error
#endif



/**
*	a base class for defining the structural elements
 */

namespace FESystemElem
{
  
  DeclareEnumClass(StructuralElemQtyEnum);
  
  DeclareEnumName(STRUCTURAL_M_MATRIX, FESystemElem::StructuralElemQtyEnum,
                  STRUCTURAL_M_MATRIX_ENUM_ID, STRUCTURAL_M_MATRIX_ENUM_NAME);
  
  DeclareEnumName(STRUCTURAL_C_MATRIX, FESystemElem::StructuralElemQtyEnum,
                  STRUCTURAL_C_MATRIX_ENUM_ID, STRUCTURAL_C_MATRIX_ENUM_NAME);
  
  DeclareEnumName(STRUCTURAL_K_MATRIX, FESystemElem::StructuralElemQtyEnum,
                  STRUCTURAL_K_MATRIX_ENUM_ID, STRUCTURAL_K_MATRIX_ENUM_NAME);
  
  DeclareEnumName(STRUCTURAL_K_G_MATRIX, FESystemElem::StructuralElemQtyEnum,
                  STRUCTURAL_K_G_MATRIX_ENUM_ID, STRUCTURAL_K_G_MATRIX_ENUM_NAME);

  DeclareEnumName(STRUCTURAL_F_T_VECTOR, FESystemElem::StructuralElemQtyEnum,
                  STRUCTURAL_F_T_VECTOR_ENUM_ID, STRUCTURAL_F_T_VECTOR_ENUM_NAME);
  
  DeclareEnumName(STRUCTURAL_F_PRESSURE_VECTOR, FESystemElem::StructuralElemQtyEnum,
                  STRUCTURAL_F_PRESSURE_VECTOR_ENUM_ID, STRUCTURAL_F_PRESSURE_VECTOR_ENUM_NAME);
  
  DeclareEnumName(STRUCTURAL_STRAIN_OPERATOR, FESystemElem::StructuralElemQtyEnum,
                  STRUCTURAL_STRAIN_OPERATOR_ENUM_ID, STRUCTURAL_STRAIN_OPERATOR_ENUM_NAME);
  
  
  class StructuralElem: public FESystemElem::FESystemElemBase
    {
public:
      StructuralElem(const unsigned int dim, 
                     const unsigned int elem_enum_ID,
                     Discipline::AnalysisDisciplineBase& discipline);
      
      virtual ~StructuralElem();
      
      
      /// @returns the number of dofs for this element
       virtual unsigned int getNDofs();

       /// @returns the transient system order for the elem equations
       virtual unsigned int getTransientSystemOrder() const;
      
      /// function will get the matrix or vector quantity for this element
      /// specified by the name in the arguement. 
      void getElementAssembledQty(const unsigned int name_enum_ID,
                                  DenseVector<double>* qty);
      
      
      /// function will get the matrix or vector quantity for this element
      /// specified by the name in the arguement. 
      void getElementAssembledQty(const unsigned int name_enum_ID,
                                  DenseMatrix<double>* qty);
      
      
      /// this function will return the sensitivity quantity by procuring the 
      /// required information by itself. It will check for the current load case and design
      /// variable and perform the operations by itself.
      void  getElementAssembledQtySensitivity(const unsigned int name_enum_ID,
                                              DenseVector<double>* qty);
      
      
      /// this function will return the sensitivity quantity by procuring the 
      /// required information by itself. It will check for the current load case and design
      /// variable and perform the operations by itself.
      void  getElementAssembledQtySensitivity(const unsigned int name_enum_ID,
                                              DenseMatrix<double>* qty);
      
      
      /// this function returns a vector of the post process quanties for this element. The quantities will 
      /// consist of the strain, stress and their sensitivities. 
      /// @param vector of load case IDs for which the post process quantities need to be calculated
      /// @param vector of DVs for which the quantities need to be calculated
      virtual std::auto_ptr<ElemPostProcessQty> 
        getElementPostProcessQty(std::vector<unsigned int> , 
                                 std::vector<DesignData::DesignParameter*>) = 0;
      
protected:
        
      /// function to direct the call to appropriate calculation routine based on the 
      /// quantity seeked. The arguements are:
      /// @param  pointer of quantity to return the data in
      /// @param quantity name
      /// @param bool specifying if sensitivity has been asked
      void calculateAssembledQty(DenseMatrix<double>* qty,
                                 const unsigned int qty_enum_ID,
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation);
      
      
      /// function to direct the call to appropriate calculation routine based on the 
      /// quantity seeked. The arguements are:
      /// @param  pointer of quantity to return the data in
      /// @param quantity name
      /// @param bool specifying if sensitivity has been asked
      void calculateAssembledQty(DenseVector<double>* qty,
                                 const unsigned int qty_enum_ID,
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation);
      
      
      /// function will return the nodal temperature vector in the first arguement of the function
      /// @param vector in which the temperatures will be returned
      /// @param bool indicating if sensitivity is being sought
      /// @param DV ID
      void extractNodalTemperatureVectorFromLoads(DenseVector<double>& temp_vec, 
                                                  bool sensitivity);
      
      
      /// calculates the transformation matrix for the structural element
      void getStructuralT_matrix(DenseMatrix<double>* matrix, 
                                 const unsigned int design_point_enum_ID,
                                 bool shape_sensitivity);
      
      
      /**
        *	function will calculate the matrix of the element,
       *	and fill the matrix passed as the arguement. The second arguement is the 
       *	name of the property, which, if passed to this function, will make the 
       *	function calculate the sensitivity of the qty wrt the property
       */
      virtual void calculate_M(DenseMatrix<double>*, 
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation)=0;
      
      /// calculates and returns the element stiffness matrix
      virtual void calculate_K(DenseMatrix<double>*, 
                               const unsigned int design_point_enum_ID,
                               bool sensitivity_calculation)=0;
      
      /// calculates and returns the element geometric stiffness matrix
      virtual void calculate_K_G(DenseMatrix<double>*, 
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation)=0;
      
      /// calculates and returns the element pressure vector
      virtual void calculate_F_Pressure(DenseVector<double>*, 
                                        const unsigned int design_point_enum_ID,
                                        bool sensitivity_calculation);
      
      
      /// calculates and returns the element thermal load vector
      virtual void calculate_F_T(DenseVector<double>*, 
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation)=0;
      
            
      /// this will calculate the body force for the element
      /// @param vector in which the force vector will be stored and returned
      /// @param bool if the sensitivity is needed or not.
      virtual void calculate_F_B(DenseVector<double>* , 
                                 const unsigned int design_point_enum_ID,
                                 bool sensitivity_calculation);
      
      
      /// initializes the property card at element local temperature, 
      /// which is taken as the average of all nodal temperatures. This is needed for temperature 
      /// dependent properties. This is an approximate way of using parameter dependent 
      /// properties, since ideally they should be calculated at each quadrature point. 
      /// That is, however, not done due to memory and CPU constraints.
      virtual void initPropertyCard();
      
      /// this clears the temperature initialization only for the property cards.
      virtual void clearPropertyCardInitialization();
      
      
      /// this is an abstract method to calculate the strain operator matrices
      /// @param the vector where all the strain operators will be returned
      /// @param vector of points, defined in the element local coordinate system where the operators will be 
      /// calculated
      //virtual void calculateStrainOperator(std::vector<DenseMatrix<double> >* , const std::vector<Point>* )=0;
    };
}


inline
unsigned int FESystemElem::StructuralElem::getNDofs()
{
  return 6 * (this->getNNodes());
}


inline
unsigned int FESystemElem::StructuralElem::getTransientSystemOrder() const
{
  return 2;
}

#endif // __fesystem_structural_elem_h__
