// $Id: PlateVonKarman.h,v 1.8.6.1 2008-02-25 04:32:54 manav Exp $

#ifndef __fesystem_plate_von_karman_h__
#define __fesystem_plate_von_karman_h__


// FESystem includes
#include "StructuralElems/structural_elem.h"


namespace Discipline
{
  class AnalysisDisciplineBase;
}


#ifndef PLATE_VON_KARMAN_NONLINEAR_STIFFNESS_COMPONENT_ENUM_ID
#define PLATE_VON_KARMAN_NONLINEAR_STIFFNESS_COMPONENT_ENUM_ID 1
#else
#error
#endif

#ifndef PLATE_VON_KARMAN_NONLINEAR_STIFFNESS_COMPONENT_ENUM_NAME
#define PLATE_VON_KARMAN_NONLINEAR_STIFFNESS_COMPONENT_ENUM_NAME "PLATE_VON_KARMAN_NONLINEAR_STIFFNESS_COMPONENT"
#else
#error
#endif



    

namespace FESystemElem
{
  DeclareEnumClass(PlateVonKarmanQtyEnum);
  
  DeclareEnumName(PLATE_VON_KARMAN_NONLINEAR_STIFFNESS_COMPONENT, PlateVonKarmanQtyEnum,
                  PLATE_VON_KARMAN_NONLINEAR_STIFFNESS_COMPONENT_ENUM_ID,
                  PLATE_VON_KARMAN_NONLINEAR_STIFFNESS_COMPONENT_ENUM_NAME);
  
  
  /// this is a nonlinear plate element, which uses a nonlinear vonKarman strain for 
  /// coupling of membrane and transverse displacements. The membrane is modeled using
  /// the iso-parametric element with lagrangian shape functions, and the plate element
  /// is modeled using the BATOZ element
  class PlateVonKarman: public FESystemElem::StructuralElem
  {
public:
    /// constructor
    PlateVonKarman(const unsigned int dim, 
                   const unsigned int elem_enum_ID,
                   Discipline::AnalysisDisciplineBase& discipline);
    
    /// destructor 
    ~PlateVonKarman();
    
      /// this function returns a vector of the post process quanties for this element. The quantities will 
      /// consist of the strain, stress and their sensitivities. 
      /// @param vector of load case IDs for which the post process quantities need to be calculated
      /// @param vector of DVs for which the quantities need to be calculated
      virtual std::auto_ptr<ElemPostProcessQty> 
      getElementPostProcessQty(std::vector<unsigned int> , std::vector<DesignData::DesignParameter*>);
    
    
    
protected:
      
      /// inserts the finite element types to be used for this element in the input 
      /// parameters
      virtual void getFETypes(std::vector<FEType>& fetypes);
    
    
    /// inserts the quadrature rules to be used for the different elements
    virtual void getQuadratureRules
      (std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures);
    
      /// calculates the mass matrix
      virtual void calculate_M(DenseMatrix<double>* matrix, 
                               const unsigned int design_point_enum_ID,
                               bool sensitivity);
    
    /// calculate stiffness matrix
    virtual void calculate_K(DenseMatrix<double>* matrix,
                             const unsigned int design_point_enum_ID,
                             bool sensitivity);
    
    /// calculate geometric stiffness matrix
    virtual void calculate_K_G(DenseMatrix<double>* matrix,
                               const unsigned int design_point_enum_ID,
                               bool sensitivity);
    
    /// calculate the thermal nodal force vector
    virtual void calculate_F_T(DenseVector<double>* matrix,
                               const unsigned int design_point_enum_ID,
                               bool sensitivity);
    
    /// calculate the linear membrane strain operator matrix
    /// @param matrix the matrix in which the operator is returned
    /// @param point the non-dimensional coordinates of the point where this needs 
    /// to be calculated.
    void calculateLinearMembraneStrainOperator
      (DenseMatrix<double>& matrix,
       const std::vector<std::vector<RealGradient> >& dphi,
       const unsigned int i_qp);
    
    /// calculate the linear plate strain operator matrix
    /// @param matrix the matrix in which the operator is returned
    /// @param point the non-dimensional coordinates of the point where this needs 
    /// to be calculated.
    void calculateLinearPlateStrainOperator
      (DenseMatrix<double>& matrix,
       const std::vector<std::vector<RealGradient> >& dphi,
       const unsigned int i_qp);

    /// calculate the nonlinear membrane strain operator matrix
    /// @param matrix the matrix in which the operator is returned
    /// @param point the non-dimensional coordinates of the point where this needs 
    /// to be calculated.
    void calculateNonlinearMembraneStrainOperator
      (DenseMatrix<double>& matrix,
       const std::vector<std::vector<RealGradient> >& dphi,
       const unsigned int i_qp,
       const DenseVector<double>& local_dofs);
    
    /// calculate the nonlinear membrane strain operator matrix sensitivity for a 
    /// given local dof index
    /// @param matrix the matrix in which the operator is returned
    /// @param point the non-dimensional coordinates of the point where this needs 
    /// to be calculated.
    void calculateNonlinearMembraneStrainOperatorDofSensitivity
      (DenseMatrix<double>& matrix,
       const std::vector<std::vector<RealGradient> >& dphi,
       const unsigned int i_qp,
       const unsigned int local_dof_index);
    
    
    /// this method calls the appropriate function for calculation of 
    /// the specified element.
    void calculatePlateVonKarmanElemQty(DenseMatrix<double>* quantity,
                                        const unsigned int qty_name,
                                        const unsigned int design_point_enum_ID,
                                        const unsigned int domain,
                                        const bool sensitivity);
    
    /// this method calculates the nonlinear part of the stiffness factor, which is 
    /// given as \f$  [K_{nl}] = \int_{\Omega} ([B_{m_{l}}]^T + [B_{m_{nl}}]^T) [D] [B_{m_{nl}}] d \Omega\f$. 
    /// If sensitivity calculate is true, and the sensitivity parameter for this element is shape
    /// sensitivity, then this element calculates the quantity at the perturbed element, 
    /// instead of the shape sensitivity. This is an important difference between the 
    /// other methods of similar name, and should be taken note of.
    void calculateNonlinearStiffnessFactor(DenseMatrix<double>* matrix, 
                                           const unsigned int design_point_enum_ID,
                                           const bool sensitivity_calculation);
    

    /// this method returns the specified quantity of this element.
    template <class QtyType >
      void getPlateVonKarmanQty(const unsigned int name,
                                QtyType* qty,
                                const unsigned int design_point_enum_ID,
                                const unsigned int domain,
                                const bool sensitivity);
      
  };
}


template <class QtyType >
void
FESystemElem::PlateVonKarman::getPlateVonKarmanQty(const unsigned int name,
                                                   QtyType* qty,
                                                   const unsigned int design_point_enum_ID,
                                                   const unsigned int domain,
                                                   const bool sensitivity)
{
  assert (qty != NULL);
  
  // the quantity here is extracted from the shape function factors obtained 
  // from the base class FESystemElemBase. A static matrix is declared for that
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  switch (qty->m() + qty->n() - (n_nodes * 12))
    {
    case 0:
      qty->resize(n_nodes * 6, n_nodes * 6);
      break;
      
    default:
      {
        // keep going
      }
    }
  
  qty->zero();
  
  // first get the factors from the FESystemElemBase
  // row and column store the information of where the block containing the
  // desired quantitiy start in the shape factor matrix.

  if (!this->local_elem_is_initialized_for_DV_map[design_point_enum_ID])
      this->initialize_element(design_point_enum_ID);
  
  switch (sensitivity)
    {
    case true:
      {
        switch (this->sensitivity_parameter)
          {
          case PROPERTY_PARAMETER_ENUM_ID:
            {
              this->calculatePlateVonKarmanElemQty(qty, name, 
                                                   design_point_enum_ID,
                                                   domain, true);
            }
            break;
            
          case SHAPE_PARAMETER_ENUM_ID:
            {
              // make sure that the design point is the base elem
              switch (design_point_enum_ID)
                {
                case BASE_ELEM_ENUM_ID:
                  {
                    // make sure that the perturbed elem is initialized
                    this->initialize_element(FESystemElem::BASE_PLUS_DELTA_ELEM::num());
                  }
                  break;
                  
                case BASE_PLUS_DELTA_ELEM_ENUM_ID:
                default:
                  {
                    Assert(false, ExcInternalError());
                  }
                }
              
              // create quantities for the perturbed and base elements
              static QtyType  base_qty;
              base_qty.zero();
              
              this->calculatePlateVonKarmanElemQty(&base_qty, name,
                                                   FESystemElem::BASE_ELEM::num(),
                                                   domain, false);
              
              this->calculatePlateVonKarmanElemQty(&base_qty, name,
                                                   FESystemElem::BASE_PLUS_DELTA_ELEM::num(),
                                                   domain, false);
              
              qty->add(-1.0, base_qty);
              qty->scale(1.0/this->perturbation);
            }
            break;
            
          default:
            Assert(false, ExcInternalError());
          }
      }
      break;
      
    case false:
    default:
      {
        this->calculatePlateVonKarmanElemQty(qty, name, design_point_enum_ID, domain, false);
      }
    }
  
}



#endif // __fesystem_plate_von_karman_h__

