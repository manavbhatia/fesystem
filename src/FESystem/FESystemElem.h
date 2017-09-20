// $Id: FESystemElem.h,v 1.16.6.5 2008-06-03 05:19:41 manav Exp $

#ifndef __fesystem_fesystem_elem_h__
#define __fesystem_fesystem_elem_h__

// C++ includes
#include <vector>
#include <string>
#include <memory>
#include <map>

// FESystem includes
#include "Database/ElementDataStorage.h"
#include "DesignData/DesignParameter.h"
#include "FESystem/FESystemElemTypeEnumHandler.h"
#include "DesignData/ShapeParameter.h"
#include "DesignData/PropertyParameter.h"

// libMesh includes
#include "geom/elem.h"
#include "fe/fe.h"
#include "quadrature/quadrature_gauss.h"
#include "numerics/dense_matrix.h"

// Forward declerations
namespace Driver
{
  class AnalysisDriver;
}

namespace Discipline
{
  class AnalysisDisciplineBase;
}

namespace Loads
{
  class VolumeLoadCombination;
  class NodalLoadCombination;
  class SurfaceLoadCombination;
}

class ElemDataCard;
class ElemPostProcessQty;

#ifndef TRANSFORM_MATRIX_ENUM_ID
#define TRANSFORM_MATRIX_ENUM_ID 1
#else
#error
#endif

#ifndef TRANSFORM_MATRIX_ENUM_NAME
#define TRANSFORM_MATRIX_ENUM_NAME "TRANSFORM_MATRIX"
#else
#error
#endif


#ifndef N_N_FACTOR_ENUM_ID
#define N_N_FACTOR_ENUM_ID 2
#else
#error
#endif

#ifndef N_N_FACTOR_ENUM_NAME
#define N_N_FACTOR_ENUM_NAME "N_N_FACTOR"
#else
#error
#endif


#ifndef N_X_N_X_FACTOR_ENUM_ID
#define N_X_N_X_FACTOR_ENUM_ID 3
#else
#error
#endif

#ifndef N_X_N_X_FACTOR_ENUM_NAME
#define N_X_N_X_FACTOR_ENUM_NAME "N_X_N_X_FACTOR"
#else
#error
#endif


#ifndef N_Y_N_Y_FACTOR_ENUM_ID
#define N_Y_N_Y_FACTOR_ENUM_ID 4
#else
#error
#endif

#ifndef N_Y_N_Y_FACTOR_ENUM_NAME
#define N_Y_N_Y_FACTOR_ENUM_NAME "N_Y_N_Y_FACTOR"
#else
#error
#endif


#ifndef N_Z_N_Z_FACTOR_ENUM_ID
#define N_Z_N_Z_FACTOR_ENUM_ID 5
#else
#error
#endif

#ifndef N_Z_N_Z_FACTOR_ENUM_NAME
#define N_Z_N_Z_FACTOR_ENUM_NAME "N_Z_N_Z_FACTOR"
#else
#error
#endif


#ifndef N_X_N_Y_FACTOR_ENUM_ID
#define N_X_N_Y_FACTOR_ENUM_ID 6
#else
#error
#endif

#ifndef N_X_N_Y_FACTOR_ENUM_NAME
#define N_X_N_Y_FACTOR_ENUM_NAME "N_X_N_Y_FACTOR"
#else
#error
#endif


#ifndef N_Y_N_Z_FACTOR_ENUM_ID
#define N_Y_N_Z_FACTOR_ENUM_ID 7
#else
#error
#endif

#ifndef N_Y_N_Z_FACTOR_ENUM_NAME
#define N_Y_N_Z_FACTOR_ENUM_NAME "N_Y_N_Z_FACTOR"
#else
#error
#endif


#ifndef N_Z_N_X_FACTOR_ENUM_ID
#define N_Z_N_X_FACTOR_ENUM_ID 8
#else
#error
#endif

#ifndef N_Z_N_X_FACTOR_ENUM_NAME
#define N_Z_N_X_FACTOR_ENUM_NAME "N_Z_N_X_FACTOR"
#else
#error
#endif


#ifndef N_X_N_FACTOR_ENUM_ID
#define N_X_N_FACTOR_ENUM_ID 9
#else
#error
#endif

#ifndef N_X_N_FACTOR_ENUM_NAME
#define N_X_N_FACTOR_ENUM_NAME "N_X_N_FACTOR"
#else
#error
#endif


#ifndef N_Y_N_FACTOR_ENUM_ID
#define N_Y_N_FACTOR_ENUM_ID 10
#else
#error
#endif

#ifndef N_Y_N_FACTOR_ENUM_NAME
#define N_Y_N_FACTOR_ENUM_NAME "N_Y_N_FACTOR"
#else
#error
#endif


#ifndef N_Z_N_FACTOR_ENUM_ID
#define N_Z_N_FACTOR_ENUM_ID 11
#else
#error
#endif

#ifndef N_Z_N_FACTOR_ENUM_NAME
#define N_Z_N_FACTOR_ENUM_NAME "N_Z_N_FACTOR"
#else
#error
#endif


#ifndef N_FACTOR_ENUM_ID
#define N_FACTOR_ENUM_ID 12
#else
#error
#endif

#ifndef N_FACTOR_ENUM_NAME
#define N_FACTOR_ENUM_NAME "N_FACTOR"
#else
#error
#endif


#ifndef SURFACE_NORMAL_ENUM_ID
#define SURFACE_NORMAL_ENUM_ID 13
#else
#error
#endif

#ifndef SURFACE_NORMAL_ENUM_NAME
#define SURFACE_NORMAL_ENUM_NAME "SURFACE_NORMAL"
#else
#error
#endif


#ifndef ELEM_VOLUME_ENUM_ID
#define ELEM_VOLUME_ENUM_ID 1
#else
#error
#endif

#ifndef ELEM_VOLUME_ENUM_NAME
#define ELEM_VOLUME_ENUM_NAME "ELEM_VOLUME"
#else
#error
#endif

#ifndef SIDE_ZERO_ENUM_ID
#define SIDE_ZERO_ENUM_ID 2
#else
#error
#endif

#ifndef SIDE_ZERO_ENUM_NAME
#define SIDE_ZERO_ENUM_NAME "SIDE_ZERO"
#else
#error
#endif


#ifndef SIDE_ONE_ENUM_ID
#define SIDE_ONE_ENUM_ID 3
#else
#error
#endif

#ifndef SIDE_ONE_ENUM_NAME
#define SIDE_ONE_ENUM_NAME "SIDE_ONE"
#else
#error
#endif


#ifndef SIDE_TWO_ENUM_ID
#define SIDE_TWO_ENUM_ID 4
#else
#error
#endif

#ifndef SIDE_TWO_ENUM_NAME
#define SIDE_TWO_ENUM_NAME "SIDE_TWO"
#else
#error
#endif


#ifndef SIDE_THREE_ENUM_ID
#define SIDE_THREE_ENUM_ID 5
#else
#error
#endif

#ifndef SIDE_THREE_ENUM_NAME
#define SIDE_THREE_ENUM_NAME "SIDE_THREE"
#else
#error
#endif


#ifndef SIDE_FOUR_ENUM_ID
#define SIDE_FOUR_ENUM_ID 6
#else
#error
#endif

#ifndef SIDE_FOUR_ENUM_NAME
#define SIDE_FOUR_ENUM_NAME "SIDE_FOUR"
#else
#error
#endif


#ifndef SIDE_FIVE_ENUM_ID
#define SIDE_FIVE_ENUM_ID 7
#else
#error
#endif

#ifndef SIDE_FIVE_ENUM_NAME
#define SIDE_FIVE_ENUM_NAME "SIDE_FIVE"
#else
#error
#endif


#ifndef BASE_ELEM_ENUM_ID
#define BASE_ELEM_ENUM_ID 1
#else
#error
#endif


#ifndef BASE_ELEM_ENUM_NAME
#define BASE_ELEM_ENUM_NAME "BASE_ELEM"
#else
#error
#endif


#ifndef BASE_PLUS_DELTA_ELEM_ENUM_ID
#define BASE_PLUS_DELTA_ELEM_ENUM_ID 2
#else
#error
#endif


#ifndef BASE_PLUS_DELTA_ELEM_ENUM_NAME
#define BASE_PLUS_DELTA_ELEM_ENUM_NAME "BASE_PLUS_DELTA_ELEM"
#else
#error
#endif


namespace FESystemElem
{
  
  /// this defines the quantity named pertinent to this class. This class provides methods
  /// to calculate these methods.
  DeclareEnumClass(FESystemElemQtyEnum);
  
  DeclareEnumName(TRANSFORM_MATRIX, FESystemElemQtyEnum,
                  TRANSFORM_MATRIX_ENUM_ID,
                  TRANSFORM_MATRIX_ENUM_NAME);
  
  DeclareEnumName( N_N_FACTOR, FESystemElemQtyEnum,
                   N_N_FACTOR_ENUM_ID,
                   N_N_FACTOR_ENUM_NAME);
  
  DeclareEnumName( N_X_N_X_FACTOR, FESystemElemQtyEnum,
                   N_X_N_X_FACTOR_ENUM_ID,
                   N_X_N_X_FACTOR_ENUM_NAME);
  
  DeclareEnumName( N_Y_N_Y_FACTOR, FESystemElemQtyEnum,
                   N_Y_N_Y_FACTOR_ENUM_ID,
                   N_Y_N_Y_FACTOR_ENUM_NAME);
  
  DeclareEnumName( N_Z_N_Z_FACTOR, FESystemElemQtyEnum,
                   N_Z_N_Z_FACTOR_ENUM_ID,
                   N_Z_N_Z_FACTOR_ENUM_NAME);
  
  DeclareEnumName( N_X_N_Y_FACTOR, FESystemElemQtyEnum,
                   N_X_N_Y_FACTOR_ENUM_ID,
                   N_X_N_Y_FACTOR_ENUM_NAME);
  
  DeclareEnumName( N_Y_N_Z_FACTOR, FESystemElemQtyEnum,
                   N_Y_N_Z_FACTOR_ENUM_ID,
                   N_Y_N_Z_FACTOR_ENUM_NAME);
  
  DeclareEnumName( N_Z_N_X_FACTOR, FESystemElemQtyEnum,
                   N_Z_N_X_FACTOR_ENUM_ID,
                   N_Z_N_X_FACTOR_ENUM_NAME);
  
  DeclareEnumName( N_X_N_FACTOR, FESystemElemQtyEnum,
                   N_X_N_FACTOR_ENUM_ID,
                   N_X_N_FACTOR_ENUM_NAME);
  
  DeclareEnumName( N_Y_N_FACTOR, FESystemElemQtyEnum,
                   N_Y_N_FACTOR_ENUM_ID,
                   N_Y_N_FACTOR_ENUM_NAME);
  
  DeclareEnumName( N_Z_N_FACTOR, FESystemElemQtyEnum,
                   N_Z_N_FACTOR_ENUM_ID,
                   N_Z_N_FACTOR_ENUM_NAME);
  
  DeclareEnumName( N_FACTOR, FESystemElemQtyEnum,
                   N_FACTOR_ENUM_ID,
                   N_FACTOR_ENUM_NAME);
  
  DeclareEnumName( SURFACE_NORMAL, FESystemElemQtyEnum,
                   SURFACE_NORMAL_ENUM_ID,
                   SURFACE_NORMAL_ENUM_NAME);
  
  
  DeclareEnumClass(IntegrationDomainEnum);
  
  DeclareEnumName( ELEM_VOLUME, IntegrationDomainEnum,
                   ELEM_VOLUME_ENUM_ID,
                   ELEM_VOLUME_ENUM_NAME);
  
  
  DeclareEnumName( SIDE_ZERO, IntegrationDomainEnum,
                   SIDE_ZERO_ENUM_ID,
                   SIDE_ZERO_ENUM_NAME);
  
  
  DeclareEnumName( SIDE_ONE, IntegrationDomainEnum,
                   SIDE_ONE_ENUM_ID,
                   SIDE_ONE_ENUM_NAME);
  
  
  DeclareEnumName( SIDE_TWO, IntegrationDomainEnum,
                   SIDE_TWO_ENUM_ID,
                   SIDE_TWO_ENUM_NAME);
  
  
  DeclareEnumName( SIDE_THREE, IntegrationDomainEnum,
                   SIDE_THREE_ENUM_ID,
                   SIDE_THREE_ENUM_NAME);
  
  
  DeclareEnumName( SIDE_FOUR, IntegrationDomainEnum,
                   SIDE_FOUR_ENUM_ID,
                   SIDE_FOUR_ENUM_NAME);
  
  
  DeclareEnumName( SIDE_FIVE, IntegrationDomainEnum,
                   SIDE_FIVE_ENUM_ID,
                   SIDE_FIVE_ENUM_NAME);
  

  
  DeclareEnumClass(DesignPointElemEnum);


  DeclareEnumName( BASE_ELEM, DesignPointElemEnum, 
		   BASE_ELEM_ENUM_ID, 
		   BASE_ELEM_ENUM_NAME);


  DeclareEnumName( BASE_PLUS_DELTA_ELEM, DesignPointElemEnum, 
		   BASE_PLUS_DELTA_ELEM_ENUM_ID, 
		   BASE_PLUS_DELTA_ELEM_ENUM_NAME);
  
  /**
    *	this is a finite element base class for elements that will know how
   *	to calculate their element quantities, like the stiffness, mass matrices, 
   *	etc. and store them in the FESystem element data storage object.
   */
  class FESystemElemBase
    {
public:
      
      
      ///	constructor
      /// @param pointer to the analysis driver that is controlling the sequence of operations
      FESystemElemBase(const unsigned int dim,
                       const unsigned int elem_enum_ID,
                       Discipline::AnalysisDisciplineBase& discipline);
      
      /**
      *	destructor
       */
      virtual ~FESystemElemBase();
      
      /// @returns the element ID
      inline unsigned int getID() const; 
      
      /// @returns the number of nodes for this elem.
      inline unsigned int getNNodes();
      
      /// @returns the number of dofs for this element
      virtual unsigned int getNDofs() = 0;
      
      /// @returns the transient system order for the elem equations
      virtual unsigned int getTransientSystemOrder() const = 0;
  
      /// @returns the element enum ID 
      inline unsigned int getFESystemElemTypeEnumID() const;
      
      /// @returns the geometric element for the unperturbed mesh
      inline const Elem* getElem(const unsigned int design_point_enum_ID) const;

      /// returns the element enum ID 
      inline const std::string getFESystemElemTypeEnumName() const;
      
      /// reinitializes an element for a new geometric element. The element specified 
      /// here is the element at base design point
      inline void reinit(const unsigned int elem_ID, 
                         Elem* elem, 
                         ElemDataCard* elem_data_card);
      
      
      /// reinitializes element for shape sensitivity. The element specified 
      /// here is the perturbed element. The forward finite dofferencing is used for this.
      inline void reinitForShapeSensitivity(const unsigned int param_ID,
                                            Elem* perturb_elem,
                                            const double perturb_value);
      
      
      /// reinitializes element for shape sensitivity
      inline void reinitForPropertySensitivity(const unsigned int param_ID);
      
      
      /// this method clears any previous initialization of the element.
      inline void clearInitialization();
      
      /// this method clears the local data that has been created for various purposes
      inline void clearElemInitialization();
      
      /// this method clears any previous initialization of the element.
      inline void clearSensitivityInitialization();
      
      
      /// this method sets the values of the dofs for the element. The vectors containing 
      /// the dof values are placed in the arguements according to their time derivative. 
      /// If a derivative order is not needed, its pointer is left to NULL.
      inline void setDofValues(std::vector<DenseVector<double>*>& dofs,
                               std::vector<DenseVector<double>*>& dof_sens);
      
      /// this method sets the values of the element loads 
      /// @param loads map of loads for this element
      /// @param load_sens sensitivities of the loads specified in the loads parameter. This
      /// is optional. For shape sensitivity, the loads are assumed for the perturned element, 
      /// for a forward difference sensitivity calculation.
      void setElementLoads
	(std::map<unsigned int, const Loads::VolumeLoadCombination*> *vol_loads,
	 std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> > *surf_loads,
	 std::map<unsigned int, std::map<unsigned int, const Loads::NodalLoadCombination*> > *nod_loads,
	 std::map<unsigned int, const Loads::VolumeLoadCombination*> *vol_load_sens,
	 std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> > *surf_load_sens,
	 std::map<unsigned int, std::map<unsigned int, const Loads::NodalLoadCombination*> > *nod_load_sens);

      /// returns the assembled quantitiy for the element
      virtual void getElementAssembledQty(const unsigned int name_enum_ID,
                                          DenseVector<double>* qty) = 0; 
      
      virtual void getElementAssembledQty(const unsigned int name_enum_ID,
                                          DenseMatrix<double>* qty) = 0;
      
      virtual void  getElementAssembledQtySensitivity(const unsigned int name_enum_ID,
                                                      DenseVector<double>* qty) = 0;
      
      virtual void  getElementAssembledQtySensitivity(const unsigned int name_enum_ID,
                                                      DenseMatrix<double>* qty) = 0;
      
      /// this method returns the quantity specified in the arguement list. This has 
      /// been templatized so that the template can be used in derived element types.
      /// This element class will provide its element quantity and calculation routines.
      /// @param name  enum ID of the quantity to be returned
      /// @param qty data structure in which the quantity will be removed
      /// @param design_point_enum_ID enum ID of the design point for which the quantity is requested
      /// @param side_num number of the side if this is a side quantity
      /// @param family enum ID of the family of FE for which the quantity is needed
      template <class QtyType> 
        void getFESystemElemQty(const unsigned int name_enum_ID,
                                QtyType* qty,
                                const unsigned int design_point_enum_ID,
                                const unsigned int domain_enum_ID,
                                const FEFamily& family);
      
      
      /// this method returns the shape sensitivity of quantity specified in the arguement list. 
      /// This has been templatized so that the template can be used in derived element types.
      /// This element class will provide its element quantity and calculation routines.
      /// @param name  name of the quantity to be returned
      /// @param qty data structure in which the quantity will be removed
      /// @param side_num number of the side if this is a side quantity
      template <class QtyType > void getFESystemElemQtyShapeSensitivity
        (const unsigned int name_enum_ID,
         QtyType* qty,
         const unsigned int domain_enum_ID,
         const FEFamily& family);

      /// calculates and returns the element post process quantity
      /// @param vector of load cases for which the quantities will be calculated
      /// @param vector of design variables for which the quantities will be calculated
      virtual std::auto_ptr<ElemPostProcessQty>
        getElementPostProcessQty(std::vector<unsigned int> load_cases,
                                 std::vector<DesignData::DesignParameter* >design_variables) = 0;
      
      
      /// calculates and returns the value of weight of the element
      double getElementMass(bool sensitivity);
            
      /// this clears the property element loads and dof values
      inline void clearLoadAndDofs();
      
      /// clears the dof values that were set for the element
      void clearDofValues();
        
protected:
        
      /// deletes the local data for the element at the specified design point
      inline void deleteElemLocalData(const unsigned int design_point_enum_ID);
        
      /// inserts the finite element types to be used for this element in the input 
      /// parameters
      virtual void getFETypes(std::vector<FEType>& fetypes) = 0;
      
      
      /// inserts the quadrature rules to be used for the different elements
      virtual void getQuadratureRules
        (std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures) = 0;
      
      
      /// initializes the set of elements to be used for this element
      void initializeFEAndQuadratureDataStructures();
      
      
      
      /// this method is an interface to the specific functions to calculate element
      /// quantities. It calls the method that performs the calculations and returns the
      /// quantity in the first arguement. 
      void calculateFESystemElemQty
        (DenseMatrix<double>* qty, 
         const unsigned int name_enum_ID,
         const unsigned int design_point_enum_ID,
         const unsigned int domain_enum_ID,
         const FEFamily& family);
      
      /// this method is an interface to the specific functions to calculate element
      /// quantities. It calls the method that performs the calculations and returns the
      /// quantity in the first arguement. 
      void calculateFESystemElemQty
        (DenseVector<double>* qty, 
         const unsigned int name_enum_ID, 
         const unsigned int design_point_enum_ID,
         const unsigned int domain_enum_ID,
         const FEFamily& family);
      
      /// this method calculates the quantities based on the shape functions 
      /// and returns it in the first arguement
      void calculateShapeFunctionFactors
        (DenseMatrix<double>* matrix,
         const unsigned int qty_name_enum_ID,
         const unsigned int design_point_enum_ID,
         const unsigned int domain_enum_ID,
         const FEFamily& family);
      
      
      /// this method calculates the quantities based on the shape functions
      /// and returns it in the first arguement 
      void calculateShapeFunctionFactors
        (DenseVector<double>* vector,
         const unsigned int qty_name_enum_ID,
         const unsigned int design_point_enum_ID,
         const unsigned int domain_enum_ID,
         const FEFamily& family);
      
      
      
      /// this method returns the side number associated with the domain specified in the 
      /// arguement list
      inline unsigned int getSideNumberFromDomainEnum(const unsigned int side) const;
      
      
      
      /// this method returns the side number associated with the domain specified in the 
      /// arguement list
      inline unsigned int getDomainEnumFromSideNumber(unsigned int  side) const;
      
      
      
      /// this function returns the name of the given quantity
      std::string getStringNameForFESystemElemQty
        ( const unsigned int name_enum_ID,
          const unsigned int design_point_enum_ID,
          const unsigned int domain_enum_ID,
          const FEFamily& family);
      
      
      
      /// this function returns the name of the sensitivity of the specified quantity
      std::string getStringNameForFESystemElemQtySensitivity
        (const unsigned int name_enum_ID,
         unsigned int DV_num,
         const unsigned int domain_enum_ID,
         const FEFamily& family);
      
      
      
      /// this function caluclates the surface normal for the element side
      /// @param vector vector in which the normal will be returned
      /// @param side  side for which the normal has to be calculated
      void calculateSurfaceNormal(DenseVector<double>* vector,
                                  const unsigned int design_point_enum_ID,
                                  const unsigned int side_enum_ID);
      
      
      
      /**
        *	function to calculate and return the transformation matrix.
       *	the arguement will be resized and filled by the matrix
       *	entries. 
       * The transformation matrix is defined by the relation Q_ij = e_local_i . e_global_j, so that
       * e_local_i = Q_ij . e_local_j
       */
      void calculate_T_matrix(DenseMatrix<double>* ,
                              const unsigned int design_point_enum_ID);
            
      /// initializes the element and the respective data structures like phi, dphi, etc.
      /// If the input arguement is NULL, the element is initialized at the quandrature points
      /// otherwise, it is initialized at the points specifies in the input parameter. It used the 
      /// _local_elem pointer to initialize this element. 
      void initialize_element(const unsigned int design_point_enum_ID,
                              const std::vector<Point>* = NULL);
      
      
      /// initialized the element side given by the side number passed in the arguement
      void initialize_element_side(unsigned int ,
				   const unsigned int design_point_enum_ID);
      
      
//      /// initializes the perturbed elem for this element
//      void initPerturbedElemForShapeParameter(DesignData::DesignParameter* param);
      
      /// initializes the property cards for local parameters
      virtual void initPropertyCard() = 0;
      
      /// this clears the temperature initialization only for the property cards.
      virtual void clearPropertyCardInitialization() = 0;
      
      /// this  function creates a local elem for this general element in 3-D space.
      /// The local element pointer is then 
      /// stored in the _local_elem pointer. This is a generic method, and accesses 
      /// the respective method for 1-D and 2-D
      /// depending on the dimension of the element
      void create_local_elem(const unsigned int design_point_enum_ID);
      
      /// this function creates  a local element in 1-D space for line elements
      void create_local_1D_elem(const unsigned int design_point_enum_ID);
      
      /// this function creates a local element in 2-D space for quad or tri elements
      void create_local_2D_elem(const unsigned int design_point_enum_ID);
      
      /// this exception is thrown when an invalid elem is called
      DeclException0(ExcInvalidElemGiven);
      
      /// the discipline that owns this element
      Discipline::AnalysisDisciplineBase& analysis_discipline;      
      
      /// dimension of the element
      const unsigned int dimension;
      
      /// enumID of the element kind
      const unsigned int elem_type_enum_ID;
      
      /// discipline of the element
      const unsigned int discipline_enumID;
      
      /// ID of the element
      unsigned int elem_ID;
            
      /// ID of the sensitivity parameter
      unsigned int sensitivity_parameter_ID;
      
      /// type of the sensitivity parameter
      unsigned int sensitivity_parameter;
      
      /// perturbation for a finite difference calculation
      double perturbation;
            
      /// boolean for keeping track of whether the FEBase and QBase data structres have
      /// been initialized or not
      bool fe_quadrature_data_initialized;

      /// boolean for keeping track of whether this element is initialized or not
      bool elem_is_initialized_for_property_sensitivity;
      
      /// boolean for keeping track of whether the property card was initialized
      bool property_card_initialized;

      //      /// the design parameter for which sensitivity is needed
      //      DesignData::DesignParameter* DV;
      
      /// card that stores the information for this element
      ElemDataCard* elem_property_card;
      
      /// database for storing element quantities
      ElementDataStorage* elem_data_storage;
      
      /// map of element loads
      std::map<unsigned int, const Loads::VolumeLoadCombination*> *volume_loads;
      std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> > *surface_loads;
      std::map<unsigned int, std::map<unsigned int, const Loads::NodalLoadCombination*> > *nodal_loads;
      
      /// map of element load sensitivity for the current design 
      /// variable for which computations are being done.
      std::map<unsigned int, const Loads::VolumeLoadCombination*> *volume_load_sens;
      std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> > *surface_load_sens;
      std::map<unsigned int, std::map<unsigned int, const Loads::NodalLoadCombination*> > *nodal_load_sens;

      /// pointer to the base (unperturbed) element. There is one elem for each design point, 
      /// including the base point and the perturbed points.
      std::map<unsigned int, Elem*> geometric_elems_for_DV_map;
      
      /// local geometric element created in the element local axis. This is stored for each DV point, 
      /// including the base design point, and all the perturbed DVs. This is essential for shape 
      /// sensitivity
      std::map<unsigned int, Elem*> local_elem_for_DV_map;
      
      /// nodes created for the local element. These are stored for each DV
      std::map<unsigned int, Node**> local_elem_nodes_for_DV_map;
      
      /// boolean for keeping track of whether this element is initialized or not. This 
      /// keeps track of each specific element for design points
      std::map<unsigned int, bool> elem_is_initialized_for_DV_map;
      
      /// boolean for keeping track of whether this element is initialized or not
      std::map<unsigned int, bool> local_elem_is_initialized_for_DV_map;
      
      
      /// boolean for keeping track of wheter this element is initialized for post processing
      /// or not
      //      bool elem_is_initialized_for_post_processing;
      
      /// keeps track of the initialized side that was initialized. This is stored for each design point.
      std::map<unsigned int, int > initialized_side_for_DV_map;

    
      /// vector of element dof values
      std::vector<DenseVector<double> > dof_values_vec;
      
      /// vector of element dof values
      std::vector<DenseVector<double> > dof_value_sensitivity_vec;

      /// finite element object for the element domain. This is stored for each design point.
      std::map<unsigned int, std::map<FEFamily, FEBase*> > fe_base_map_for_DV;
      
      /// finite element object for the side. This is stored for each design point.
      std::map<unsigned int, std::map<FEFamily, FEBase*> > fe_side_map_for_DV;
      
      /// quadrature rule for the element volume integration. This is stored for each design point.
      std::map<unsigned int, std::map<FEFamily, QBase*> > qbase_map_for_DV;
      
      /// quadrature rule for the side integration. This is stored for each design point.
      std::map<unsigned int, std::map<FEFamily, QBase*> > qbase_side_map_for_DV;
      
//      /// vector of JxW for all quadrature points on the volume of this element. 
//      /// This is stored for each design point.
//      std::map<unsigned int, std::map<FEFamily, std::vector<Real>*> > JxW_map_for_DV;
//      
//      /// vector of shape functions for this element at all quadrature points.
//      /// This is stored for each design point.
//      std::map<unsigned int, std::map<FEFamily, std::vector<std::vector<Real> >*> > phi_map_for_DV;
//      
//      /// vector of gradient of shape functions at all quadrature points. 
//      /// This is stored for each design point.
//      std::map<unsigned int, std::map<FEFamily, std::vector<std::vector<RealGradient> >* > > dphi_map_for_DV;
//      
//      /// vector of shape functions on the element side. This is stored for each design point.
//      std::map<unsigned int, std::map<FEFamily, std::vector<std::vector<Real> >* > > phi_face_map_for_DV;
//      
//      /// vector of JxW on the element side. This is stored for each design point.
//      std::map<unsigned int, std::map<FEFamily, std::vector<Real>*> > JxW_face_map_for_DV;
    };
  
  
  std::auto_ptr<FESystemElem::FESystemElemBase> 
    createFESystemElem(const unsigned int elem_kind_enum_ID, 
                       Discipline::AnalysisDisciplineBase& discipline_base);
}




inline
unsigned int
FESystemElem::FESystemElemBase::getID() const
{
  Assert(this->elem_is_initialized_for_DV_map.find(FESystemElem::BASE_ELEM::num())->second,
         ExcInvalidState());
  
  return this->elem_ID;
}



inline 
unsigned int 
FESystemElem::FESystemElemBase::getNNodes()
{
  return this->geometric_elems_for_DV_map[FESystemElem::BASE_ELEM::num()]->n_nodes();
}



inline
unsigned int
FESystemElem::FESystemElemBase::getFESystemElemTypeEnumID() const
{
  return this->elem_type_enum_ID;
}




inline
const std::string
FESystemElem::FESystemElemBase::getFESystemElemTypeEnumName() const
{
  return FESystemElem::FESystemElemTypeEnum::enumName(this->elem_type_enum_ID);
}



inline
const Elem*
FESystemElem::FESystemElemBase::getElem(const unsigned int design_point_enum_ID) const
{
  const Elem* return_elem = NULL;

  switch(this->elem_is_initialized_for_DV_map.find(design_point_enum_ID)->second)
    {
    case true:
      return_elem = this->geometric_elems_for_DV_map.find(design_point_enum_ID)->second;
      break;

    case false:
    default:
      Assert(false, ExcInvalidState());
    }
  
  return return_elem;
}




inline 
unsigned int 
FESystemElem::FESystemElemBase::getSideNumberFromDomainEnum(const unsigned int side) const
{
  switch (side)
    {
    case ELEM_VOLUME_ENUM_ID:
      return this->geometric_elems_for_DV_map.find(FESystemElem::BASE_ELEM::num())->second->n_sides();
      break;
      
    case SIDE_ZERO_ENUM_ID:
      return 0;
      break;
      
    case SIDE_ONE_ENUM_ID:
      return 1;
      break;
      
    case SIDE_TWO_ENUM_ID:
      return 2;
      break;
      
    case SIDE_THREE_ENUM_ID:
      return 3;
      break;
      
    case SIDE_FOUR_ENUM_ID:
      return 4;
      break;
      
    case SIDE_FIVE_ENUM_ID:
      return 5;
      break;
      
    default:
      abort();
      break;
    }
}





inline 
unsigned int
FESystemElem::FESystemElemBase::getDomainEnumFromSideNumber(unsigned int side) const
{
  if (side == this->geometric_elems_for_DV_map.find(FESystemElem::BASE_ELEM::num())->second->n_sides())
    return FESystemElem::ELEM_VOLUME::num();
  
  switch (side)
    {
    case 0:
      return FESystemElem::SIDE_ZERO::num();
      break;
      
    case 1:
      return FESystemElem::SIDE_ONE::num();
      break;
      
    case 2:
      return FESystemElem::SIDE_TWO::num();
      break;
      
    case 3:
      return FESystemElem::SIDE_THREE::num();
      break;
      
    case 4:
      return FESystemElem::SIDE_FOUR::num();
      break;
      
    case 5:
      return FESystemElem::SIDE_FIVE::num();
      break;
      
    default:
      abort();
      break;
    }
  
}





template <class QtyType >
void 
FESystemElem::FESystemElemBase::getFESystemElemQty(const unsigned int name,
                                                   QtyType* qty,
                                                   const unsigned int design_point_enum_ID,
                                                   const unsigned int domain,
                                                   const FEFamily& family)
{
  assert (qty != NULL);
  
  // get the string name for the quantity
  static std::string qty_string_name;
  qty_string_name.clear();
  qty_string_name = this->getStringNameForFESystemElemQty(name, 
                                                          design_point_enum_ID,
                                                          domain, family);
  
  // check if the quantity exists in the map element database or not
  // if not, calculate and store it
  if (this->elem_data_storage->checkIfQtyExists
      (this->discipline_enumID, this->elem_ID, qty_string_name))
    {
    // get the factor the database
    this->elem_data_storage->get_element_data(*qty, this->discipline_enumID, 
                                              this->elem_ID, qty_string_name);
    }
  else
    {
    this->calculateFESystemElemQty(qty, name, design_point_enum_ID,
                                   domain, family);
    
    this->elem_data_storage->add_element_data(*qty, this->discipline_enumID,
                                              this->elem_ID, qty_string_name);
    }
}





template <class QtyType >
void 
FESystemElem::FESystemElemBase::getFESystemElemQtyShapeSensitivity(const unsigned int name,
                                                                   QtyType* qty,
                                                                   const unsigned int domain,
                                                                   const FEFamily& family)
{
  Assert(this->elem_is_initialized_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()],
         ExcInvalidState());
  assert (qty != NULL);
  
  // create quantities for the perturbed and base elements
  static QtyType  base_qty;
  
  // get the string name for the quantity
  static std::string qty_string_name, qty_sensitivity_string_name;
  qty_string_name.clear();  qty_sensitivity_string_name.clear();
  
  qty_string_name = this->getStringNameForFESystemElemQty(name, 
                                                          FESystemElem::BASE_ELEM::num(),
                                                          domain, family);
  qty_sensitivity_string_name = 
    this->getStringNameForFESystemElemQtySensitivity(name,
                                                     this->sensitivity_parameter_ID,
                                                     domain,
                                                     family);
  
  // check if the quantity exists in the element database or not
  // if not, calculate and store it
  if (this->elem_data_storage->checkIfQtyExists
      (this->discipline_enumID, this->elem_ID, qty_sensitivity_string_name))
    {
    // get the factor from the database
    this->elem_data_storage->get_element_data(*qty, 
                                              this->discipline_enumID, this->elem_ID, 
                                              qty_sensitivity_string_name);
    }
  else
    {
    // next, obtain the quantity for the perturbed element
    this->calculateFESystemElemQty(qty, name,
                                   FESystemElem::BASE_PLUS_DELTA_ELEM::num(), 
                                   domain, family);
    
    // also, obtain the quantity for the base element. The base quantity also needs to
    // be evaluated at the same load case, dof_value for sensitivity calculation
    base_qty.zero();
    this->getFESystemElemQty(name, &base_qty,  FESystemElem::BASE_ELEM::num(),
                             domain, family);
    
    qty->add(-1.0, base_qty);
    qty->scale(1.0/this->perturbation);
    
    this->elem_data_storage->add_element_data(*qty,
                                              this->discipline_enumID, this->elem_ID, 
                                              qty_sensitivity_string_name);
    }
}


inline
void 
FESystemElem::FESystemElemBase::reinit(const unsigned int ID, 
                                       Elem* elem, 
                                       ElemDataCard* data_card)
{
  // first make sure that this element has been cleared of any previous 
  // initializations
  Assert(!this->elem_is_initialized_for_DV_map[FESystemElem::BASE_ELEM::num()], ExcInvalidState());
  Assert(!this->elem_is_initialized_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()], 
         ExcInvalidState());
  Assert(!this->elem_is_initialized_for_property_sensitivity, ExcInvalidState());
  
  // now make sure that the pointers are not NULL
  Assert(elem != NULL, ExcEmptyObject());
  Assert(data_card != NULL, ExcEmptyObject());
  
  // make sure that the provided element is of the same kind as this elem will handle
  Assert(elem->type() == FESystemElem::FESystemElemTypeEnum::elemType(this->elem_type_enum_ID),
         ExcInternalError());
   
  // finally, store the data provided
  this->elem_ID = ID;
  this->geometric_elems_for_DV_map[FESystemElem::BASE_ELEM::num()] = elem;
  this->elem_property_card = data_card;

  // the fe and quad data structures need to be initialized once for this finite element
  if (!this->fe_quadrature_data_initialized)
    this->initializeFEAndQuadratureDataStructures();

  
  // also, resize the dof vectors if that has not already been done
  static const unsigned int n_dofs = this->getNDofs();
  
  if (this->dof_values_vec.size() != (this->getTransientSystemOrder()+1))
    {
    this->dof_values_vec.clear();
    this->dof_value_sensitivity_vec.clear();
    
    for (unsigned int i=0; i <= this->getTransientSystemOrder(); i++)
      {
      this->dof_values_vec.push_back(DenseVector<double>(n_dofs));
      this->dof_values_vec[i].zero();
      this->dof_value_sensitivity_vec.push_back(DenseVector<double>(n_dofs));
      this->dof_value_sensitivity_vec[i].zero();
      }
    }
    
  this->elem_is_initialized_for_DV_map[FESystemElem::BASE_ELEM::num()] = true;
}




inline
void 
FESystemElem::FESystemElemBase::reinitForShapeSensitivity(const unsigned int param_ID,
                                                          Elem* perturb_elem,
                                                          const double perturb_value)
{
  // first make sure that this element has been cleared of any previous 
  // initializations
  Assert(this->elem_is_initialized_for_DV_map[FESystemElem::BASE_ELEM::num()], 
	 ExcInvalidState());
  Assert(!this->elem_is_initialized_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()], 
	 ExcInvalidState());
  Assert(!this->elem_is_initialized_for_property_sensitivity, ExcInvalidState());
  
  // now make sure that the pointers are not NULL
  Assert(perturb_elem != NULL, ExcEmptyObject());
  Assert(fabs(perturb_value) > FESystemNumbers::Epsilon, ExcInternalError());
  
  // make sure that the provided element is of the same kind as this elem will handle
  Assert(perturb_elem->type() == 
         FESystemElem::FESystemElemTypeEnum::elemType(this->elem_type_enum_ID),
         ExcInternalError());
  
  
  // finally, store the data provided
  this->sensitivity_parameter_ID = param_ID;
  this->geometric_elems_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = perturb_elem;
  this->perturbation = perturb_value;
  
  this->sensitivity_parameter = DesignData::SHAPE_PARAMETER::num();
  
  this->elem_is_initialized_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = true;
}



inline
void
FESystemElem::FESystemElemBase::reinitForPropertySensitivity(const unsigned int param_ID)
{
  // first make sure that this element has been cleared of any previous 
  // initializations
  Assert(this->elem_is_initialized_for_DV_map[FESystemElem::BASE_ELEM::num()], 
	 ExcInvalidState());
  Assert(!this->elem_is_initialized_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()], 
	 ExcInvalidState());
  Assert(!this->elem_is_initialized_for_property_sensitivity, ExcInvalidState());
  
  
  // now make sure that the pointers are not NULL
  this->sensitivity_parameter_ID = param_ID;
  
  this->sensitivity_parameter = DesignData::PROPERTY_PARAMETER::num();
  
  this->elem_is_initialized_for_property_sensitivity = true;
}


  

inline
void 
FESystemElem::FESystemElemBase::deleteElemLocalData(const unsigned int design_point_enum_ID)
{
  unsigned int n_nodes;
  
  // clear the base element data
  // first release the local elem pointer in case it is 3D, since for that case
  //  the local element is just a pointer to the base elem
  
  switch (this->local_elem_for_DV_map[design_point_enum_ID] == NULL)
    {
    case true:
      {
        // nothing to be done here, since the element was not created.
      }
      break;
      
    case false:
    default:
      {
        switch (this->local_elem_for_DV_map[design_point_enum_ID]->dim())
          {
          case 3:
            this->local_elem_for_DV_map[design_point_enum_ID] = NULL; 
            break;
            
          default:
            {
              // delete the element
              n_nodes = this->local_elem_for_DV_map[design_point_enum_ID]->n_nodes();
              delete this->local_elem_for_DV_map[design_point_enum_ID];
              this->local_elem_for_DV_map[design_point_enum_ID] = NULL;
              
              // and then delete the nodes. it is assumed that if the local elem was created, 
              // the nodes too will be present, and need to be deleted
              for (unsigned int i = 0; i< n_nodes; i++)
                delete this->local_elem_nodes_for_DV_map[design_point_enum_ID][i];
              
              // finally, delete the array that was created to store the nodes.
              delete[] this->local_elem_nodes_for_DV_map[design_point_enum_ID];
              this->local_elem_nodes_for_DV_map[design_point_enum_ID] = NULL;
            }
            break;
          }
      }
    }
  
  // next, set the initialized bools to false
  this->local_elem_is_initialized_for_DV_map[design_point_enum_ID] = false;
}



inline
void 
FESystemElem::FESystemElemBase::clearInitialization()
{
  // clear the property card initialization
  this->clearPropertyCardInitialization();

  // first clear any sensitivity data
  this->clearSensitivityInitialization();

  // clear the local data
  this->clearElemInitialization();

  // clear element load and dof data
  this->clearLoadAndDofs();
}





inline
void 
FESystemElem::FESystemElemBase::clearLoadAndDofs()
{
  for (unsigned int i=0; i <= this->getTransientSystemOrder(); i++)
    {
    this->dof_values_vec[i].zero();
    this->dof_value_sensitivity_vec[i].zero();
    }
    
  this->volume_loads = NULL;
  this->surface_loads = NULL;
  this->nodal_loads = NULL;

  this->volume_load_sens = NULL;
  this->surface_load_sens = NULL;
  this->nodal_load_sens = NULL;
}



inline
void 
FESystemElem::FESystemElemBase::clearDofValues()
{
  for (unsigned int i=0; i <= this->getTransientSystemOrder(); i++)
    {
    this->dof_values_vec[i].zero();
    this->dof_value_sensitivity_vec[i].zero();
    }
}



inline
void 
FESystemElem::FESystemElemBase::clearElemInitialization()
{
  // the necessary flags are set, so that the element will have to be reinitialized 
  // for any sort of calculation
  this->elem_ID = FESystemNumbers::InvalidID;
  this->geometric_elems_for_DV_map[FESystemElem::BASE_ELEM::num()] = NULL;
  this->elem_property_card = NULL;

  // next, set the initialization flags to false
  this->elem_is_initialized_for_DV_map[FESystemElem::BASE_ELEM::num()] = false;
  this->local_elem_is_initialized_for_DV_map[FESystemElem::BASE_ELEM::num()] = false;
  this->initialized_side_for_DV_map[FESystemElem::BASE_ELEM::num()] = -1;

  //  this->elem_is_initialized_for_post_processing = false;
}



inline
void 
FESystemElem::FESystemElemBase::clearSensitivityInitialization()
{
  this->sensitivity_parameter_ID = FESystemNumbers::InvalidID;
  this->perturbation = 0.0;

  this->sensitivity_parameter = FESystemNumbers::InvalidID;
  // this is needed for shape parameters only, however, they are still set to false, 
  // since it does not harm anything.
  this->geometric_elems_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = NULL;
  this->elem_is_initialized_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = false;
  this->local_elem_is_initialized_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = false;
  this->initialized_side_for_DV_map[FESystemElem::BASE_PLUS_DELTA_ELEM::num()] = -1;
  
  // next, set the initialized bools to false
  this->elem_is_initialized_for_property_sensitivity = false;
}



inline
void
FESystemElem::FESystemElemBase::setDofValues(std::vector<DenseVector<double>*>& dofs, 
                                             std::vector<DenseVector<double>*>& dof_sens)
{
  Assert(this->dof_values_vec.size() == (this->getTransientSystemOrder()+1), ExcInternalError());
  Assert(this->dof_value_sensitivity_vec.size() == (this->getTransientSystemOrder()+1), 
         ExcInternalError());

  Assert(dofs.size() == (this->getTransientSystemOrder()+1), ExcInternalError());
  Assert(dof_sens.size() == (this->getTransientSystemOrder()+1), ExcInternalError());
  Assert(dofs[0] != NULL, ExcInternalError());
  
  for (unsigned int i=0; i <= this->getTransientSystemOrder(); i++)
    {
    if (dofs[i] != NULL)
      this->dof_values_vec[i] = *dofs[i];
    else
      this->dof_values_vec[i].zero();
    
    if (dof_sens[i] != NULL)
      this->dof_value_sensitivity_vec[i] = *dof_sens[i];
    else
      this->dof_value_sensitivity_vec[i].zero();
    }
}



inline
void 
FESystemElem::FESystemElemBase::setElementLoads
(std::map<unsigned int, const Loads::VolumeLoadCombination*> *vol_loads,
 std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> > *surf_loads,
 std::map<unsigned int, std::map<unsigned int, const Loads::NodalLoadCombination*> > *nod_loads,
 std::map<unsigned int, const Loads::VolumeLoadCombination*> *vol_load_sens,
 std::map<unsigned int, std::map<unsigned int, const Loads::SurfaceLoadCombination*> > *surf_load_sens,
 std::map<unsigned int, std::map<unsigned int, const Loads::NodalLoadCombination*> > *nod_load_sens)
{
  Assert(vol_loads != NULL, ExcInternalError());
  Assert(surf_loads != NULL, ExcInternalError());
  Assert(nod_loads != NULL, ExcInternalError());

  this->volume_loads = vol_loads;
  this->surface_loads = surf_loads;
  this->nodal_loads = nod_loads;

  if (vol_load_sens != NULL || surf_load_sens != NULL || nod_load_sens != NULL)
    {
      Assert(vol_load_sens != NULL, ExcInternalError());
      Assert(surf_load_sens != NULL, ExcInternalError());
      Assert(nod_load_sens != NULL, ExcInternalError());
      
      this->volume_load_sens = vol_load_sens;
      this->surface_load_sens = surf_load_sens;
      this->nodal_load_sens = nod_load_sens;
    }
}

#endif // __fesystem_fesystem_elem_h__
