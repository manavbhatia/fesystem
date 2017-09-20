// $Id: ThermalAnalysis.h,v 1.11.6.6 2008-06-03 05:19:33 manav Exp $

#ifndef __fesystem_thermal_analysis_h__
#define __fesystem_thermal_analysis_h__

// C++ includes 
#include <set>

// FESystem includes 
#include "Discipline/AnalysisDisciplineBase.h"
#include "Utilities/NameEnumHandler.h"
#include "Loads/load.h"


// libMesh includes
#include "libmesh_common.h"

// Forward Decleration
namespace FESystem
{
  class FESystemController;
  class FESystemElem;
}

class RadiationCavityAnalysis;


#ifndef THERMAL_DISCIPLINE_ENUM_ID
#define THERMAL_DISCIPLINE_ENUM_ID 2
#else
#error
#endif

#ifndef THERMAL_DISCIPLINE_ENUM_NAME
#define THERMAL_DISCIPLINE_ENUM_NAME "THERMAL_DISCIPLINE"
#else
#error
#endif



#ifndef TEMPERATURE_BOUNDARY_CONDITION_ENUM_ID
#define TEMPERATURE_BOUNDARY_CONDITION_ENUM_ID 1
#else
#error
#endif

#ifndef TEMPERATURE_BOUNDARY_CONDITION_ENUM_NAME
#define TEMPERATURE_BOUNDARY_CONDITION_ENUM_NAME "TEMPERATURE_BOUNDARY_CONDITION"
#else
#error
#endif


#ifndef NODAL_HEAT_LOAD_ENUM_ID
#define NODAL_HEAT_LOAD_ENUM_ID 2
#else
#error
#endif

#ifndef NODAL_HEAT_LOAD_ENUM_NAME
#define NODAL_HEAT_LOAD_ENUM_NAME "NODAL_HEAT_LOAD"
#else
#error
#endif


#ifndef VOLUME_HEAT_LOAD_ENUM_ID
#define VOLUME_HEAT_LOAD_ENUM_ID 3
#else
#error
#endif

#ifndef VOLUME_HEAT_LOAD_ENUM_NAME
#define VOLUME_HEAT_LOAD_ENUM_NAME "VOLUME_HEAT_LOAD"
#else
#error
#endif


#ifndef SURFACE_HEAT_LOAD_ENUM_ID
#define SURFACE_HEAT_LOAD_ENUM_ID 4
#else
#error
#endif

#ifndef SURFACE_HEAT_LOAD_ENUM_NAME
#define SURFACE_HEAT_LOAD_ENUM_NAME "SURFACE_HEAT_LOAD"
#else
#error
#endif


#ifndef SURFACE_RADIATION_HEAT_LOAD_ENUM_ID
#define SURFACE_RADIATION_HEAT_LOAD_ENUM_ID 5
#else
#error
#endif

#ifndef SURFACE_RADIATION_HEAT_LOAD_ENUM_NAME
#define SURFACE_RADIATION_HEAT_LOAD_ENUM_NAME "SURFACE_RADIATION_HEAT_LOAD"
#else
#error
#endif


#ifndef SURFACE_CONVECTION_HEAT_LOAD_ENUM_ID
#define SURFACE_CONVECTION_HEAT_LOAD_ENUM_ID 6
#else
#error
#endif

#ifndef SURFACE_CONVECTION_HEAT_LOAD_ENUM_NAME
#define SURFACE_CONVECTION_HEAT_LOAD_ENUM_NAME "SURFACE_CONVECTION_HEAT_LOAD"
#else
#error
#endif




DeclareEnumName(TEMPERATURE_BOUNDARY_CONDITION, LoadNameEnum,
                TEMPERATURE_BOUNDARY_CONDITION_ENUM_ID,
                TEMPERATURE_BOUNDARY_CONDITION_ENUM_NAME);


DeclareEnumName(NODAL_HEAT_LOAD, LoadNameEnum,
		NODAL_HEAT_LOAD_ENUM_ID,
		NODAL_HEAT_LOAD_ENUM_NAME);
  
DeclareEnumName(VOLUME_HEAT_LOAD, LoadNameEnum,
		VOLUME_HEAT_LOAD_ENUM_ID,
		VOLUME_HEAT_LOAD_ENUM_NAME);
  
DeclareEnumName(SURFACE_HEAT_LOAD, LoadNameEnum,
		SURFACE_HEAT_LOAD_ENUM_ID,
		SURFACE_HEAT_LOAD_ENUM_NAME);
  
DeclareEnumName(SURFACE_RADIATION_HEAT_LOAD, LoadNameEnum,
		SURFACE_RADIATION_HEAT_LOAD_ENUM_ID,
		SURFACE_RADIATION_HEAT_LOAD_ENUM_NAME);

DeclareEnumName(SURFACE_CONVECTION_HEAT_LOAD, LoadNameEnum,
		SURFACE_CONVECTION_HEAT_LOAD_ENUM_ID,
		SURFACE_CONVECTION_HEAT_LOAD_ENUM_NAME);
  


namespace Discipline
{
  
  class ThermalDisciplineInfo;
  
  DeclareEnumName(THERMAL_DISCIPLINE, Discipline::AnalysisDisciplineEnum,
                  THERMAL_DISCIPLINE_ENUM_ID,
                  THERMAL_DISCIPLINE_ENUM_NAME);
  
  
  
 
  
  /// This class provides methods to calculate the matrices and vectors related to the 
  /// thermal discipline. 
  class ThermalAnalysis: public Discipline::AnalysisDisciplineBase
    {
public:
      
      /// constructor. It takes two arguements,
      /// the arguement is a pointer to the AnalysisDriver that owns this object
      ThermalAnalysis(FESystem::FESystemController& controller,
                      const Discipline::ThermalDisciplineInfo& info);
      
      /// destructor
      ~ThermalAnalysis();
      
      /// returns if same system matrix can be used for different load cases
      virtual inline bool sameSystemMatrixForLoadCases();
      
      /// @returns the eigenproblem kind for the discipline
      virtual unsigned int getEigenProblemKindEnumID();
      
      /// @returns the order of the transient nature of the problem
      virtual unsigned int getTransientSystemOrder();

      /// returns the enum ID of the boundary condition for this disicipline
      virtual std::auto_ptr<Loads::DirichletBoundaryConditionDataInfo> 
        getBoundaryConditionLoadInfo() const;
      
      /// @returns true if the specified matrix has exists for the discipline
      virtual bool disciplineHasMatrix(const unsigned int qty_enum_ID) const;

      /// @returns a FESystemDatabase::DataInfo pointer, which defines the quantity 
      /// passed in the parameter, for the current analysis
      virtual std::auto_ptr<FESystemDatabase::DataInfoBase> getDataInfoForQty
	(const unsigned int qty_enum_ID);

protected:
        
        /// function to add analysis variables
        virtual void addAnalysisVariables();

        
      /// function will fill the input vector with the element quantities required for the calculation
      /// of the specified global quantity. This is an abstract function and will be implemented in the
      /// derived class
      /// @param global quantitiy
      /// @param element pointer
      /// @param vector in which the quantities will be added
      void getElemQty(FESystemElem::FESystemElemBase* elem,
                      std::map<unsigned int, double>& real_elem_data,
                      std::map<unsigned int, DenseMatrix<double> >& mat_elem_data,
                      std::map<unsigned int, DenseVector<double> >& vec_elem_data,
                      const bool sensitivity, 
                      const unsigned int DV_ID);
      
      
      /// function to calculate the conductance matrix. First arguement is the matrix to store 
      ///  the result in, and the second arguement is a bool to tell the function if sensitivity
      /// is desired
      virtual bool calculateTransientMatrix(const unsigned int order, bool = false);

      /// function to calculate the conductance matrix. First arguement is the matrix to store 
      ///  the result in, and the second arguement is a bool to tell the function if sensitivity
      /// is desired
      virtual bool calculateK(bool = false);
      
      /// function to calculate the Jacobian matrix. First arguement is the matrix to store 
      ///  the result in
      virtual bool calculateJac();
      
      /// function to calculate the A matrix for eigen problem analysis
      virtual bool calculateEigenProblemAMatrix(bool sensitivity = false);
      
      /// function to calculate the B matrix for eigen problem analysis
      virtual bool calculateEigenProblemBMatrix(bool sensitivity = false);

      /// function to calculate the force vector. First arguement is the vector to store 
      ///  the result in, and the second arguement is a bool to tell the function if sensitivity
      /// is desired
      virtual bool calculateF(bool = false);
      
      /// this method will perform some discipline specific operations once the 
      /// global quantities have been calculated
      virtual void performAdditionalOperationsOnGlobalQuantities
        (std::map<unsigned int, SparseMatrix<double>* >& matrix_map,
         std::map<unsigned int, NumericVector<double>* >& vector_map);
      
      /// this method fills the vector with the solution vector necessary for filling the quantity
      /// for the current analysis. It checks the current status of the discipline, and 
      /// apropriately loads the solution vector
      virtual const NumericVector<double>& getSolutionVector(const unsigned int transient_order,
                                                             const bool sensitivity);

      /// vector of radiation cavity analysis 
      std::vector <RadiationCavityAnalysis*> radiation_cavity_analysis_vector;
      
    };
}


inline
unsigned int
Discipline::ThermalAnalysis::getTransientSystemOrder()
{
  return 1;
}



inline 
bool
Discipline::ThermalAnalysis::sameSystemMatrixForLoadCases()
{
  switch (this->analysis_type_enum_ID)
    {
    case LINEAR_ANALYSIS_ENUM_ID:
      return false;
      break;
      
    case NONLINEAR_ANALYSIS_ENUM_ID:
      return false;
      break;
      
    default:
      abort();
      break;
    }
}



#endif // __fesystem_thermal_analysis_h__
