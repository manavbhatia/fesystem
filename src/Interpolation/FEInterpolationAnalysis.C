// $Id: FEInterpolationAnalysis.C,v 1.5 2006-09-05 20:41:48 manav Exp $

// C++ includes
#include <sstream>

// FESystem includes
#include "Interpolation/FEInterpolationAnalysis.h"
#include "FESystem/FESystemController.h"
#include "FESystem/AnalysisCase.h"
#include "Mesh/MeshList.h"
#include "Mesh/FEMeshData.h"
#include "Database/GlobalDataStorage.h"
#include "Loads/LoadDatabase.h"
#include "AnalysisDriver/AnalysisDriver.h"


#include "Interpolation/FEInterpolationElem.h"
#include "Interpolation/InterpolationBar2.h"
#include "Interpolation/InterpolationQuad4.h"
#include "Interpolation/InterpolationHex8.h"

// libMesh includes
#include "numerics/sparse_matrix.h"
#include "numerics/numeric_vector.h"
#include "numerics/dense_matrix.h"
#include "numerics/dense_vector.h"




FEInterpolationAnalysis::FEInterpolationAnalysis(FESystem::FESystemController& controller,
                                                 Discipline::InterpolationDisciplineInfo& info):
Discipline::AnalysisDisciplineBase(controller, info)
{

}



FEInterpolationAnalysis::~FEInterpolationAnalysis()
{
	
}


void FEInterpolationAnalysis::addVariablesForAnalysis()
{
  FEType fetype(FIRST, LAGRANGE);
  // add a temperature variable to the analysis driver
  this->addVariable("Variable", fetype);
}



void FEInterpolationAnalysis::initMatricesAndVectors()
{
  // the dof map used here will be the one for the secondary mesh, which is the one
  // on which the variable has to be interpolated

  // hence, from the analysis case, get the interpolaton analysis mesh ID, and 
  // from the analysis driver, get the analysis dof map.

  
  Assert(this->analysis_driver != NULL,
         ExcEmptyObject());
  
  unsigned int n_dofs = this->dof_map->n_dofs();
	

  // initialize the system matrix to the dof_map connectivity
  SparseMatrix<double>& system_matrix = this->analysis_driver->getMatrix("SystemMatrix");
	
  system_matrix.attach_dof_map(*(this->dof_map));
  system_matrix.init();
	
  NumericVector<double>& sol_vec = this->analysis_driver->getVector("Solution");
  NumericVector<double>& rhs_vec = this->analysis_driver->getVector("RHS");
	
  sol_vec.init(n_dofs);
  rhs_vec.init(n_dofs);

}






std::string 
FEInterpolationAnalysis::getStringNameForQty(const unsigned int qty_enum_ID,
                                           bool sensitivity)
{
  std::string name;
	
  switch (qty_enum_ID)
    {
    case SYSTEM_MATRIX_ENUM_ID:
      name = "Interpolation_K_matrix";
      break;
			
    case FORCE_VECTOR_ENUM_ID:
      name = "Interpolation_F_vector";
      break;
			
    case JACOBIAN_MATRIX_ENUM_ID:
    default:
      abort();
      break;
    }

  if (!sensitivity)
    {
    return name;
    }
  if (sensitivity)
    {
    // if sensitivity is being calculated, the DV ID will be appended to the name of the quantity
    std::string sensitivity_name = "d";
      sensitivity_name += name;
      sensitivity_name += "_dDV_";
      std::ostringstream dv_ID;
      dv_ID << this->analysis_driver->getCurrentDesignParameter().getID();
      sensitivity_name += dv_ID.str();
      return sensitivity_name;
    }
}






void 
FEInterpolationAnalysis::calculateK(SparseMatrix<double>& matrix, bool sensitivity)
{
  // get the name by which the quantity will be referred
  std::string name = this->getStringNameForQty(Driver::SYSTEM_MATRIX::num(), sensitivity);
  
  // check if the K matrix exists in the database
  // this matrix is independent of load case
  bool K_matrix_present = 
    this->fesystem_controller.global_data_storage->checkIfMatrixExists(name);
	
  // if not, calculate and save it
  if (! K_matrix_present)
    {
    unsigned int qty_enum_ID;
		
    if (!sensitivity)	
      qty_enum_ID = Driver::SYSTEM_MATRIX::num();
    else
      qty_enum_ID = Driver::SYSTEM_MATRIX_SENSITIVITY::num();
		
      // calculate the matrix, save it and return
      this->calculateGlobalQty(matrix, qty_enum_ID);
		
      this->fesystem_controller.global_data_storage->storeMatrix(name, matrix);
    }
  else 
    {
    // get the matrix from the database, and return
    this->fesystem_controller.global_data_storage->fillMatrix(name, matrix);
    }
}





void FEInterpolationAnalysis::calculateJac(SparseMatrix<double>& matrix)
{
  // nothing to be done here since the analysis is not nonlinear
  abort();
}





void 
FEInterpolationAnalysis::calculateF(NumericVector<double>& vector, bool sensitivity)
{
  unsigned int load_case = this->analysis_driver->getCurrentLoadCase();
	
  // get the name by which the quantity will be referred
  std::string name = this->getStringNameForQty(Driver::FORCE_VECTOR::num(), sensitivity);
	
  // check if the vector exists in the database
  bool Force_vector_present = 
    this->fesystem_controller.global_data_storage->checkIfVectorExists(load_case, name);
	
  // if not, calculate and save it
  if (! Force_vector_present)
    {
    unsigned int qty_enum_ID;
		
    if (!sensitivity)
      qty_enum_ID = Driver::FORCE_VECTOR::num();
    else
      qty_enum_ID = Driver::FORCE_VECTOR_SENSITIVITY::num();
    
		
    // calculate the matrix, save it and return
    this->calculateGlobalQty(vector, qty_enum_ID);

    this->fesystem_controller.global_data_storage->storeVector(load_case,name,vector);
    }
  else 
    {
      // get the matrix from the database, and return
      this->fesystem_controller.global_data_storage->fillVector(load_case,name,vector);
    }
}






//std::auto_ptr<FESystemElem>
//FEInterpolationAnalysis::createElem(const unsigned int elemID, 
//				  const Elem* elem)
//{
//  std::auto_ptr<FEInterpolationElem> interpolation_elem(NULL);
//
//  // get the element tag for this element
//  const std::string& elem_tag = this->mesh_data->getElemTag(elem);
//	
//  switch(elem->dim())
//    {
//    case 1:
//      {
//	if (elem_tag == "INTERPOLATION_BAR2")
//	  {
//	    interpolation_elem.reset(new InterpolationBar2(elemID,elem,this->analysis_driver));
//	  }
//	else
//	  abort();
//			
//      }
//      break;
//			
//    case 2:
//      {
//	if (elem_tag == "INTERPOLATION_QUAD4")
//	  {
//	    interpolation_elem.reset(new InterpolationQuad4(elemID,elem,this->analysis_driver));
//	  }
//	else 
//	  abort();
//      }
//      break;
//			
//    case 3:
//      {
//	if (elem_tag == "INTERPOLATION_HEX8")
//	  {
//	    interpolation_elem.reset(new InterpolationHex8(elemID,elem,this->analysis_driver));
//	  }
//	else 
//	  abort();
//      }
//      break;
//			
//    default:
//      abort();
//      break;
//    }
//	
//  // now create an auto_ptr with a dynamic case to FESystemElem
//	
//  std::auto_ptr<FESystemElem> return_qty(dynamic_cast<FESystemElem*>(interpolation_elem.release()));
//	
//  return return_qty;
//}



void
  FEInterpolationAnalysis::getElemQty(const unsigned int qty,
                                    FESystemElem::FESystemElemBase* elem,
                                    std::vector<DenseMatrix<double>* >* elem_data)
{
  // make sure that valid pointer have been given 
  assert (elem != NULL);
  assert (elem_data != NULL);
	
  static DenseMatrix<double> matrix;
  
  elem_data->clear();
	
  switch(qty)
    {
    case SYSTEM_MATRIX_ENUM_ID:
      {
        // get K matrix
        elem->getElementAssembledQty(FE_INTERPOLATION_ELEM_K_MATRIX::num(),
                                     &matrix);
        elem_data->push_back(&matrix);
      }
      break;
			
			
    case SYSTEM_MATRIX_SENSITIVITY_ENUM_ID:
      {
        // get K matrix sensitivity
        elem->getElementAssembledQtySensitivity(FE_INTERPOLATION_ELEM_K_MATRIX::num(),
                                                &matrix);
        elem_data->push_back(&matrix);
      }
      break;
			
			
    case JACOBIAN_MATRIX_ENUM_ID:			
    case FORCE_VECTOR_ENUM_ID:	
    case FORCE_VECTOR_SENSITIVITY_ENUM_ID:	
    default:
      abort();
      break;
    }
}



void FEInterpolationAnalysis::getElemQty(const unsigned int qty,
                                       FESystemElem::FESystemElemBase* elem,
                                       std::vector<DenseVector<double>* >* elem_data)
{
  // make sure that valid pointer have been given 
  assert (elem != NULL);
  assert (elem_data != NULL);

  static DenseVector<double> vector;

  elem_data->clear();
	
  switch(qty)
    {
		
    case FORCE_VECTOR_ENUM_ID:
      {
        // get F
        elem->getElementAssembledQty(FE_INTERPOLATION_ELEM_F_VECTOR::num(),
                                     &vector);
        elem_data->push_back(&vector);
      }
      break;
			
    case FORCE_VECTOR_SENSITIVITY_ENUM_ID:
      {
        // get F  sensitivity
        elem->getElementAssembledQtySensitivity(FE_INTERPOLATION_ELEM_F_VECTOR::num(),
                                                &vector);
        elem_data->push_back(&vector);
      }
      break;
      
      
    case SYSTEM_MATRIX_ENUM_ID:
    case JACOBIAN_MATRIX_ENUM_ID:
    case SYSTEM_MATRIX_SENSITIVITY_ENUM_ID:
    default:
      abort();
      break;
    }
	
}

