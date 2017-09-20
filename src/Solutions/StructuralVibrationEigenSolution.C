// $Id: StructuralVibrationEigenSolution.C,v 1.5.6.4 2008-08-25 04:39:11 manav Exp $

// FESystem includes
#include "Solutions/StructuralVibrationEigenSolution.h"
#include "FESystem/FESystemNumbers.h"
#include "Utilities/InputOutputUtility.h"
#include "Discipline/StructuralAnalysis.h"
#include "Discipline/ThermalAnalysis.h"
#include "Solvers/FESystemSolverBase.h"
#include "DesignData/DesignDatabase.h"
#include "Database/DataInfo.h"


Solution::StructuralVibrationEigenSolution::StructuralVibrationEigenSolution():
Solution::SolutionBase(Solution::STRUCTURAL_VIBRATION_EIGEN_SOLUTION::num())
{
  // add the participating discipline to this solution
  this->addParticipatingDiscipline(Discipline::STRUCTURAL_DISCIPLINE::num());
}



Solution::StructuralVibrationEigenSolution::~StructuralVibrationEigenSolution()
{
  
}




unsigned int 
Solution::StructuralVibrationEigenSolution::getDisciplineTransientNatureEnumID(const unsigned int discipline) const 
{
  unsigned int enum_ID = 0;

  switch (discipline)
    {
    case STRUCTURAL_DISCIPLINE_ENUM_ID:
      enum_ID = Discipline::STEADY_STATE_SOLUTION::num();
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
  
  return enum_ID;
}



void Solution::StructuralVibrationEigenSolution::solve()
{
  // solve the linear system if the geometric effects are needed
  if (include_geometric_effects)
    this->linearSolve(Discipline::STRUCTURAL_DISCIPLINE::num());
  
  // solve the eigen problem
  this->eigenSolve(Discipline::STRUCTURAL_DISCIPLINE::num());
}





std::istream&
Solution::StructuralVibrationEigenSolution::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int num = 0;
  
  FESystemIOUtility::readFromInput(input, Solution::STRUCTURAL_VIBRATION_EIGEN_SOLUTION::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  
  tag.clear();
  FESystemIOUtility::readFromInput(input, "GEOMETRIC_EFFECTS", this->include_geometric_effects);
  
  FESystemIOUtility::readFromInput(input, "COUPLED_THERMAL_STRUCTURAL", 
                                   this->coupled_thermoelastic);
  if (this->coupled_thermoelastic)
    {
      this->addParticipatingDiscipline(Discipline::THERMAL_DISCIPLINE::num());
      FESystemIOUtility::readFromInput(input, "INTERPOLATION_CASE_ID", 
                                       this->interpolation_case_ID);
    }
  
  // Read in the list of load cases
  unsigned int case_ID =0;
  FESystemIOUtility::readFromInput(input, "N_ANALYSIS_LOAD_CASES", num);
  for (unsigned int case_it =0; case_it < num; case_it++)
    {
    case_ID = 0;
    input >> case_ID;
		
    this->load_case_IDs.push_back(case_ID);
    }
  
  FESystemIOUtility::readFromInput(input, "N_SOLVER_INFO_ID", num);
  unsigned int solver_class = 0, info_id = 0;
  for (unsigned int i=0; i < num; i++)
    {
    tag.clear();
    input >> tag;
    solver_class = Solver::SolverClassEnum::enumID(tag);
    input >> info_id;
    
    this->addSolverID(solver_class, info_id);
    }
  
  FESystemIOUtility::readFromInput(input, Solution::STRUCTURAL_VIBRATION_EIGEN_SOLUTION::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}



FESystemUtility::AutoPtrVector<FESystemDatabase::DataInfoBase>
Solution::StructuralVibrationEigenSolution::getTimeIndependentNodalSolutionDataInfo
(const unsigned int discipline_enum_ID)
{
  FESystemUtility::AutoPtrVector<FESystemDatabase::DataInfoBase> data_info_vec;
  
  switch (discipline_enum_ID)
    {
    case STRUCTURAL_DISCIPLINE_ENUM_ID:
      {
        std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
        dv_vector(this->fesystem_controller->design_database->getParameters().release());
        
        FESystemDatabase::TimeIndependentDataInfo eig_val_data_info;
                  
        FESystemDatabase::DataInfoBase *data_info  = NULL;
        
        DenseVector<double> eig_value_real;

        // iterate over each load case and DV and add the solutions
        for (unsigned int i=0; i < this->load_case_IDs.size(); i++)
          {
          // first get the number of converged eigen values for this load case
          eig_val_data_info.clear();
          eig_val_data_info.setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
          eig_val_data_info.setName("EigenValuesReal");
          eig_val_data_info.setLoadCase(this->load_case_IDs[i]);
          this->fesystem_controller->global_data_storage->fillVector(eig_val_data_info, eig_value_real);

          for (unsigned int j=0; j < eig_value_real.size(); j++)
            {
            data_info = new FESystemDatabase::TimeIndependentEigenDataInfo;
            
            data_info->setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
            data_info->setName("EigenVector_Real");
            data_info->setLoadCase(this->load_case_IDs[i]);
            dynamic_cast<FESystemDatabase::TimeIndependentEigenDataInfo*>(data_info)->setModeInfo(j+1);
            
            data_info_vec.push_back(data_info);
            data_info = NULL;
            }
          }
        
        for (unsigned int i=0; i < this->load_case_IDs.size(); i++)
          for (unsigned int j=0; j < dv_vector->size(); j++)
            {
            // first get the number of converged eigen values for this load case
            eig_val_data_info.clear();
            eig_val_data_info.setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
            eig_val_data_info.setName("EigenValuesReal");
            eig_val_data_info.setLoadCase(this->load_case_IDs[i]);
            eig_val_data_info.setDVID((*dv_vector)[j]->getID());
            this->fesystem_controller->global_data_storage->fillVector(eig_val_data_info,
                                                                       eig_value_real);
            
            for (unsigned int k=0; k < eig_value_real.size(); k++)
              {
              data_info = new FESystemDatabase::TimeIndependentEigenDataInfo;
              
              data_info->setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
              data_info->setName("EigenVector_Real");
              data_info->setLoadCase(this->load_case_IDs[i]);
              data_info->setDVID((*dv_vector)[j]->getID());
              dynamic_cast<FESystemDatabase::TimeIndependentEigenDataInfo*>(data_info)->setModeInfo(k+1);
              
              data_info_vec.push_back(data_info);
              data_info = NULL;
              }
            }
      }
        break;
        
        default:
          Assert(false, ExcInternalError());
    }
      
      return data_info_vec;
}


FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet>
Solution::StructuralVibrationEigenSolution::getTimeDependentNodalSolutionDataInfo
(const unsigned int discipline_enum_ID)
{
  // unused parameter
  (void) discipline_enum_ID;
  
  // this solution does not have any transients. Hence, return an empty vector
  FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet> vec;
  return vec;
}



