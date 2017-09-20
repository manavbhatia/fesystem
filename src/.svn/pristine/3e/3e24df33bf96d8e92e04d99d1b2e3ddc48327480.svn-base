// $Id $

// FESystem includes
#include "Solutions/LinearizedBucklingEigenSolution.h"
#include "FESystem/AnalysisCase.h"
#include "FESystem/FESystemNumbers.h"
#include "Utilities/InputOutputUtility.h"
#include "Discipline/StructuralAnalysis.h"
#include "Discipline/ThermalAnalysis.h"
#include "Discipline/DisciplineInfo.h"
#include "Solvers/FESystemSolverBase.h"
#include "DesignData/DesignDatabase.h"
#include "Database/DataInfo.h"
#include "Interpolation/DirectInterpolation.h"



Solution::LinearizedBucklingEigenSolution::LinearizedBucklingEigenSolution():
Solution::SolutionBase(Solution::LINEARIZED_BUCKLING_EIGEN_SOLUTION::num()),
coupled_thermoelastic(false),
interpolation_case_ID(FESystemNumbers::InvalidID)
{
  // add the participating discipline to this solution
  this->addParticipatingDiscipline(Discipline::STRUCTURAL_DISCIPLINE::num());
}



Solution::LinearizedBucklingEigenSolution::~LinearizedBucklingEigenSolution()
{

}



unsigned int 
Solution::LinearizedBucklingEigenSolution::getDisciplineTransientNatureEnumID(const unsigned int discipline) const 
{
  unsigned int enum_ID = 0;

  switch (discipline)
    {
    case STRUCTURAL_DISCIPLINE_ENUM_ID:
    case THERMAL_DISCIPLINE_ENUM_ID:
      enum_ID = Discipline::STEADY_STATE_SOLUTION::num();
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
  
  return enum_ID;
}




void Solution::LinearizedBucklingEigenSolution::solve()
{
  // if a coupled thermal-structural response is requested, then
  // first perform a thermal analysis, import the loads to structural and
  // then perform the structural solution
  if (this->coupled_thermoelastic)
    {
    const Discipline::DisciplineInfo& thermal_info = 
    this->fesystem_controller->analysis_case->getAnalysisDisciplineInfo
    (Discipline::THERMAL_DISCIPLINE::num());
    
    switch (thermal_info.getAnalysisTypeEnumID())
      {
      case LINEAR_ANALYSIS_ENUM_ID:
        this->linearSolve(Discipline::THERMAL_DISCIPLINE::num());
        break;
        
      case NONLINEAR_ANALYSIS_ENUM_ID:
        this->nonlinearSolve(Discipline::THERMAL_DISCIPLINE::num());
        break;
        
      default:
        Assert(false, ExcInternalError());
        break;
      }
    
    // get the interpolation case from the analysis case
    const InterpolationCase& interpolation_case = 
      this->fesystem_controller->analysis_case->getInterpolationCase(this->interpolation_case_ID);
    std::auto_ptr<DirectInterpolation> interpolation
      (new DirectInterpolation(*(this->fesystem_controller), interpolation_case));
    interpolation->interpolate();
    }

  this->linearSolve(Discipline::STRUCTURAL_DISCIPLINE::num());
  
  // solve the eigen problem
  this->eigenSolve(Discipline::STRUCTURAL_DISCIPLINE::num());
}





std::istream&
Solution::LinearizedBucklingEigenSolution::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int num = 0;
  
  FESystemIOUtility::readFromInput(input, Solution::LINEARIZED_BUCKLING_EIGEN_SOLUTION::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "ID", this->ID);
  
  
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
  
  FESystemIOUtility::readFromInput(input, Solution::LINEARIZED_BUCKLING_EIGEN_SOLUTION::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}




FESystemUtility::AutoPtrVector<FESystemDatabase::DataInfoBase>
Solution::LinearizedBucklingEigenSolution::getTimeIndependentNodalSolutionDataInfo
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
          // set the data info for the static solution
          data_info = new FESystemDatabase::TimeIndependentDataInfo;
          data_info->setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
          data_info->setName("Solution");
          data_info->setLoadCase(this->load_case_IDs[i]);
          
          data_info_vec.push_back(data_info);
          data_info = NULL;
          
          // now the eigen vectors
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
            // sensitivity of static solutions
            data_info = new FESystemDatabase::TimeIndependentDataInfo;
            
            data_info->setDisciplineEnumID(Discipline::STRUCTURAL_DISCIPLINE::num());
            data_info->setName("Solution");
            data_info->setLoadCase(this->load_case_IDs[i]);
            data_info->setDVID((*dv_vector)[j]->getID());
            data_info_vec.push_back(data_info);
            data_info = NULL;
            
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
        
        case THERMAL_DISCIPLINE_ENUM_ID:
          {
            std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
            dv_vector(this->fesystem_controller->design_database->getParameters().release());
            
            data_info_vec.resize(this->load_case_IDs.size() * (1 + dv_vector->size()));
            
            FESystemDatabase::DataInfoBase *data_info  = NULL;
            
            unsigned int location_id = 0;
            
            // iterate over each load case and DV and add the solutions
            for (unsigned int i=0; i < this->load_case_IDs.size(); i++)
              {
              data_info = new FESystemDatabase::TimeIndependentDataInfo;
              
              data_info->setDisciplineEnumID(discipline_enum_ID);
              data_info->setName("Solution");
              data_info->setLoadCase(this->load_case_IDs[i]);
              
              data_info_vec.reset(location_id, data_info);
              location_id++;
              data_info = NULL;
              }
            
            for (unsigned int i=0; i < this->load_case_IDs.size(); i++)
              for (unsigned int j=0; j < dv_vector->size(); j++)
                {
                data_info = new FESystemDatabase::TimeIndependentDataInfo;
                
                data_info->setDisciplineEnumID(discipline_enum_ID);
                data_info->setName("Solution");
                data_info->setLoadCase(this->load_case_IDs[i]);
                data_info->setDVID((*dv_vector)[j]->getID());
                
                data_info_vec.reset(location_id, data_info);
                location_id++;
                data_info = NULL;
                }
          }
            break;

            
            default:
              Assert(false, ExcInternalError());
    }
          
          return data_info_vec;
}



FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet>
Solution::LinearizedBucklingEigenSolution::getTimeDependentNodalSolutionDataInfo
(const unsigned int discipline_enum_ID)
{
  // unused parameter
  (void) discipline_enum_ID;

  // this solution does not have any transients. Hence, return an empty vector
  FESystemUtility::AutoPtrVector<FESystemDatabase::TransientDataInfoSet> vec;
  return vec;
}

