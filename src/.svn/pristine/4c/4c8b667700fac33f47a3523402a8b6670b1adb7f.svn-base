// $Id: LoadDatabase.C,v 1.10.6.2 2007-06-13 14:58:05 manav Exp $

// C++ inlcudes

// FESystem includes
#include "Loads/LoadDatabase.h"
#include "Utilities/InputOutputUtility.h"
#include "FESystem/FESystemNumbers.h"
#include "Loads/load.h"
#include "Loads/LoadSet.h"
#include "Loads/LoadCase.h"
#include "Loads/LoadCombination.h"
#include "Loads/LoadDataInfo.h"


LoadDatabase::LoadDatabase(FESystem::FESystemController& controller):
fesystem_controller(controller)
{
  
}



LoadDatabase::~LoadDatabase()
{
  this->clear();
}




void LoadDatabase::clear()
{
  // iterate over all elements in the map, and delete the pointers
  std::map<unsigned int, LoadSetBase*>::iterator load_set_it = this->load_set_map.begin();
  std::map<unsigned int, LoadSetBase*>::const_iterator load_set_end = this->load_set_map.end();
  
  for (; load_set_it != load_set_end; load_set_it++)
    {
      delete load_set_it->second; 
      load_set_it->second = NULL; 
    }
  this->load_set_map.clear();
  
  
  std::map<unsigned int, LoadCaseBase*>::iterator load_case_it = 
  this->load_case_map.begin();
  std::map<unsigned int, LoadCaseBase*>::const_iterator load_case_end = 
  this->load_case_map.end();
  
  for (; load_case_it != load_case_end; load_case_it++)
    {
      delete load_case_it->second; 
      load_case_it->second = NULL;
    }
  this->load_case_map.clear();
  
}	




void LoadDatabase::addLoadSet(LoadSetBase* load_set_ptr)
{
  // make sure that the pointer is not NULL and that ID does not already exist
  assert (load_set_ptr != NULL);
  
  unsigned int load_set_ID = load_set_ptr->getLoadSetID();
  
  assert (this->load_set_map.find(load_set_ID) == 
          this->load_set_map.end());
  
  std::pair<std::map<unsigned int, LoadSetBase*>::iterator, bool > insert_return_value = 
  this->load_set_map.insert(std::map<unsigned int, LoadSetBase*>::value_type(load_set_ID, 
                                                                             load_set_ptr));
	
  assert (insert_return_value.second == true);	
}





void LoadDatabase::addLoadCase( LoadCaseBase* load_case_ptr)
{
  // make sure that the pointer is not NULL and that ID does not already exist
  assert (load_case_ptr != NULL);
  
  unsigned int load_case_ID = load_case_ptr->getLoadCaseID();
  
  // make sure that this load case does not already exist
  AssertThrow(this->load_case_map.count(load_case_ID) == 0, 
              ExcInternalError());
  
  std::pair<std::map<unsigned int, LoadCaseBase*>::iterator, bool > insert_return_value = 
  this->load_case_map.insert(std::map<unsigned int, LoadCaseBase*>::value_type
                             (load_case_ID, load_case_ptr));
  
  assert (insert_return_value.second);
}





std::istream& LoadDatabase::readFromInputStream(std::istream& input)
{
  unsigned int num_load_set=0, num_load_case=0, enum_ID=0;
  std::string tag;
	
  // the data should begin with tag BEGIN_LOAD_SET_DATA 
  // which should be followed by an 
  // integer specifying the number of load sets
  FESystemIOUtility::readFromInput(input, "LOAD_DATABASE");
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "LOAD_SET_DATA");
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "N_LOAD_SETS", num_load_set);
  
  std::auto_ptr<LoadSetBase> load_set;
  
  for (unsigned int load_set_incr=0; load_set_incr < num_load_set; load_set_incr++)
    {
      tag.clear();
      FESystemIOUtility::peekFromInput(input, tag);
      enum_ID = LoadSetKindEnum::enumID(tag);
      
      load_set.reset(createLoadSet(enum_ID).release());
      
      load_set->readFromInputStream(input);
      
      this->addLoadSet(load_set.release());
    }
  
  // the load set data should end with the END tag
  FESystemIOUtility::readFromInput(input, "LOAD_SET_DATA");
  FESystemIOUtility::readFromInput(input, "END");
	
	
  // this should be followed by load case data, which will begin with the 
  // tag BEGIN_LOAD_CASE_DATA, followed by an integer specifying the number
  // of load cases
  FESystemIOUtility::readFromInput(input, "LOAD_CASE_DATA");
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  FESystemIOUtility::readFromInput(input, "N_LOAD_CASES", num_load_case);
  
  unsigned int load_case_type_enum_ID=0;
  std::auto_ptr<LoadCaseBase> load_case_ptr;
  
  for (unsigned int load_case_incr=0; load_case_incr < num_load_case; load_case_incr++)
    {
      // check which kind of load case is provided
      tag.clear();
      FESystemIOUtility::peekFromInput(input, tag);
      load_case_type_enum_ID = LoadCaseKindEnum::enumID(tag);
      
      load_case_ptr.reset(createLoadCase(load_case_type_enum_ID, *this).release());
      
      load_case_ptr->readFromInputStream(input);
      
      this->addLoadCase(load_case_ptr.release());
    }
	
	
  // the load case data will end with an END tag
  FESystemIOUtility::readFromInput(input, "LOAD_CASE_DATA");
  FESystemIOUtility::readFromInput(input, "END");
	
  FESystemIOUtility::readFromInput(input, "LOAD_DATABASE");
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}




const LoadCaseBase&
LoadDatabase::getLoadCaseFromID(const unsigned int ID) const
{
  // make sure that the load case with this ID exists
  std::map<unsigned int, LoadCaseBase* >::const_iterator load_case_it = 
  this->load_case_map.find(ID);
  
  assert (load_case_it != this->load_case_map.end());
  
  return *(load_case_it->second);
}






const LoadSetBase& 
LoadDatabase::getLoadSetFromID(const unsigned int ID) const
{	
  // make sure that the load set with this ID exists
  // make sure that the load case with this ID exists
  std::map<unsigned int, LoadSetBase* >::const_iterator load_set_it = 
  this->load_set_map.find(ID);
  
  assert (load_set_it != this->load_set_map.end());
  
  return *(load_set_it->second);
}




unsigned int
LoadDatabase::getAvailableLoadSetID() const
{	
  unsigned int load_set_ID = 0;
  
  if (this->load_set_map.size() == 0)
    load_set_ID = 1;
  else
    load_set_ID = (this->load_set_map.rbegin()->first) + 1;
  
  return load_set_ID;
}






template < >
void 
LoadDatabase::getLoadCombination<Loads::VolumeLoadDataInfo, Loads::VolumeLoadCombination>
(const Loads::VolumeLoadDataInfo& load_info,
 Loads::VolumeLoadCombination& load_combo)
{
  load_combo.clear();
  
  static std::vector<const VolumeLoad*> load_vec;
  load_vec.clear();
  
  static std::vector<std::pair<unsigned int, double> > load_pair;
  load_pair.clear();
  
  const LoadCaseBase& load_case = this->getLoadCaseFromID(load_info.getLoadCaseID());
  
  load_case.getLoadFactorLoadSetPair(load_info, load_pair);
  
  // iterate over the load pairs, and add the loads to the load vector
  for (unsigned int i=0; i < load_pair.size(); i++)
    {
      const VolumeLoadSet& load_set = 
      dynamic_cast<const VolumeLoadSet&>(this->getLoadSetFromID(load_pair[i].first));
      
      // make sure that the loads are of the same kind
      AssertThrow(load_info.getLoadNameEnumID() == load_set.getLoadNameEnumID(), ExcInternalError());
      
      load_vec.clear();
      // finally, ask the load set to return the load vector
      load_set.getLoadsForElement(load_info.getElemID(), load_vec);
      
      for (unsigned int j=0; j < load_vec.size(); j++)
        load_combo.addLoad(load_pair[i].second, load_vec[j]);
    }
}




template < >
void 
LoadDatabase::getLoadCombination<Loads::SurfaceLoadDataInfo, Loads::SurfaceLoadCombination>
(const Loads::SurfaceLoadDataInfo& load_info,
 Loads::SurfaceLoadCombination& load_combo)
{
  load_combo.clear();
  
  static std::vector<const SurfaceLoad*> load_vec;
  load_vec.clear();
  
  static std::vector<std::pair<unsigned int, double> > load_pair;
  load_pair.clear();
  
  const LoadCaseBase& load_case = this->getLoadCaseFromID(load_info.getLoadCaseID());
  
  load_case.getLoadFactorLoadSetPair(load_info, load_pair);
  
  // iterate over the load pairs, and add the loads to the load vector
  for (unsigned int i=0; i < load_pair.size(); i++)
    {
      const SurfaceLoadSet& load_set = 
      dynamic_cast<const SurfaceLoadSet&>(this->getLoadSetFromID(load_pair[i].first));
      
      // make sure that the loads are of the same kind
      AssertThrow(load_info.getLoadNameEnumID() == load_set.getLoadNameEnumID(), ExcInternalError());
      
      load_vec.clear();
      // finally, ask the load set to return the load vector
      load_set.getLoadsForElementSide(load_info.getElemID(), load_info.getSurfaceID(),
                                      load_vec);
      
      for (unsigned int j=0; j < load_vec.size(); j++)
        load_combo.addLoad(load_pair[i].second, load_vec[j]);
    }
}





template < >
void 
LoadDatabase::getLoadCombination<Loads::NodalLoadDataInfo, Loads::NodalLoadCombination>
(const Loads::NodalLoadDataInfo& load_info,
 Loads::NodalLoadCombination& load_combo)
{
  load_combo.clear();
  
  static std::vector<const NodalLoad*> load_vec;
  load_vec.clear();
  
  static std::vector<std::pair<unsigned int, double> > load_pair;
  load_pair.clear();
  
  const LoadCaseBase& load_case = this->getLoadCaseFromID(load_info.getLoadCaseID());
  
  load_case.getLoadFactorLoadSetPair(load_info, load_pair);
  
  // iterate over the load pairs, and add the loads to the load vector
  for (unsigned int i=0; i < load_pair.size(); i++)
    {
      const NodalLoadSet& load_set = 
      dynamic_cast<const NodalLoadSet&>(this->getLoadSetFromID(load_pair[i].first));
      
      // make sure that the loads are of the same kind
      AssertThrow(load_info.getLoadNameEnumID() == load_set.getLoadNameEnumID(), ExcInternalError());
      
      load_vec.clear();
      // finally, ask the load set to return the load vector
      load_set.getLoadsForNode(load_info.getNodeID(), load_vec);
      
      for (unsigned int j=0; j < load_vec.size(); j++)
        load_combo.addLoad(load_pair[i].second, load_vec[j]);
    }
}





template <>
FESystemUtility::AutoPtrVector<Loads::NodalLoadCombination>
LoadDatabase::getAllLoadCombinations<Loads::NodalLoadDataInfo , Loads::NodalLoadCombination>
(const Loads::NodalLoadDataInfo& load_info)
{
  static std::vector<const NodalLoad*> load_vec;
  load_vec.clear();
  
  static std::vector<std::pair<unsigned int, double> > load_pair;
  load_pair.clear();
  
  const LoadCaseBase& load_case = this->getLoadCaseFromID(load_info.getLoadCaseID());
  
  load_case.getLoadFactorLoadSetPair(load_info, load_pair);
  
  // create a map of load combos for the nodes
  static std::map<unsigned int, Loads::NodalLoadCombination*> load_combo_map;
  static std::map<unsigned int, Loads::NodalLoadCombination*>::iterator it;
  load_combo_map.clear();
  
  // iterate over the load pairs, and add the loads to the load vector
  for (unsigned int i=0; i < load_pair.size(); i++)
    {
      const NodalLoadSet& load_set = 
      dynamic_cast<const NodalLoadSet&>(this->getLoadSetFromID(load_pair[i].first));
      
      // make sure that the loads are of the same kind
      AssertThrow(load_info.getLoadNameEnumID() == load_set.getLoadNameEnumID(), ExcInternalError());
      
      load_vec.clear();
      // finally, ask the load set to return the load vector
      load_set.getAllLoads(load_vec);
      
      for (unsigned int j=0; j < load_vec.size(); j++)
        {
          it = load_combo_map.find(load_vec[j]->getNodeID());
          
          // if the node does not exist in the map, add it 
          if (it == load_combo_map.end())
            it = load_combo_map.insert
            (std::map<unsigned int, Loads::NodalLoadCombination*>::value_type
             (load_vec[j]->getNodeID(), Loads::createNodalLoadCombination
              (load_info.getLoadClassEnumID()).release())).first;
          
          it->second->addLoad(load_pair[i].second, load_vec[j]);
        }
    }
  
  // finally, transfer the created load combo pointers to the return vector
  FESystemUtility::AutoPtrVector<Loads::NodalLoadCombination> return_load_vec(load_combo_map.size());
  it = load_combo_map.begin();
  for (unsigned int i=0 ; i < load_combo_map.size(); i++)
    {
      return_load_vec.reset(i, it->second);
      it++;
    }
  
  load_combo_map.clear();
  
  return return_load_vec; 
}




template <>
FESystemUtility::AutoPtrVector<Loads::DirichletBoundaryConditionCombination>
LoadDatabase::getAllLoadCombinations<Loads::DirichletBoundaryConditionDataInfo,
Loads::DirichletBoundaryConditionCombination>
(const Loads::DirichletBoundaryConditionDataInfo& load_info)
{
  static std::vector<const DirichletBoundaryCondition*> load_vec;
  load_vec.clear();
  
  static std::vector<std::pair<unsigned int, double> > load_pair;
  load_pair.clear();
  
  const LoadCaseBase& load_case = this->getLoadCaseFromID(load_info.getLoadCaseID());
  
  load_case.getLoadFactorLoadSetPair(load_info, load_pair);
  
  // create a map of load combos for the nodes
  static std::map<unsigned int, std::map<unsigned int, 
  Loads::DirichletBoundaryConditionCombination*> > load_combo_map;
  load_combo_map.clear();
  
  static std::map<unsigned int, std::map<unsigned int, 
  Loads::DirichletBoundaryConditionCombination*> >::iterator it;
  static std::map<unsigned int, Loads::DirichletBoundaryConditionCombination*>::iterator dof_it;
  
  unsigned int n_bcs = 0;
  
  // iterate over the load pairs, and add the loads to the load vector
  for (unsigned int i=0; i < load_pair.size(); i++)
    {
      const BoundaryConditionLoadSet& load_set = 
      dynamic_cast<const BoundaryConditionLoadSet&>(this->getLoadSetFromID(load_pair[i].first));
      
      // make sure that the loads are of the same kind
      AssertThrow(load_info.getLoadNameEnumID() == load_set.getLoadNameEnumID(), ExcInternalError());
      
      load_vec.clear();
      // finally, ask the load set to return the load vector
      load_set.getAllLoads(load_vec);
      
      for (unsigned int j=0; j < load_vec.size(); j++)
        {
          it = load_combo_map.find(load_vec[j]->getNodeID());
          
          // if the node does not exist in the map, add it 
          if (it == load_combo_map.end())
            it = load_combo_map.insert
            (std::map<unsigned int, std::map<unsigned int, 
             Loads::DirichletBoundaryConditionCombination*> >::value_type
             (load_vec[j]->getNodeID(), 
              std::map<unsigned int, Loads::DirichletBoundaryConditionCombination*>())).first;
          
          // now, look for the dof number for this boundary condition
          dof_it = it->second.find(load_vec[j]->getDofNumber());
          
          if (dof_it == it->second.end())
            {
              dof_it = it->second.insert
              (std::map<unsigned int, Loads::DirichletBoundaryConditionCombination*>::value_type
               (load_vec[j]->getDofNumber(), Loads::createDirichletBoundaryConditionCombination
                (load_info.getLoadClassEnumID()).release())).first;
              n_bcs++;
            }
          
          dof_it->second->addLoad(load_pair[i].second, load_vec[j]);
        }
    }
  
  // finally, transfer the created load combo pointers to the return vector
  FESystemUtility::AutoPtrVector<Loads::DirichletBoundaryConditionCombination> 
  return_load_vec(n_bcs);
  it = load_combo_map.begin();
  unsigned int bc_it = 0;
  for ( ; it != load_combo_map.end(); it++)
    {
      dof_it = it->second.begin();
      for ( ; dof_it != it->second.end(); dof_it++)
        {
          return_load_vec.reset(bc_it, dof_it->second);
          bc_it++;
        }
    }
  
  load_combo_map.clear();
  
  return return_load_vec; 
}



