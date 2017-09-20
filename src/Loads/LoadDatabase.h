// $Id: LoadDatabase.h,v 1.9.6.2 2007-06-13 14:58:05 manav Exp $

#ifndef __fesystem_load_database_h__
#define __fesystem_load_database_h__

// C++ includes
#include <iostream>
#include <map>
#include <vector>
#include <string>

// FESystem includes
#include "Utilities/AutoptrVector.h"

// Forward declerations
namespace FESystem
{
  class FESystemController;
}

class LoadSetBase;
class LoadCaseBase;
class TimeIndependentLoadCase;
class TimeDependentLoadCase;


/// This class will store all the loads and will provide methods to obtain loads 
/// based on certain criterion, which could be load case, time instant, element ID,
/// load tag, etc. The results for all these will be returned as a vector of Load.
/// The idea of this class is similar to that of a database where you store all the 
/// data, perform a search on it, and return the results. 

class LoadDatabase
{
 public:	

  LoadDatabase(FESystem::FESystemController& controller);
	
  ~LoadDatabase();

  /// clear the data stored in the database
  void clear();
	
  /// @returns the fesystem controller
  FESystem::FESystemController& getFESystemController() const;


  /// function to read the load set and load case definition from a text file
  std::istream& readFromInputStream(std::istream& );
		
  /// this will return the loads acting on an element
  /// from the specified load case and time instant
  template <typename LoadInfoType, typename LoadComboType>
    void getLoadCombination(const LoadInfoType& load_info,
                            LoadComboType& load_combo);
  
  /// this will return all the loads from a load case for a specified time instant
  /// and tag. this is specifically for loads that are not element related, like 
  /// nodal loads
  template <typename LoadInfoType, typename LoadComboType>
    FESystemUtility::AutoPtrVector<LoadComboType> 
    getAllLoadCombinations(const LoadInfoType& load_info);
  
  

	
	
  /// returns the load set ID for the load given by the name in the load
  /// case
  /// @param load_case number of the load case
  /// @param time_instant the time instant for which load is needed
  /// @param load_name name of the load for which set ID is being seeked
  
  template <typename LoadInfoType>
    unsigned int getLoadSetID(const LoadInfoType& load_info);
  

  const LoadCaseBase& getLoadCaseFromID(const unsigned int ) const;
  
  const LoadSetBase& getLoadSetFromID(const unsigned int ) const;
		
  void addLoadSet(LoadSetBase* load_set);
	
  void addLoadCase(LoadCaseBase* load_case);

  unsigned int getAvailableLoadSetID() const;

 protected:
  
  template <typename LoadInfoType>
    void getLoadFactorAndSetIDPair(const LoadInfoType& load_info,
				   std::vector<std::pair<unsigned int, double> >& load_pair);
  
  template <typename LoadType, typename LoadSetType, typename LoadComboType>
    void getElemVolumeLoadCombination(const unsigned int load_case_ID,
                                      const unsigned int elem_ID,
                                      const unsigned int load_enum_ID,
                                      const bool if_time_dep, const double time_instant,
                                      const bool if_sensitivity, const unsigned int DV_ID,
                                      LoadComboType& load_combo);
    
  
  template <typename LoadType, typename LoadSetType, typename LoadComboType>
    void getNodalLoadCombination(const unsigned int load_case_ID,
                                 const unsigned int node_ID,
                                 const unsigned int load_enum_ID,
                                 const bool if_time_dep, const double time_instant,
                                 const bool if_sensitivity, const unsigned int DV_ID,
                                 LoadComboType& load_combo);
    
  
  template <typename LoadType, typename LoadSetType, typename LoadComboType>
    void getElemSurfaceLoadCombinationForSide(const unsigned int load_case_ID,
                                              const unsigned int elem_ID,
                                              const unsigned int side_ID,
                                              const unsigned int load_enum_ID,
                                              const bool if_time_dep, const double time_instant,
                                              const bool if_sensitivity, const unsigned int DV_ID,
                                              LoadComboType& load_combo);
    
    
  template <typename LoadType, typename LoadSetType, typename LoadComboType>
    FESystemUtility::AutoPtrVector<LoadComboType> getAllNodalLoadCombination
    (const unsigned int load_case_ID,
     const unsigned int load_enum_ID,
     const bool if_time_dep, const double time_instant,
     const bool if_sensitivity, const unsigned int DV_ID);

  
  FESystem::FESystemController& fesystem_controller;

  
  std::map<unsigned int, LoadSetBase*> load_set_map;
  
  
  std::map<unsigned int, LoadCaseBase*> load_case_map;
	
};



inline
FESystem::FESystemController&
LoadDatabase::getFESystemController() const
{
  return this->fesystem_controller;
}


#endif // __fesystem_load_database_h__
