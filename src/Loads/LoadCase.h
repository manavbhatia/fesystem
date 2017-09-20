// $Id: LoadCase.h,v 1.8.6.2 2007-06-13 14:58:05 manav Exp $

#ifndef __fesystem_load_case_base_h__
#define __fesystem_load_case_base_h__

// C++ includes
#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <memory>

// FESystem includes
#include "FESystem/FESystemNumbers.h"
#include "FESystem/FESystemExceptions.h"
#include "Utilities/NameEnumHandler.h"


// Forward declerations
class LoadDatabase;

namespace Loads
{
  class LoadDataInfoBase;
}


/// this class stores the ID of the loads acting on the finite element model at each time
/// instant. The IDs are stored against the load name, so that when a finite element analysis
/// problem is defined, only the names of the loads acting on it need be declared, and the
/// loadcase will store all the data.

#ifndef STATIC_LOAD_CASE_ENUM_ID 
#define STATIC_LOAD_CASE_ENUM_ID 1
#else
#error
#endif

#ifndef STATIC_LOAD_CASE_ENUM_NAME 
#define STATIC_LOAD_CASE_ENUM_NAME "STATIC_LOAD_CASE"
#else
#error
#endif


#ifndef TRANSIENT_LOAD_CASE_ENUM_ID 
#define TRANSIENT_LOAD_CASE_ENUM_ID 2
#else
#error
#endif

#ifndef TRANSIENT_LOAD_CASE_ENUM_NAME 
#define TRANSIENT_LOAD_CASE_ENUM_NAME "TRANSIENT_LOAD_CASE"
#else
#error
#endif

DeclareEnumClass(LoadCaseKindEnum);

DeclareEnumName(STATIC_LOAD_CASE, LoadCaseKindEnum, 
		STATIC_LOAD_CASE_ENUM_ID, 
		STATIC_LOAD_CASE_ENUM_NAME);


DeclareEnumName(TRANSIENT_LOAD_CASE, LoadCaseKindEnum, 
		TRANSIENT_LOAD_CASE_ENUM_ID,
		TRANSIENT_LOAD_CASE_ENUM_NAME);


class LoadCaseBase
{
 public:
  /// constructor
  LoadCaseBase(LoadDatabase& database);
	
  /// destructor
  virtual ~LoadCaseBase();
	
  /// clears the data structures for this object
  virtual void clear();

  /// @returns the load case kind
  virtual unsigned int getLoadCaseKindEnumID() const = 0;
  
  /// fills the vector with the load set ID and its load factor
  virtual void getLoadFactorLoadSetPair(const Loads::LoadDataInfoBase& load_info,
					std::vector<std::pair<unsigned int, double> >& load_pair) const =0;
  

  /// @returns the ID of this load case
  unsigned int getLoadCaseID() const;
  
  /// this method reads the load case data from the input stream
  virtual std::istream& readFromInputStream(std::istream& ) = 0;
  
  
 protected:
  
  /// load database reference
  LoadDatabase& load_database;

  /// load case ID for this LoadCase object
  unsigned int load_case_ID;
		
};



typedef std::map<double, std::pair<unsigned int, double> > TimeToLoadFactorValuePairMap;
typedef std::map<unsigned int, TimeToLoadFactorValuePairMap> DVToTimeLoadFactorValueMap;
typedef std::map<unsigned int, DVToTimeLoadFactorValueMap> LoadEnumToDVTimeMap;
typedef std::map<unsigned int, std::pair<unsigned int, unsigned int> > DVToLoadFactorIDPairMap;
typedef std::map<unsigned int, DVToLoadFactorIDPairMap> LoadEnumToLoadFactorIDMap;
typedef std::map<unsigned int, std::pair<unsigned int, double> > DVToLoadFactorValuePairMap;
typedef std::map<unsigned int, DVToLoadFactorValuePairMap> LoadEnumIDToDVLoadFactorValueMap;

class TimeIndependentLoadCase: public LoadCaseBase
{
 public:
  /// constructor
  TimeIndependentLoadCase(LoadDatabase& database);
	
  /// destructor
  virtual ~TimeIndependentLoadCase();
	
  /// clears the data structures for this object
  virtual void clear();
  
  /// @returns the load case kind
  virtual unsigned int getLoadCaseKindEnumID() const;

  /// fills the vector with the load set ID and its load factor
  virtual void getLoadFactorLoadSetPair(const Loads::LoadDataInfoBase& load_info,
					std::vector<std::pair<unsigned int, double> >& load_pair) const;
  
  
  /// this method reads the load case data from the input stream
  virtual std::istream& readFromInputStream(std::istream& );
  
  
  /// this will add a loadset with its reference tag.
  void addLoadSetID(const unsigned int load_enum_ID,
                    const unsigned int load_set_ID,
                    const double load_factor,
                    const unsigned int DV_ID);

protected:
  
  /// fills the vector with the load set ID, load factor value pair
  void getLoadPair(const unsigned int load_enum_ID,
		   const unsigned int DV_ID,
		   std::vector<std::pair<unsigned int, double> >& load_pair) const;
  
  
	

	
  /// this will store the LoadSetIDs for each time instant. At each time instant, 
  /// multiple kinds of loads can exist. The LoadSetIDMap stores these as a tag 
  /// to LoadSetID map
  LoadEnumIDToDVLoadFactorValueMap load_map;	
};




class TimeDependentLoadCase: public LoadCaseBase
{
 public:
  /// constructor
  TimeDependentLoadCase(LoadDatabase& database);
	
  /// destructor
  virtual ~TimeDependentLoadCase();
	
  /// clears the data structures for this object
  virtual void clear();
  
  /// @returns the load case kind
  virtual unsigned int getLoadCaseKindEnumID() const;

  /// @returns true if the specified load enum ID is load factor dependent
  bool ifLoadFactorDependent(unsigned int load_enum_ID);

  /// fills the vector with the load set ID and its load factor
  virtual void getLoadFactorLoadSetPair(const Loads::LoadDataInfoBase& load_info,
					std::vector<std::pair<unsigned int, double> >& load_pair) const;
  
  
  /// this method reads the load case data from the input stream
  virtual std::istream& readFromInputStream(std::istream& );
  
  /// this will add a loadset with its reference tag to a specified
  /// time instant. If the specified time instant does not exist,
  /// this method will create it and add the load
  void addLoadSetID(const double time,
                    const unsigned int load_enum_ID,
                    const unsigned int load_set_ID,
                    const double load_factor,
                    const unsigned int DV_ID );
  
  
  /// this will add a loadset with its reference tag to a specified
  /// time instant. If the specified time instant does not exist,
  /// this method will create it and add the load
  void addLoadSetID(const unsigned int load_enum_ID,
                    const unsigned int load_set_ID,
                    const unsigned int load_factor_ID,
                    const unsigned int DV_ID);
  
 protected:
 
  /// fills the vector with the load pair
  void getLoadPair(const double time, 
		   const unsigned int load_enum_ID,
		   const unsigned int DV_ID,
		   std::vector<std::pair<unsigned int, double> >& load_pair) const;


  /// this stores the load set IDs for loads associated with a time dependent load factor
  LoadEnumToLoadFactorIDMap load_factor_dependent_load_ID_map;
  
  /// this will store the LoadSetIDs for each time instant. At each time instant, 
  /// multiple kinds of loads can exist. The LoadSetIDMap stores these as a tag 
  /// to LoadSetID map
  LoadEnumToDVTimeMap time_load_ID_map;
  
};

std::auto_ptr<LoadCaseBase> createLoadCase(const unsigned int enum_ID, LoadDatabase& database);



inline 
unsigned int
LoadCaseBase::getLoadCaseID() const
{
  return this->load_case_ID;
}





/* inline */
/* void */
/* TimeIndependentLoadCase::getLoadPair(const unsigned int load_enum_ID, */
/* 				     const unsigned int DV_ID, */
/* 				     std::vector<std::pair<unsigned int, double> >& load_pair) const */
/* { */
/*   load_pair.clear(); */

/*   // get the iterator for the specified time instant */
/*   LoadEnumIDToDVLoadFactorValueMap::const_iterator load_enum_it, load_enum_end; */
/*   load_enum_it = this->load_map.find(load_enum_ID); */
/*   load_enum_end = this->load_map.end(); */
  
/*   if(load_enum_it == load_enum_end) */
/*     return; */
  

/*   // get the iterators for the load enum from this */
/*   DVToLoadFactorValuePairMap::const_iterator DV_it, DV_end; */
/*   DV_it = load_enum_it->second.find(DV_ID); */
/*   DV_end = load_enum_it->second.end(); */
  
/*   if (DV_it == DV_end) */
/*     return; */

/*   // If the code gets here, then a valid iterator has been found */
/*   load_pair.push_back(std::pair<unsigned int, double>(DV_it->second.first, DV_it->second.second)); */
/* } */




#endif // __fesystem_load_case_base_h__
