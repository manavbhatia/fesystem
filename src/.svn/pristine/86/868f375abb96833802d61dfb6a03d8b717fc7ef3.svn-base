// $Id: LoadCase.C,v 1.9.4.4 2008-02-25 04:26:22 manav Exp $

// C++ includes


// FESystem includes
#include "Loads/LoadCase.h"
#include "Loads/load.h"
#include "Loads/LoadDataInfo.h"
#include "Loads/LoadDatabase.h"
#include "Utilities/InputOutputUtility.h"
#include "FESystem/FESystemNumbers.h"
#include "FESystem/FESystemController.h"
#include "Database/FunctionDatabase.h"
#include "Discipline/AnalysisDisciplineBase.h"


LoadCaseBase::LoadCaseBase(LoadDatabase& database):
  load_database(database)
{
  
}



LoadCaseBase::~LoadCaseBase()
{
  
}



void
LoadCaseBase::clear()
{
  this->load_case_ID = FESystemNumbers::InvalidID;
}


TimeIndependentLoadCase::TimeIndependentLoadCase(LoadDatabase& database):
LoadCaseBase(database)
{
  
}


TimeIndependentLoadCase::~TimeIndependentLoadCase()
{
  this->clear();
  LoadCaseBase::clear();
}



void
TimeIndependentLoadCase::clear()
{
  this->load_map.clear();
  LoadCaseBase::clear();
}


unsigned int 
TimeIndependentLoadCase::getLoadCaseKindEnumID() const
{
  return STATIC_LOAD_CASE::num();
}



void
TimeIndependentLoadCase::addLoadSetID(const unsigned int load_enum_ID,
				      const unsigned int load_set_ID, 
				      const double load_factor,
				      const unsigned int DV_ID)
{
  // get the iterator for the specified time instant
  LoadEnumIDToDVLoadFactorValueMap::iterator load_enum_it, load_enum_end;
  load_enum_it = this->load_map.find(load_enum_ID);
  load_enum_end = this->load_map.end();
  
  if(load_enum_it == load_enum_end)
    load_enum_it = this->load_map.insert
      (LoadEnumIDToDVLoadFactorValueMap::value_type(load_enum_ID, DVToLoadFactorValuePairMap::map())).first;
  

  // get the iterators for the load enum from this
  DVToLoadFactorValuePairMap::iterator DV_it, DV_end;
  DV_it = load_enum_it->second.find(DV_ID);
  DV_end = load_enum_it->second.end();
  
  AssertThrow(DV_it == DV_end,  FESystemExceptions::ExcDuplicateID("Design Variable for Load Set",
								   DV_ID));
  load_enum_it->second.insert
    (DVToLoadFactorValuePairMap::value_type(DV_ID, std::pair<unsigned int, double>
					    (load_set_ID, load_factor)));
}



void
TimeIndependentLoadCase::getLoadFactorLoadSetPair
(const Loads::LoadDataInfoBase& load_info,
 std::vector<std::pair<unsigned int, double> >& load_pair) const
{
  Assert(!load_info.ifTimeDependent(), ExcInternalError());
  
  load_pair.clear();
  
  unsigned int DV_ID = FESystemNumbers::InvalidID;
  if (load_info.ifSensitivity())
    DV_ID = load_info.getDVID();
  
  this->getLoadPair(load_info.getLoadNameEnumID(), DV_ID,
                    load_pair);
}



void
TimeIndependentLoadCase::getLoadPair
(const unsigned int load_enum_ID,
 const unsigned int DV_ID,
 std::vector<std::pair<unsigned int, double> >& load_pair) const
{
  load_pair.clear();
  
  // find where the load exists: in the time loads, or the load factor loads
  // first check in the load factor based load
  // get the iterator for the specified time instant
  LoadEnumIDToDVLoadFactorValueMap::const_iterator load_enum_it, load_enum_end;
  load_enum_it = this->load_map.find(load_enum_ID);
  load_enum_end = this->load_map.end();
  
  if (load_enum_it != load_enum_end) // implies that the load exists in this map
    {
    // check for the DV
    DVToLoadFactorValuePairMap::const_iterator DV_it, DV_end;
    DV_it = load_enum_it->second.find(DV_ID);
    DV_end = load_enum_it->second.end();
    
    if (DV_it == DV_end)
      return;
    
    // If the code gets here, then a valid iterator has been found
    load_pair.push_back(std::pair<unsigned int, double>(DV_it->second.first, 
                                                        DV_it->second.second));
    }  
}




std::istream& TimeIndependentLoadCase::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int  n_load_sets=0, n_DV=0, load_set_ID;
  double load_factor = 0.0;
  
  // the input should begin with BEGIN
  FESystemIOUtility::readFromInput(input, STATIC_LOAD_CASE::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  // read in the load case ID
  FESystemIOUtility::readFromInput(input, "ID", this->load_case_ID);
  
  // read in the number of time instants and the number of loadsets defined for 
  // each time instant. The same loads will have to be redefined at each time instant
  
  FESystemIOUtility::readFromInput(input, "N_DESIGN_VARIABLES", n_DV);
  
  
  // for each time instant, read in the time value, the loadset tags and IDs
  unsigned int load_enum_ID = 0, DV_ID =FESystemNumbers::InvalidID;
  for (unsigned int dv_it= 0; dv_it <= n_DV; dv_it++)
    {

    // apart from the first data set, all other will be for a specific design variable
    if (dv_it > 0)
      {
	FESystemIOUtility::readFromInput(input, "DESIGN_VARIABLE", DV_ID);
	// no given IDs should equal the invalid ID of FESystem
	AssertThrow(DV_ID != FESystemNumbers::InvalidID, ExcInternalError());
      }
    
      FESystemIOUtility::readFromInput(input, "N_LOAD_SETS", n_load_sets);
      
      for (unsigned int load_set_incr=0; load_set_incr < n_load_sets; load_set_incr++)
        {
        input >> tag; //load set name
        load_enum_ID = LoadNameEnum::enumID(tag);
	input >> load_factor;
        input >> load_set_ID; // load set ID
        
	this->addLoadSetID(load_enum_ID, load_set_ID, load_factor, DV_ID);
        
        tag.clear();
        }
    }
  
  // the input should end with END
  FESystemIOUtility::readFromInput(input, STATIC_LOAD_CASE::name());
  FESystemIOUtility::readFromInput(input, "END");
  return input;
}





TimeDependentLoadCase::TimeDependentLoadCase(LoadDatabase& database):
LoadCaseBase(database)
{
  
}


TimeDependentLoadCase::~TimeDependentLoadCase()
{
  
}



void
TimeDependentLoadCase::clear()
{
  this->load_factor_dependent_load_ID_map.clear();
  this->time_load_ID_map.clear();
  LoadCaseBase::clear();
}



unsigned int 
TimeDependentLoadCase::getLoadCaseKindEnumID() const
{
  return TRANSIENT_LOAD_CASE::num();
}





void
TimeDependentLoadCase::addLoadSetID(const double time_value,
				    const unsigned int load_enum_ID,
				    const unsigned int load_set_ID, 
				    const double load_factor,
				    const unsigned int DV_ID)
{
  // get the iterator to the specified time instant
  LoadEnumToDVTimeMap::iterator load_enum_it, load_enum_end;
  load_enum_it = this->time_load_ID_map.find(load_enum_ID);
  load_enum_end = this->time_load_ID_map.end();
  
  if (load_enum_it == load_enum_end)
    load_enum_it = this->time_load_ID_map.insert
      (LoadEnumToDVTimeMap::value_type(load_enum_ID, DVToTimeLoadFactorValueMap::map())).first;
  
  // get the iterator for the specified time instant
  DVToTimeLoadFactorValueMap::iterator DV_it, DV_end;
  DV_it = load_enum_it->second.find(DV_ID);
  DV_end = load_enum_it->second.end();
  
  if(DV_it == DV_end)
    DV_it = load_enum_it->second.insert
      (DVToTimeLoadFactorValueMap::value_type(DV_ID, TimeToLoadFactorValuePairMap::map())).first;
  

  // get the iterators for the load enum from this
  TimeToLoadFactorValuePairMap::iterator time_it, time_end;
  time_it = DV_it->second.find(time_value);
  time_end = DV_it->second.end();
  
  // duplicate time value specifies
  AssertThrow(time_it == time_end,  ExcInternalError());
  DV_it->second.insert
    (TimeToLoadFactorValuePairMap::value_type(time_value, std::pair<unsigned int, double>
					      (load_set_ID, load_factor)));
}





void
TimeDependentLoadCase::addLoadSetID(const unsigned int load_enum_ID,
				    const unsigned int load_set_ID, 
				    const unsigned int load_factor_ID,
				    const unsigned int DV_ID)
{
  // get the iterator for the specified time instant
  LoadEnumToLoadFactorIDMap::iterator load_enum_it, load_enum_end;
  load_enum_it = this->load_factor_dependent_load_ID_map.find(load_enum_ID);
  load_enum_end = this->load_factor_dependent_load_ID_map.end();
  
  if(load_enum_it == load_enum_end)
    load_enum_it = this->load_factor_dependent_load_ID_map.insert
      (LoadEnumToLoadFactorIDMap::value_type(load_enum_ID, DVToLoadFactorIDPairMap::map())).first;
  

  // get the iterators for the load enum from this
  DVToLoadFactorIDPairMap::iterator DV_it, DV_end;
  DV_it = load_enum_it->second.find(DV_ID);
  DV_end = load_enum_it->second.end();
  
  AssertThrow(DV_it == DV_end,  FESystemExceptions::ExcDuplicateID("Design Variable for Load Set",
								   DV_ID));
  load_enum_it->second.insert
    (DVToLoadFactorIDPairMap::value_type(DV_ID, std::pair<unsigned int, unsigned int>
					 (load_set_ID, load_factor_ID)));
}





std::istream& 
TimeDependentLoadCase::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int n_time_instants=0, n_load_sets=0, n_DV=0, load_set_ID,
    n_load_factor_vals =0, load_factor_ID;
  double load_factor = 0.0;

  // the input should begin with BEGIN
  FESystemIOUtility::readFromInput(input, TRANSIENT_LOAD_CASE::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  // read in the load case ID
  FESystemIOUtility::readFromInput(input, "ID", this->load_case_ID);
	
  // read in the number of time instants and the number of loadsets defined for 
  // each time instant. The same loads will have to be redefined at each time instant
  
  FESystemIOUtility::readFromInput(input, "N_TRANSIENT_LOAD_FACTOR_DEPENDENT_LOAD_SETS", n_load_factor_vals);  
  FESystemIOUtility::readFromInput(input, "N_TIME_INSTANTS", n_time_instants);
  FESystemIOUtility::readFromInput(input, "N_DESIGN_VARIABLES", n_DV);
	

  // for each time instant, read in the time value, the loadset tags and IDs
  double time=0.0;
  unsigned int load_enum_ID = 0, DV_ID = FESystemNumbers::InvalidID;
  for (unsigned int dv_it= 0; dv_it <= n_DV; dv_it++)
    {

    // apart from the first data set, all other will be for a specific design variable
      if (dv_it > 0)
	{
	  FESystemIOUtility::readFromInput(input, "DESIGN_VARIABLE", DV_ID);
	  // no given IDs should equal the invalid ID of FESystem
	  AssertThrow(DV_ID != FESystemNumbers::InvalidID, ExcInternalError());	
	  FESystemIOUtility::readFromInput(input, "N_LOAD_SETS", n_load_factor_vals);
	}

      // iterate over all load factor dependent load sets and store them
      for (unsigned int nl_sets =0; nl_sets < n_load_factor_vals; nl_sets++)
	{
	  input >> tag;
	  load_enum_ID = LoadNameEnum::enumID(tag);
	  input >> load_factor_ID;
	  input >> load_set_ID;

	  this->addLoadSetID(load_enum_ID, load_set_ID, load_factor_ID, DV_ID);
	}
      
      for (unsigned int time_incr=0; time_incr < n_time_instants; time_incr++)
	{
	  FESystemIOUtility::readFromInput(input, "TIME_INSTANT", time);
	  FESystemIOUtility::readFromInput(input, "N_LOAD_SETS", n_load_sets);
	  
	  for (unsigned int load_set_incr=0; load_set_incr < n_load_sets; load_set_incr++)
	    {
	      input >> tag; //load set name
	      load_enum_ID = LoadNameEnum::enumID(tag);
	      input >> load_factor;
	      input >> load_set_ID; // load set ID
	      
	      this->addLoadSetID(time, load_enum_ID, load_set_ID, load_factor, DV_ID);
	      
	      tag.clear();
	    }
	}
    } 
  
  // the input should end with END
  FESystemIOUtility::readFromInput(input, TRANSIENT_LOAD_CASE::name());
  FESystemIOUtility::readFromInput(input, "END");
  return input;
}




void
TimeDependentLoadCase::getLoadFactorLoadSetPair
(const Loads::LoadDataInfoBase& load_info,
 std::vector<std::pair<unsigned int, double> >& load_pair) const
{
  Assert(load_info.ifTimeDependent(), ExcInternalError());
  
  load_pair.clear();

  unsigned int DV_ID = FESystemNumbers::InvalidID;
  if (load_info.ifSensitivity())
    DV_ID = load_info.getDVID();
  
  this->getLoadPair(load_info.getTime(),
                    load_info.getLoadNameEnumID(),
                    DV_ID,
                    load_pair);
}



void
TimeDependentLoadCase::getLoadPair(const double time_value, 
				   const unsigned int load_enum_ID,
				   const unsigned int DV_ID,
				   std::vector<std::pair<unsigned int, double> >& load_pair) const
{
  load_pair.clear();
  
  // find where the load exists: in the time loads, or the load factor loads
  // first check in the load factor based load
  // get the iterator for the specified time instant
  LoadEnumToLoadFactorIDMap::const_iterator load_enum_it, load_enum_end;
  load_enum_it = this->load_factor_dependent_load_ID_map.find(load_enum_ID);
  load_enum_end = this->load_factor_dependent_load_ID_map.end();
  
  if (load_enum_it != load_enum_end) // implies that the load exists in this map
    {
      // check for the DV
      DVToLoadFactorIDPairMap::const_iterator DV_it, DV_end;
      DV_it = load_enum_it->second.find(DV_ID);
      DV_end = load_enum_it->second.end();
      
      if (DV_it == DV_end)
	return;
      
      // get the function from the function database, and add the value to the pair
      double value = 0.0;
      this->load_database.getFESystemController().function_database->getScalarFunctionValue
	(DV_it->second.second, time_value, value);

      // If the code gets here, then a valid iterator has been found
      load_pair.push_back(std::pair<unsigned int, double>(DV_it->second.first, value));
    }
  // check for the load in the time dependent load. If the specified time instant does not exist in 
  // the map, then the loads are interpolated
  // get the iterator to the specified time instant
  {
    LoadEnumToDVTimeMap::const_iterator load_enum_it, load_enum_end;
    load_enum_it = this->time_load_ID_map.find(load_enum_ID);
    load_enum_end = this->time_load_ID_map.end();
  
    if (load_enum_it == load_enum_end)
      return;
  
    // get the iterator for the specified time instant
    DVToTimeLoadFactorValueMap::const_iterator DV_it, DV_end;
    DV_it = load_enum_it->second.find(DV_ID);
    DV_end = load_enum_it->second.end();
  
    if(DV_it == DV_end)
      return;
    
    // get the iterators for the load enum from this
    TimeToLoadFactorValuePairMap::const_iterator time_it, time_end;
    time_it = DV_it->second.find(time_value);
    time_end = DV_it->second.end();
    
    // if the time value does not exist, then the load needs to be interpolated. 
    if (time_it == time_end)
      Assert(false, ExcInternalError()); // needs to be implemented
      
  }  
  
}



std::auto_ptr<LoadCaseBase> 
createLoadCase(const unsigned int enum_ID, LoadDatabase& load_database)
{
  std::auto_ptr<LoadCaseBase> load_case;

  switch(enum_ID)
    {
    case STATIC_LOAD_CASE_ENUM_ID:
      load_case.reset(new TimeIndependentLoadCase(load_database));
      break;

    case TRANSIENT_LOAD_CASE_ENUM_ID:
      load_case.reset(new TimeDependentLoadCase(load_database));
      break;

    default:
      Assert(false, ExcInternalError());
    }
  
  return load_case;
}

