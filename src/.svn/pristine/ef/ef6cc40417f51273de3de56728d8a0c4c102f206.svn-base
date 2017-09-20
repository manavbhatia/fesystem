// $Id: DesignDatabase.h,v 1.2 2006-09-05 20:41:41 manav Exp $

#ifndef __fesystem_design_database_h__
#define __fesystem_design_database_h__

// C++ includes
#include <iostream>
#include <vector>
#include <map>
#include <memory>

// FESystem includes 
#include "FESystem/FESystemExceptions.h"
#include "DesignData/PropertyParameter.h"
#include "DesignData/ShapeParameter.h"



namespace DesignData
{
class DesignDatabase
{
public:
  DesignDatabase();
  
  ~DesignDatabase();
  
  /// returns the vector of all parameters
  inline std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
    getParameters() const;

  /// returns the vector of property parameters
  inline std::auto_ptr<std::vector<DesignData::PropertyParameter*> > 
    getPropertyParameters() const;
  
  /// returns the vector of shape parameters
  inline std::auto_ptr<std::vector<DesignData::ShapeParameter*> > 
    getShapeParameters() const;
  
  /// returns a writable reference to parameter for the given ID
  inline DesignParameter& getParameterForID(const unsigned int ID) const;
  
  /// returns the number of parameters
  inline const unsigned int getNParams() const;
  
  /// overloaded operator to read from the input stream
  friend std::istream& operator>>(std::istream& input, 
                                  DesignData::DesignDatabase& database);
  
  // reads the parameter database from an input stream
  std::istream& readFromInputStream(std::istream& input);

protected:
  
  
  /// local type definitions
  typedef std::map<unsigned int, DesignData::DesignParameter*> IDToParameterMap; 
    
  /// vector of DVs
  std::map<unsigned int, DesignData::DesignParameter*> ID_to_design_parameter_map;
  
  /// this expception is thrown if the requested name does not exist in the map
  DeclException1(ExcFESystemIDDoesNotExist, unsigned int, 
                 << "Requested Parameter ID: " << arg1 << " does not exist in map.");
  
  DeclException1(ExcDuplicateParameterID, unsigned int, 
                 << "Parameter ID: " << arg1 << " already exists in map.");

private:
  /// adds a parameter to the map. This method should be used only by this class
  inline void addParameterToDatabase(DesignData::DesignParameter *param);
};
}


inline 
std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
DesignData::DesignDatabase::getParameters() const
{
  // iterate over all the params in the map, add it to the vector
  std::auto_ptr<std::vector<DesignData::DesignParameter*> > 
  vector(new std::vector<DesignData::DesignParameter*>());
	
  DesignData::DesignDatabase::IDToParameterMap::const_iterator param_it, param_end;
  param_it = this->ID_to_design_parameter_map.begin();
  param_end = this->ID_to_design_parameter_map.end();
	
  for (; param_it != param_end; param_it++)
    {
    vector->push_back(param_it->second);
    }
	
  return vector;
}



inline 
std::auto_ptr<std::vector<DesignData::PropertyParameter*> > 
DesignData::DesignDatabase::getPropertyParameters() const
{
  // iterate over all the DVs in the map, if type == property DV, add it to the vector
  std::auto_ptr<std::vector<DesignData::PropertyParameter*> > 
  vector(new std::vector<DesignData::PropertyParameter*>());
	
  DesignData::DesignDatabase::IDToParameterMap::const_iterator param_it, param_end;
  param_it = this->ID_to_design_parameter_map.begin();
  param_end = this->ID_to_design_parameter_map.end();
	
  for (; param_it != param_end; param_it++)
    {
    if (param_it->second->getParameterTypeEnumID() == 
        DesignData::PROPERTY_PARAMETER::num())
      {
      DesignData::PropertyParameter* ptr = 
      dynamic_cast<DesignData::PropertyParameter*>(param_it->second);
      vector->push_back(ptr);
      }
    }
	
  return vector;  
}



inline 
std::auto_ptr<std::vector<DesignData::ShapeParameter*> > 
DesignData::DesignDatabase::getShapeParameters() const
{
  // iterate over all the DVs in the map, if type == property DV, add it to the vector
  std::auto_ptr<std::vector<DesignData::ShapeParameter*> > 
  vector(new std::vector<DesignData::ShapeParameter*>());
	
  DesignData::DesignDatabase::IDToParameterMap::const_iterator param_it, param_end;
  param_it = this->ID_to_design_parameter_map.begin();
  param_end = this->ID_to_design_parameter_map.end();
	
  for (; param_it != param_end; param_it++)
    {
    if (param_it->second->getParameterTypeEnumID() == 
        DesignData::SHAPE_PARAMETER::num())
      {
      DesignData::ShapeParameter* ptr = 
      dynamic_cast<DesignData::ShapeParameter*>(param_it->second);
      vector->push_back(ptr);
      }
    }
	
  return vector;  
}



inline 
DesignData::DesignParameter&
DesignData::DesignDatabase::getParameterForID(const unsigned int ID) const
{
  // search the map for the existence of this ID, if found, return the param reference
  // pointer
  DesignData::DesignDatabase::IDToParameterMap::const_iterator param_it = 
  this->ID_to_design_parameter_map.find(ID);
	
  Assert(param_it != this->ID_to_design_parameter_map.end(),
         DesignData::DesignDatabase::ExcFESystemIDDoesNotExist(ID));
	
  return *(param_it->second);  
}



inline 
const unsigned int 
DesignData::DesignDatabase::getNParams() const
{
  return this->ID_to_design_parameter_map.size();
}




inline 
void 
DesignData::DesignDatabase::addParameterToDatabase(DesignData::DesignParameter *param)
{
  Assert(param != NULL, 
         FESystemExceptions::ExcNullPointer());
  
  // make sure that the ID does not already exists
  unsigned int ID = param->getID();
  if (this->ID_to_design_parameter_map.find(ID) != 
      this->ID_to_design_parameter_map.end())
    {
    delete param; param = NULL;
    Assert(false, DesignData::DesignDatabase::ExcDuplicateParameterID(ID));
    }
		
  bool insert_success = 
  this->ID_to_design_parameter_map.insert(DesignData::DesignDatabase::IDToParameterMap::
                                         value_type(ID, param)).second;
		
  Assert(insert_success, FESystemExceptions::ExcInsertUnsuccessful());
}





#endif // __fesystem_design_database_h__
