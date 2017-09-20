// $Id: DisciplineInfo.h,v 1.4 2006-09-05 20:09:54 manav Exp $

#ifndef __fesystem_discipline_info_h__
#define __fesystem_discipline_info_h__

// C++ includes
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <memory>

// FESystem include
#include "Discipline/AnalysisDisciplineBase.h"
#include "Utilities/NameEnumHandler.h"

namespace Discipline
{

  DeclareEnumClass(DisciplineInfoTypeEnum);
  
  /// this forms a base class for storing the information of disciplinary analysis. 
  /// the basic information for all finite element disciplinary analysis is the mesh, analysis
  /// kind, parameter dependence, etc.
  class DisciplineInfo
  {
public:
    
    /// constructor
    DisciplineInfo();

    /// destructor
    virtual ~DisciplineInfo();
    
    inline unsigned int getID() const;
    
    /// @returns mesh ID of for the disciplinary analysis
    inline unsigned int getMeshID() const;

    /// @returns the discipline enumeration ID
    virtual unsigned int getDisciplineEnumID() const =0;
    
    /// @returns the discipline enumeration name
    virtual const std::string getDisciplineEnumName() const =0;
    
    /// @returns the discipline analysis type enumeration ID,
    /// which can be linear or nonlinear.
    inline unsigned int getAnalysisTypeEnumID() const;
    
    /// @returns the discipline analysis type enumeration name,
    /// which can be linear or nonlinear.
    inline const std::string getAnalysisTypeEnumName() const;

    
    /// @returns set of local parameters on which this analysis depends
    inline const std::set<unsigned int>& getLocalParameterSet() const;
    
    /// @returns true if this discipline is dependent on this local parameter
    /// @param param_enum_ID parameter enum ID to be checked
    inline bool checkLocalParameterDependence(const unsigned int param_enum_ID) const;

    /// reads from the specified input stream
    virtual std::istream& readFromInputStream(std::istream& input) = 0;
    
    /// overloaded operator to read from input stream
    friend std::istream& operator>>(std::istream& input, Discipline::DisciplineInfo& info);

protected:
      
      /// ID
      unsigned int ID;

    /// mesh ID
    unsigned int mesh_ID;
    
    /// analysis type
    unsigned int analysis_type;
    
    /// dependence of the analysis on local parameters
    std::set<unsigned int> local_parameters;
  };
  

  /// @returns a pointer to a newly created DisciplineInfo
  std::auto_ptr<Discipline::DisciplineInfo> 
    createDisciplineInfo(const unsigned int info_enum_id);
  
}


inline 
unsigned int
Discipline::DisciplineInfo::getID() const
{
  return this->ID;
}


inline 
unsigned int
Discipline::DisciplineInfo::getMeshID() const
{
  return this->mesh_ID;
}




inline 
unsigned int
Discipline::DisciplineInfo::getAnalysisTypeEnumID() const
{
  return this->analysis_type;
}





inline 
const std::string
Discipline::DisciplineInfo::getAnalysisTypeEnumName() const
{
  return Discipline::DisciplineAnalysisTypeEnum::enumName(this->analysis_type);
}





inline
const std::set<unsigned int>& 
Discipline::DisciplineInfo::getLocalParameterSet() const
{
  return this->local_parameters;
}




inline
bool
Discipline::DisciplineInfo::checkLocalParameterDependence(const unsigned int param_enum_ID) const
{
  if (this->local_parameters.count(param_enum_ID) == 0)
    return false;
  else
    return true;
}



#endif //  __fesystem_analysis_discipline_info_h__
