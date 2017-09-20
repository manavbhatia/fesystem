// $Id: LoadSet.h,v 1.4.6.1 2007-05-08 05:18:55 manav Exp $

#ifndef __fesystem_load_set_h__
#define __fesystem_load_set_h__

// C++ includes
#include <memory>
#include <string>
#include <vector>

// FESystem includes
#include "Utilities/NameEnumHandler.h"


// Forward declerations
class LoadBase;
class NodalLoad;
class SurfaceLoad;
class VolumeLoad;
class DirichletBoundaryCondition;


#ifndef VOLUME_LOAD_SET_ENUM_ID
#define VOLUME_LOAD_SET_ENUM_ID 1
#else
#error
#endif



#ifndef VOLUME_LOAD_SET_ENUM_NAME
#define VOLUME_LOAD_SET_ENUM_NAME "VOLUME_LOAD_SET"
#else
#error
#endif


#ifndef SURFACE_LOAD_SET_ENUM_ID
#define SURFACE_LOAD_SET_ENUM_ID 2
#else
#error
#endif



#ifndef SURFACE_LOAD_SET_ENUM_NAME
#define SURFACE_LOAD_SET_ENUM_NAME "SURFACE_LOAD_SET"
#else
#error
#endif


#ifndef NODAL_LOAD_SET_ENUM_ID
#define NODAL_LOAD_SET_ENUM_ID 3
#else
#error
#endif



#ifndef NODAL_LOAD_SET_ENUM_NAME
#define NODAL_LOAD_SET_ENUM_NAME "NODAL_LOAD_SET"
#else
#error
#endif



#ifndef BOUNDARY_CONDITION_LOAD_SET_ENUM_ID
#define BOUNDARY_CONDITION_LOAD_SET_ENUM_ID 4
#else
#error
#endif



#ifndef BOUNDARY_CONDITION_LOAD_SET_ENUM_NAME
#define BOUNDARY_CONDITION_LOAD_SET_ENUM_NAME "BOUNDARY_CONDITION_LOAD_SET"
#else
#error
#endif



#ifndef AEROELASTIC_LOAD_SET_ENUM_ID
#define AEROELASTIC_LOAD_SET_ENUM_ID 5
#else
#error
#endif



#ifndef AEROELASTIC_LOAD_SET_ENUM_NAME
#define AEROELASTIC_LOAD_SET_ENUM_NAME "AEROELASTIC_LOAD_SET"
#else
#error
#endif



DeclareEnumClass(LoadSetKindEnum);

DeclareEnumName(VOLUME_LOAD_SET, LoadSetKindEnum,
		VOLUME_LOAD_SET_ENUM_ID,
		VOLUME_LOAD_SET_ENUM_NAME);

DeclareEnumName(SURFACE_LOAD_SET, LoadSetKindEnum,
		SURFACE_LOAD_SET_ENUM_ID,
		SURFACE_LOAD_SET_ENUM_NAME);

DeclareEnumName(NODAL_LOAD_SET, LoadSetKindEnum,
		NODAL_LOAD_SET_ENUM_ID,
		NODAL_LOAD_SET_ENUM_NAME);


DeclareEnumName(BOUNDARY_CONDITION_LOAD_SET, LoadSetKindEnum,
		BOUNDARY_CONDITION_LOAD_SET_ENUM_ID,
		BOUNDARY_CONDITION_LOAD_SET_ENUM_NAME);

DeclareEnumName(AEROELASTIC_LOAD_SET, LoadSetKindEnum,
                AEROELASTIC_LOAD_SET_ENUM_ID,
                AEROELASTIC_LOAD_SET_ENUM_NAME);

class LoadSetBase
{
 public:	

  LoadSetBase();
  
  virtual ~LoadSetBase();

  virtual unsigned int getLoadSetKindEnumID() const = 0;
  
  unsigned int getLoadSetID() const;
  
  virtual std::istream& readFromInputStream(std::istream& input) = 0;
  
  unsigned int getLoadNameEnumID() const;

 protected: 

  unsigned int load_set_ID;

  unsigned int load_kind_enum_ID;
  
  std::string set_tag;
};





class VolumeLoadSet: public LoadSetBase
{
 public: 

  VolumeLoadSet();

  virtual ~VolumeLoadSet();

  virtual unsigned int getLoadSetKindEnumID() const;

  /// This will search the loads in this set and add create a vector containing loads 
  /// that act on an element
  void getLoadsForElement(const unsigned int,
                          std::vector<const VolumeLoad*> &loads) const;
  
  
  void getAllLoads(std::vector<const VolumeLoad*>& loads) const;
  
  virtual std::istream& readFromInputStream(std::istream& input);

 protected:
  
  std::multimap<unsigned int, VolumeLoad*> elem_ID_to_load_map;

};





class SurfaceLoadSet: public LoadSetBase
{
 public:
  
  SurfaceLoadSet();


  ~SurfaceLoadSet();
  
  virtual unsigned int getLoadSetKindEnumID() const;
  
  /// This will search the loads in this set and add create a vector containing loads 
  /// that act on an element
  void getLoadsForElement(const unsigned int,
                          std::vector<const SurfaceLoad*> &loads) const;
  
  /// This will find loads for the specified side of an element
  void getLoadsForElementSide(const unsigned int, 
                              const unsigned int,
                              std::vector<const SurfaceLoad*>& loads ) const;
  
  void getAllLoads(std::vector<const SurfaceLoad*>& loads) const;
  
  
  virtual std::istream& readFromInputStream(std::istream& input);

 protected:

  std::multimap<unsigned int, SurfaceLoad*> elem_ID_to_load_map;
  
};







class NodalLoadSet: public LoadSetBase
{
 public:
  NodalLoadSet();
  
  ~NodalLoadSet();
  
  virtual unsigned int getLoadSetKindEnumID() const;

  /// This will search the loads in this set and add create a vector containing loads 
  /// that act on a node
  void getLoadsForNode(const unsigned int ,
                       std::vector<const NodalLoad*>& loads) const;
  
  void getAllLoads(std::vector<const NodalLoad*>& loads) const;
  
  virtual std::istream& readFromInputStream(std::istream& );
  
 protected:
  
  std::multimap<unsigned int, NodalLoad*> node_ID_to_load_map;
};





class BoundaryConditionLoadSet: public LoadSetBase
{
 public:
  BoundaryConditionLoadSet();
  
  ~BoundaryConditionLoadSet();
  
  virtual unsigned int getLoadSetKindEnumID() const;

  /// This will search the loads in this set and add create a vector containing loads 
  /// that act on a node
  void getLoadsForNode(const unsigned int node_id,
                       std::vector<const DirichletBoundaryCondition*>& loads) const;
  
  void getAllLoads(std::vector<const DirichletBoundaryCondition*>& loads) const;
  
  virtual std::istream& readFromInputStream(std::istream& );
  
 protected:
  
  std::multimap<unsigned int, DirichletBoundaryCondition*> node_ID_to_load_map;
};


/// this class defines the operational conditions for an aeroelastic solution. The 
/// user defines the operational altitude and the mach numbers at which solution has 
/// to be performed
class AeroelasticLoadSet: public LoadSetBase
  {
  public:
    AeroelasticLoadSet();
    
    ~AeroelasticLoadSet();

    /// @returns an integer defining the enumeration ID of this object
    virtual unsigned int getLoadSetKindEnumID() const;
    
    /// @returns a reference to the vector of mach numbers for which analysis is to be performed
    const std::vector<double>& getMachNumbers() const;

    /// @returns a reference to the vector of dynamic pressures for which analysis is to be performed
    const std::vector<double>& getDynamicPressures() const;

    /// reads from the input stream and initializes itself
    virtual std::istream& readFromInputStream(std::istream& );
    
  protected:
    
    /// the altitude of operation, which will be used to calculate the
    /// density, pressure and sonic velocity from standard atmosphere
    double altitude;
    
    /// list of mach numbers for which the analysis is supposed to be performed
    std::vector<double> mach_numbers;
    
    /// list of dynamic pressures corrsponding to the mach numbers requested for analysis
    std::vector<double> dynamic_pressures;
    
  };


std::auto_ptr<LoadSetBase> createLoadSet(const unsigned int enum_ID);


#endif // __fesystem_load_set_h__

