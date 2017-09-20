// $Id: LoadSet.C,v 1.7.4.2 2007-05-08 05:18:55 manav Exp $

// C++ include


// FESystem include
#include "Loads/LoadSet.h"
#include "Loads/load.h"
#include "Utilities/InputOutputUtility.h"


LoadSetBase::LoadSetBase():
load_set_ID(FESystemNumbers::InvalidID),
load_kind_enum_ID(FESystemNumbers::InvalidID)
{
  
}





LoadSetBase::~LoadSetBase()
{
  
}





unsigned int 
LoadSetBase::getLoadSetID() const
{
  Assert(this->load_set_ID != FESystemNumbers::InvalidID, ExcInternalError());

  return this->load_set_ID;	
}




unsigned int 
LoadSetBase::getLoadNameEnumID() const
{
  Assert(this->load_kind_enum_ID != FESystemNumbers::InvalidID, ExcInternalError());
  
  return this->load_kind_enum_ID;	
}






VolumeLoadSet::VolumeLoadSet():
  LoadSetBase()
{
  
}



VolumeLoadSet::~VolumeLoadSet()
{
  // iterate over all loads and delete them
  std::multimap<unsigned int, VolumeLoad*>::iterator it, end;
  it = this->elem_ID_to_load_map.begin();
  end = this->elem_ID_to_load_map.end();
  
  for (; it != end; it++)
    {
      delete it->second;
      it->second = NULL;
    }

  this->elem_ID_to_load_map.clear();
}



unsigned int 
VolumeLoadSet::getLoadSetKindEnumID() const
{
  return VOLUME_LOAD_SET::num();	
}




void
VolumeLoadSet::getLoadsForElement(const unsigned int elem_ID,
				  std::vector<const VolumeLoad*>& loads) const
{
  std::multimap<unsigned int, VolumeLoad*>::const_iterator it, end;
  it = this->elem_ID_to_load_map.lower_bound(elem_ID);
  end = this->elem_ID_to_load_map.upper_bound(elem_ID);
  for ( ; it != end; it++)
    loads.push_back(it->second);
}




void
VolumeLoadSet::getAllLoads(std::vector<const VolumeLoad*>& loads) const
{
  std::multimap<unsigned int, VolumeLoad*>::const_iterator it, end;
  it = this->elem_ID_to_load_map.begin();
  end = this->elem_ID_to_load_map.end();
  loads.resize(this->elem_ID_to_load_map.size());
  unsigned int i=0;
  for (; it != end; it++)
    {
      loads[i] = it->second;
      i++;
    }
}



std::istream& 
VolumeLoadSet::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int n_loads=0, elem_ID=0, load_class_enum_ID = 0;
	
  // the load set should begin with the BEGIN tag, and followed by 
  // the number of loads in the set
  FESystemIOUtility::readFromInput(input, VOLUME_LOAD_SET::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
	
	
  // next, it should contain the name, ID and kind, and number of DOFs of load
  FESystemIOUtility::readFromInput(input, "NAME", this->set_tag);
  FESystemIOUtility::readFromInput(input, "ID", this->load_set_ID);
  
  tag.clear();
  FESystemIOUtility::readFromInput(input, "KIND", tag);
  this->load_kind_enum_ID = LoadNameEnum::enumID(tag);
  load_class_enum_ID = getLoadClassEnumIDForLoadName(this->load_kind_enum_ID);

  FESystemIOUtility::readFromInput(input, "N_LOADS", n_loads);
	
  
  std::auto_ptr<VolumeLoad> load;
  for (unsigned int load_incr=0; load_incr < n_loads; load_incr++)
    {
      load.reset(createVolumeLoad(load_class_enum_ID).release());
      load->readFromInputStream(input);
      
      elem_ID = load->getLoadID();
      this->elem_ID_to_load_map.insert
	(std::multimap<unsigned int, VolumeLoad*>::value_type
	 (elem_ID, load.release()));
    }
	
	
  FESystemIOUtility::readFromInput(input, VOLUME_LOAD_SET::name());
  FESystemIOUtility::readFromInput(input, "END");
	
  return input;
}





NodalLoadSet::NodalLoadSet():
  LoadSetBase()
{
  
}



NodalLoadSet::~NodalLoadSet()
{
  // iterate over all loads and delete them
  std::multimap<unsigned int, NodalLoad*>::iterator it, end;
  it = this->node_ID_to_load_map.begin();
  end = this->node_ID_to_load_map.end();
  
  for (; it != end; it++)
    {
      delete it->second;
      it->second = NULL;
    }

  this->node_ID_to_load_map.clear();
}



unsigned int 
NodalLoadSet::getLoadSetKindEnumID() const
{
  return NODAL_LOAD_SET::num();	
}




void
NodalLoadSet::getLoadsForNode(const unsigned int elem_ID,
			      std::vector<const NodalLoad*>& loads) const
{
  std::multimap<unsigned int, NodalLoad*>::const_iterator it, end;
  it = this->node_ID_to_load_map.lower_bound(elem_ID);
  end = this->node_ID_to_load_map.upper_bound(elem_ID);
  for ( ; it != end; it++)
    loads.push_back(it->second);
}




void
NodalLoadSet::getAllLoads(std::vector<const NodalLoad*>& loads) const
{
  std::multimap<unsigned int, NodalLoad*>::const_iterator it, end;
  it = this->node_ID_to_load_map.begin();
  end = this->node_ID_to_load_map.end();
  loads.resize(this->node_ID_to_load_map.size());
  unsigned int i=0;
  for (; it != end; it++)
    {
      loads[i] = it->second;
      i++;
    }
}



std::istream& 
NodalLoadSet::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int n_loads=0, n_dofs=0, node_ID = 0, load_class_enum_ID = 0;

  // the load set should begin with the BEGIN tag, and followed by 
  // the number of loads in the set
  FESystemIOUtility::readFromInput(input, NODAL_LOAD_SET::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  
  // next, it should contain the name, ID and kind, and number of DOFs of load
  FESystemIOUtility::readFromInput(input, "NAME", this->set_tag);
  FESystemIOUtility::readFromInput(input, "ID", this->load_set_ID);

  tag.clear();
  FESystemIOUtility::readFromInput(input, "KIND", tag);
  this->load_kind_enum_ID = LoadNameEnum::enumID(tag);
  load_class_enum_ID = getLoadClassEnumIDForLoadName(this->load_kind_enum_ID);

  FESystemIOUtility::readFromInput(input, "N_DOFS", n_dofs);
  
  FESystemIOUtility::readFromInput(input, "N_LOADS", n_loads);
  
  
  std::auto_ptr<NodalLoad> load;

  for (unsigned int load_incr=0; load_incr < n_loads; load_incr++)
    {
      load.reset(createNodalLoad(load_class_enum_ID, n_dofs).release());
      load->readFromInputStream(input);
      
      node_ID = load->getNodeID();
      this->node_ID_to_load_map.insert
	(std::multimap<unsigned int, NodalLoad*>::value_type
	 (node_ID, load.release()));
    }
  
  
  FESystemIOUtility::readFromInput(input, NODAL_LOAD_SET::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}






BoundaryConditionLoadSet::BoundaryConditionLoadSet():
  LoadSetBase()
{
  
}



BoundaryConditionLoadSet::~BoundaryConditionLoadSet()
{
  // iterate over all loads and delete them
  std::multimap<unsigned int, DirichletBoundaryCondition*>::iterator it, end;
  it = this->node_ID_to_load_map.begin();
  end = this->node_ID_to_load_map.end();
  
  for (; it != end; it++)
    {
      delete it->second;
      it->second = NULL;
    }

  this->node_ID_to_load_map.clear();
}



unsigned int 
BoundaryConditionLoadSet::getLoadSetKindEnumID() const
{
  return BOUNDARY_CONDITION_LOAD_SET::num();	
}




void
BoundaryConditionLoadSet::getLoadsForNode(const unsigned int elem_ID,
					  std::vector<const DirichletBoundaryCondition*>& loads) const
{
  std::multimap<unsigned int, DirichletBoundaryCondition*>::const_iterator it, end;
  it = this->node_ID_to_load_map.lower_bound(elem_ID);
  end = this->node_ID_to_load_map.upper_bound(elem_ID);
  for ( ; it != end; it++)
    loads.push_back(it->second);
}




void
BoundaryConditionLoadSet::getAllLoads(std::vector<const DirichletBoundaryCondition*>& loads) const
{
  std::multimap<unsigned int, DirichletBoundaryCondition*>::const_iterator it, end;
  it = this->node_ID_to_load_map.begin();
  end = this->node_ID_to_load_map.end();
  loads.resize(this->node_ID_to_load_map.size());
  unsigned int i=0;
  for (; it != end; it++)
    {
      loads[i] = it->second;
      i++;
    }
}



std::istream& 
BoundaryConditionLoadSet::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int n_loads=0, node_ID = 0, load_class_enum_ID = 0;

  // the load set should begin with the BEGIN tag, and followed by 
  // the number of loads in the set
  FESystemIOUtility::readFromInput(input, BOUNDARY_CONDITION_LOAD_SET::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
  
  
  // next, it should contain the name, ID and kind, and number of DOFs of load
  FESystemIOUtility::readFromInput(input, "NAME", this->set_tag);
  FESystemIOUtility::readFromInput(input, "ID", this->load_set_ID);
  
  tag.clear();
  FESystemIOUtility::readFromInput(input, "KIND", tag);
  this->load_kind_enum_ID = LoadNameEnum::enumID(tag);
  load_class_enum_ID = getLoadClassEnumIDForLoadName(this->load_kind_enum_ID);

  FESystemIOUtility::readFromInput(input, "N_LOADS", n_loads);
  
  
  std::auto_ptr<DirichletBoundaryCondition> load;
  
  for (unsigned int load_incr=0; load_incr < n_loads; load_incr++)
    {
      load.reset(createBoundaryConditionLoad(load_class_enum_ID).release());
      load->readFromInputStream(input);
      
      node_ID = load->getNodeID();
      this->node_ID_to_load_map.insert
        (std::multimap<unsigned int, DirichletBoundaryCondition*>::value_type
         (node_ID, load.release()));
    }
  
  
  FESystemIOUtility::readFromInput(input, BOUNDARY_CONDITION_LOAD_SET::name());
  FESystemIOUtility::readFromInput(input, "END");
  
  return input;
}





SurfaceLoadSet::SurfaceLoadSet():
  LoadSetBase()
{
  
}



SurfaceLoadSet::~SurfaceLoadSet()
{
  // iterate over all loads and delete them
  std::multimap<unsigned int, SurfaceLoad*>::iterator it, end;
  it = this->elem_ID_to_load_map.begin();
  end = this->elem_ID_to_load_map.end();
  
  for (; it != end; it++)
    {
      delete it->second;
      it->second = NULL;
    }

  this->elem_ID_to_load_map.clear();
}



unsigned int 
SurfaceLoadSet::getLoadSetKindEnumID() const
{
  return SURFACE_LOAD_SET::num();	
}




void
SurfaceLoadSet::getLoadsForElement(const unsigned int elem_ID,
				   std::vector<const SurfaceLoad*>& loads) const
{
  std::pair<std::multimap<unsigned int, SurfaceLoad*>::const_iterator, 
    std::multimap<unsigned int, SurfaceLoad*>::const_iterator> it;
  it = this->elem_ID_to_load_map.equal_range(elem_ID);

  for (; it.first != it.second; it.first++)
    loads.push_back(it.first->second);
}



void
SurfaceLoadSet::getLoadsForElementSide(const unsigned int elem_ID,
				       const unsigned int side_num,
				       std::vector<const SurfaceLoad*>& loads) const
{
  std::pair<std::multimap<unsigned int, SurfaceLoad*>::const_iterator, 
    std::multimap<unsigned int, SurfaceLoad*>::const_iterator> it;
  it = this->elem_ID_to_load_map.equal_range(elem_ID);

  for (; it.first != it.second; it.first++)
    {
      if (it.first->second->getSurfaceID() == side_num)
	loads.push_back(it.first->second);
    }
}




void
SurfaceLoadSet::getAllLoads(std::vector<const SurfaceLoad*>& loads) const
{
  std::multimap<unsigned int, SurfaceLoad*>::const_iterator it, end;
  it = this->elem_ID_to_load_map.begin();
  end = this->elem_ID_to_load_map.end();
  loads.resize(this->elem_ID_to_load_map.size());
  unsigned int i=0;
  for (; it != end; it++)
    {
      loads[i] = it->second;
      i++;
    }
}



std::istream& 
SurfaceLoadSet::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int n_loads=0,  elem_ID = 0, load_class_enum_ID = 0;
	
  // the load set should begin with the BEGIN tag, and followed by 
  // the number of loads in the set
  FESystemIOUtility::readFromInput(input, SURFACE_LOAD_SET::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
	
	
  // next, it should contain the name, ID and kind, and number of DOFs of load
  FESystemIOUtility::readFromInput(input, "NAME", this->set_tag);
  FESystemIOUtility::readFromInput(input, "ID", this->load_set_ID);
  
  tag.clear();
  FESystemIOUtility::readFromInput(input, "KIND", tag);
  this->load_kind_enum_ID = LoadNameEnum::enumID(tag);
  load_class_enum_ID = getLoadClassEnumIDForLoadName(this->load_kind_enum_ID);
  
  FESystemIOUtility::readFromInput(input, "N_LOADS", n_loads);
	
  
  std::auto_ptr<SurfaceLoad> load;
  
  for (unsigned int load_incr=0; load_incr < n_loads; load_incr++)
    {
      load.reset(createSurfaceLoad(load_class_enum_ID).release());
      load->readFromInputStream(input);
      elem_ID = load->getElemID();
      this->elem_ID_to_load_map.insert
	(std::multimap<unsigned int, SurfaceLoad*>::value_type
	 (elem_ID, load.release()));
    }
	
	
  FESystemIOUtility::readFromInput(input, SURFACE_LOAD_SET::name());
  FESystemIOUtility::readFromInput(input, "END");
	
  return input;
}




AeroelasticLoadSet::AeroelasticLoadSet():
LoadSetBase(),
altitude(0.0)
{
  
}


AeroelasticLoadSet::~AeroelasticLoadSet()
{
  
}


unsigned int
AeroelasticLoadSet::getLoadSetKindEnumID() const
{
  return AEROELASTIC_LOAD_SET::num();
}


const std::vector<double>&
AeroelasticLoadSet::getMachNumbers() const
{
  Assert(this->mach_numbers.size() > 0, ExcInternalError());
  return this->mach_numbers;
}


const std::vector<double>& 
AeroelasticLoadSet::getDynamicPressures() const
{
  Assert(this->dynamic_pressures.size() > 0, ExcInternalError());
  return this->dynamic_pressures;
}


std::istream& 
AeroelasticLoadSet::readFromInputStream(std::istream& input)
{
  unsigned int n_mach_numbers=0;
	
  // the load set should begin with the BEGIN tag, and followed by 
  // the number of loads in the set
  FESystemIOUtility::readFromInput(input, AEROELASTIC_LOAD_SET::name());
  FESystemIOUtility::readFromInput(input, "BEGIN");
	
	
  // next, it should contain the name, ID and kind, and number of DOFs of load
  FESystemIOUtility::readFromInput(input, "NAME", this->set_tag);
  FESystemIOUtility::readFromInput(input, "ID", this->load_set_ID);
  
  FESystemIOUtility::readFromInput(input, "ALTITUDE", this->altitude);
  
  FESystemIOUtility::readFromInput(input, "N_MACH_NUMBERS", n_mach_numbers);
	
  this->mach_numbers.resize(n_mach_numbers);
  
  for (unsigned int i=0; i<n_mach_numbers; i++)
    input >> mach_numbers[i];
  
  FESystemIOUtility::readFromInput(input, AEROELASTIC_LOAD_SET::name());
  FESystemIOUtility::readFromInput(input, "END");
	
  return input;  
}




std::auto_ptr<LoadSetBase> 
createLoadSet(const unsigned int enum_ID)
{
  
  std::auto_ptr<LoadSetBase>  return_ptr;
  
  switch (enum_ID)
    {
    case VOLUME_LOAD_SET_ENUM_ID:
      return_ptr.reset(new VolumeLoadSet);
      break;
      
    case SURFACE_LOAD_SET_ENUM_ID:
      return_ptr.reset(new SurfaceLoadSet);
      break;

    case NODAL_LOAD_SET_ENUM_ID:
      return_ptr.reset(new NodalLoadSet);
      break;

    case BOUNDARY_CONDITION_LOAD_SET_ENUM_ID:
      return_ptr.reset(new BoundaryConditionLoadSet);
      break;
      
    case AEROELASTIC_LOAD_SET_ENUM_ID:
      return_ptr.reset(new AeroelasticLoadSet);
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
  
  return return_ptr;
}

