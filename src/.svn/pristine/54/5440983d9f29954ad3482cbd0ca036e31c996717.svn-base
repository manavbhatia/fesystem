// $Id: LoadDataInfo.C,v 1.1.2.1 2007-05-15 20:38:52 manav Exp $

// C++ includes


// FESystem includes
#include "Loads/LoadDataInfo.h"
#include "FESystem/FESystemNumbers.h"
#include "FESystem/FESystemExceptions.h"


Loads::LoadDataInfoBase::LoadDataInfoBase():
qty_scope(Loads::INVALID_SCOPE),
load_case_ID(FESystemNumbers::InvalidID),
load_name_enum_ID(FESystemNumbers::InvalidID),
load_type_enum_ID(FESystemNumbers::InvalidID),
if_sensitivity(false),
DV_ID(FESystemNumbers::InvalidID),
if_time_dependent(false),
time_value(0.0)
{
  
}
      
  
Loads::LoadDataInfoBase::~LoadDataInfoBase()
{
  this->clear();
}

      

void
Loads::LoadDataInfoBase::clear()
{
  this->qty_scope = Loads::INVALID_SCOPE;
  this->load_case_ID = FESystemNumbers::InvalidID;
  this->load_name_enum_ID = FESystemNumbers::InvalidID;
  this->load_type_enum_ID = FESystemNumbers::InvalidID;
  this->if_sensitivity = false;
  this->DV_ID = FESystemNumbers::InvalidID;
  this->if_time_dependent = false;
  this->time_value = 0.0;
}


      
void 
Loads::LoadDataInfoBase::setQtyScope(const Loads::LoadDataInfoScope enum_ID)
{
  Assert(enum_ID != Loads::INVALID_SCOPE, ExcInternalError());
  Assert(this->qty_scope == Loads::INVALID_SCOPE, ExcInternalError());
  this->qty_scope = enum_ID;
}

      

Loads::LoadDataInfoScope 
Loads::LoadDataInfoBase::getQtyScope() const
{
  Assert(this->qty_scope != Loads::INVALID_SCOPE, ExcInternalError());
  return this->qty_scope;
}

      

void
Loads::LoadDataInfoBase::setLoadCaseID(const unsigned int ID)
{
  Assert(ID != FESystemNumbers::InvalidID, ExcInternalError());
  Assert(this->load_case_ID == FESystemNumbers::InvalidID, ExcInternalError());
  this->load_case_ID = ID;
}
      


unsigned int
Loads::LoadDataInfoBase::getLoadCaseID() const
{
  Assert(this->load_case_ID != FESystemNumbers::InvalidID, ExcInternalError());
  return this->load_case_ID;
}

      

void
Loads::LoadDataInfoBase::setLoadClassEnumID(const unsigned int enum_ID)
{
  Assert(enum_ID != FESystemNumbers::InvalidID, ExcInternalError());
  Assert(this->load_type_enum_ID == FESystemNumbers::InvalidID, ExcInternalError());
  this->load_type_enum_ID = enum_ID;
}
      

unsigned int
Loads::LoadDataInfoBase::getLoadClassEnumID() const
{
  Assert(this->load_type_enum_ID != FESystemNumbers::InvalidID, ExcInternalError());
  return this->load_type_enum_ID;  
}
      


void
Loads::LoadDataInfoBase::setLoadNameEnumID(const unsigned int enum_ID)
{
  Assert(enum_ID != FESystemNumbers::InvalidID, ExcInternalError());
  Assert(this->load_name_enum_ID == FESystemNumbers::InvalidID, ExcInternalError());
  this->load_name_enum_ID = enum_ID;
}
      


unsigned int
Loads::LoadDataInfoBase::getLoadNameEnumID() const
{
  Assert(this->load_name_enum_ID != FESystemNumbers::InvalidID, ExcInternalError());
  return this->load_name_enum_ID;
}
      


bool
Loads::LoadDataInfoBase::ifSensitivity() const
{
  return this->if_sensitivity;
}
      


void
Loads::LoadDataInfoBase::setDVID(const unsigned int ID)
{
  Assert(ID != FESystemNumbers::InvalidID, ExcInternalError());
  Assert(this->DV_ID == FESystemNumbers::InvalidID, ExcInternalError());
  Assert(this->if_sensitivity == false, ExcInternalError());
  this->DV_ID = ID;
  this->if_sensitivity = true;
}
      


unsigned int
Loads::LoadDataInfoBase::getDVID() const
{
  Assert(this->DV_ID != FESystemNumbers::InvalidID, ExcInternalError());
  return this->DV_ID;
}
      


bool
Loads::LoadDataInfoBase::ifTimeDependent() const
{
  return this->if_time_dependent;
}
      


void
Loads::LoadDataInfoBase::setTime(const double t_val)
{
  Assert(t_val >= 0.0, ExcInternalError());
  Assert(this->time_value == 0.0, ExcInternalError());
  Assert(!this->if_time_dependent, ExcInternalError());
  
  this->time_value = t_val;
  this->if_time_dependent = true;
}
      


double
Loads::LoadDataInfoBase::getTime() const
{
  Assert(this->if_time_dependent, ExcInternalError());
  return this->time_value;
}
      




Loads::ElemLoadDataInfo::ElemLoadDataInfo():
Loads::LoadDataInfoBase(),
elem_ID(FESystemNumbers::InvalidID)
{
  
}
 


Loads::ElemLoadDataInfo::~ElemLoadDataInfo()
{
  this->clear();
}
      



void
Loads::ElemLoadDataInfo::clear()
{
  this->elem_ID = FESystemNumbers::InvalidID;
  Loads::LoadDataInfoBase::clear();
}
      


void
Loads::ElemLoadDataInfo::setElemID(const unsigned int ID)
{
  Assert(ID  != FESystemNumbers::InvalidID, ExcInternalError());
  Assert(this->elem_ID  == FESystemNumbers::InvalidID, ExcInternalError());
  this->elem_ID = ID;
}
      


unsigned int
Loads::ElemLoadDataInfo::getElemID() const
{
  Assert(this->elem_ID  != FESystemNumbers::InvalidID, ExcInternalError());
  return this->elem_ID;
}
      




Loads::VolumeLoadDataInfo::VolumeLoadDataInfo():
Loads::ElemLoadDataInfo()
{
  
}
 



Loads::VolumeLoadDataInfo::~VolumeLoadDataInfo()
{
  
}
      



Loads::SurfaceLoadDataInfo::SurfaceLoadDataInfo():
Loads::ElemLoadDataInfo(),
surface_ID(FESystemNumbers::InvalidID)
{
  
}


      

Loads::SurfaceLoadDataInfo::~SurfaceLoadDataInfo()
{
  this->clear();
}
 



void
Loads::SurfaceLoadDataInfo::clear()
{
  this->surface_ID = FESystemNumbers::InvalidID;
  Loads::ElemLoadDataInfo::clear();
}
      


void
Loads::SurfaceLoadDataInfo::setSurfaceID(const unsigned int ID)
{
  Assert(this->surface_ID == FESystemNumbers::InvalidID, ExcInternalError());
  this->surface_ID = ID;
}
      

unsigned int
Loads::SurfaceLoadDataInfo::getSurfaceID() const
{
  return this->surface_ID;
}
      



Loads::NodalLoadDataInfo::NodalLoadDataInfo():
Loads::LoadDataInfoBase(),
n_dofs(FESystemNumbers::InvalidID),
node_ID(FESystemNumbers::InvalidID)
{
  
}
      


Loads::NodalLoadDataInfo::~NodalLoadDataInfo()
{
  this->clear();
}
      


void
Loads::NodalLoadDataInfo::clear()
{
  this->n_dofs = FESystemNumbers::InvalidID;
  this->node_ID = FESystemNumbers::InvalidID;
  Loads::LoadDataInfoBase::clear();
}
      


void
Loads::NodalLoadDataInfo::setNodeID(const unsigned int ID)
{
  Assert(ID != FESystemNumbers::InvalidID, ExcInternalError());
  Assert(this->node_ID == FESystemNumbers::InvalidID, ExcInternalError());
  this->node_ID = ID;  
}
      


unsigned int
Loads::NodalLoadDataInfo::getNodeID() const
{
  Assert(this->node_ID != FESystemNumbers::InvalidID, ExcInternalError());
  return this->node_ID;  
}
      


void
Loads::NodalLoadDataInfo::setNDofs(const unsigned int ID)
{
  Assert(ID != FESystemNumbers::InvalidID, ExcInternalError());
  Assert(this->n_dofs == FESystemNumbers::InvalidID, ExcInternalError());
  this->n_dofs = ID;  
}
      


unsigned int
Loads::NodalLoadDataInfo::getNDofs() const
{
  Assert(this->n_dofs != FESystemNumbers::InvalidID, ExcInternalError());
  return this->n_dofs;
}
      



Loads::DirichletBoundaryConditionDataInfo::DirichletBoundaryConditionDataInfo():
Loads::LoadDataInfoBase(),
node_ID(FESystemNumbers::InvalidID)
{
  
}
 



Loads::DirichletBoundaryConditionDataInfo::~DirichletBoundaryConditionDataInfo()
{
  this->clear();
}



void
Loads::DirichletBoundaryConditionDataInfo::clear()
{
  this->node_ID = FESystemNumbers::InvalidID;
  Loads::LoadDataInfoBase::clear();
}
      


void
Loads::DirichletBoundaryConditionDataInfo::setNodeID(const unsigned int ID)
{
  Assert(ID != FESystemNumbers::InvalidID, ExcInternalError());
  Assert(this->node_ID == FESystemNumbers::InvalidID, ExcInternalError());
  this->node_ID = ID;  
}
      

unsigned int
Loads::DirichletBoundaryConditionDataInfo::getNodeID() const
{
  Assert(this->node_ID != FESystemNumbers::InvalidID, ExcInternalError());
  return this->node_ID;  
}
      

std::auto_ptr<Loads::VolumeLoadDataInfo>
Loads::createVolumeLoadDataInfo(const unsigned int load_type_enum_ID)
{
  std::auto_ptr<Loads::VolumeLoadDataInfo> return_ptr(new Loads::VolumeLoadDataInfo);
  return_ptr->setLoadClassEnumID(load_type_enum_ID);
  return return_ptr;
}
  

std::auto_ptr<Loads::NodalLoadDataInfo>
Loads::createNodalLoadDataInfo(const unsigned int load_type_enum_ID)
{
  std::auto_ptr<Loads::NodalLoadDataInfo> return_ptr(new Loads::NodalLoadDataInfo);
  return_ptr->setLoadClassEnumID(load_type_enum_ID);
  return return_ptr;
}
  

std::auto_ptr<Loads::SurfaceLoadDataInfo>
Loads::createSurfaceLoadDataInfo(const unsigned int load_type_enum_ID)
{
  std::auto_ptr<Loads::SurfaceLoadDataInfo> return_ptr(new Loads::SurfaceLoadDataInfo);
  return_ptr->setLoadClassEnumID(load_type_enum_ID);
  return return_ptr;
}

