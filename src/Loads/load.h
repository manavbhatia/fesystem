// $Id: load.h,v 1.8.4.3 2008-08-21 00:37:02 manav Exp $

#ifndef __fesystem_load_h__
#define __fesystem_load_h__

// C++ includes
#include <iostream>
#include <set>
#include <vector>
#include <memory>

// FESystem includes
#include "Utilities/NameEnumHandler.h"
#include "FESystem/FESystemNumbers.h"
#include "FESystem/FESystemExceptions.h"

//#ifndef NODAL_LOAD_ENUM_ID 
//#define NODAL_LOAD_ENUM_ID 1
//#else
//#error
//#endif
//
//#ifndef NODAL_LOAD_ENUM_NAME 
//#define NODAL_LOAD_ENUM_NAME "NODAL_LOAD"
//#else
//#error
//#endif
//
//#ifndef SURFACE_LOAD_ENUM_ID 
//#define SURFACE_LOAD_ENUM_ID 2
//#else
//#error
//#endif
//
//#ifndef SURFACE_LOAD_ENUM_NAME 
//#define SURFACE_LOAD_ENUM_NAME "SURFACE_LOAD"
//#else
//#error
//#endif
//
//
//#ifndef VOLUME_LOAD_ENUM_ID 
//#define VOLUME_LOAD_ENUM_ID 3
//#else
//#error
//#endif
//
//#ifndef VOLUME_LOAD_ENUM_NAME 
//#define VOLUME_LOAD_ENUM_NAME "VOLUME_LOAD"
//#else
//#error
//#endif


#ifndef NODAL_POINT_LOAD_ENUM_ID 
#define NODAL_POINT_LOAD_ENUM_ID 1
#else
#error
#endif

#ifndef NODAL_POINT_LOAD_ENUM_NAME 
#define NODAL_POINT_LOAD_ENUM_NAME "NODAL_POINT_LOAD"
#else
#error
#endif


#ifndef DIRICHLET_BOUNDARY_CONDITION_ENUM_ID 
#define DIRICHLET_BOUNDARY_CONDITION_ENUM_ID 2
#else
#error
#endif

#ifndef DIRICHLET_BOUNDARY_CONDITION_ENUM_NAME 
#define DIRICHLET_BOUNDARY_CONDITION_ENUM_NAME "DIRICHLET_BOUNDARY_CONDITION"
#else
#error
#endif


#ifndef SCALAR_VOLUME_LOAD_ENUM_ID 
#define SCALAR_VOLUME_LOAD_ENUM_ID 3
#else
#error
#endif

#ifndef SCALAR_VOLUME_LOAD_ENUM_NAME 
#define SCALAR_VOLUME_LOAD_ENUM_NAME "SCALAR_VOLUME_LOAD"
#else
#error
#endif


#ifndef SURFACE_RADIATION_LOAD_ENUM_ID 
#define SURFACE_RADIATION_LOAD_ENUM_ID 4
#else
#error
#endif

#ifndef SURFACE_RADIATION_LOAD_ENUM_NAME 
#define SURFACE_RADIATION_LOAD_ENUM_NAME "SURFACE_RADIATION_LOAD"
#else
#error
#endif


#ifndef SCALAR_SURFACE_LOAD_ENUM_ID 
#define SCALAR_SURFACE_LOAD_ENUM_ID 5
#else
#error
#endif

#ifndef SCALAR_SURFACE_LOAD_ENUM_NAME 
#define SCALAR_SURFACE_LOAD_ENUM_NAME "SCALAR_SURFACE_LOAD"
#else
#error
#endif



#ifndef SURFACE_CONVECTION_LOAD_ENUM_ID 
#define SURFACE_CONVECTION_LOAD_ENUM_ID 6
#else
#error
#endif

#ifndef SURFACE_CONVECTION_LOAD_ENUM_NAME 
#define SURFACE_CONVECTION_LOAD_ENUM_NAME "SURFACE_CONVECTION_LOAD"
#else
#error
#endif


#ifndef VECTOR_VOLUME_LOAD_ENUM_ID 
#define VECTOR_VOLUME_LOAD_ENUM_ID 7
#else
#error
#endif

#ifndef VECTOR_VOLUME_LOAD_ENUM_NAME 
#define VECTOR_VOLUME_LOAD_ENUM_NAME "VECTOR_VOLUME_LOAD"
#else
#error
#endif



DeclareEnumClass(LoadClassEnum);
DeclareEnumClass(LoadNameEnum);

//DeclareEnumName(NODAL_LOAD, LoadClassEnum,
//		NODAL_LOAD_ENUM_ID,
//		NODAL_LOAD_ENUM_NAME);
//
//DeclareEnumName(SURFACE_LOAD, LoadClassEnum,
//		SURFACE_LOAD_ENUM_ID,
//		SURFACE_LOAD_ENUM_NAME);
//
//DeclareEnumName(VOLUME_LOAD, LoadClassEnum,
//		VOLUME_LOAD_ENUM_ID,
//		VOLUME_LOAD_ENUM_NAME);

DeclareEnumName(NODAL_POINT_LOAD, LoadClassEnum,
		NODAL_POINT_LOAD_ENUM_ID,
		NODAL_POINT_LOAD_ENUM_NAME);


DeclareEnumName(DIRICHLET_BOUNDARY_CONDITION, LoadClassEnum,
		DIRICHLET_BOUNDARY_CONDITION_ENUM_ID,
		DIRICHLET_BOUNDARY_CONDITION_ENUM_NAME);

DeclareEnumName(SCALAR_VOLUME_LOAD, LoadClassEnum,
		SCALAR_VOLUME_LOAD_ENUM_ID,
		SCALAR_VOLUME_LOAD_ENUM_NAME);

DeclareEnumName(VECTOR_VOLUME_LOAD, LoadClassEnum,
                VECTOR_VOLUME_LOAD_ENUM_ID,
                VECTOR_VOLUME_LOAD_ENUM_NAME);

DeclareEnumName(SURFACE_RADIATION_LOAD, LoadClassEnum,
		SURFACE_RADIATION_LOAD_ENUM_ID,
		SURFACE_RADIATION_LOAD_ENUM_NAME);

DeclareEnumName(SCALAR_SURFACE_LOAD, LoadClassEnum,
		SCALAR_SURFACE_LOAD_ENUM_ID,
		SCALAR_SURFACE_LOAD_ENUM_NAME);


DeclareEnumName(SURFACE_CONVECTION_LOAD, LoadClassEnum,
		SURFACE_CONVECTION_LOAD_ENUM_ID,
		SURFACE_CONVECTION_LOAD_ENUM_NAME);


class LoadBase 
{
 public:
  
  LoadBase():
    load_ID(FESystemNumbers::InvalidID)
    { }

  virtual ~LoadBase()
    { }
	
  virtual unsigned int getLoadClassEnumID() const = 0;

//  virtual unsigned int getLoadNameEnumID() const = 0;

  /// return the ID of the load 
  unsigned int getLoadID() const
    { 
      Assert (this->load_ID != FESystemNumbers::InvalidID, ExcInternalError()); 
      return this->load_ID;
    }	
  
  virtual std::istream& readFromInputStream(std::istream& input) = 0;
  
 protected:
  
  unsigned int load_ID;
};



class NodalLoad: public LoadBase
{
 public:
  NodalLoad():
    LoadBase(),
    node_ID(FESystemNumbers::InvalidID)
    { }

  virtual ~NodalLoad()
    { }

//  virtual unsigned int getLoadClassEnumID() const 
//    {
//      return NODAL_LOAD::num();
//    }


//  void setNodeID(const unsigned int ID)
//    {
//      Assert(ID != FESystemNumbers::InvalidID, ExcInternalError());
//      this->node_ID = ID;
//    }

  unsigned int getNodeID() const 
    {
      Assert (this->node_ID != FESystemNumbers::InvalidID, ExcInternalError());
      return this->node_ID;
    }
  
   
 protected:
  
  unsigned int node_ID;
};




class NodalPointLoad: public NodalLoad
{
 public: 
  NodalPointLoad():
  NodalLoad(),
  n_dofs(FESystemNumbers::InvalidID)
  {
  }
  
  virtual ~NodalPointLoad()
    { }

  virtual unsigned int getLoadClassEnumID() const 
    {
      return NODAL_POINT_LOAD::num();
    }

  unsigned int getNDofs() const
  {
    Assert(this->n_dofs > 0, ExcInternalError());
    return this->n_dofs;
  }
  
  double getValue(const unsigned int i) const
    {
      Assert( i < this->n_dofs, ExcInternalError());
      return this->values[i];
    }
  
  virtual std::istream& readFromInputStream(std::istream& input) 
    {
      input >> this->load_ID 
	    >> this->node_ID;
      for (unsigned int i=0; i < this->n_dofs;  i++)
        input >> this->values[i];
      
      AssertThrow(this->load_ID != FESystemNumbers::InvalidID, ExcInternalError());
      AssertThrow(this->node_ID != FESystemNumbers::InvalidID, ExcInternalError());
      
      return input;
    }

  friend std::auto_ptr<NodalLoad> 
  createNodalLoad(const unsigned int nodal_load_kind_enum_ID, const unsigned int dim);

 protected:
  void setNDofs(const unsigned int n)
  {
    // make sure that this value has not already been set
    Assert(this->n_dofs == FESystemNumbers::InvalidID, ExcInternalError());
    Assert(n > 0, ExcInternalError());
    this->n_dofs = n;
    this->values.resize(this->n_dofs);
    for (unsigned int i=0; i < this->n_dofs; i++)
      this->values[i] = 0.0;
  }
  
  unsigned int n_dofs;
  
  std::vector<double> values;
};




class ElemLoad: public LoadBase
{
 public: 
  ElemLoad():
    LoadBase(),
    elem_ID(FESystemNumbers::InvalidID)
    {
      
    }

  virtual ~ElemLoad()
    {
      
    }

  unsigned int getElemID() const
    {
      Assert(this->elem_ID != FESystemNumbers::InvalidID, ExcInternalError());
      return this->elem_ID;
    }
  
 protected:

  unsigned int elem_ID;
};



class SurfaceLoad: public ElemLoad
{
 public: 
  SurfaceLoad():
    ElemLoad(),
    surface_id(FESystemNumbers::InvalidID)
    {

    }

  virtual ~SurfaceLoad()
    {
      
    }

//  virtual unsigned int getLoadClassEnumID() const 
//    {
//      return SURFACE_LOAD::num();
//    }

  unsigned int getSurfaceID() const
    {
      return this->surface_id;
    }

 protected:
  
  unsigned int surface_id;  
};



class ScalarSurfaceLoad: public SurfaceLoad
{
 public:
  ScalarSurfaceLoad(): 
    SurfaceLoad(),
    value(0.0)
    {
      
    }

  virtual ~ScalarSurfaceLoad()
    {
      
    }

  virtual unsigned int getLoadClassEnumID() const 
    {
      return SCALAR_SURFACE_LOAD::num();
    }


  double getValue() const
    {
      return this->value;
    }

  virtual std::istream& readFromInputStream(std::istream& input) 
    {
      input >> this->load_ID 
	    >> this->elem_ID
	    >> this->surface_id
	    >> this->value;
      
      AssertThrow(this->load_ID != FESystemNumbers::InvalidID, ExcInternalError());
      AssertThrow(this->elem_ID != FESystemNumbers::InvalidID, ExcInternalError());
      
      return input;
    }
  
 protected:

  double value;
};




class SurfaceConvectionLoad: public SurfaceLoad
{
 public:
  SurfaceConvectionLoad():
    SurfaceLoad(),
    temperature(0.0),
    convection_coeff(0.0)
    {
      
    }

  virtual unsigned int getLoadClassEnumID() const 
    {
      return SURFACE_CONVECTION_LOAD::num();
    }

  double getAmbientTemperature() const
    {
      return this->temperature;
    }
  

  double getConvectionCoeff() const
    {
      return this->convection_coeff;
    }


  virtual std::istream& readFromInputStream(std::istream& input) 
    {
      input >> this->load_ID 
	    >> this->elem_ID
	    >> this->surface_id
	    >> this->convection_coeff
	    >> this->temperature;

      AssertThrow(this->load_ID != FESystemNumbers::InvalidID, ExcInternalError());
      AssertThrow(this->elem_ID != FESystemNumbers::InvalidID, ExcInternalError());
      
      return input;
    }


 protected:

  double temperature;

  double convection_coeff;
};






class SurfaceRadiationLoad: public SurfaceLoad
{
 public: 
  SurfaceRadiationLoad():
    SurfaceLoad(),
    temperature(0.0)
    {
      
    }

  virtual ~SurfaceRadiationLoad()
    {
      
    }

  virtual unsigned int getLoadClassEnumID() const 
    {
      return SURFACE_RADIATION_LOAD::num();
    }


  double getTemperature() const
    {
      return this->temperature;
    }

  virtual std::istream& readFromInputStream(std::istream& input) 
    {
      input >> this->load_ID 
	    >> this->elem_ID
	    >> this->surface_id
	    >> this->temperature;

      AssertThrow(this->load_ID != FESystemNumbers::InvalidID, ExcInternalError());
      AssertThrow(this->elem_ID != FESystemNumbers::InvalidID, ExcInternalError());

      return input;
    }

 protected:
  
  double temperature;
};







class VolumeLoad: public ElemLoad
{
 public: 
  VolumeLoad():
    ElemLoad()
    {
      
    }
  
  virtual ~VolumeLoad()
    {
      
    }

//  virtual unsigned int getLoadClassEnumID() const 
//    {
//      return VOLUME_LOAD::num();
//    }

 protected:

};



class ScalarVolumeLoad: public VolumeLoad
{
 public: 
  ScalarVolumeLoad():
    VolumeLoad(),
    value(0.0)
    {
      
    }


  virtual ~ScalarVolumeLoad()
    {
      
    }
  

  virtual unsigned int getLoadClassEnumID() const 
    {
      return SCALAR_VOLUME_LOAD::num();
    }


  double getValue() const
    {
      return this->value;
    }

  virtual std::istream& readFromInputStream(std::istream& input) 
    {
      input >> this->load_ID 
	    >> this->elem_ID
	    >> this->value;

      AssertThrow(this->load_ID != FESystemNumbers::InvalidID, ExcInternalError());
      AssertThrow(this->elem_ID != FESystemNumbers::InvalidID, ExcInternalError());

      return input;
    }

 protected:
  
  double value;
};




//class VectorVolumeLoad: public VolumeLoad
//  {
//  public: 
//    VectorVolumeLoad():
//    VolumeLoad(),
//    ndofs(0)
//    {
//      
//    }
//    
//    
//    virtual ~VectorVolumeLoad()
//    {
//      
//    }
//    
//    
//    virtual unsigned int getLoadClassEnumID() const 
//    {
//      return VECTOR_VOLUME_LOAD::num();
//    }
//    
//    
//    double getValue() const
//    {
//      return this->value;
//    }
//    
//    virtual std::istream& readFromInputStream(std::istream& input) 
//    {
//      input >> this->load_ID 
//	    >> this->elem_ID
//	    >> this->value;
//      
//      AssertThrow(this->load_ID != FESystemNumbers::InvalidID, ExcInternalError());
//      AssertThrow(this->elem_ID != FESystemNumbers::InvalidID, ExcInternalError());
//      
//      return input;
//    }
//    
//  protected:
//    
//    unsigned int n_dofs;
//    std::vector<double> values;
//  };








class DirichletBoundaryCondition: public NodalLoad
{
 public:  
  DirichletBoundaryCondition():
    NodalLoad(),
    dof_number(FESystemNumbers::InvalidID),
    value(0.0)
    {

    }


  virtual ~DirichletBoundaryCondition()
    {

    }

  virtual unsigned int getLoadClassEnumID() const 
    {
      return DIRICHLET_BOUNDARY_CONDITION::num();
    }

  double getValue() const
    {
      return this->value;
    }

  unsigned int getDofNumber() const
    {
      Assert(this->dof_number != FESystemNumbers::InvalidID, ExcInternalError());
      return this->dof_number;
    }

  virtual std::istream& readFromInputStream(std::istream& input) 
    {
      input >> this->load_ID
	    >> this->node_ID
	    >> this->dof_number
	    >> this->value;

      AssertThrow(this->load_ID != FESystemNumbers::InvalidID, ExcInternalError());
      AssertThrow(this->node_ID != FESystemNumbers::InvalidID, ExcInternalError());
      AssertThrow(this->dof_number != FESystemNumbers::InvalidID, ExcInternalError());
      
      return input;
    }

 protected:

  unsigned int dof_number;

  double value;
};



std::auto_ptr<SurfaceLoad> createSurfaceLoad(const unsigned int load_kind_enum_ID);
std::auto_ptr<NodalLoad> createNodalLoad(const unsigned int nodal_load_kind_enum_ID,
                                         const unsigned int dim);
std::auto_ptr<VolumeLoad> createVolumeLoad(const unsigned int volume_load_kind_enum_ID);
std::auto_ptr<DirichletBoundaryCondition> 
createBoundaryConditionLoad(const unsigned int bc_load_kind_enum_ID);

unsigned int 
getLoadClassEnumIDForLoadName(unsigned int load_name_enum_ID);


#endif //__fesystem_load_h__
