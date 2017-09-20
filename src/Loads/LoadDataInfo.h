// $Id: LoadDataInfo.h,v 1.1.2.1 2007-05-15 20:38:52 manav Exp $

#ifndef __fesystem_load_data_info_h__
#define __fesystem_load_data_info_h__


// C++ includes
#include <vector>
#include <memory>

// FESystem includes


namespace Loads
{
  enum LoadDataInfoScope
    {
      SINGLE,
      ALL,
      INVALID_SCOPE
    };
  
  
  class LoadDataInfoBase
    {
    public:
      LoadDataInfoBase();
      
      virtual ~LoadDataInfoBase();

      virtual void clear();
      
      void setQtyScope(const Loads::LoadDataInfoScope enum_ID);
      
      Loads::LoadDataInfoScope getQtyScope() const;
      
      void setLoadCaseID(const unsigned int ID);

      unsigned int getLoadCaseID() const;

      void setLoadClassEnumID(const unsigned int enum_ID);

      unsigned int getLoadClassEnumID() const;

      void setLoadNameEnumID(const unsigned int enum_ID);

      unsigned int getLoadNameEnumID() const;
      
      bool ifSensitivity() const;

      void setDVID(const unsigned int ID);

      unsigned int getDVID() const;

      bool ifTimeDependent() const;

      void setTime(const double t_val);

      double getTime() const;

    protected:
      
      /// stated if the quantity is for a single load or the entire set of loads
      Loads::LoadDataInfoScope qty_scope;
      
      /// load case ID 
      unsigned int load_case_ID;

      /// enum ID of the load for the discipline
      unsigned int load_name_enum_ID;
      
      /// enum ID of the load type that has been specified, for example, SurfaceRadiationLoad, etc.
      unsigned int load_type_enum_ID;
      
      /// boolean to check if the load specification is for sensitivity
      bool if_sensitivity;
      
      /// design variable of ID
      unsigned int DV_ID;
      
      /// boolean to check if the load is time dependent
      bool if_time_dependent;
      
      /// time value
      double time_value;
    };
  
  
  class ElemLoadDataInfo: public Loads::LoadDataInfoBase
    {
    public:
      ElemLoadDataInfo();

      virtual ~ElemLoadDataInfo();

      virtual void clear();

      void setElemID(const unsigned int ID);

      unsigned int getElemID() const;
      
    protected:

      
      unsigned int elem_ID;
    };


  class VolumeLoadDataInfo: public Loads::ElemLoadDataInfo
    {
    public:
      VolumeLoadDataInfo();

      virtual ~VolumeLoadDataInfo();
      
    protected:
      

    };


  class SurfaceLoadDataInfo: public Loads::ElemLoadDataInfo
    {
    public:
      SurfaceLoadDataInfo();

      virtual ~SurfaceLoadDataInfo();

      virtual void clear();

      void setSurfaceID(const unsigned int ID);

      unsigned int getSurfaceID() const;
      
    protected:

      unsigned int surface_ID;
    };


  class NodalLoadDataInfo: public Loads::LoadDataInfoBase
    {
    public:
      NodalLoadDataInfo();

      virtual ~NodalLoadDataInfo();

      virtual void clear();

      void setNodeID(const unsigned int ID);

      unsigned int getNodeID() const;
      
      void setNDofs(const unsigned int ID);

      unsigned int getNDofs() const;
      
    protected:
      
      unsigned int n_dofs;
      
      unsigned int node_ID;
    };
  
  
  class DirichletBoundaryConditionDataInfo: public Loads::LoadDataInfoBase
    {
    public:
      DirichletBoundaryConditionDataInfo();

      virtual ~DirichletBoundaryConditionDataInfo();
      
      virtual void clear();

      void setNodeID(const unsigned int ID);

      unsigned int getNodeID() const;
      
    protected:
      
      unsigned int node_ID;
    };


  std::auto_ptr<Loads::VolumeLoadDataInfo>
    createVolumeLoadDataInfo(const unsigned int load_type_enum_ID);

  std::auto_ptr<Loads::NodalLoadDataInfo>
    createNodalLoadDataInfo(const unsigned int load_type_enum_ID);
  
  std::auto_ptr<Loads::SurfaceLoadDataInfo>
    createSurfaceLoadDataInfo(const unsigned int load_type_enum_ID);
  
}


#endif // __fesystem_load_data_info_h__
