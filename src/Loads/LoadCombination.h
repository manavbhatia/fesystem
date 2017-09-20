// $Id: LoadCombination.h,v 1.1.2.2 2007-06-13 14:58:05 manav Exp $

#ifndef __fesystem_load_combination_h__
#define __fesystem_load_combination_h__

// C++ includes
#include <iostream>
#include <set>
#include <vector>


// FESystem includes
#include "Loads/load.h"


namespace Loads
{
  
  class LoadCombinationBase 
    {
    public:
      
      LoadCombinationBase()
      { }
      
      virtual ~LoadCombinationBase()
      { }
      
      virtual void clear() = 0;
      
      virtual unsigned int nLoads() const = 0;
      
      virtual unsigned int getLoadClassEnumID() const = 0;
      
      //      virtual unsigned int getLoadNameEnumID() const = 0;
      
    protected:
      
    };
  
  
  
  class NodalLoadCombination: public LoadCombinationBase
    {
    public:
      NodalLoadCombination():
      LoadCombinationBase()
      { }
      
      virtual ~NodalLoadCombination()
      { }
      
      //      virtual unsigned int getLoadClassEnumID() const 
      //	{
      //	  return NODAL_LOAD::num();
      //	}
      
      virtual unsigned int nLoads() const
      {
        return this->load_set.size();
      }
      
      unsigned int getNodeID() const 
      {
        Assert (this->nLoads() > 0, ExcInternalError());
        return (*this->load_set.begin())->getNodeID();
      }
      
      virtual void addLoad(const double factor, const NodalLoad* load) = 0;
      
    protected:
      
      std::set<const NodalLoad*> load_set;
    };
  
  
  
  
  class NodalPointLoadCombination: public NodalLoadCombination
    {
    public: 
      NodalPointLoadCombination():
      NodalLoadCombination()
      {
      }
      
      virtual ~NodalPointLoadCombination()
      { }
      
      virtual unsigned int getLoadClassEnumID() const 
      {
        return NODAL_POINT_LOAD::num();
      }
      
      virtual void clear()
      {
        this->load_set.clear();
        this->values.clear();
      }
      
      unsigned int getNDofs() const
      {
        Assert(this->nLoads() > 0, ExcInternalError());
        return dynamic_cast<const NodalPointLoad*>(*this->load_set.begin())->getNDofs();
      }
      
      double getValue(const unsigned int i) const
      {
        Assert( i < this->getNDofs(), ExcInternalError());
        return this->values[i];
      }
      
      virtual void addLoad(const double factor, const NodalLoad* load_ptr) 
      {
        Assert(load_ptr != NULL, ExcInternalError());
        
        Assert(load_ptr->getLoadClassEnumID() == NODAL_POINT_LOAD::num(), 
               ExcInternalError());
        
        const NodalPointLoad* ptr = dynamic_cast<const NodalPointLoad*>(load_ptr);
        
        // set the node ID
        if (this->load_set.size() == 0)
          {
            this->values.resize(ptr->getNDofs());
            for (unsigned int i=0; i < ptr->getNDofs(); i++)
              this->values[i] = 0.0;
          }
        else
          Assert(this->getNodeID() == ptr->getNodeID() &&
                 this->getNDofs() == ptr->getNDofs(), ExcInternalError());
        
        bool insert_success = this->load_set.insert(load_ptr).second;
        // add this pointer to the set, and adjust the value of the load
        Assert(insert_success, ExcInternalError());
        
        for (unsigned int i=0; i < this->getNDofs(); i++)
          this->values[i] += factor * ptr->getValue(i);
      }
      
    protected:
      
      std::vector<double> values;
    };
  
  
  
  
  
  
  class SurfaceLoadCombination: public LoadCombinationBase
    {
    public: 
      SurfaceLoadCombination():
      LoadCombinationBase()
      {
        
      }
      
      virtual ~SurfaceLoadCombination()
      {
        
      }
      
      //      virtual unsigned int getLoadClassEnumID() const 
      //	{
      //	  return SURFACE_LOAD::num();
      //	}
      
      
      virtual unsigned int nLoads() const
      {
        return this->load_set.size();
      }
      
      unsigned int getElemID() const
      {
        Assert(this->nLoads() > 0, ExcInternalError());
        return (*this->load_set.begin())->getElemID();
      }
      
      unsigned int getSurfaceID() const
      {
        Assert(this->nLoads() > 0, ExcInternalError());
        return (*this->load_set.begin())->getSurfaceID();
      }
      
      virtual void addLoad(const double factor, const SurfaceLoad* load_ptr) =0;
      
    protected:
      
      std::set<const SurfaceLoad*> load_set;
    };
  
  
  
  
  
  class ScalarSurfaceLoadCombination: public SurfaceLoadCombination
    {
    public:
      ScalarSurfaceLoadCombination(): 
      SurfaceLoadCombination(),
      value(0.0)
      {
        
      }
      
      virtual ~ScalarSurfaceLoadCombination()
      {
        
      }
      
      
      virtual void clear() 
      {
        this->value = 0.0;
        this->load_set.clear();
      }
      
      virtual unsigned int getLoadClassEnumID() const 
      {
        return SCALAR_SURFACE_LOAD::num();
      }
      
      
      double getValue() const
      {
        return this->value;
      }
      
      
      virtual void addLoad(const double factor, const SurfaceLoad* load_ptr)
      {
        Assert(load_ptr != NULL, ExcInternalError());
        
        Assert(load_ptr->getLoadClassEnumID() == SCALAR_SURFACE_LOAD::num(),
               ExcInternalError());
        
        const ScalarSurfaceLoad* ptr = dynamic_cast<const ScalarSurfaceLoad*>(load_ptr);
        
        // set the node ID
        if (this->nLoads() > 0)
          Assert(this->getElemID() == ptr->getElemID() &&
                 this->getSurfaceID() == ptr->getSurfaceID(),
                 ExcInternalError());
        
        bool insert_success = this->load_set.insert(load_ptr).second;
        // add this pointer to the set, and adjust the value of the load
        Assert(insert_success, ExcInternalError());
        
        this->value += factor * ptr->getValue();
      }
      
    protected:
      
      double value;
    };
  
  
  
  
  
  class SurfaceConvectionLoadCombination: public SurfaceLoadCombination
    {
    public:
      SurfaceConvectionLoadCombination():
      SurfaceLoadCombination(),
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
      
      virtual void clear()
      {
        this->temperature = 0.0;
        this->convection_coeff = 0.0;
        this->load_set.clear();
      }
      
      virtual void addLoad(const double factor, const SurfaceLoad* load_ptr)
      {
        Assert(load_ptr != NULL, ExcInternalError());
        
        Assert(load_ptr->getLoadClassEnumID() == SURFACE_CONVECTION_LOAD::num(),
               ExcInternalError());
        
        const SurfaceConvectionLoad* ptr = dynamic_cast<const SurfaceConvectionLoad*>(load_ptr);
        
        // set the node ID
        if (this->load_set.size() > 0)
          Assert(this->getElemID() == ptr->getElemID() &&
                 this->getSurfaceID() == ptr->getSurfaceID(),
                 ExcInternalError());
        
        bool insert_success = this->load_set.insert(load_ptr).second;
        // add this pointer to the set, and adjust the value of the load
        Assert(insert_success, ExcInternalError());
        
        this->temperature += factor * ptr->getAmbientTemperature();
        this->convection_coeff += factor * ptr->getConvectionCoeff();
      }
      
      
    protected:
      
      double temperature;
      
      double convection_coeff;
    };
  
  
  
  
  
  
  class SurfaceRadiationLoadCombination: public SurfaceLoadCombination
    {
    public: 
      SurfaceRadiationLoadCombination():
      SurfaceLoadCombination(),
      temperature(0.0)
      {
        
      }
      
      virtual ~SurfaceRadiationLoadCombination()
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
      
      virtual void clear()
      {
        this->temperature = 0.0;
        this->load_set.clear();
      }
      
      virtual void addLoad(const double factor, const SurfaceLoad* load_ptr)
      {
        Assert(load_ptr != NULL, ExcInternalError());
        
        Assert(load_ptr->getLoadClassEnumID() == SURFACE_RADIATION_LOAD::num(),
               ExcInternalError());
        
        const SurfaceRadiationLoad* ptr = dynamic_cast<const SurfaceRadiationLoad*>(load_ptr);
        
        // set the node ID
        if (this->load_set.size() > 0)
          Assert(this->getElemID() == ptr->getElemID() &&
                 this->getSurfaceID() == ptr->getSurfaceID(),
                 ExcInternalError());
        
        bool insert_success = this->load_set.insert(load_ptr).second;
        // add this pointer to the set, and adjust the value of the load
        Assert(insert_success, ExcInternalError());
        
        this->temperature += factor * ptr->getTemperature();
      }
      
      
    protected:
      
      double temperature;
    };
  
  
  
  
  
  class VolumeLoadCombination: public LoadCombinationBase
    {
    public: 
      VolumeLoadCombination():
      LoadCombinationBase()
      {
        
      }
      
      virtual ~VolumeLoadCombination()
      {
        
      }
      
      //      virtual unsigned int getLoadClassEnumID() const 
      //	{
      //	  return VOLUME_LOAD::num();
      //	}
      
      virtual unsigned int nLoads() const
      {
        return this->load_set.size();
      }
      
      unsigned int getElemID() const
      {
        Assert(this->nLoads() > 0, ExcInternalError());
        return (*this->load_set.begin())->getElemID();
      }
      
      virtual void addLoad(const double factor, const VolumeLoad* load_ptr) =0;
      
    protected:
      
      std::set<const VolumeLoad*> load_set;
    };
  
  
  
  
  class ScalarVolumeLoadCombination: public VolumeLoadCombination
    {
    public: 
      ScalarVolumeLoadCombination():
      VolumeLoadCombination(),
      value(0.0)
      {
        
      }
      
      
      virtual ~ScalarVolumeLoadCombination()
      {
        
      }
      
      virtual void clear() 
      {
        this->value = 0.0;
        this->load_set.clear();
      }
      
      
      virtual unsigned int getLoadClassEnumID() const 
      {
        return SCALAR_VOLUME_LOAD::num();
      }
      
      
      double getValue() const
      {
        return this->value;
      }
      
      virtual void addLoad(const double factor, const VolumeLoad* load_ptr)
      {
        Assert(load_ptr != NULL, ExcInternalError());
        
        Assert(load_ptr->getLoadClassEnumID() == SCALAR_VOLUME_LOAD::num(),
               ExcInternalError());
        
        const ScalarVolumeLoad* ptr = dynamic_cast<const ScalarVolumeLoad*>(load_ptr);
        
        // set the node ID
        if (this->nLoads() > 0)
          Assert(this->getElemID() == ptr->getElemID(), ExcInternalError());
        
        bool insert_success = this->load_set.insert(load_ptr).second;
        // add this pointer to the set, and adjust the value of the load
        Assert(insert_success, ExcInternalError());
        
        this->value += factor * ptr->getValue();
      }
      
    protected:
      
      double value;
    };
  
  
  
  
  
  
  class DirichletBoundaryConditionCombination: public NodalLoadCombination
    {
    public:  
      DirichletBoundaryConditionCombination():
      NodalLoadCombination(),
      value(0.0)
      {
        
      }
      
      
      virtual ~DirichletBoundaryConditionCombination()
      {
        
      }
      
      virtual unsigned int getLoadClassEnumID() const 
      {
        return DIRICHLET_BOUNDARY_CONDITION::num();
      }
      
      void clear()
      {
        this->value = 0.0;
      }  
      
      double getValue() const
      {
        return this->value;
      }
      
      unsigned int getDofNumber() const
      {
        Assert(this->nLoads() > 0, ExcInternalError());
        return dynamic_cast<const DirichletBoundaryCondition*>(*this->load_set.begin())->getDofNumber();
      }
      
      void addLoad(const double factor, const NodalLoad* load_ptr)
      {
        Assert(load_ptr != NULL, ExcInternalError());
        
        Assert(load_ptr->getLoadClassEnumID() == DIRICHLET_BOUNDARY_CONDITION::num(),
               ExcInternalError());
        
        const DirichletBoundaryCondition* ptr = 
        dynamic_cast<const DirichletBoundaryCondition*>(load_ptr);
        
        // set the node ID
        if (this->nLoads() > 0)
          Assert(this->getNodeID() == ptr->getNodeID() &&
                 this->getDofNumber() == ptr->getDofNumber(), ExcInternalError());
        
        bool insert_success = this->load_set.insert(load_ptr).second;
        // add this pointer to the set, and adjust the value of the load
        Assert(insert_success, ExcInternalError());
        
        this->value += factor * ptr->getValue();
      }
      
    protected:
      
      double value;
    };
  
  
  std::auto_ptr<Loads::VolumeLoadCombination>
  createVolumeLoadCombination(const unsigned int load_type_enum_ID);
  
  std::auto_ptr<Loads::NodalLoadCombination>
  createNodalLoadCombination(const unsigned int load_type_enum_ID);
  
  std::auto_ptr<Loads::DirichletBoundaryConditionCombination>
  createDirichletBoundaryConditionCombination(const unsigned int load_type_enum_ID);
  
  std::auto_ptr<Loads::SurfaceLoadCombination>
  createSurfaceLoadCombination(const unsigned int load_type_enum_ID);
  
}


#endif //__fesystem_load_combination_h__
