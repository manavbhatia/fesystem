// $Id: LoadedSolution.h,v 1.2.4.3 2007-05-09 18:45:36 manav Exp $

#ifndef __fesystem_loaded_solution_h__
#define __fesystem_loaded_solution_h__

// C++ includes
#include <string>
#include <memory>

// FESystem includes
#include "FESystem/FESystemNumbers.h"
#include "FESystem/FESystemExceptions.h"

// libMesh includes
#include "numerics/numeric_vector.h"

// forward declerations
namespace FESystemDatabase
{
  class GlobalDataStorage;
  class DataInfoBase;
}


namespace FESystemUtility
{
  /// this is a data structure to store the vector and its relevant data
  class LoadedVectorData
  {
public:
    /// constructor
    LoadedVectorData():
    vec(NULL)
    {
      this->vec = NumericVector<double>::build().release();
    }

    
    /// destructor
    ~LoadedVectorData()
    {
      delete this->vec;
    }

    
    void clear()
    {
      if (this->vec->initialized())
        this->vec->zero();
      this->data_info.reset();
    }

    
    /// @returns true if the data info has been set
    inline bool initialized(){
      bool init = false;
      if (this->data_info.get() != NULL)
        init = true;
      return init;}
    
    
    /// @returns the load case for this vector
    inline const FESystemDatabase::DataInfoBase& getDataInfo(){
      Assert(this->initialized(), ExcInternalError());
      return *(this->data_info.get());}
    
    // this stores a copy of the data info
    void setDataInfo(const FESystemDatabase::DataInfoBase& data_info_obj);
    
    /// @returns a reference to the vector
    inline NumericVector<double>& getVector(){
      return *(this->vec);}

    
protected:

    // data info object to store info about this vector
    std::auto_ptr<FESystemDatabase::DataInfoBase> data_info;
    
    // vector 
    NumericVector<double>* vec;
  };
  

  
  
  /// this class is a data structure that stores loaded solution from the database
  class LoadedSolution
  {
public: 
    /// constructor
    LoadedSolution();
    
    /// destructor 
    ~LoadedSolution();
    
    /// clears the data
    void clear();
    
    /// loads the solution from the database. The method first checks if the
    /// loaded solution is the same as the given solution. If it is, it simply returns, 
    /// otherwise it loads the solution
    const NumericVector<double>& loadSolutionFromDatabase
      (FESystemDatabase::GlobalDataStorage& database,
       const FESystemDatabase::DataInfoBase& data_info);
    
protected:

    /// defines the maximum number of vectors that can be loaded
    const unsigned int maximum_vectors;
      
    /// location of the previous vector that was loaded. this is used to decide
    /// which vector to overwrite when a new vector has to be loaded
    unsigned int previous_vec;
    
    /// vector of loaded data
    std::vector<FESystemUtility::LoadedVectorData*> loaded_vectors;
  };
}

#endif // __fesystem_loaded_solution_h__
