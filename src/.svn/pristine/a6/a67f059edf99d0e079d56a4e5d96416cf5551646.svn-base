// $Id: FESystemController.h,v 1.20.4.4 2007-06-13 14:55:51 manav Exp $

#ifndef __fesystem_controller_h__
#define __fesystem_controller_h__

// C++ includes
#include <memory>

// FESystem includes
#include "Utilities/NameEnumHandler.h"

// MPI includes
#include "mpi.h"

// libMesh includes
#include "parameters.h"

#ifdef WITH_MATHEMATICA
// mathlink include
#include "mathlink.h"
#endif

#ifndef PETSC_PACKAGE_ENUM_ID 
#define PETSC_PACKAGE_ENUM_ID  1
#else
#error
#endif


#ifndef PETSC_PACKAGE_ENUM_NAME
#define PETSC_PACKAGE_ENUM_NAME "PETSC_PACKAGE_ENUM_NAME"
#else
#error
#endif


// Forward declerations

class LoadDatabase;
class ElementDataStorage;
class FunctionDatabase;
class OutputProcessor;
class PostProcessQtyDatabase;
class Log;
class ElemSetList;
class LibMeshInit;

namespace FESystemDatabase
{
  class GlobalDataStorage;
}

namespace  FESystemUtility
{
  class TimeLogs;
}

namespace FESystem
{
  class AnalysisCase;
}

namespace Driver
{
  class AnalysisDriver;
}

namespace Discipline
{
  class AnalysisDisciplineBase;
  class DisciplineInfo;
}

namespace MeshDS
{
  class MeshList;
}

namespace Property
{
  class PropertyDatabase;
}

namespace DesignData
{
  class DesignDatabase;
}

namespace Solver
{
  class FESystemSolverBase;
}


namespace FESystem
{
  DeclareEnumClass(PackageNameEnum);
  
  DeclareEnumName(PETSC_PACKAGE, FESystem::PackageNameEnum,
                  PETSC_PACKAGE_ENUM_ID,
                  PETSC_PACKAGE_ENUM_NAME);
  

  
  /// world MPI communicator for all processors
  extern MPI_Comm COMM_WORLD;
  
  /// MPI info
  extern MPI_Info MPI_INFO;
  
  /// total number of processors in this communicator
  extern unsigned int total_processors;
  
  /// local processor ID
  extern unsigned int local_processor;
  
#ifdef WITH_MATHEMATICA
  /// math link data
  extern MLENV mathlink_environment;
  
  /// mathlink link to the kernerl
  extern MLINK mathlink_link;
#endif
  
  
  class FESystemController
    {
public:
      
      FESystemController(int , char**);
      ~FESystemController();
      
      // FESystem data objects
      
      /// data structure to store the mesh used for FE analysis
      std::auto_ptr<MeshDS::MeshList> mesh_list;
      
      /// stores the loads used in the FE analysis
      std::auto_ptr<LoadDatabase> load_database;
      
      /// stores element level data
      std::auto_ptr<ElementDataStorage> element_data_storage;
      
      /// stores global matrices and vectors
      std::auto_ptr<FESystemDatabase::GlobalDataStorage> global_data_storage;
      
      /// stores the post process quantities
      std::auto_ptr<PostProcessQtyDatabase> post_process_qty_database;
      
      /// stores property cards
      std::auto_ptr<Property::PropertyDatabase> property_database;

      /// stores design data
      std::auto_ptr<DesignData::DesignDatabase> design_database;

      /// stores function data
      std::auto_ptr<FunctionDatabase> function_database;

      /// element set list
      std::auto_ptr<ElemSetList> elem_set_list;
      
      /// stores analysis details 
      std::auto_ptr<FESystem::AnalysisCase> analysis_case;
      
      /// output processor
      std::auto_ptr<OutputProcessor> output_processor;
      
      /// object to log messages
      std::auto_ptr<Log> log_file;
            
      /// Performance logging 
      std::auto_ptr<FESystemUtility::TimeLogs> performance_logging;
      
      /// this function will read from the input file, and initialize for the
      /// solution process
      void readAndInitialize();
      
      /// this function will solve the finite element system defined in the input 
      /// files
      void solve();
      
      /// this will write the outputs for a plotting package
      void writeOutput();
      
      // this is to be used to get the analysis discipline from the fesystem controller
      Discipline::AnalysisDisciplineBase* 
        getAnalysisDiscipline(const unsigned int discipline_enum_ID);
      
      /// @returns the analysis driver for the discipline
      Driver::AnalysisDriver* getAnalysisDriver(const unsigned int driver_enum_ID);
            
      /// initialized the property cards with the global parameters
      void initPropertyCardsForGlobalParameters();

protected:	

      /// initializes the object for computations. This method also initializes the 
      /// associated libraries
      void initialize(int argc , char** argv);
      
      /// initializes the object for computations. This method also initializes the 
      /// associated libraries
      void finalize();

      /// name of the input file for this analysis
      std::string input_file_name;
      
      /// the libmesh initialization object
      std::auto_ptr<LibMeshInit> libmesh_init;
      
      /// analysis driver map for their kind enum IDs
      std::map<unsigned int, Driver::AnalysisDriver*> analysis_driver_map;
      
      /// analysis discipline map for their kind enum IDs
      std::map<unsigned int, Discipline::AnalysisDisciplineBase*> analysis_discipline_map;
      
//      /// solver map for their kind enum IDs
//      std::map<unsigned int, Solver::FESystemSolverBase*> solver_map;
    };
  
}


#endif // __fesystem_controller_h__
