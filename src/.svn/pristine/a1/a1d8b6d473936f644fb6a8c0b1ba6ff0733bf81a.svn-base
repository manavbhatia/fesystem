// $Id: FESystemController.C,v 1.27.4.7 2007-06-13 14:55:41 manav Exp $

// C++ includes
#include <fstream>

// FESystem inlcudes
#include "FESystem/FESystemController.h"
#include "FESystem/FESystemCommon.h"
#include "FESystem/AnalysisCase.h"
#include "Mesh/MeshList.h"
#include "Loads/LoadDatabase.h"
#include "Properties/PropertyDatabase.h"
#include "Database/ElementDataStorage.h"
#include "Database/GlobalDataStorage.h"
#include "Database/FunctionDatabase.h"
#include "OutputProcessor/OutputProcessor.h"
#include "PostProcess/PostProcessQtyDatabase.h"
#include "DesignData/DesignDatabase.h"
#include "Utilities/Log.h"
#include "Utilities/ElemSetList.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "AnalysisDriver/LinearAnalysisDriver.h"
#include "AnalysisDriver/NonLinearAnalysisDriver.h"
#include "AnalysisDriver/EigenProblemAnalysisDriver.h"
#include "AnalysisDriver/TransientAnalysisDriver.h"
#include "AnalysisDriver/RogerApproximationAeroelasticity.h"
#include "Discipline/ThermalDisciplineInfo.h"
#include "Discipline/StructuralDisciplineInfo.h"
#include "Discipline/PistonTheoryInfo.h"
#include "Solvers/FESystemSolverBase.h"
#include "Solutions/SolutionBase.h"
#include "Utilities/TimeLogs.h"
#include "Utilities/ParallelUtility.h"

// external library includes
#include "petsc.h"
#include "slepc.h"
#include "hdf5.h"


void
FESystem::FESystemController::initialize(int argc , char** argv)
{
  
  // first init the libraries
  // MPI needs to be initialized first
  {
    MPI_Init(&argc, &argv);
    MPI_Comm_dup (MPI_COMM_WORLD, &FESystem::COMM_WORLD);
    FESystem::MPI_INFO = MPI_INFO_NULL;
    
    int num = 0;
    // now get the total number of processors in the rank
    MPI_Comm_size(FESystem::COMM_WORLD, &num);
    FESystem::total_processors = static_cast<unsigned int>(num);
    // now get the local rank of the processor 
    num = 0;
    MPI_Comm_rank(FESystem::COMM_WORLD, &num);
    FESystem::local_processor = static_cast<unsigned int>(num);
  }
    
  // then libMesh
  this->libmesh_init.reset(new LibMeshInit(argc, argv, FESystem::COMM_WORLD));
    
  // then petsc
  PetscErrorCode ierr = 0;
  PETSC_COMM_WORLD = FESystem::COMM_WORLD;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // then slepc
  ierr = SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // then HDF5
  ierr = H5open();
  assert (ierr >= 0);
  
  // initialize the link to mathematica
#ifdef WITH_MATHEMATICA
  long err;
  int mathlink_argc = 4;
  char * mathlink_argv[5] = {"-linkmode", "launch", "-linkname", "/Applications/Mathematica\\ 5.2.app/Contents/MacOS/MathKernel -mathlink", NULL};

	mathlink_environment =  MLInitialize( (MLParametersPointer)0);
	if( mathlink_environment == (MLENV)0) exit(1);
  
	mathlink_link = MLOpenArgv( mathlink_environment, mathlink_argv, mathlink_argv + mathlink_argc, &err);
	if(mathlink_link == (MLINK)0) exit(2);
	
  // next, read the notebook
  MLPutFunction( FESystem::mathlink_link, "EvaluatePacket", 1L);
  MLPutFunction(FESystem::mathlink_link, "ToExpression",1L);
  MLPutString(FESystem::mathlink_link, "<</Users/manav/Documents/codes/fem/ShapeFactors/ff.m;");
  MLEndPacket( FESystem::mathlink_link);
  MLFlush(FESystem::mathlink_link);
    
  int pkt;
  
  while((pkt = MLNextPacket( FESystem::mathlink_link),pkt) && pkt != RETURNPKT) 
    {
    MLNewPacket( FESystem::mathlink_link);
    if( MLError( FESystem::mathlink_link))
      fprintf( stderr, "Error detected by MathLink: %s.\n",
               MLErrorMessage(FESystem::mathlink_link));
    }
  MLNewPacket(FESystem::mathlink_link);
#endif
}




void
FESystem::FESystemController::finalize()
{
  this->mesh_list.reset();
  this->load_database.reset();
  this->element_data_storage.reset();
  this->global_data_storage.reset();
  this->post_process_qty_database.reset();
  this->property_database.reset();
  this->design_database.reset();
  this->function_database.reset();
  this->elem_set_list.reset();
  this->analysis_case.reset();
  this->output_processor.reset();
  this->log_file.reset();
  this->performance_logging.reset();

  
  // now close all the libraries in the reverse order of their initialization
#ifdef WITH_MATHEMATICA
  // mathlink
  MLPutFunction( FESystem::mathlink_link, "Exit", 0);
  if( FESystem::mathlink_link) MLClose( FESystem::mathlink_link);
  if( FESystem::mathlink_environment) MLDeinitialize( FESystem::mathlink_environment);
#endif  
  
  // HDF5 first
  PetscErrorCode ierr = 0;
  ierr = H5close();
  assert (ierr >= 0);
  
  // then slepc 
  ierr = SlepcFinalize();
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // then petsc
  ierr = PetscFinalize();
  CHKERRABORT(FESystem::COMM_WORLD, ierr);
  
  // then libmesh
  this->libmesh_init.reset();
  
  // then MPI
  MPI_Finalize();  
}


FESystem::FESystemController::FESystemController(int argc , char** argv):
mesh_list(NULL),
load_database(NULL),
element_data_storage(NULL),
global_data_storage(NULL),
post_process_qty_database(NULL),
property_database(NULL),
design_database(NULL),
function_database(NULL),
elem_set_list(NULL),
analysis_case(NULL),
output_processor(NULL),
log_file(NULL),
performance_logging(NULL),
input_file_name()
{
  assert (argc > 1);
  this->input_file_name = argv[1];
  
  this->initialize(argc, argv);
    
  // start logging the runtime from this point
  this->performance_logging.reset(new FESystemUtility::TimeLogs);
}






FESystem::FESystemController::~FESystemController()
{	  
  // iterate over the analysis driver map and delete them
  {
    std::map<unsigned int, Driver::AnalysisDriver*>::iterator it, end;
    it = this->analysis_driver_map.begin();
    end = this->analysis_driver_map.end();
    
    for (; it != end; it++)
      {
      delete it->second;
      it->second = NULL;
      }
    this->analysis_driver_map.clear();
  }
  
  // iterate over the analysis discipline map and delete them
  {
    std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it, end;
    it = this->analysis_discipline_map.begin();
    end = this->analysis_discipline_map.end();
    
    for (; it != end; it++)
      {
      delete it->second;
      it->second = NULL;
      }
    
    this->analysis_discipline_map.clear();
  }
  
  
//  // iterate over the solver map and delete them
//  {
//    std::map<unsigned int, Solver::FESystemSolverBase*>::iterator it, end;
//    it = this->solver_map.begin();
//    end = this->solver_map.end();
//    
//    for (; it != end; it++)
//      {
//      delete it->second;
//      it->second = NULL;
//      }
//    
//    this->solver_map.clear();
//  }
  
  this->finalize();
}






void 
FESystem::FESystemController::readAndInitialize()
{
  
  this->performance_logging->setEvent("FESystemController::readAndInitialize()", 
                                         "FESystem Run Time");

  // create the hdf5 database file name
  std::string hdf_name = this->input_file_name;
  hdf_name += ".hdf5";
  
  // init all the independent data structures
  {
    this->mesh_list.reset(new MeshDS::MeshList);
    this->load_database.reset(new LoadDatabase(*this));
    this->element_data_storage.reset(new ElementDataStorage);
    this->global_data_storage.reset(new FESystemDatabase::GlobalDataStorage(*this, hdf_name));
    this->post_process_qty_database.reset(new PostProcessQtyDatabase);
    this->property_database.reset(new Property::PropertyDatabase(*this));
    this->design_database.reset(new DesignData::DesignDatabase);
    this->function_database.reset(new FunctionDatabase);
    this->elem_set_list.reset(new ElemSetList);
    this->analysis_case.reset(new FESystem::AnalysisCase);
    // the log file is created only on processor = 0
    if (FESystem::local_processor == 0)
      this->log_file.reset(new Log);
  }
  
  // create an input stream
  std::fstream  input_file;
  input_file.open(this->input_file_name.c_str(), std::fstream::in);
	
  // tell each system to read from the file one by one
  this->function_database->readFromInputStream(input_file);
  this->property_database->readFromInputStream(input_file);
  this->mesh_list->readFromInputStream(input_file);
  this->elem_set_list->readFromInputStream(input_file);
  this->load_database->readFromInputStream(input_file);
  this->analysis_case->readFromInputStream(input_file);
  this->design_database->readFromInputStream(input_file);
  
  const OutputInfo& output_info = this->analysis_case->getOutputInfo();
  if (FESystem::local_processor == 0)
    this->output_processor.reset(OutputProcessor::createOutputProcessor
                                 (*this,output_info).release());
  	
  // close the input stream
  input_file.close();

  this->performance_logging->unsetEvent("FESystemController::readAndInitialize()",
                                        "FESystem Run Time");
}






void 
FESystem::FESystemController::solve()
{	
  this->performance_logging->setEvent("FESystemController::solve()", 
                                         "FESystem Run Time");

  // get the solutions from the analysis case, and ask them to solve
  std::vector<Solution::SolutionBase*> solutions = 
  this->analysis_case->getSolutions();
  
  std::vector<Solution::SolutionBase*>::iterator it, end;
  it = solutions.begin();
  end = solutions.end();
  
  for ( ; it != end; it++)
    {
    (*it)->attachFESystemController(*this);
    (*it)->solve();
    }

  this->performance_logging->unsetEvent("FESystemController::solve()", 
                                        "FESystem Run Time");
}




Driver::AnalysisDriver*
FESystem::FESystemController::getAnalysisDriver(const unsigned int driver_enum_ID)
{
  std::map<unsigned int, Driver::AnalysisDriver*>::iterator it;
  it = this->analysis_driver_map.find(driver_enum_ID);
  
  // if it does not exist, create the driver
  if (it == this->analysis_driver_map.end())
    {
    
    Driver::AnalysisDriver* driver = NULL;
    
    switch (driver_enum_ID)
      {
      case LINEAR_ANALYSIS_DRIVER_ENUM_ID:
        // for now, the ID is being kept as the driver type enum ID
        driver = new Driver::LinearAnalysisDriver(driver_enum_ID,
                                                  *this);
        break;

      case NONLINEAR_ANALYSIS_DRIVER_ENUM_ID:
        // for now, the ID is being kept as the driver type enum ID
        driver = new Driver::NonLinearAnalysisDriver(driver_enum_ID,
                                                     *this);
        break;
        
      case EIGENPROBLEM_ANALYSIS_DRIVER_ENUM_ID:
        // for now, the ID is being kept as the driver type enum ID
        driver = new Driver::EigenProblemAnalysisDriver(driver_enum_ID,
                                                        *this);
        break;

      case LINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_ID:
        // for now, the ID is being kept as the driver type enum ID
        driver = new Driver::LinearTransientAnalysisDriver(driver_enum_ID,
                                                        *this);
        break;

      case NONLINEAR_TRANSIENT_ANALYSIS_DRIVER_ENUM_ID:
        // for now, the ID is being kept as the driver type enum ID
        driver = new Driver::NonlinearTransientAnalysisDriver(driver_enum_ID,
                                                           *this);
        break;

        case ROGER_APPROXIMATION_AEROELASTICITY_DRIVER_ENUM_ID:
          // for now, the ID is being kept as the driver type enum ID
          driver = new Driver::RogerApproximationAeroelasticityDriver(driver_enum_ID,
                                                                *this);
          break;

        default:
        Assert(false, 
               FESystemExceptions::ExcEnumCaseNotImplemented
               (Driver::AnalysisDriverTypeEnum::enumName(driver_enum_ID)));
      }
  
    bool insert_success = this->analysis_driver_map.insert
      (std::map<unsigned int, Driver::AnalysisDriver*>::value_type
       (driver_enum_ID, driver)).second;
  
    Assert(insert_success, ExcInternalError());
    
    return driver;
    }
  else
    return it->second;
}



Discipline::AnalysisDisciplineBase*
FESystem::FESystemController::getAnalysisDiscipline(const unsigned int discipline_enum_ID)
{
  const Discipline::DisciplineInfo& info = 
  this->analysis_case->getAnalysisDisciplineInfo(discipline_enum_ID);

  std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::iterator it;
  it = this->analysis_discipline_map.find(discipline_enum_ID);  
  
  // if it does not exist, create the driver
  if (it == this->analysis_discipline_map.end())
    {
    
    Discipline::AnalysisDisciplineBase* discipline = NULL;
    
    switch (discipline_enum_ID)
      {
      case THERMAL_DISCIPLINE_ENUM_ID:
        {
          const Discipline::ThermalDisciplineInfo& thermal_info =
          dynamic_cast<const Discipline::ThermalDisciplineInfo&>(info);
          discipline = new Discipline::ThermalAnalysis(*this, thermal_info);
        }
        break;
        
      case STRUCTURAL_DISCIPLINE_ENUM_ID:
        {
          const Discipline::StructuralDisciplineInfo& structural_info =
          dynamic_cast<const Discipline::StructuralDisciplineInfo&>(info);
          discipline = new Discipline::StructuralAnalysis(*this, structural_info);
        }
        break;
        
        case PISTON_THEORY_ENUM_ID:
        {
          const Discipline::PistonTheoryInfo& piston_theory_info =
          dynamic_cast<const Discipline::PistonTheoryInfo&>(info);
          discipline = new Discipline::PistonTheory(*this, piston_theory_info);
        }
          break;

        default:
        Assert(false, 
               FESystemExceptions::ExcEnumCaseNotImplemented
               (Discipline::AnalysisDisciplineEnum::enumName(discipline_enum_ID)));
      }
    
    bool insert_success = this->analysis_discipline_map.insert
      (std::map<unsigned int, Discipline::AnalysisDisciplineBase*>::value_type
       (discipline_enum_ID, discipline)).second;
    
    Assert(insert_success, ExcInternalError());
    
    return discipline;
    }
  else
    return it->second;
}




void 
FESystem::FESystemController::writeOutput()
{
  this->performance_logging->setEvent("FESystemController::writeOutput()", 
                                      "FESystem Run Time");

  if (FESystem::local_processor == 0)
    this->output_processor->writeData();

  this->performance_logging->unsetEvent("FESystemController::writeOutput()", 
                                        "FESystem Run Time");
}


void 
FESystem::FESystemController::initPropertyCardsForGlobalParameters()
{
  // create a map of all values of global property parameters and then 
  // ask the property database to initialize them
  std::map<unsigned int, double> value_map;
  std::auto_ptr<std::vector<DesignData::PropertyParameter*> > 
  property_params(this->design_database->getPropertyParameters().release());

  std::vector<DesignData::PropertyParameter*>::const_iterator it, end;
  it = property_params->begin();
  end = property_params->end();
  
  bool insert_success = false;
  
  for (; it != end; it++)
    {
    insert_success = value_map.insert(std::map<unsigned int, double>::value_type
                                      ((*it)->getID(), (*it)->getValue())).second;
    Assert(insert_success == true, ExcInternalError());
    }
  
  this->property_database->initAllCardsForGlobalParameters(value_map);
}



