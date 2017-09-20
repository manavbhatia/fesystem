// $Id: fem_main.C,v 1.13 2006-09-05 20:41:37 manav Exp $


// C++ include files
#include <iostream>

// Basic include files needed for mesh functionality
#include "base/libmesh.h"
#include "FESystem/FESystemController.h"
#include "Mesh/MeshList.h"
#include "Mesh/FEMeshData.h"
#include "Utilities/ParallelUtility.h"

//#include <>

int main (int argc, char* argv[])
{  
    
  unsigned int val = 0;
  while (val == 1)
    {
      pid_t pid, ppid;
      
      /* get the process id */
      if ((pid = getpid()) < 0) 
        std::cout << "unable to get pid" << std::endl;
      else 
        std::cout << "Process# : " << pid << std::endl;
      
      /* get the parent process id */
      if ((ppid = getppid()) < 0) 
        std::cout << "unable to get ppid" << std::endl;
      else 
        std::cout << "Parent#:" << ppid << std::endl;
      
      std::cout << "sleeping for 10 seconds... " << std::endl;
      sleep(10);
    }
  
  FESystem::FESystemController system_to_analyze(argc, argv);
  system_to_analyze.readAndInitialize();
  system_to_analyze.solve();
  system_to_analyze.writeOutput();
}

