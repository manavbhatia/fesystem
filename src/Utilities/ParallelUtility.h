//$Id:$
/*
 *  Parallel.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 11/30/08.
 *  Copyright 2008 Virginia Tech. All rights reserved.
 *
 */
#ifndef __fesystem_utility_parallel_h__
#define __fesystem_utility_parallel_h__

// C++ includes
#include <sstream>
#include <string>
#include <iostream> 
#include <vector>

// mpi includes
#include "mpi.h"

// FESystem includes
#include "FESystem/FESystemExceptions.h"


namespace FESystemUtility
{
  namespace Parallel 
  {
  
  /// this method will print the \p data to the screen with the number
  /// of the process from which it is coming. The << operator should be
  /// supported by data
  template <typename T> void PrintAll(MPI_Comm& comm, const T& data)
  {
    int n_proc=0, proc_id=-1, err=-1, num=0; 
    
    err = MPI_Comm_size(comm, &n_proc);
    AssertThrow(err == MPI_SUCCESS, ExcInternalError());
    err = MPI_Comm_rank(comm, &proc_id);
    AssertThrow(err == MPI_SUCCESS, ExcInternalError());
    
    if (proc_id != 0)
      {
        // send a message
        std::ostringstream oss;
        oss << data;
        std::string str = oss.str();
        num = str.size();
        
        err = MPI_Send(&num, 1, MPI_INT, 0, 1, comm);
        AssertThrow(err == MPI_SUCCESS, ExcInternalError());
        
        err = MPI_Send(&(str[0]), num, MPI_CHAR, 0, 1, comm);
        AssertThrow(err == MPI_SUCCESS, ExcInternalError());
      }
    else
      {
        // receive the message from each processor
        std::cout << "Message from Processor #: " << proc_id << " :: "
        << data << std::endl;
        
        for (int i=1; i<n_proc; i++)
          {
            // receive the size of the text
            err = MPI_Recv(&num, 1, MPI_INT, i, 1, comm, MPI_STATUS_IGNORE);
            AssertThrow(err == MPI_SUCCESS, ExcInternalError());
            
            std::string str;
            str.resize(num);
            
            // receive the text
            err = MPI_Recv(&(str[0]), num, MPI_CHAR, i, 1, comm, MPI_STATUS_IGNORE);
            AssertThrow(err == MPI_SUCCESS, ExcInternalError());
            
            std::cout << "Message from Processor #: " << i << " :: "
            << str << std::endl;
          }
        
      }
  }
  }
}


#endif // __fesystem_utility_parallel_h__



