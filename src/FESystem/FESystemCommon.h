// $Id: FESystemCommon.h,v 1.4 2006-11-13 00:08:46 manav Exp $

#ifndef __fesystem_common_h__
#define __fesystem_common_h__

#ifdef WITH_MATHEMATICA
#include "mathlink.h"
#endif

namespace FESystem
{
  /// world MPI communicator for all processors
  MPI_Comm COMM_WORLD;
  
  /// MPI info
  MPI_Info MPI_INFO;
  
  /// total number of processors in this communicator
  unsigned int total_processors = 0;
  
  /// local processor ID
  unsigned int local_processor = 0;
  
#ifdef WITH_MATHEMATICA
  /// mathlink environment
  MLENV mathlink_environment = (MLENV)0;
  
  /// mathlink link to the kernerl
  MLINK mathlink_link = (MLINK)0;
#endif
}  
  

#endif //__fesystem_common_h__
