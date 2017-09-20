// $Id: FESystemExceptions.h,v 1.2 2006-09-05 20:17:36 manav Exp $

#ifndef __fesystem_exceptions_h__
#define __fesystem_exceptions_h__

// C++ incldes
#include <string>

// deal.II includes
#include "base/exceptions.h"


// the following exceptions have been defined on top of deal.II exceptions to define some local 
// checks in the code
namespace FESystemExceptions
{
  /// this is an exception for a NULL pointer, when it should have had a valid value
  DeclException0(ExcNullPointer);

  /// this is an exception when the insert operation on maps, etc turns out to be
  /// unsuccessful
  DeclException0(ExcInsertUnsuccessful);
  
  /// this is an exception when an invalid tag is found in the input. The first input param
  /// is the tag found, and the second param is the tag that was expected
  DeclException2(ExcIOBadTag, std::string, std::string, 
                 << "Error in input! Expected: " << arg2 << "   Found: " << arg1 );
  

  /// this exception indicates that an invalid ID was found. Invalid ID is defined in the 
  /// \p FESystemNmbers namespace.
  DeclException1(ExcInvalidID, unsigned int, 
                 << "Invalid ID: " << arg1 );
  
  
  /// this is an exception when an enumeration case was found, and was not 
  /// handled in the switch block.
  DeclException1(ExcEnumCaseNotImplemented, std::string,
                 << "Enumeration Case:  " << arg1 << "  not handled."  );
  
  /// this is an exception when a duplicate ID is specified in the input
  DeclException2(ExcDuplicateID, std::string, unsigned int,
                 << "Duplicate " <<  arg1 << " ID: " << arg2 << "  specified.");

  /// this is an exception when a requested ID does not exist
  DeclException2(ExcIDDoesNotExist, std::string, unsigned int,
                 << "Requested " <<  arg1 << " ID: " << arg2 << "  does not exist.");
  
}



#endif // __fesystem_exceptions_h__
