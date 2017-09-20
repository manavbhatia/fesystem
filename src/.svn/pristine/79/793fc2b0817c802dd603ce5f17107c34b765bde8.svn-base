// $Id: Log.h,v 1.3.6.1 2007-05-08 05:19:13 manav Exp $

#ifndef __fesystem_log_h__
#define __fesystem_log_h__


// C++ includes
#include <fstream>
#include <string>


// FESystem includes



// libMesh includes


/// This is a utility class to log the progress of an analysis in a file for debug purposes and also to record solutions 
/// as they are calculated.

class Log
{
 public:
  /// constuctor
  /// param name of input file
  Log( );
	
	
  /// destructor
  ~Log();
	
	
  /// write to the log file
  template<class T>
    void write(T& data)
    {
      log_file << data;
      log_file << std::endl << std::endl;
    }


  /// write to the log file
  template<class T>
    void write2(T& data)
    {
      data.write(log_file);
      log_file << std::endl << std::endl;
    }
	
	
 protected:
		
  /// log file to which all the entries are written
  std::fstream log_file;
};


#endif // __fesystem_log_h__
