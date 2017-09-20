// $Id: Log.C,v 1.2 2006-09-05 20:41:44 manav Exp $

// C++ includes


// FESystem includes
#include "Log.h"


Log::Log()
{
	// open the log file 
	this->log_file.open("log_file", std::fstream::out);
}


Log::~Log()
{
	// close the log file
	this->log_file.close();
}
