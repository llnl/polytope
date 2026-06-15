#ifndef POLYTOPE_ErrorHandler_HH
#define POLYTOPE_ErrorHandler_HH

#include <string>

namespace polytope
{

//! Issues an error with the given status. By default, an error issues a 
//! message to cout and exits the program with the status.
void error(const std::string& message, int status = -1);

//! Sets the error handler for the polytope library.
void setErrorHandler(void (*handler)(const std::string&, int));

}

#endif
