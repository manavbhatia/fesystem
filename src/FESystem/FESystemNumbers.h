// $Id: FESystemNumbers.h,v 1.4 2006-12-07 01:40:48 manav Exp $

#ifndef __fesystem_numbers_h__
#define __fesystem_numbers_h__


namespace FESystemNumbers
{
  /// an ID that will always be considered to be invalid in the code.
  static const unsigned int InvalidID = 0;
  
  /// the smallest size limit for division in the code
  static const double Epsilon = 1.0e-13;
  
  static const double Pi = 3.1415926535897932385;
}


#endif // __fesystem_numbers_h__
