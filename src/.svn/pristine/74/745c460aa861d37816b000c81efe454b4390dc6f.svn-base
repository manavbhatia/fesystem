/*
 *  PistonTheoryInfo.h
 *  fesystem_xcode
 *
 *  Created by Manav Bhatia on 9/21/08.
 *  Copyright 2008 Manav Bhatia. All rights reserved.
 *
 */

// $Id: PistonTheoryInfo.h,v 1.3 2006/09/05 20:09:54 manav Exp $

#ifndef __fesystem_piston_theory_info_h__
#define __fesystem_piston_theory_info_h__

// FESystem includes
#include "Discipline/DisciplineInfo.h"
#include "Discipline/PistonTheory.h"


#ifndef PISTON_THEORY_INFO_ENUM_ID
#define PISTON_THEORY_INFO_ENUM_ID 3
#else
#error
#endif

#ifndef PISTON_THEORY_INFO_ENUM_NAME
#define PISTON_THEORY_INFO_ENUM_NAME "PISTON_THEORY_INFO"
#else
#error
#endif



namespace Discipline
{
  DeclareEnumName(PISTON_THEORY_INFO, Discipline::DisciplineInfoTypeEnum,
                  PISTON_THEORY_INFO_ENUM_ID, PISTON_THEORY_INFO_ENUM_NAME);
  
  
  class PistonTheoryInfo: public DisciplineInfo
    {
    public:
      
      /// constructor
      PistonTheoryInfo();
      
      /// destructor
      virtual ~PistonTheoryInfo();
      
      /// @returns the discipline enum ID
      virtual inline unsigned int getDisciplineEnumID() const;
      
      /// @returns the discipline enum name
      virtual inline const std::string getDisciplineEnumName() const;
      
      /// reads from the input stream
      virtual std::istream& readFromInputStream(std::istream& input);
            
    protected:
      
      /// the order of this theory
      unsigned int theory_order;
      
      /// the set of element IDs on which the loads will be calculated
      unsigned int element_set_ID;
      
    };
}




inline 
unsigned int
Discipline::PistonTheoryInfo::getDisciplineEnumID() const
{
  return Discipline::PISTON_THEORY::num();
}





inline 
const std::string
Discipline::PistonTheoryInfo::getDisciplineEnumName() const
{
  return Discipline::PISTON_THEORY::name();
}



#endif //  __fesystem_piston_theory_info_h__
