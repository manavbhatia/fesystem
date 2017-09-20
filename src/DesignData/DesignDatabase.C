// $Id: DesignDatabase.C,v 1.3 2006-09-05 20:41:41 manav Exp $

// FESystem includes
#include "DesignData/DesignDatabase.h"
#include "Utilities/InputOutputUtility.h"

DesignData::DesignDatabase::DesignDatabase()
{
  
}




DesignData::DesignDatabase::~DesignDatabase()
{
  // iterate over all the parameters and delete them, since they were created using
  // the new operator
  DesignData::DesignDatabase::IDToParameterMap::iterator param_it, param_end;
  param_it = this->ID_to_design_parameter_map.begin();
  param_end = this->ID_to_design_parameter_map.end();
  
  for (; param_it != param_end; param_it++)
    {
    delete param_it->second;
    param_it->second = NULL;
    }
}



std::istream& 
DesignData::DesignDatabase::readFromInputStream(std::istream& input)
{
  std::string tag;
  unsigned int num = 0, param_kind =0;
  DesignParameter *param = NULL;

  FESystemIOUtility::readFromInput(input, "DESIGN_DATABASE");
  FESystemIOUtility::readFromInput(input, "BEGIN");

  // read in the list of design variables
  FESystemIOUtility::readFromInput(input, "DESIGN_PARAMETERS");
  FESystemIOUtility::readFromInput(input, "BEGIN");
  FESystemIOUtility::readFromInput(input, "N_PARAMETERS", num);

  for (unsigned int param_it=0; param_it < num; param_it++)
    {
    tag.clear();
    FESystemIOUtility::peekFromInput(input, tag);
    param_kind = DesignParameterTypeEnum::enumID(tag);
    switch(param_kind)
      {
      case PROPERTY_PARAMETER_ENUM_ID:
        {
          PropertyParameter *prop_param = new DesignData::PropertyParameter;
          input >> (*prop_param);
          param = prop_param;
        }
        break;

      case SHAPE_PARAMETER_ENUM_ID:
        {
          ShapeParameter *shape_param = new DesignData::ShapeParameter;
          input >> (*shape_param);
          param = shape_param;
        }
        break;
      }
    
    this->addParameterToDatabase(param);
    }
  
  FESystemIOUtility::readFromInput(input, "DESIGN_PARAMETERS");
  FESystemIOUtility::readFromInput(input, "END");
  
  FESystemIOUtility::readFromInput(input, "DESIGN_DATABASE");
  FESystemIOUtility::readFromInput(input, "END");

  return input;
}


std::istream& 
DesignData::operator>>(std::istream& input, DesignData::DesignDatabase& database)
{
  database.readFromInputStream(input);
  
  return input;
}

