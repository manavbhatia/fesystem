// $Id: DataInfo.C,v 1.1.2.5 2007-05-14 16:41:36 manav Exp $

// C++ include
#include <string>
#include <sstream>

// FESystem includes
#include "Database/DataInfo.h"
#include "FESystem/FESystemNumbers.h"
#include "FESystem/FESystemExceptions.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "Utilities/AutoptrVector.h"

FESystemDatabase::GenericDataInfoBase::GenericDataInfoBase()
{
  this->clear();
}
    


FESystemDatabase::GenericDataInfoBase::GenericDataInfoBase
(const FESystemDatabase::GenericDataInfoBase& data_info)
{
  this->discipline_enum_ID = data_info.discipline_enum_ID;
  this->if_load_case_dependent = data_info.if_load_case_dependent;
  this->load_case = data_info.load_case;
  this->if_sensitivity = data_info.if_sensitivity;
  this->DV_ID = data_info.DV_ID;
  this->name = data_info.name;
}



FESystemDatabase::GenericDataInfoBase::~GenericDataInfoBase()
{
  
}
    



void
FESystemDatabase::GenericDataInfoBase::clear()
{
  this->discipline_enum_ID = FESystemNumbers::InvalidID;
  this->if_load_case_dependent = false;
  this->load_case = FESystemNumbers::InvalidID;
  this->if_sensitivity = false;
  this->DV_ID = FESystemNumbers::InvalidID;
  this->name.clear();
}




const std::string&
FESystemDatabase::GenericDataInfoBase::getName() const
{
  Assert(this->name != "", ExcInternalError());
  return this->name;
}



bool
FESystemDatabase::GenericDataInfoBase::ifLoadCaseDependent() const
{
  return this->if_load_case_dependent;
}




unsigned int
FESystemDatabase::GenericDataInfoBase::getLoadCase() const
{
  // make sure that this is load case dependent. Else, this should not be called
  Assert(this->if_load_case_dependent, ExcInternalError());
  
  return this->load_case;
}




void
FESystemDatabase::GenericDataInfoBase::setLoadCase(const unsigned int lcase)
{
  // make sure that this object does not have any other load case stored. If it does, 
  // the user should have called clear before calling this method
  Assert(! this->if_load_case_dependent, ExcInternalError());
  Assert(lcase != FESystemNumbers::InvalidID, ExcInternalError());
  
  this->load_case = lcase;
  this->if_load_case_dependent = true;
}




bool
FESystemDatabase::GenericDataInfoBase::ifSensitivityData() const
{
  return this->if_sensitivity;
}




unsigned int
FESystemDatabase::GenericDataInfoBase::getDVID() const
{
  // make sure that this is load case dependent. Else, this should not be called
  Assert(this->if_sensitivity, ExcInternalError());
  
  return this->DV_ID;
}



void
FESystemDatabase::GenericDataInfoBase::setDVID(const unsigned int ID)
{
  // make sure that this object does not have any other DV_ID stored. If it does, 
  // the user should have called clear before calling this method
  Assert(! this->if_sensitivity, ExcInternalError());
  Assert(ID != FESystemNumbers::InvalidID, ExcInternalError());
  
  this->DV_ID = ID;
  this->if_sensitivity = true;
}



void
FESystemDatabase::GenericDataInfoBase::setName(const std::string& qty_name)
{
  Assert(this->name == "", ExcInternalError());
  Assert(qty_name != "", ExcInternalError());
  this->name = qty_name;
}




unsigned int 
FESystemDatabase::GenericDataInfoBase::getDisciplineEnumID() const
{
  Assert(this->discipline_enum_ID != FESystemNumbers::InvalidID, 
	 ExcInternalError());

  return this->discipline_enum_ID;
}



const std::string& 
FESystemDatabase::GenericDataInfoBase::getDisciplineEnumName() const
{
  Assert(this->discipline_enum_ID != FESystemNumbers::InvalidID, 
	 ExcInternalError());
  
  return Discipline::AnalysisDisciplineEnum::enumName(this->discipline_enum_ID);
}



void 
FESystemDatabase::GenericDataInfoBase::setDisciplineEnumID(const unsigned int enum_ID)
{
  // make sure that this object has been cleared
  Assert(this->discipline_enum_ID == FESystemNumbers::InvalidID,
	 ExcInternalError());

  // make sure that the number is valid
  Assert(enum_ID != FESystemNumbers::InvalidID, ExcInternalError());
  
  this->discipline_enum_ID = enum_ID;
}


std::ostream&
FESystemDatabase::GenericDataInfoBase::write(std::ostream& out)
{
  out << "Discipline :   " << this->getDisciplineEnumName() << std::endl;
  if (this->ifLoadCaseDependent())
    out << "Load Case Number :   " << this->getLoadCase() << std::endl;
  if (this->ifSensitivityData())
    out << "Sensitivity for DV ID :   " << this->getDVID() << std::endl;
  out << "Qty Name :   " << this->getName() << std::endl;
  
  return out;
}



bool 
FESystemDatabase::GenericDataInfoBase::isEqual
(const FESystemDatabase::GenericDataInfoBase& data_info) const
{
  // check each data and return the result
  if (this->discipline_enum_ID != data_info.discipline_enum_ID ||
      this->if_load_case_dependent != data_info.if_load_case_dependent ||
      this->load_case != data_info.load_case ||
      this->if_sensitivity != data_info.if_sensitivity || 
      this->DV_ID  != data_info.DV_ID || 
      this->name != data_info.name)
    return false;
  
  // if it gets here, means that this the two objects are same
  return true;
}




FESystemDatabase::DataInfoBase::DataInfoBase():
  FESystemDatabase::GenericDataInfoBase()
{
  
}



FESystemDatabase::DataInfoBase::DataInfoBase
(const DataInfoBase& data_info):
  FESystemDatabase::GenericDataInfoBase(data_info)
{
  
}



FESystemDatabase::DataInfoBase::~DataInfoBase()
{
  
}



bool 
FESystemDatabase::DataInfoBase::operator== 
  (const FESystemDatabase::DataInfoBase& data_info) const
{
  return (FESystemDatabase::GenericDataInfoBase::isEqual(data_info) &&
	  this->getDataInfoKindEnumID() == data_info.getDataInfoKindEnumID());
}



FESystemDatabase::TimeIndependentDataInfo::TimeIndependentDataInfo():
FESystemDatabase::DataInfoBase()
{
  
}



FESystemDatabase::TimeIndependentDataInfo::TimeIndependentDataInfo
(const TimeIndependentDataInfo& data_info):
  FESystemDatabase::DataInfoBase(data_info)
{
  
}



FESystemDatabase::TimeIndependentDataInfo::~TimeIndependentDataInfo()
{
  
}

    


unsigned int
FESystemDatabase::TimeIndependentDataInfo::getDataInfoKindEnumID() const
{
  return FESystemDatabase::TIME_INDEPENDENT_DATA::num();
}




const std::string
FESystemDatabase::TimeIndependentDataInfo::getPath() const
{
  std::string full_path_name = this->getDisciplineEnumName();
  full_path_name +="_";
  full_path_name += this->getName();
  
  if (this->ifSensitivityData())
    {
    std::ostringstream dv_id_oss;
    dv_id_oss << this->getDVID();
    
    full_path_name += "_SENSITIVITY_DV_";
    full_path_name += dv_id_oss.str();
    }
  
  if (this->ifLoadCaseDependent())
    {
    std::ostringstream lc_oss;
    lc_oss << this->getLoadCase();

    full_path_name += "_LC_";
    full_path_name += lc_oss.str();
    }

  return full_path_name;
}



const std::string
FESystemDatabase::TimeIndependentDataInfo::getDisplayName() const
{
  std::string full_path_name = this->getName();
  
  if (this->ifSensitivityData())
    {
    std::ostringstream dv_id_oss;
    dv_id_oss << this->getDVID();
    
    full_path_name += "_SENS_DV_";
    full_path_name += dv_id_oss.str();
    }
  
  if (this->ifLoadCaseDependent())
    {
    std::ostringstream lc_oss;
    lc_oss << this->getLoadCase();

    full_path_name += "_LC_";
    full_path_name += lc_oss.str();
    }

  return full_path_name;
}



bool 
FESystemDatabase::TimeIndependentDataInfo::operator== 
  (const FESystemDatabase::DataInfoBase& data_info) const
{
  // check each data and return the result
  if (! FESystemDatabase::DataInfoBase::operator==(data_info))
    return false;
  
  // if it gets here, means that this the two objects are same
  return true;
}



FESystemDatabase::TimeDependentDataInfo::TimeDependentDataInfo():
FESystemDatabase::DataInfoBase()
{
  this->clear();
}
    



FESystemDatabase::TimeDependentDataInfo::TimeDependentDataInfo
(const TimeDependentDataInfo& data_info):
  FESystemDatabase::DataInfoBase(data_info)
{
  this->transient_iteration_number = data_info.transient_iteration_number;
  this->time_instant = data_info.time_instant;
  this->derivative_order = data_info.derivative_order;
  this->set_order = data_info.set_order;
}




FESystemDatabase::TimeDependentDataInfo::~TimeDependentDataInfo()
{
  
}





void
FESystemDatabase::TimeDependentDataInfo::clear()
{
  this->transient_iteration_number = FESystemNumbers::InvalidID;
  this->time_instant = 0.0;
  this->derivative_order = FESystemNumbers::InvalidID;
  this->set_order = false;
  
  FESystemDatabase::DataInfoBase::clear();
}
    


unsigned int
FESystemDatabase::TimeDependentDataInfo::getTransientIterationNumber() const
{
  Assert(this->transient_iteration_number != FESystemNumbers::InvalidID,
         ExcInternalError());
  return this->transient_iteration_number;
}



double
FESystemDatabase::TimeDependentDataInfo::getTimeValue() const
{
  Assert(this->transient_iteration_number != FESystemNumbers::InvalidID,
         ExcInternalError());
  return this->time_instant; 
}



void
FESystemDatabase::TimeDependentDataInfo::setOrder(unsigned int i)
{
  Assert(!this->set_order, ExcInternalError());
  this->derivative_order = i; 
  this->set_order = true;
}


unsigned int
FESystemDatabase::TimeDependentDataInfo::getOrder() const
{
  Assert(this->set_order, ExcInternalError());
  return this->derivative_order;
}



void
FESystemDatabase::TimeDependentDataInfo::setTransientIterationInfo
(const unsigned int iter_num, const double time_val)
{
  // make sure that this object was cleared before setting new values
  Assert(this->transient_iteration_number == FESystemNumbers::InvalidID,
         ExcInternalError());
  // the given number should be valid
  Assert(iter_num > 0, ExcInternalError());
  Assert(time_val >= 0.0, ExcInternalError());
  
  this->transient_iteration_number = iter_num;
  this->time_instant = time_val;
}



unsigned int
FESystemDatabase::TimeDependentDataInfo::getDataInfoKindEnumID() const
{
  return FESystemDatabase::TIME_DEPENDENT_DATA::num();
}



const std::string
FESystemDatabase::TimeDependentDataInfo::getPath() const
{
  std::string full_path_name = this->getDisciplineEnumName();
  full_path_name +="_";
  full_path_name += this->getName();

  // append the derivative order
  std::ostringstream order_iter_oss;
  order_iter_oss << this->getOrder();
  
  full_path_name += "_ORDER_";
  full_path_name += order_iter_oss.str();

  // append the transient information to the name
  std::ostringstream time_iter_oss;
  time_iter_oss << this->getTransientIterationNumber();
  
  full_path_name += "_TIME_ITER_";
  full_path_name += time_iter_oss.str();
  
  // append the sensitivity information to the name
  if (this->ifSensitivityData())
    {
    std::ostringstream dv_id_oss;
    dv_id_oss << this->getDVID();
    
    full_path_name += "_SENSITIVITY_DV_";
    full_path_name += dv_id_oss.str();
    }
  
  // append the load case information to the name
  if (this->ifLoadCaseDependent())
    {
    std::ostringstream lc_oss;
    lc_oss << this->getLoadCase();
    
    full_path_name += "_LC_";
    full_path_name += lc_oss.str();
    }
  
  return full_path_name;
}



const std::string
FESystemDatabase::TimeDependentDataInfo::getDisplayName() const
{
  std::string full_path_name = this->getName();

  // append the derivative order
  std::ostringstream order_iter_oss;
  order_iter_oss << this->getOrder();
  
  full_path_name += "_OR_";
  full_path_name += order_iter_oss.str();

  // append the transient information to the name
  std::ostringstream time_iter_oss;
  time_iter_oss << this->getTransientIterationNumber();
  
  full_path_name += "_TIME_";
  full_path_name += time_iter_oss.str();
  
  // append the sensitivity information to the name
  if (this->ifSensitivityData())
    {
    std::ostringstream dv_id_oss;
    dv_id_oss << this->getDVID();
    
    full_path_name += "_SENS_DV_";
    full_path_name += dv_id_oss.str();
    }
  
  // append the load case information to the name
  if (this->ifLoadCaseDependent())
    {
    std::ostringstream lc_oss;
    lc_oss << this->getLoadCase();
    
    full_path_name += "_LC_";
    full_path_name += lc_oss.str();
    }
  
  return full_path_name;
}



std::ostream&
FESystemDatabase::TimeDependentDataInfo::write(std::ostream& out)
{
  FESystemDatabase::GenericDataInfoBase::write(out);

  out << "Time Derivative Order :   " << this->getOrder() << std::endl;
  out << "Transient Iteration Number :   " << this->getTransientIterationNumber() << std::endl;
  out << "Time Value :   " << this->getTimeValue() << std::endl;
  
  return out;
}





bool 
FESystemDatabase::TimeDependentDataInfo::operator== 
  (const FESystemDatabase::DataInfoBase& data_info) const
{
  if (! FESystemDatabase::DataInfoBase::operator==(data_info))
    return false;
  else
    {
      const FESystemDatabase::TimeDependentDataInfo &time_data_info = 
	dynamic_cast<const FESystemDatabase::TimeDependentDataInfo&>(data_info);
      
      // check each data and return the result
      if (this->transient_iteration_number != time_data_info.transient_iteration_number ||
          this->time_instant != time_data_info.time_instant ||
          this->derivative_order != time_data_info.derivative_order)
	return false;
    }
  
  // if it gets here, means that this the two objects are same
  return true;
}



FESystemDatabase::TimeIndependentEigenDataInfo::TimeIndependentEigenDataInfo():
FESystemDatabase::TimeIndependentDataInfo()
{
  this->clear();
}




FESystemDatabase::TimeIndependentEigenDataInfo::TimeIndependentEigenDataInfo
(const TimeIndependentEigenDataInfo& data_info):
  FESystemDatabase::TimeIndependentDataInfo(data_info)
{
  this->mode_number = data_info.mode_number;
}



FESystemDatabase::TimeIndependentEigenDataInfo::~TimeIndependentEigenDataInfo()
{
  
}




void
FESystemDatabase::TimeIndependentEigenDataInfo::clear()
{
  this->mode_number = FESystemNumbers::InvalidID;
  
  FESystemDatabase::TimeIndependentDataInfo::clear();
}




unsigned int
FESystemDatabase::TimeIndependentEigenDataInfo::getModeNumber() const
{
  // make sure that this was set by the user
  Assert(this->mode_number != FESystemNumbers::InvalidID, ExcInternalError());
  return this->mode_number;
}


    


void
FESystemDatabase::TimeIndependentEigenDataInfo::setModeInfo(const unsigned int mode_num)
{
  // make sure that this object was cleared
  Assert(this->mode_number == FESystemNumbers::InvalidID, ExcInternalError());
  // make sure that the number provided is correct
  Assert(mode_num != FESystemNumbers::InvalidID, ExcInternalError());
  
  this->mode_number = mode_num;
}






unsigned int
FESystemDatabase::TimeIndependentEigenDataInfo::getDataInfoKindEnumID() const
{
  return FESystemDatabase::TIME_INDEPENDENT_EIGEN_DATA::num();
}

    



const std::string 
FESystemDatabase::TimeIndependentEigenDataInfo::getPath() const
{
  std::string full_path_name = this->getDisciplineEnumName();
  full_path_name +="_";
  full_path_name += this->getName();
  
  // append the mode information to the name
  std::ostringstream mode_num_oss;
  mode_num_oss << this->getModeNumber();
  
  full_path_name += "_MODE_NUM_";
  full_path_name += mode_num_oss.str();
  
  // append the sensitivity information to the name
  if (this->ifSensitivityData())
    {
    std::ostringstream dv_id_oss;
    dv_id_oss << this->getDVID();
    
    full_path_name += "_SENSITIVITY_DV_";
    full_path_name += dv_id_oss.str();
    }
  
  // append the load case information to the name
  if (this->ifLoadCaseDependent())
    {
    std::ostringstream lc_oss;
    lc_oss << this->getLoadCase();
    
    full_path_name += "_LC_";
    full_path_name += lc_oss.str();
    }
  
  return full_path_name;  
}
    




const std::string 
FESystemDatabase::TimeIndependentEigenDataInfo::getDisplayName() const
{
  std::string full_path_name = this->getName();
  
  // append the mode information to the name
  std::ostringstream mode_num_oss;
  mode_num_oss << this->getModeNumber();
  
  full_path_name += "_MODE_";
  full_path_name += mode_num_oss.str();
  
  // append the sensitivity information to the name
  if (this->ifSensitivityData())
    {
    std::ostringstream dv_id_oss;
    dv_id_oss << this->getDVID();
    
    full_path_name += "_SENS_DV_";
    full_path_name += dv_id_oss.str();
    }
  
  // append the load case information to the name
  if (this->ifLoadCaseDependent())
    {
    std::ostringstream lc_oss;
    lc_oss << this->getLoadCase();
    
    full_path_name += "_LC_";
    full_path_name += lc_oss.str();
    }
  
  return full_path_name;  
}



std::ostream&
FESystemDatabase::TimeIndependentEigenDataInfo::write(std::ostream& out)
{
  FESystemDatabase::GenericDataInfoBase::write(out);
  
  out << "Mode Number :   " << this->getModeNumber() << std::endl;

  return out;
}



bool 
FESystemDatabase::TimeIndependentEigenDataInfo::operator== 
  (const FESystemDatabase::DataInfoBase& data_info) const
{
  if (! FESystemDatabase::TimeIndependentDataInfo::operator==(data_info))
    return false;
  else
    {
      const FESystemDatabase::TimeIndependentEigenDataInfo &eigen_data_info = 
	dynamic_cast<const FESystemDatabase::TimeIndependentEigenDataInfo&>(data_info);
      
      // check each data and return the result
      if (this->mode_number != eigen_data_info.mode_number)
	return false;
    }
  
  // if it gets here, means that this the two objects are same
  return true;
}





FESystemDatabase::TimeDependentEigenDataInfo::TimeDependentEigenDataInfo():
FESystemDatabase::TimeDependentDataInfo()
{
  this->clear();
}




FESystemDatabase::TimeDependentEigenDataInfo::TimeDependentEigenDataInfo
(const TimeDependentEigenDataInfo& data_info):
  FESystemDatabase::TimeDependentDataInfo(data_info)
{
  this->mode_number = data_info.mode_number;
}




FESystemDatabase::TimeDependentEigenDataInfo::~TimeDependentEigenDataInfo()
{
  
}




void
FESystemDatabase::TimeDependentEigenDataInfo::clear()
{
  this->mode_number = FESystemNumbers::InvalidID;
  
  FESystemDatabase::TimeDependentDataInfo::clear();  
}




unsigned int
FESystemDatabase::TimeDependentEigenDataInfo::getModeNumber() const
{
  // make sure that this was set by the user
  Assert(this->mode_number != FESystemNumbers::InvalidID, ExcInternalError());
  return this->mode_number;  
}





void
FESystemDatabase::TimeDependentEigenDataInfo::setModeInfo(const unsigned int mode_num)
{
  // make sure that this object was cleared
  Assert(this->mode_number == FESystemNumbers::InvalidID, ExcInternalError());
  // make sure that the number provided is correct
  Assert(mode_num != FESystemNumbers::InvalidID, ExcInternalError());
  
  this->mode_number = mode_num;  
}






unsigned int
FESystemDatabase::TimeDependentEigenDataInfo::getDataInfoKindEnumID() const
{
  return FESystemDatabase::TIME_INDEPENDENT_EIGEN_DATA::num();
}





const std::string 
FESystemDatabase::TimeDependentEigenDataInfo::getPath() const
{
  std::string full_path_name = this->getDisciplineEnumName();
  full_path_name +="_";
  full_path_name += this->getName();
  
  // append the mode information to the name
  std::ostringstream mode_num_oss;
  mode_num_oss << this->getModeNumber();
  
  full_path_name += "_MODE_NUM_";
  full_path_name += mode_num_oss.str();

  // append the transient information to the name
  std::ostringstream time_iter_oss;
  time_iter_oss << this->getTransientIterationNumber();
  
  full_path_name += "_TIME_ITER_";
  full_path_name += time_iter_oss.str();
  
  // append the sensitivity information to the name
  if (this->ifSensitivityData())
    {
    std::ostringstream dv_id_oss;
    dv_id_oss << this->getDVID();
    
    full_path_name += "_SENSITIVITY_DV_";
    full_path_name += dv_id_oss.str();
    }
  
  // append the load case information to the name
  if (this->ifLoadCaseDependent())
    {
    std::ostringstream lc_oss;
    lc_oss << this->getLoadCase();
    
    full_path_name += "_LC_";
    full_path_name += lc_oss.str();
    }
  
  return full_path_name;  
}




const std::string 
FESystemDatabase::TimeDependentEigenDataInfo::getDisplayName() const
{
  std::string full_path_name = this->getName();
  
  // append the mode information to the name
  std::ostringstream mode_num_oss;
  mode_num_oss << this->getModeNumber();
  
  full_path_name += "_MODE_";
  full_path_name += mode_num_oss.str();

  // append the transient information to the name
  std::ostringstream time_iter_oss;
  time_iter_oss << this->getTransientIterationNumber();
  
  full_path_name += "_TIME_";
  full_path_name += time_iter_oss.str();
  
  // append the sensitivity information to the name
  if (this->ifSensitivityData())
    {
    std::ostringstream dv_id_oss;
    dv_id_oss << this->getDVID();
    
    full_path_name += "_SENS_DV_";
    full_path_name += dv_id_oss.str();
    }
  
  // append the load case information to the name
  if (this->ifLoadCaseDependent())
    {
    std::ostringstream lc_oss;
    lc_oss << this->getLoadCase();
    
    full_path_name += "_LC_";
    full_path_name += lc_oss.str();
    }
  
  return full_path_name;  
}


std::ostream&
FESystemDatabase::TimeDependentEigenDataInfo::write(std::ostream& out)
{
  FESystemDatabase::TimeDependentDataInfo::write(out);
  
  out << "Mode Number :   " << this->getModeNumber() << std::endl;
  
  return out;
}



bool 
FESystemDatabase::TimeDependentEigenDataInfo::operator== 
		  (const FESystemDatabase::DataInfoBase& data_info) const
{
  if (! FESystemDatabase::TimeDependentDataInfo::operator==(data_info))
    return false;
  else
    {
      const FESystemDatabase::TimeDependentEigenDataInfo &eigen_data_info = 
	dynamic_cast<const FESystemDatabase::TimeDependentEigenDataInfo&>(data_info);
      
      // check each data and return the result
      if (this->mode_number != eigen_data_info.mode_number)
	return false;
    }
  
  // if it gets here, means that this the two objects are same
  return true;
}




FESystemDatabase::TransientDataInfoSet::TransientDataInfoSet():
  FESystemDatabase::GenericDataInfoBase(),
derivative_order(0)
{
  this->clear();
}



FESystemDatabase::TransientDataInfoSet::~TransientDataInfoSet()
{
  this->clear();
}



void
FESystemDatabase::TransientDataInfoSet::clear()
{
  this->derivative_order = 0;
  this->set_order = false;
  this->time_values.clear();
  FESystemDatabase::GenericDataInfoBase::clear();
}



void
FESystemDatabase::TransientDataInfoSet::setOrder(unsigned int i)
{
  Assert(!this->set_order, ExcInternalError());
  this->derivative_order = i; 
  this->set_order = true;
}


unsigned int
FESystemDatabase::TransientDataInfoSet::getOrder() const
{
  Assert(this->set_order, ExcInternalError());
  return this->derivative_order;
}


unsigned int 
FESystemDatabase::TransientDataInfoSet::getNDataInfo() const
{
  Assert(this->time_values.size() > 0, ExcInternalError());
  return this->time_values.size();
}




FESystemUtility::AutoPtrVector<FESystemDatabase::TimeDependentDataInfo> 
FESystemDatabase::TransientDataInfoSet::getDataInfoForAllTimeSteps() const
{
  unsigned int n = this->getNDataInfo();

  FESystemUtility::AutoPtrVector<FESystemDatabase::TimeDependentDataInfo> 
    data_info_vec(n);

  FESystemDatabase::TimeDependentDataInfo *info=NULL;
  
  for (unsigned int i=0; i < n; i++)
    {
      info = new FESystemDatabase::TimeDependentDataInfo();
      info->setDisciplineEnumID(this->getDisciplineEnumID());
      info->setName(this->getName());
      if (this->ifLoadCaseDependent())
        info->setLoadCase(this->getLoadCase());
      if (this->ifSensitivityData())
        info->setDVID(this->getDVID());
      info->setTransientIterationInfo(i+1, this->time_values[i]);
      info->setOrder(this->getOrder());
      
      data_info_vec.reset(i, info);
      info = NULL;
    }

  return data_info_vec;
}




void 
FESystemDatabase::TransientDataInfoSet::setTransientIterationTimeVector
(const std::vector<double>& time_vals)
{
  Assert(this->time_values.size() == 0, ExcInternalError());
  Assert(time_vals.size() > 0, ExcInternalError());
  this->time_values = time_vals;
}



const std::string 
FESystemDatabase::TransientDataInfoSet::getDisplayName() const
{
  std::string full_path_name = this->getName();
  
  // append the sensitivity information to the name
  if (this->ifSensitivityData())
    {
    std::ostringstream dv_id_oss;
    dv_id_oss << this->getDVID();
    
    full_path_name += "_SENS_DV_";
    full_path_name += dv_id_oss.str();
    }
  
  // append the load case information to the name
  if (this->ifLoadCaseDependent())
    {
    std::ostringstream lc_oss;
    lc_oss << this->getLoadCase();
    
    full_path_name += "_LC_";
    full_path_name += lc_oss.str();
    }
  
  return full_path_name;  
}



std::auto_ptr<FESystemDatabase::DataInfoBase> 
FESystemDatabase::createDataInfoCopy(const FESystemDatabase::DataInfoBase& data_info)
{
  std::auto_ptr<FESystemDatabase::DataInfoBase> copy;
  
  switch(data_info.getDataInfoKindEnumID())
    {
    case  TIME_INDEPENDENT_DATA_ENUM_ID:
      copy.reset(new FESystemDatabase::TimeIndependentDataInfo
		 (dynamic_cast<const FESystemDatabase::TimeIndependentDataInfo&>(data_info)));
      break;

    case TIME_DEPENDENT_DATA_ENUM_ID:
      copy.reset(new FESystemDatabase::TimeDependentDataInfo
		 (dynamic_cast<const FESystemDatabase::TimeDependentDataInfo&>(data_info)));
      break;

    case  TIME_INDEPENDENT_EIGEN_DATA_ENUM_ID:
      copy.reset(new FESystemDatabase::TimeIndependentEigenDataInfo
		 (dynamic_cast<const FESystemDatabase::TimeIndependentEigenDataInfo&>(data_info)));
      break;

    case  TIME_DEPENDENT_EIGEN_DATA_ENUM_ID:
      copy.reset(new FESystemDatabase::TimeDependentEigenDataInfo
		 (dynamic_cast<const FESystemDatabase::TimeDependentEigenDataInfo&>(data_info)));      
      break;

    default:
      Assert(false, ExcInternalError());
    }

  return copy;
}
