// $Id: DataInfo.h,v 1.1.2.4 2007-05-14 16:41:36 manav Exp $

#ifndef __fesystem_data_info_h__
#define __fesystem_data_info_h__

// C++ include
#include <string>
#include <memory>
#include <vector>
#include <ostream>

// FESystem includes
#include "Utilities/NameEnumHandler.h"


namespace FESystemUtility
{
  template<typename T> class AutoPtrVector;
}


#ifndef TIME_INDEPENDENT_DATA_ENUM_ID
#define TIME_INDEPENDENT_DATA_ENUM_ID 1
#else
#error
#endif

#ifndef TIME_INDEPENDENT_DATA_ENUM_NAME
#define TIME_INDEPENDENT_DATA_ENUM_NAME "TIME_INDEPENDENT_DATA"
#else
#error
#endif



#ifndef TIME_DEPENDENT_DATA_ENUM_ID
#define TIME_DEPENDENT_DATA_ENUM_ID 2
#else
#error
#endif

#ifndef TIME_DEPENDENT_DATA_ENUM_NAME
#define TIME_DEPENDENT_DATA_ENUM_NAME "TIME_DEPENDENT_DATA"
#else
#error
#endif



#ifndef TIME_INDEPENDENT_EIGEN_DATA_ENUM_ID
#define TIME_INDEPENDENT_EIGEN_DATA_ENUM_ID 3
#else
#error
#endif

#ifndef TIME_INDEPENDENT_EIGEN_DATA_ENUM_NAME
#define TIME_INDEPENDENT_EIGEN_DATA_ENUM_NAME "TIME_INDEPENDENT_EIGEN_DATA"
#else
#error
#endif


#ifndef TIME_DEPENDENT_EIGEN_DATA_ENUM_ID
#define TIME_DEPENDENT_EIGEN_DATA_ENUM_ID 4
#else
#error
#endif

#ifndef TIME_DEPENDENT_EIGEN_DATA_ENUM_NAME
#define TIME_DEPENDENT_EIGEN_DATA_ENUM_NAME "TIME_DEPENDENT_EIGEN_DATA"
#else
#error
#endif



namespace FESystemDatabase
{
  
  DeclareEnumClass(DataInfoKindEnum);
  
  DeclareEnumName(TIME_INDEPENDENT_DATA, FESystemDatabase::DataInfoKindEnum,
                  TIME_INDEPENDENT_DATA_ENUM_ID, TIME_INDEPENDENT_DATA_ENUM_NAME);
  
  DeclareEnumName(TIME_DEPENDENT_DATA, FESystemDatabase::DataInfoKindEnum,
                  TIME_DEPENDENT_DATA_ENUM_ID, TIME_DEPENDENT_DATA_ENUM_NAME);
  
  DeclareEnumName(TIME_INDEPENDENT_EIGEN_DATA, FESystemDatabase::DataInfoKindEnum,
                  TIME_INDEPENDENT_EIGEN_DATA_ENUM_ID, TIME_INDEPENDENT_EIGEN_DATA_ENUM_NAME);
  
  DeclareEnumName(TIME_DEPENDENT_EIGEN_DATA, FESystemDatabase::DataInfoKindEnum,
                  TIME_DEPENDENT_EIGEN_DATA_ENUM_ID, TIME_DEPENDENT_EIGEN_DATA_ENUM_NAME);

  
  /// this id the base class for all data info objects. The primary purpose of the data info 
  /// objects is to give information about a numeric object (matrix/vector) about its kind, 
  /// the time iteration and time value, nonlinear iteration, etc. and all such details that 
  /// are necessary to property store the object. 
  class GenericDataInfoBase
  {
  public:
    /// constructor
    GenericDataInfoBase();
    
    /// copy constructor
    GenericDataInfoBase(const GenericDataInfoBase& data_info);
    
    /// Destructor
    virtual ~GenericDataInfoBase();
    
    /// clears the data structures of this object
    virtual void clear();
    
    /// @returns the name of the quantity with which it is referred
    const std::string& getName() const;

    // @returns a string that will be assciated with this qty if it is written to output files. 
    // This is essentially a shorter form of the path, since some softwares have a restriction
    // on the number of characters to be used for displaying the qty name. The total length of 
    // this will be 256 characters
    virtual const std::string getDisplayName() const=0;

    /// @returns true if the quantity is load case dependent, else false
    bool ifLoadCaseDependent() const;

    /// @returns the load case of the data
    unsigned int getLoadCase() const;
    
    /// sets the load case value, and also sets to true the boolean value of if_load_case_dependent
    void setLoadCase(const unsigned int load_case);
    
    /// @returns true if the quantity is a sensitivity quantity, else false
    bool ifSensitivityData() const;
    
    /// @returns the DV ID if this is a sensitivity quantity
    unsigned int getDVID() const;

    /// sets the DV_ID, and also sets to true the boolean value of if_sensitivity
    void setDVID(const unsigned int ID);
    
    /// sets the name of the quantity
    void setName(const std::string& qty_name);

    /// @returns the name of the discipline to which this quantity belongs
    const std::string& getDisciplineEnumName() const;

    /// @returns the enum ID of the discipline to which this quantity belongs
    unsigned int getDisciplineEnumID() const;

    /// sets the discipline enum ID
    void setDisciplineEnumID(const unsigned int enum_ID);

    /// @returns true if the given data info object and this, point to the same data
    virtual bool isEqual (const FESystemDatabase::GenericDataInfoBase& data_info) const;

    /// writes the data info to an output stream
    virtual std::ostream& write(std::ostream& out);

  protected:
      
    /// enum ID of the discipline
    unsigned int discipline_enum_ID;

    /// boolean to tell if this quantity is dependent on a load case or not
    bool if_load_case_dependent;
      
    /// load case value of the object
    unsigned int load_case;
    
    /// boolean to tell if this is a sensitivity quantity
    bool if_sensitivity;
    
    /// DV ID for this quantitiy if it is a sensitivity quantity
    unsigned int DV_ID;
      
    /// string name of the quantity
    std::string name;
  };



  /// this id the base class for all data info objects. The primary purpose of the data info 
  /// objects is to give information about a numeric object (matrix/vector) about its kind, 
  /// the time iteration and time value, nonlinear iteration, etc. and all such details that 
  /// are necessary to property store the object. 
  class DataInfoBase: public FESystemDatabase::GenericDataInfoBase
  {
  public:
    /// constructor
    DataInfoBase();

    /// copy constructor
    DataInfoBase(const DataInfoBase& data_info);
    
    /// Destructor
    virtual ~DataInfoBase();
    
    /// @returns the enumeration ID of the kind of the data 
    virtual unsigned int getDataInfoKindEnumID() const = 0;
    
    /// @returns the exact path to the location of the data for this object in the database
    virtual const std::string getPath() const = 0;
    
    /// @returns true if the given data info object and this, point to the same data
    virtual bool operator== (const FESystemDatabase::DataInfoBase& data_info) const;

  protected:
    
  };
  
  

  /// this defines a time independent quantity
  class TimeIndependentDataInfo: public FESystemDatabase::DataInfoBase
  {
  public:
    /// constructor
    TimeIndependentDataInfo();
    
    /// copy constructor
    TimeIndependentDataInfo(const TimeIndependentDataInfo& data_info);

    /// destructor
    ~TimeIndependentDataInfo();
    
    /// @returns the enumeration ID of the kind of the data 
    virtual unsigned int getDataInfoKindEnumID() const;
    
    /// @returns the exact path to the location of the data for this object in the database
    virtual const std::string getPath() const;

    // @returns a string that will be assciated with this qty if it is written to output files. 
    // This is essentially a shorter form of the path, since some softwares have a restriction
    // on the number of characters to be used for displaying the qty name. The total length of 
    // this will be 64 characters
    virtual const std::string getDisplayName() const;

    /// @returns true if the given data info object and this, point to the same data
    virtual bool operator== (const FESystemDatabase::DataInfoBase& data_info) const;

  protected:
    
  };



  /// this defines a time dependent quantity
  class TimeDependentDataInfo: public FESystemDatabase::DataInfoBase
  {
  public:
    /// constructor
    TimeDependentDataInfo();
    
    /// copy constructor
    TimeDependentDataInfo(const TimeDependentDataInfo& data_info);

    /// destructor
    ~TimeDependentDataInfo();
    
    /// clears the data structures of this object
    virtual void clear();
    
    /// @returns the time iteration number
    unsigned int getTransientIterationNumber() const;
    
    /// @returns the time value at which this quantity is calculated
    double getTimeValue() const;
    
    /// sets the transient iteration info
    void setTransientIterationInfo(const unsigned int iter_num, 
                                   const double time_val);
    
    /// @returns the enumeration ID of the kind of the data 
    virtual unsigned int getDataInfoKindEnumID() const;

    /// sets the time derivative of the solution vector
    void setOrder(const unsigned int i);
    
    /// @returns the time derivative of the solution vector
    unsigned int getOrder() const;
    
    /// @returns the exact path to the location of the data for this object in the database
    virtual const std::string getPath() const;
    
    // @returns a string that will be assciated with this qty if it is written to output files. 
    // This is essentially a shorter form of the path, since some softwares have a restriction
    // on the number of characters to be used for displaying the qty name. The total length of 
    // this will be 256 characters
    virtual const std::string getDisplayName() const;

    /// this method writes the info to an output stream
    virtual std::ostream& write(std::ostream& out);

    /// @returns true if the given data info object and this, point to the same data
    virtual bool operator== (const FESystemDatabase::DataInfoBase& data_info) const;

  protected:
      
    /// transient iteration
    unsigned int transient_iteration_number;
    
    /// time value of this quantity
    double time_instant;
    
    /// time derivative order of the quantity
    unsigned int derivative_order;
    
    /// boolean to keep track of whether the order was set
    bool set_order;
  };
  
  
  
  /// this defines a time independent eigen quantity
  class TimeIndependentEigenDataInfo: 
  public FESystemDatabase::TimeIndependentDataInfo
  {
  public:
    /// constructor
    TimeIndependentEigenDataInfo();
    
    /// copy constructor
    TimeIndependentEigenDataInfo(const TimeIndependentEigenDataInfo& data_info);

    /// destructor
    ~TimeIndependentEigenDataInfo();
    
    /// clears the data structures of this object
    virtual void clear();
    
    /// @returns the time iteration number
    unsigned int getModeNumber() const;
    
    /// sets the modal information for this quantity
    void setModeInfo(const unsigned int mode_num);

    /// @returns the enumeration ID of the kind of the data 
    virtual unsigned int getDataInfoKindEnumID() const;
    
    /// @returns the exact path to the location of the data for this object in the database
    virtual const std::string getPath() const;

    /// this method writes the info to an output stream
    virtual std::ostream& write(std::ostream& out);

    // @returns a string that will be assciated with this qty if it is written to output files. 
    // This is essentially a shorter form of the path, since some softwares have a restriction
    // on the number of characters to be used for displaying the qty name. The total length of 
    // this will be 256 characters
    virtual const std::string getDisplayName() const;
    
    /// @returns true if the given data info object and this, point to the same data
    virtual bool operator== (const FESystemDatabase::DataInfoBase& data_info) const;

  protected:
    
    /// transient iteration
    unsigned int mode_number;
  };




  /// this defines a time independent eigen quantity
  class TimeDependentEigenDataInfo: 
  public FESystemDatabase::TimeDependentDataInfo
  {
  public:
    /// constructor
    TimeDependentEigenDataInfo();
      
    /// copy constructor
    TimeDependentEigenDataInfo(const TimeDependentEigenDataInfo& data_info);
      
    /// destructor
    ~TimeDependentEigenDataInfo();
      
    /// clears the data structures of this object
    virtual void clear();
      
    /// @returns the time iteration number
    unsigned int getModeNumber() const;
      
    /// sets the modal information for this quantity
    void setModeInfo(const unsigned int mode_num);

    /// @returns the enumeration ID of the kind of the data 
    virtual unsigned int getDataInfoKindEnumID() const;
      
    /// this method writes the info to an output stream
    virtual std::ostream& write(std::ostream& out);

    /// @returns the exact path to the location of the data for this object in the database
    virtual const std::string getPath() const;
      
    // @returns a string that will be assciated with this qty if it is written to output files. 
    // This is essentially a shorter form of the path, since some softwares have a restriction
    // on the number of characters to be used for displaying the qty name. The total length of 
    // this will be 256 characters
    virtual const std::string getDisplayName() const;

    /// @returns true if the given data info object and this, point to the same data
    virtual bool operator== (const FESystemDatabase::DataInfoBase& data_info) const;

  protected:
        
    /// transient iteration
    unsigned int mode_number;
  };


  /// this is a class to define a set of solutions in a transient. It will create and return the
  /// DataInfoBase objects for the solutions in the set
  class TransientDataInfoSet: public FESystemDatabase::GenericDataInfoBase
    {
    public:
      TransientDataInfoSet();

      ~TransientDataInfoSet();

      void clear();
      
      virtual const std::string getDisplayName() const;

      unsigned int getNDataInfo() const;

      void setOrder(const unsigned int order);

      unsigned int getOrder() const;

      FESystemUtility::AutoPtrVector<FESystemDatabase::TimeDependentDataInfo> 
	getDataInfoForAllTimeSteps() const;
      
      void setTransientIterationTimeVector(const std::vector<double>& time_vals);

    protected:
      
      /// derivative order of the solution
      unsigned int derivative_order;
      
      /// boolean to track if the derivative order was set
      bool set_order;
      
      /// vector of time values
      std::vector<double> time_values;
    };

  
  /// this creates a copy of data info and returns it
  std::auto_ptr<FESystemDatabase::DataInfoBase> createDataInfoCopy
    (const FESystemDatabase::DataInfoBase& data_info);
}


#endif // __fesystem_data_info_h__

