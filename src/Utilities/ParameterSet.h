// $Id: ParameterSet.h,v 1.1 2006-09-05 21:11:42 manav Exp $

#ifndef __parameter_set_h__
#define __parameter_set_h__


// C++ includes
#include <string>
#include <map>



/// this is an abstract base class for parameter sets, so that pointers to the 
/// derived class can be cast to this class for storage in a container
class ParameterSetBase
{
public:
  /// defalut constructor
  ParameterSetBase();
  
  
  /// destructor 
  virtual ~ParameterSetBase();
  
  
  /// an abstract method for clearing the object
  virtual void clear() = 0;

  
  /// @returns the ID of the parameter set
  unsigned int getID() const;
  
private:
  
  /// integer id of the class
  unsigned int ID;
};



ParameterSetBase::ParameterSetBase():
ID(FESystemNumber::InvalidID)
{}



ParameterSetBase::~ParameterSetBase()
{}



unsigned int ParameterSetBase::getID()
{
  assert (this->ID != FESystemNumber::InvalidID);
  return this->ID;
}


/// this class provides the necessary functionality for storing a 
/// set of values by a key. These values could be integer, real or string 
/// in nature.
template <class KeyType>
class ParameterSet: public ParameterSetBase
{
 public:
  
  /// this enumeration provides options for storing the parameter values
  enum StorageOption
  {
    ADD_VALUES,           /// if the parameter is already present, add the two values
                          /// and store it. In case of strings, the new value is 
                          /// appended to the old one
    REPLACE_VALUE,        /// replaces value if already present
    STORE_IF_ABSENT       /// if the parameter is already present, do nothing, otherwise
                          /// add it
  };
  
  /// default constructor
  ParameterSet();

  /// copy constructor
  ParameterSet(ParameterSet& param);

  
  /// destructor
  ~ParameterSet();
  
  /// clears the parameters stored in the map
  inline void clear();
  
  /// add integer valued parameter
  /// @param name string name of the parameter to be added
  /// @param value integer value of the parameter
  /// @param option storage options for the parameter
  void addIntegerParameter(const KeyType& key,
                           unsigned int value,
                           ParameterSet::StorageOption option =
                           ParameterSet::REPLACE_VALUE);
  
  /// add real valued parameter
  /// @param name string name of the parameter to be added
  /// @param value float value of the parameter
  /// @param option storage options for the parameter
  void addFloatParameter(const KeyType& name,
                         double value,
                         ParameterSet::StorageOption option =
                         ParameterSet::REPLACE_VALUE);
  
  /// add string valued parameter
  /// @param name string name of the parameter to be added
  /// @param value string value of the parameter
  /// @param option storage options for the parameter
  void addStringParameter(const KeyType& name,
                          const std::string& value,
                          ParameterSet::StorageOption option =
                          ParameterSet::REPLACE_VALUE);
  

  
  /// add string valued parameter
  /// @param name string name of the parameter to be added
  /// @param value string value of the parameter
  /// @param option storage options for the parameter
  void addBoolParameter(const KeyType& name,
                          const bool value,
                          ParameterSet::StorageOption option =
                          ParameterSet::REPLACE_VALUE);
  
  /// add string valued parameter
  /// @param name string name of the parameter to be added
  /// @param value string value of the parameter
  /// @param option storage options for the parameter
  void addTensorParameter(const KeyType& name,
                          const TensorBase& value,
                          ParameterSet::StorageOption option =
                          ParameterSet::REPLACE_VALUE);
  
  
  /// @returns the parameter value 
  /// @param name string name of the parameter to be returned
  inline unsigned int getIntegerParameter(const KeyType& name) const;

  
  /// @returns the parameter value 
  /// @param name string name of the parameter to be returned
  inline double getFloatParameter(const KeyType& name) const;

  
  /// @returns the parameter value 
  /// @param name string name of the parameter to be returned
  inline std::string getStringParameter(const KeyType& name) const;
  
  /// @returns the parameter value 
  /// @param name string name of the parameter to be returned
  inline bool getBoolParameter(const KeyType& name) const;

  /// @returns the parameter value 
  /// @param name string name of the parameter to be returned
  inline TensorBase& getTensorParameter(const KeyType& name) const;
  
protected:

  typedef std::map<KeyType, unsigned int> IntParamMap;  
  typedef std::map<KeyType, double> FloatParamMap;  
  typedef std::map<KeyType, std::string> StringParamMap;  
  typedef std::map<KeyType, bool> BoolParamMap;  
  typedef std::map<KeyType, TensorBase*> TensorParamMap;  
  
  /// this method provides a generic interface for obtaining the 
  /// parameter values from a map
  /// @param name a string name of the parameter whose value is to be retuned
  /// @param map_ref a reference to the map from which value has to be obtained
  /// @returns value of the parameter of type specified in the first template arguement
  template <class ValueType, class MapType>
    ValueType getParameter(const KeyType& name,
                           const MapType& map_ref) const;
  
  /// this method provides a generic interface for storing the 
  /// parameter value to a map
  /// @param name string name of the parameter whose value is to be added
  /// @param value value of the parameter to be added
  /// @param map_ref reference to the map from which value has to be obtained
  /// @option option about what to do if the parameter already exists
  template <class ValueType, class MapType>
    void addParameter(const KeyType& name,
                      const ValueType value,
                      MapType& map_ref,
                      ParameterSet::StorageOption option);
  
  
  /// integer parameters 
  IntParamMap integer_params;
  
  /// float paramters
  FloatParamMap float_params;
  
  /// string parameters
  StringParamMap string_params;

  /// boolean parameters
  BoolParamMap boolean_params;
  
  /// tensor parameters
  TensorParamMap tensor_params;
};


template <class KeyType>
inline 
void ParameterSet<KeyType>::clear()
{
  this->integer_params.clear();
  this->float_params.clear();
  this->string_params.clear();
  this->boolean_params.clear();
  // the tensors need to be deleted, since they were stored as a local quantity
  // created by the new operator
  {
    TensorParamMap::iterator it, end;
    it = this->tensor_map.begin();
    end = this->tensor_map.end();
    for (; it != end; it++)
      delete it->second;
    this->tensor_params.clear();
  }
}




template <class KeyType>
template<class ValueType, class MapType>
ValueType ParameterSet<KeyType>::getParameter(const KeyType& name,
                                              const MapType& map_ref) const
{
  // get an iterator to the parameter name in the map
  typename MapType::const_iterator it, end;
  it = map_ref.find(name);
  end = map_ref.end();
  
  // make sure that the parameter exists in the map
  assert (it != end);
  
  // next return the parameter value
  return it->second;
}



template<class KeyType>
template <class ValueType, class MapType>
void ParameterSet<KeyType>::addParameter(const KeyType& name,
                                         const ValueType value,
                                         MapType& map_ref,
                                         ParameterSet::StorageOption option)
{
  // get an iterator to the parameter name in the map
  typename MapType::iterator it, end;
  it = map_ref.find(name);
  end = map_ref.end();
  
  // if the parameter does not exist, add the value
  if (it == end)
    {
    // insert the pair
    std::pair<typename MapType::iterator, bool> insert_return_pair = 
    map_ref.insert(typename MapType::value_type(name, value));
    
    // make sure that the insert was successful
    assert (insert_return_pair.second == true);
    }
  else
    {
    switch (option)
      {
      case ParameterSet::ADD_VALUES:
        {
          it->second += value;
        }
        break;
        
      case ParameterSet::REPLACE_VALUE:
        {
          it->second = value;
        }
        break;
        
      case ParameterSet::STORE_IF_ABSENT:
        // nothing to be done here, since the value is already present
        break;
        
      default:
        abort();
      }
    }
}



template<class KeyType>
inline 
unsigned int ParameterSet<KeyType>::getIntegerParameter(const KeyType& name) const
{
  return this->getParameter<unsigned int, IntParamMap>(name, 
                                                       this->integer_params);
}



template<class KeyType>
inline 
double ParameterSet<KeyType>::getFloatParameter(const KeyType& name) const
{
  return this->getParameter<double, FloatParamMap>(name, 
                                                   this->float_params);
}



template<class KeyType>
inline 
std::string ParameterSet<KeyType>::getStringParameter(const KeyType& name) const
{
  return this->getParameter<std::string, StringParamMap>(name, 
                                                         this->string_params);
}

template<class KeyType>
inline 
const std::string ParameterSet<KeyType>::getBoolParameter(const KeyType& name) const
{
  return this->getParameter<bool, BoolParamMap>(name, 
                                                this->boolean_params);
}


template<class KeyType>
inline 
const TensorBase& ParameterSet<KeyType>::getTensorParameter(const KeyType& name) const
{
  TensorBase* tbase = this->getParameter<TensorBase*, TensorParamMap>(name, 
                                                                      this->tensor_params);
  return *tbase;
}


template<class KeyType>
ParameterSet<KeyType>::ParameterSet():
ParameterSetBase()
{
  
}



template<class KeyType>
ParameterSet<KeyType>::ParameterSet(ParameterSet& param_set):
ParameterSetBase()
{
  // copy the maps
  // integer
  {
    typename IntParamMap::const_iterator it, end;
    it = param_set.integer_params.begin();
    end = param_set.integer_params.end();
    for (; it != end; it++)
      {
      std::pair<typename IntParamMap::iterator,bool> insert_return_pair = 
      this->integer_params.insert(IntParamMap::value_type(it->first, it->second));
      assert (insert_return_pair.second == true);
      }
  }
  
  // float
  {
    typename FloatParamMap::const_iterator it, end;
    it = param_set.float_params.begin();
    end = param_set.float_params.end();
    for (; it != end; it++)
      {
      std::pair<typename FloatParamMap::iterator,bool> insert_return_pair = 
      this->float_params.insert(FloatParamMap::value_type(it->first, it->second));
      assert (insert_return_pair.second == true);
      }
  }
  
  // string parameters
  {
    typename StringParamMap::const_iterator it, end;
    it = param_set.string_params.begin();
    end = param_set.string_params.end();
    for (; it != end; it++)
      {
      std::pair<typename StringParamMap::iterator,bool> insert_return_pair = 
      this->string_params.insert(StringParamMap::value_type(it->first, it->second));
      assert (insert_return_pair.second == true);
      }
  }
  
}



template<class KeyType>
ParameterSet<KeyType>::~ParameterSet()
{
  
}



template<class KeyType>
void ParameterSet<KeyType>::addIntegerParameter(const KeyType& name,
                                           unsigned int value,
                                           ParameterSet::StorageOption option)
{
  this->addParameter<unsigned int, IntParamMap>(name, 
                                                value,
                                                this->integer_params,
                                                option);
}


template<class KeyType>
void ParameterSet<KeyType>::addFloatParameter(const KeyType& name,
                                         double value,
                                         ParameterSet::StorageOption option)
{
  this->addParameter<double, FloatParamMap>(name, 
                                            value,
                                            this->float_params,
                                            option);
}



template<class KeyType>
void ParameterSet<KeyType>::addStringParameter(const KeyType& name,
                                          const std::string& value,
                                          ParameterSet::StorageOption option)
{
  this->addParameter<std::string, StringParamMap>(name, 
                                                  value,
                                                  this->string_params,
                                                  option);
}


template<class KeyType>
void ParameterSet<KeyType>::addBoolParameter(const KeyType& name,
                                        const bool value,
                                        ParameterSet::StorageOption option)
{
  this->addParameter<bool, BoolParamMap>(name, 
                                         value,
                                         this->boolean_params,
                                         option);
}

template<class KeyType>
void ParameterSet<KeyType>::addTensorParameter(const KeyType& name,
                                        const TensorBase& tensor,
                                        ParameterSet::StorageOption option)
{
  TensorBase *local_tensor = TensorBase::return_copy(tensor).release();
  
  this->addParameter<TensorBase*, TensorParamMap>(name, 
                                                  value,
                                                  this->boolean_params,
                                                  option);
}

#endif // __parameter_set_h__
