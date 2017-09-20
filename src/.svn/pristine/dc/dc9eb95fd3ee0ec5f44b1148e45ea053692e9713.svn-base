// $Id: LoadedSolution.C,v 1.2.4.4 2008-02-25 04:33:47 manav Exp $

// C++ includes


// FESystem includes
#include "Utilities/LoadedSolution.h"
#include "Database/GlobalDataStorage.h"
#include "Database/DataInfo.h"
#include "FESystem/FESystemExceptions.h"



void 
FESystemUtility::LoadedVectorData::setDataInfo(const FESystemDatabase::DataInfoBase& data_info_obj)
{
  Assert(!this->initialized(), ExcInternalError());
  this->data_info.reset(FESystemDatabase::createDataInfoCopy(data_info_obj).release());
}


FESystemUtility::LoadedSolution::LoadedSolution():
maximum_vectors(2),
previous_vec(0)
{
  for (unsigned int i=0; i < this->maximum_vectors; i++)
    this->loaded_vectors.push_back(new FESystemUtility::LoadedVectorData);
}



FESystemUtility::LoadedSolution::~LoadedSolution()
{
  std::vector<FESystemUtility::LoadedVectorData*>::iterator it, end;
  it = this->loaded_vectors.begin();
  end = this->loaded_vectors.end();
  
  for (; it != end; it++)
    delete (*it);
}



void 
FESystemUtility::LoadedSolution::clear()
{
  std::vector<FESystemUtility::LoadedVectorData*>::iterator it, end;
  it = this->loaded_vectors.begin();
  end = this->loaded_vectors.end();
  
  for (; it != end; it++)
    (*it)->clear();
}




const NumericVector<double>&
FESystemUtility::LoadedSolution::loadSolutionFromDatabase
(FESystemDatabase::GlobalDataStorage& database,
 const FESystemDatabase::DataInfoBase& data_info)
{
  // if the vector already exists, return it, else load it
  std::vector<FESystemUtility::LoadedVectorData*>::iterator it, end;
  it = this->loaded_vectors.begin();
  end = this->loaded_vectors.end();
  
  for (; it != end; it++)
    {
      if ((*it)->initialized() && (*it)->getDataInfo() == data_info)
        return (*it)->getVector();
    }
  
  // if the code gets here, that means that the vector was not loaded. In this case, 
  // load the specified vector.
  unsigned int next_vector = 0;
  if (previous_vec == (this->loaded_vectors.size() - 1))
    next_vector = 0;
  else
    next_vector = previous_vec + 1;
  
  // if the vector has been used, then clear it before setting the next data
  if (this->loaded_vectors[next_vector]->initialized())
    this->loaded_vectors[next_vector]->clear();
  
  this->loaded_vectors[next_vector]->setDataInfo(data_info);

  Assert (database.checkIfQtyExists(data_info), ExcInternalError());
  database.fillVector(data_info, this->loaded_vectors[next_vector]->getVector());
  
  this->previous_vec = next_vector;
  
  return this->loaded_vectors[next_vector]->getVector();
}
