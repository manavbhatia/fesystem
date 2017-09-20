// $Id: PostProcessQtyDatabase.C,v 1.4 2006-09-05 20:41:50 manav Exp $

// C++ includes
#include <cassert>

// FESystem includes
#include "PostProcess/PostProcessQtyDatabase.h"

// libMesh includes


PostProcessQtyDatabase::PostProcessQtyDatabase()
{
  
}



PostProcessQtyDatabase::~PostProcessQtyDatabase()
{
  // iterate over all elements in the map and delete the pointers to the quantities
  AnalysisPostProcessMap::iterator an_it, an_end;
  an_it = this->analysis_post_process_map.begin();
  an_end = this->analysis_post_process_map.end();
  
  // iterate over the post process map for each analysis, and clear those
  for (; an_it != an_end; an_it++)
    {
    PostProcessQtyMap::iterator it, end;
    it = an_it->second->begin();
    end = an_it->second->end();
    
    for (; it != end; it++)
      {
      delete it->second;
      }
    
    delete an_it->second;
    an_it->second = NULL;
    }
}





PostProcessQty& 
PostProcessQtyDatabase::getElementPostProcessQty(const unsigned int discipline,
                                                 unsigned int elem_ID)
{
  // find the discipline in the iterator
  AnalysisPostProcessMap::iterator it = 
  this->analysis_post_process_map.find(discipline);
  
  assert (it != this->analysis_post_process_map.end());
  
  // find the element in the iterator
  PostProcessQtyMap::iterator qty_it = it->second->find(elem_ID);
  
  // make sure that this does not already exist
  assert (qty_it != it->second->end());
  
  return (*qty_it->second);
}


void PostProcessQtyDatabase::addElementPostProcessQty(const unsigned int discipline,
                                                      unsigned int elem_ID, 
                                                      std::auto_ptr<ElemPostProcessQty> qty )
{
  // find the discipline in the iterator
  AnalysisPostProcessMap::iterator it = 
  this->analysis_post_process_map.find(discipline);
  
  if (it == this->analysis_post_process_map.end())
    {
    std::pair<AnalysisPostProcessMap::iterator, bool> insert_pair =
    this->analysis_post_process_map.insert 
    (AnalysisPostProcessMap::value_type(discipline, new PostProcessQtyMap));
    
    assert (insert_pair.second == true);
    it = insert_pair.first;
    }
  
  // find the element in the iterator
  PostProcessQtyMap::iterator qty_it = it->second->find(elem_ID);
  
  // make sure that this does not already exist
  assert (qty_it == it->second->end());
  
  bool insert = 
    it->second->insert(std::map<unsigned int, PostProcessQty*>::value_type
                       (elem_ID, qty.release())).second;
  
  assert (insert == true);
}

