// $Id: ElementDataStorage.C,v 1.12.6.2 2008-02-25 04:20:25 manav Exp $

// C++ includes
#include <cassert>


// FESystem includes
#include "Database/ElementDataStorage.h"
#include "FESystem/FESystemExceptions.h"

// libMesh Includes


ElementDataStorage::ElementDataStorage()
{
	
}



ElementDataStorage::~ElementDataStorage()
{
  this->clear();
}


void ElementDataStorage::clear()
{
  // first iterate over all the elements, then over the tags, and then over load case IDs	
  {
    ElementDataStorage::DisciplineToElemIDMatrixMap::iterator
    discipline_tag_it, discipline_tag_end;
    discipline_tag_it = this->matrix_data.begin();
    discipline_tag_end = this->matrix_data.end();
    
    for (; discipline_tag_it != discipline_tag_end; discipline_tag_it++)
      {    
      ElementDataStorage::ElemIDToTagMatrixMap::iterator  elem_tag_it, elem_tag_end;
      elem_tag_it = discipline_tag_it->second.begin();
      elem_tag_end = discipline_tag_it->second.end();
      
      for (; elem_tag_it != elem_tag_end; elem_tag_it++)
        {
        ElementDataStorage::TagToMatrixMap::iterator tag_caseID_it, tag_caseID_end;
        tag_caseID_it= elem_tag_it->second.begin();
        tag_caseID_end= elem_tag_it->second.end();
        
        for (; tag_caseID_it != tag_caseID_end; tag_caseID_it++)
          {
          ElementDataStorage::LoadCaseIDToMatrixMap::iterator  caseID_matrix_it, caseID_matrix_end;
          caseID_matrix_it= tag_caseID_it->second.begin();
          caseID_matrix_end= tag_caseID_it->second.end();
          
          for (; caseID_matrix_it != caseID_matrix_end; caseID_matrix_it++)
            {
            delete caseID_matrix_it->second;
            }
          
          tag_caseID_it->second.clear();
          }
        
        elem_tag_it->second.clear();
        }
      
      discipline_tag_it->second.clear();
      }  
    
    this->matrix_data.clear();
  }
  
  // repeat the above procedure for vector data
  {
    ElementDataStorage::DisciplineToElemIDVectorMap::iterator
    discipline_tag_it, discipline_tag_end;
    discipline_tag_it = this->vector_data.begin();
    discipline_tag_end = this->vector_data.end();
    
    for (; discipline_tag_it != discipline_tag_end; discipline_tag_it++)
      {    
      ElementDataStorage::ElemIDToTagVectorMap::iterator  elem_tag_it, elem_tag_end;
      elem_tag_it = discipline_tag_it->second.begin();
      elem_tag_end = discipline_tag_it->second.end();
      
      for (; elem_tag_it != elem_tag_end; elem_tag_it++)
        {
        ElementDataStorage::TagToVectorMap::iterator tag_caseID_it, tag_caseID_end;
        tag_caseID_it= elem_tag_it->second.begin();
        tag_caseID_end= elem_tag_it->second.end();
        
        for (; tag_caseID_it != tag_caseID_end; tag_caseID_it++)
          {
          ElementDataStorage::LoadCaseIDToVectorMap::iterator  caseID_vector_it, caseID_vector_end;
          caseID_vector_it= tag_caseID_it->second.begin();
          caseID_vector_end= tag_caseID_it->second.end();
          
          for (; caseID_vector_it != caseID_vector_end; caseID_vector_it++)
            {
            delete caseID_vector_it->second;
            }
          
          tag_caseID_it->second.clear();
          }
        
        elem_tag_it->second.clear();
        }
      
      discipline_tag_it->second.clear();
      }  
    
    this->vector_data.clear();
  }	
}



void ElementDataStorage::add_element_data(const DenseMatrix<double>& matrix,
                                          const unsigned int discipline_enum_ID,
                                          const unsigned int elem, 
                                          const std::string& quantity,
                                          const unsigned int load_case_ID)
{	
//  // first, check if the element has been set to null. If it has been, abort
//  if (this->check_if_null(discipline_enum_ID, elem,quantity,load_case_ID) == true)
//    abort();
	
  
  // check if the discipline ID exists in the map. If not, create it.
  ElementDataStorage::DisciplineToElemIDMatrixMap::iterator discipline_tag_it = 
  this->matrix_data.find(discipline_enum_ID);
	
  if (discipline_tag_it == this->matrix_data.end())
    {
    std::pair<ElementDataStorage::DisciplineToElemIDMatrixMap::iterator, bool> 
    insert_return_value = 
    this->matrix_data.insert(ElementDataStorage::DisciplineToElemIDMatrixMap::value_type
                             (discipline_enum_ID, ElementDataStorage::ElemIDToTagMatrixMap::map()));
		
    assert (insert_return_value.second == true);
		
    discipline_tag_it = insert_return_value.first;
    }
  
  
  // check if the elem ID exists, if not, creat it in the map
  ElementDataStorage::ElemIDToTagMatrixMap::iterator elem_tag_it = 
    discipline_tag_it->second.find(elem);
	
  if (elem_tag_it == discipline_tag_it->second.end())
    {
    std::pair<ElementDataStorage::ElemIDToTagMatrixMap::iterator, bool> insert_return_value = 
    discipline_tag_it->second.insert(ElementDataStorage::ElemIDToTagMatrixMap::value_type
                                     (elem, TagToMatrixMap::map()));
		
    assert (insert_return_value.second == true);
		
    elem_tag_it = insert_return_value.first;
    }
	
	
  // check if the tag exists, if not, create
  TagToMatrixMap::iterator tag_caseID_it = elem_tag_it->second.find(quantity);
	
  if (tag_caseID_it == elem_tag_it->second.end())
    {
    std::pair<TagToMatrixMap::iterator, bool> insert_return_value = 
    elem_tag_it->second.insert(TagToMatrixMap::value_type (quantity, LoadCaseIDToMatrixMap::map()));
		
    assert (insert_return_value.second == true);
		
    tag_caseID_it = insert_return_value.first;
    }
	
  // check if the loadcase ID exists, 
  // if it does, delete the matrix entry, and create a new one, if not, add it
  LoadCaseIDToMatrixMap::iterator caseID_matrix_it = tag_caseID_it->second.find(load_case_ID);
  if (caseID_matrix_it != tag_caseID_it->second.end())
    {
    delete caseID_matrix_it->second;
    tag_caseID_it->second.erase(caseID_matrix_it);
    }
	
  DenseMatrix<double>* local_matrix = new DenseMatrix<double>(matrix);
		
  std::pair<LoadCaseIDToMatrixMap::iterator,bool> insert_return_value = 
    tag_caseID_it->second.insert(LoadCaseIDToMatrixMap::value_type(load_case_ID, local_matrix));
	
  assert (insert_return_value.second == true);
}




void ElementDataStorage::add_element_data(const DenseVector<double>& vector,
                                          const unsigned int discipline_enum_ID,
                                          const unsigned int elem, 
                                          const std::string& quantity,
                                          const unsigned int load_case_ID)
{
//  // first, check if the element has been set to null. If it has been, abort
//  if (this->check_if_null(discipline_enum_ID, elem,quantity,load_case_ID) == true)
//    abort();
	
  
  // check if the discipline ID exists in the map. If not, create it.
  ElementDataStorage::DisciplineToElemIDVectorMap::iterator discipline_tag_it = 
  this->vector_data.find(discipline_enum_ID);
	
  if (discipline_tag_it == this->vector_data.end())
    {
    std::pair<ElementDataStorage::DisciplineToElemIDVectorMap::iterator, bool> 
    insert_return_value = 
    this->vector_data.insert(ElementDataStorage::DisciplineToElemIDVectorMap::value_type
                             (discipline_enum_ID, ElementDataStorage::ElemIDToTagVectorMap::map()));
		
    assert (insert_return_value.second == true);
		
    discipline_tag_it = insert_return_value.first;
    }
  
  
  // check if the elem ID exists, if not, creat it in the map
  ElementDataStorage::ElemIDToTagVectorMap::iterator elem_tag_it = 
    discipline_tag_it->second.find(elem);
	
  if (elem_tag_it == discipline_tag_it->second.end())
    {
    std::pair<ElementDataStorage::ElemIDToTagVectorMap::iterator, bool> insert_return_value = 
    discipline_tag_it->second.insert(ElementDataStorage::ElemIDToTagVectorMap::value_type
                                     (elem, TagToVectorMap::map()));
		
    assert (insert_return_value.second == true);
		
    elem_tag_it = insert_return_value.first;
    }
	
	
  // check if the tag exists, if not, create
  TagToVectorMap::iterator tag_caseID_it = elem_tag_it->second.find(quantity);
	
  if (tag_caseID_it == elem_tag_it->second.end())
    {
    std::pair<TagToVectorMap::iterator, bool> insert_return_value = 
    elem_tag_it->second.insert(TagToVectorMap::value_type (quantity, LoadCaseIDToVectorMap::map()));
		
    assert (insert_return_value.second == true);
		
    tag_caseID_it = insert_return_value.first;
    }
	
  // check if the loadcase ID exists, 
  // if it does, delete the vector entry, and create a new one, if not, add it
  LoadCaseIDToVectorMap::iterator caseID_vector_it = tag_caseID_it->second.find(load_case_ID);
  if (caseID_vector_it != tag_caseID_it->second.end())
    {
    delete caseID_vector_it->second;
    tag_caseID_it->second.erase(caseID_vector_it);
    }
	
  DenseVector<double>* local_vector = new DenseVector<double>(vector);
		
  std::pair<LoadCaseIDToVectorMap::iterator,bool> insert_return_value = 
    tag_caseID_it->second.insert(LoadCaseIDToVectorMap::value_type(load_case_ID, local_vector));
	
  assert (insert_return_value.second == true);
}




bool ElementDataStorage::checkIfQtyExists(const unsigned int discipline_enum_ID,
                                          const unsigned int elem,
                                          const std::string& quantity,
                                          const unsigned int load_case_ID)
{
    bool in_matrix_qtys = this->checkIfMatrixQtyExists(discipline_enum_ID, elem, quantity, load_case_ID);
    bool in_vector_qtys = this->checkIfVectorQtyExists(discipline_enum_ID, elem, quantity, load_case_ID);

    // it is assumed that the qty exists either in vector or in matrix map.
    return (in_matrix_qtys + in_vector_qtys);
    
}




bool ElementDataStorage::checkIfMatrixQtyExists(const unsigned int discipline_enum_ID,
                                                const unsigned int elem,
                                                const std::string& quantity,
                                                const unsigned int load_case_ID)
{
//  // first, check if the element has been set to null. If it has been, abort
//  if (this->check_if_null(discipline_enum_ID, elem,quantity,load_case_ID))
//    return true;
		
  // check if the element exists in the map
  ElementDataStorage::DisciplineToElemIDMatrixMap::iterator discipline_tag_it = 
  this->matrix_data.find(discipline_enum_ID);
		
  // if the element doesn't exist in the map, return false, else check further
  if (discipline_tag_it == this->matrix_data.end())
    return false; 
  else
    {
    ElementDataStorage::ElemIDToTagMatrixMap::iterator elem_tag_it = 
    discipline_tag_it->second.find(elem);
    
    if (elem_tag_it == discipline_tag_it->second.end())
      return false;
    else 
      {
      // get the iterator to the load to matrix map
      ElementDataStorage::TagToMatrixMap::iterator tag_caseID_it = 
      elem_tag_it->second.find(quantity);
      
      if (tag_caseID_it == elem_tag_it->second.end())
        return false;
      else 
        {
        // get the iterator to the load case to qty map
        ElementDataStorage::LoadCaseIDToMatrixMap::iterator caseID_matrix_it = 
        tag_caseID_it->second.find(load_case_ID);
        
        if (caseID_matrix_it == tag_caseID_it->second.end())
          return false;
        else 
          return true;
        }
      }      
    }
}	



bool ElementDataStorage::checkIfVectorQtyExists(const unsigned int discipline_enum_ID,
                                                const unsigned int elem,
                                                const std::string& quantity,
                                                const unsigned int load_case_ID)
{
//  // first, check if the element has been set to null. If it has been, abort
//  if (this->check_if_null(discipline_enum_ID, elem,quantity,load_case_ID))
//    return true;
		
  // check if the element exists in the map
  ElementDataStorage::DisciplineToElemIDVectorMap::iterator discipline_tag_it = 
  this->vector_data.find(discipline_enum_ID);
		
  // if the element doesn't exist in the map, return false, else check further
  if (discipline_tag_it == this->vector_data.end())
    return false; 
  else
    {
    ElementDataStorage::ElemIDToTagVectorMap::iterator elem_tag_it = 
    discipline_tag_it->second.find(elem);
    
    if (elem_tag_it == discipline_tag_it->second.end())
      return false;
    else 
      {
      // get the iterator to the load to vector map
      ElementDataStorage::TagToVectorMap::iterator tag_caseID_it = 
      elem_tag_it->second.find(quantity);
      
      if (tag_caseID_it == elem_tag_it->second.end())
        return false;
      else 
        {
        // get the iterator to the load case to qty map
        ElementDataStorage::LoadCaseIDToVectorMap::iterator caseID_vector_it = 
        tag_caseID_it->second.find(load_case_ID);
        
        if (caseID_vector_it == tag_caseID_it->second.end())
          return false;
        else 
          return true;
        }
      }      
    }
}



void 
ElementDataStorage::get_element_data(DenseMatrix<double>& qty,
                                     const unsigned int discipline_enum_ID,
                                     const unsigned int elem,
                                     const std::string& quantity,
                                     const unsigned int load_case_ID)
{
//  // first, check if the element has been set to null. If it has been, abort
//  if (this->check_if_null(discipline_enum_ID, elem,quantity,load_case_ID) == true)
//    abort();
	
  // check if the element exists in the map
  ElementDataStorage::DisciplineToElemIDMatrixMap::iterator discipline_tag_it = 
  this->matrix_data.find(discipline_enum_ID);
  
  Assert(discipline_tag_it != this->matrix_data.end(), ExcInternalError());
  
  ElementDataStorage::ElemIDToTagMatrixMap::iterator elem_tag_it = 
  discipline_tag_it->second.find(elem);

  Assert(elem_tag_it != discipline_tag_it->second.end(), ExcInternalError());

  // get the iterator to the load to matrix map
  ElementDataStorage::TagToMatrixMap::iterator tag_caseID_it = 
    elem_tag_it->second.find(quantity);

  Assert(tag_caseID_it != elem_tag_it->second.end(), ExcInternalError());

  // get the iterator to the load case to qty map
  ElementDataStorage::LoadCaseIDToMatrixMap::iterator caseID_matrix_it = 
    tag_caseID_it->second.find(load_case_ID);

  Assert(caseID_matrix_it != tag_caseID_it->second.end(), ExcInternalError());

  qty = *(caseID_matrix_it->second);
}



void 
ElementDataStorage::get_element_data(DenseVector<double>& qty,
                                     const unsigned int discipline_enum_ID,
                                     const unsigned int elem,
                                     const std::string& quantity,
                                     const unsigned int load_case_ID)
{
//  // first, check if the element has been set to null. If it has been, abort
//  if (this->check_if_null(discipline_enum_ID, elem,quantity,load_case_ID) == true)
//    abort();
	
  // check if the element exists in the map
  ElementDataStorage::DisciplineToElemIDVectorMap::iterator discipline_tag_it = 
  this->vector_data.find(discipline_enum_ID);
  
  Assert(discipline_tag_it != this->vector_data.end(), ExcInternalError());

  ElementDataStorage::ElemIDToTagVectorMap::iterator elem_tag_it = 
  discipline_tag_it->second.find(elem);
  
  Assert(elem_tag_it != discipline_tag_it->second.end(), ExcInternalError());

  // get the iterator to the load to vector map
  ElementDataStorage::TagToVectorMap::iterator tag_caseID_it = 
  elem_tag_it->second.find(quantity);
  
  Assert(tag_caseID_it != elem_tag_it->second.end(), ExcInternalError());

  // get the iterator to the load case to qty map
  ElementDataStorage::LoadCaseIDToVectorMap::iterator caseID_vector_it = 
  tag_caseID_it->second.find(load_case_ID);

  Assert(caseID_vector_it != tag_caseID_it->second.end(), ExcInternalError());
  
  qty = *(caseID_vector_it->second);
}






//void ElementDataStorage::set_to_null(const unsigned int discipline_enum_ID,
//                                     const unsigned int elem_ID,
//                                     const std::string& tag,
//                                     const unsigned int load_case_ID)
//{
//  // check if the elem ID exists in the map or not. If it does not, 
//  // create it
//  std::map<unsigned int, std::multimap<std::string, unsigned int> >::iterator elem_tag_it =
//  this->null_qty_map.find(elem_ID);
//	
//  if (elem_tag_it == this->null_qty_map.end())
//    {
//    std::pair<std::map<unsigned int, 
//    std::multimap<std::string, unsigned int> >::iterator, bool> insert_return_value = 
//    this->null_qty_map.insert(std::map<unsigned int, std::multimap<std::string, unsigned int> >::
//                              value_type(elem_ID, std::multimap<std::string, unsigned int>::multimap()));
//		
//    assert (insert_return_value.second == true);
//		
//    elem_tag_it = insert_return_value.first;
//    }
//	
//  // insert the tag, case ID pair into the multimap
//  elem_tag_it->second.insert(std::multimap<std::string, unsigned int>::value_type
//                             (tag, load_case_ID));
//}






//bool ElementDataStorage::check_if_null(const unsigned int discipline_enum_ID,
//                                       const unsigned int elem_ID, 
//                                       const std::string& tag,
//                                       const unsigned int load_case_ID)
//{	
//  // check if the elem ID exists, then the tag, then the caseID, 
//  // if all three exist, return true, else flase
//  
//  std::map<unsigned int, std::multimap<std::string, unsigned int> >::const_iterator elem_tag_it =
//  this->null_qty_map.find(elem_ID);
//	
//  if (elem_tag_it == this->null_qty_map.end())
//    {
//    return false;
//    }
//  else
//    {
//    std::pair<	std::multimap<std::string, unsigned int>::const_iterator,
//    std::multimap<std::string, unsigned int>::const_iterator> range_check = 
//    elem_tag_it->second.equal_range(tag);
//		
//    if (range_check.first == elem_tag_it->second.end())
//      {
//      return false;
//      }
//    else
//      {
//      std::multimap<std::string, unsigned int>::const_iterator case_it = range_check.first;
//      std::multimap<std::string, unsigned int>::const_iterator case_end = range_check.second;
//      for (;case_it != case_end; case_it++)
//        {
//	      if (case_it->second == load_case_ID) 
//          return true;
//        }
//      
//      // did not find anything, hence, return false
//      return  false;
//      }
//    }
//}
