// $Id: ElementDataStorage.h,v 1.11.6.1 2008-02-25 04:20:26 manav Exp $

#ifndef __fesystem_element_data_storage_h__
#define __fesystem_element_data_storage_h__

// C++ includes
#include <map>
#include <string>
#include <memory>

// FESystem includes


// libMesh includes
#include "numerics/dense_vector.h"
#include "numerics/dense_matrix.h"



/**
*	this class will store all the calculated element quantities 
 *  for use in the assembly and sensitivity analysis process
 */

class ElementDataStorage
{
public:
	
  ElementDataStorage();
  
  ~ElementDataStorage();
		
  /**
    *	function to clear the data stored
   */
  void clear();
	
  /**
    *	function to add the element data. This function has been 
   *	overloaded to allow usage for both vector and matrix data.
   *	the DenseMatrix<double> and DenseVector<double> should be used
   *	as the data
   *	The data is stored for elem, tag, and load_case
   */
	
	
  void add_element_data(const DenseMatrix<double>& ,
                        const unsigned int discipline_enum_ID,
                        const unsigned int elem_ID, 
                        const std::string& tag, 
                        const unsigned int case_ID = 0);
	
  void add_element_data(const DenseVector<double>& ,
                        const unsigned int discipline_enum_ID,
                        const unsigned int elem_ID, 
                        const std::string& tag, 
                        const unsigned int case_ID = 0);
		
    /// checks if the specified quantitiy exists in the database or not
    bool checkIfQtyExists(const unsigned int discipline_enum_ID,
                                const unsigned int elem_ID,
                                const std::string& tag, 
                                const unsigned int case_ID=0);
	
//  /// this function sets the data of this quantity to null, so that 
//  /// instead of storing a zero matrix or vector, it could be set to null, 
//  /// to save space. Also, while retireving the quantity, existence of 
//  /// a null will be checked for it
//  void set_to_null(const unsigned int discipline_enum_ID,
//                   const unsigned int elem_ID,
//                   const std::string& tag, 
//                   const unsigned int case_ID=0);
//	
//  bool check_if_null(const unsigned int discipline_enum_ID,
//                     const unsigned int elem_ID,
//                     const std::string& tag, 
//                     const unsigned int case_ID=0);
//	
  /**
    *	function to get the element quantity from the stored data
   */
  void get_element_data(DenseMatrix<double>& qty,
                        const unsigned int discipline_enum_ID,
                        const unsigned int elem_ID, 
                        const std::string& tag,
                        const unsigned int case_ID=0);
	
  /**
    *	function to get the element quantity from the stored data
   */
  void get_element_data(DenseVector<double>& qty,
                        const unsigned int discipline_enum_ID,
                        const unsigned int elem_ID, 
                        const std::string& tag,
                        const unsigned int case_ID=0);

private:

    /// checks if the specified quantitiy exists in the database or not
    bool checkIfMatrixQtyExists(const unsigned int discipline_enum_ID,
                                const unsigned int elem_ID,
                                const std::string& tag, 
                                const unsigned int case_ID=0);
    
    
    /// checks if the specified quantitiy exists in the database or not
    bool checkIfVectorQtyExists(const unsigned int discipline_enum_ID,
                                const unsigned int elem_ID,
                                const std::string& tag, 
                                const unsigned int case_ID=0);
    
//    /// this map will strore names of all quantities that have been set to null	
//    std::map<unsigned int, std::multimap<std::string, unsigned int> > null_qty_map;
		
  
  typedef std::map<unsigned int, DenseMatrix<double>* > LoadCaseIDToMatrixMap;
  typedef std::map<unsigned int, DenseVector<double>* > LoadCaseIDToVectorMap;
  typedef std::map<const std::string , LoadCaseIDToMatrixMap> TagToMatrixMap;
  typedef std::map<const std::string , LoadCaseIDToVectorMap> TagToVectorMap;
  typedef std::map<unsigned int, TagToMatrixMap> ElemIDToTagMatrixMap;
  typedef std::map<unsigned int, TagToVectorMap> ElemIDToTagVectorMap;
  typedef std::map<unsigned int, ElemIDToTagMatrixMap> DisciplineToElemIDMatrixMap;
  typedef std::map<unsigned int, ElemIDToTagVectorMap> DisciplineToElemIDVectorMap;
  
  
  /**
    *	the map for matrix quantities
   */
  DisciplineToElemIDMatrixMap matrix_data;
		
  
  
  /**
    *	the map for vector quantities
   */
  DisciplineToElemIDVectorMap vector_data;
	
};


#endif // __fesystem_element_data_storage_h__
