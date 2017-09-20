// $Id: GlobalDataStorage.h,v 1.10.4.3 2007-05-08 05:18:36 manav Exp $

#ifndef __fesystem_global_data_storage_h__
#define __fesystem_global_data_storage_h__

// C++ includes
#include <set>
#include <string>
#include <memory>
#include <fstream>

// FESystem includes
#include "Database/DataInfo.h"
#include "FESystem/FESystemExceptions.h"

// libMesh includes
#include "numerics/numeric_vector.h"
#include "numerics/sparse_matrix.h"


namespace FESystem
{
  class FESystemController;
}



namespace FESystemDatabase
{
  // forward decleration of data info classes
  class DataInfoBase;
  
  /**
  *	this class will store all the calculated element quantities 
   *  for use in the assembly and sensitivity analysis process
   */
  
  class GlobalDataStorage
  {
public:
    
    GlobalDataStorage(FESystem::FESystemController& fesystem,
                      const std::string& name);
    
    ~GlobalDataStorage();
    
    /**
      *	function to clear the data stored
     */
    void clear();
		
    /// stores the matrix in a PETSC binary format. the arguements are
    /// @param load case number
    /// @param matrix name
    /// @param matrix to be stored
    template<class MatrixType>
      void storeMatrix( const FESystemDatabase::DataInfoBase& data_info,
                        MatrixType& );
    
    /// stores the vector in an ascii format. the arguements are
    /// @param load case number
    /// @param matrix name
    /// @param matrix to be stored
    template<class VectorType>
      void storeVector( const FESystemDatabase::DataInfoBase& data_info,
                        VectorType& );
    
    /// stores the matrix in an ascii format. the arguements are
    /// @param load case number
    /// @param matrix name
    /// @param matrix to be stored
    template<class DataType>
      void asciiStore( const FESystemDatabase::DataInfoBase& data_info,
                       DataType& );
    
    /**
      *	function to get a quantity from the stored data
     */
    template<class MatrixType>
      void fillMatrix(const FESystemDatabase::DataInfoBase& data_info,
                      MatrixType& );
    
    /**
      *	function to get a quantity from the stored data
     */
    template<class VectorType>
      void fillVector(const FESystemDatabase::DataInfoBase& data_info,
                      VectorType& );
    
    
    /// checks if the matrix exists, and returns it as a boolean.
    /// the input arguements are 
    /// @param load case number, 
    /// @param matrix name
    bool checkIfQtyExists(const FESystemDatabase::DataInfoBase& data_info);
    
private:
      
      /// this method stores a given array to the HDF5 database file. This can 
      /// be used for both vectors and matrices. For vectors, specify the 
      /// second dimension as 1.
      void storeMatrixToHDF5(const std::string& file_name, 
                             const unsigned int m,
                             const unsigned int n,
                             const double *mat_data,
                             const bool if_overwrite);
    
    
    /// this method stores a given array to the HDF5 database file. This can 
    /// be used for both vectors and matrices. For vectors, specify the 
    /// second dimension as 1.
    /// The user should now have any data pointed to by mat_data.
    /// Upon return, the matrix pointer \p mat_data will have a newly created
    /// array of data of the same dimensions as mxn. This should be deleted by 
    /// user.
    void fillMatrixFromHDF5(const std::string& file_name, 
                            unsigned int& m,
                            unsigned int& n,
                            double **mat_data);
    
    
    /// reference to the FESystemController object
    FESystem::FESystemController& fesystem_controller;
    
    /// name of the HDF5 database name
    std::string database_name;
    
    /// HDF5 database file id
    unsigned int hdf5_file_id;
		
    /// map to keep track of quantities that have been stored in the database
    std::set<std::string> stored_quantity;
    
  };
  
  
  
  
  
  template<class DataType>
    void GlobalDataStorage::asciiStore(const FESystemDatabase::DataInfoBase& data_info,
                                       DataType& matrix)
    {	
      std::string file_name = data_info.getPath();
      
      std::fstream output;
      output.open(file_name.c_str(), std::fstream::out);
      
      // write the matrix elements to the file
      output << matrix;
      
      // close the file stream
      output.close();
    }
  
  
    
  
  template< typename VectorType>
    void GlobalDataStorage::storeVector(const FESystemDatabase::DataInfoBase& data_info,
                                        VectorType& vector)
    {		
      std::string file_name = data_info.getPath();
      
      // create an array and copy the values in it
      unsigned int n_vals = vector.size();
      double* vals = new double[n_vals];
      
      for (unsigned int i=0; i<n_vals; i++)
        vals[i] = vector.el(i);
      
      bool if_overwrite = false;
      if (this->checkIfQtyExists(data_info))
        if_overwrite = true;

      this->storeMatrixToHDF5(file_name, n_vals, 1, &(vals[0]), if_overwrite);

      delete[] vals;
      
      // record the name in the map, if it already does not exist
      if ( ! this->checkIfQtyExists(data_info))
        this->stored_quantity.insert(file_name);
    }
  
  
  
  template< typename VectorType>
    void GlobalDataStorage::fillVector(const FESystemDatabase::DataInfoBase& data_info,
                                       VectorType& vector)
    {
      // make sure that the data has been stored earlier
      Assert(this->checkIfQtyExists(data_info), ExcInternalError());

      std::string file_name = data_info.getPath();
      
      // create an array and copy the values in it
      unsigned int n_vals, tmp;
      double *vals = NULL;
      
      this->fillMatrixFromHDF5(file_name, n_vals, tmp, &vals);
      
      Assert(tmp == 1, ExcInternalError());
      Assert(vals != NULL, ExcInternalError());
      
      if (vector.size() != n_vals)
        vector.resize(n_vals);
      
      for (unsigned int i=0; i<n_vals; i++)
        vector.set(i, *(vals+i));
      
      delete[] vals;
    }
  
  
  
  template< typename MatrixType>
    void GlobalDataStorage::storeMatrix(const FESystemDatabase::DataInfoBase& data_info,
                                        MatrixType& matrix)
    {	
      std::string file_name = data_info.getPath();
      
      // create an array and copy the values in it
      unsigned int m = matrix.m(), n = matrix.n();
      double* vals = new double[m*n];
      
      for (unsigned int i=0; i<m; i++)
        for (unsigned int j=0; j<n; j++)
          vals[i*n+j] = matrix.el(i,j);
      
      bool if_overwrite = false;
      if (this->checkIfQtyExists(data_info))
        if_overwrite = true;

      this->storeMatrixToHDF5(file_name, m, n, &(vals[0]), if_overwrite);
      
      delete[] vals;
      
      // record the name in the map, if it already does not exist
      if ( ! this->checkIfQtyExists(data_info))
        this->stored_quantity.insert(file_name);
    }
  
  
  
  
  template<typename MatrixType >
    void GlobalDataStorage::fillMatrix(const FESystemDatabase::DataInfoBase& data_info,
                                       MatrixType& matrix)
    {
      // make sure that the data has been stored earlier
      Assert(this->checkIfQtyExists(data_info), ExcInternalError());

      std::string file_name = data_info.getPath();
      
      // create an array and copy the values in it
      unsigned int m, n;
      double* vals = NULL;
      
      this->fillMatrixFromHDF5(file_name, m, n, &vals);
      
      Assert(vals != NULL, ExcInternalError());
      
      if (matrix.m() != m || matrix.n() != n)
        matrix.resize(m,n);
      
      for (unsigned int i=0; i<m; i++)
        for (unsigned int j=0; j<n; j++)
          matrix.set(i, j, *(vals+(i*n+j)));
      
      delete[] vals;
    }
  
  
  
  
  
  template< >
    inline
    void GlobalDataStorage::storeMatrix(const FESystemDatabase::DataInfoBase& data_info,
                                        SparseMatrix<double>& matrix)
    {	
      std::string file_name = data_info.getPath();
      
      matrix.write_to_binary_file(file_name);
      
      // record the name in the map, if it already does not exist
      if ( ! this->checkIfQtyExists(data_info))
        this->stored_quantity.insert(file_name);
    }
  
  
  
  
  template< >
    inline
    void GlobalDataStorage::storeVector(const FESystemDatabase::DataInfoBase& data_info,
                                        NumericVector<double>& vector)
    {	
      std::string file_name = data_info.getPath();
      
      vector.write_to_binary_file(file_name);	
      
      
      // record the name in the map, if it already does not exist
      if ( ! this->checkIfQtyExists(data_info))
        this->stored_quantity.insert(file_name);
    }
  
  
  
  template< >
    inline
    void GlobalDataStorage::fillMatrix(const FESystemDatabase::DataInfoBase& data_info,
                                       SparseMatrix<double>& matrix)
    {
      // make sure that the data has been stored earlier
      Assert(this->checkIfQtyExists(data_info), ExcInternalError());

      std::string file_name = data_info.getPath();
      
      // make sure that the file exists
      
      matrix.load_from_binary_file(file_name);
    }
  
  
  
  template< >
    inline
    void GlobalDataStorage::fillVector(const FESystemDatabase::DataInfoBase& data_info,
                                       NumericVector<double>& vector)
    {
      // make sure that the data has been stored earlier
      Assert(this->checkIfQtyExists(data_info), ExcInternalError());

      std::string file_name = data_info.getPath();
      
      // make sure that the file exists
      
      
      vector.load_from_binary_file(file_name);
    }
  
}


#endif // __fesystem_global_data_storage_h__
