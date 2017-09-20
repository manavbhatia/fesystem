// $Id: GlobalDataStorage.C,v 1.11.4.4 2008-08-21 00:53:20 manav Exp $

// C++ includes
#include <sstream>
#include <fstream>

// FESystem includes
#include "Database/GlobalDataStorage.h"
#include "FESystem/FESystemExceptions.h"
#include "FESystem/FESystemController.h"
#include "FESystem/FESystemNumbers.h"
#include "Numerics/VectorBase.h"
#include "Numerics/MatrixBase.h"

// libMesh Includes
#include "numerics/dense_matrix.h"
#include "numerics/dense_vector.h"

// HDF5 include
#include "hdf5.h"


FESystemDatabase::GlobalDataStorage::GlobalDataStorage
(FESystem::FESystemController& fesystem,
 const std::string& name):
fesystem_controller(fesystem),
database_name(name),
hdf5_file_id(FESystemNumbers::InvalidID)
{
  // first create a CREATION property list
  int creation_property_list_id = H5Pcreate(H5P_FILE_CREATE);
  
  // and an ACCESS property list
  int access_property_list_id = H5Pcreate(H5P_FILE_ACCESS);
  
  // if only one processor is being used, then a local accessible file can be created, 
  // otherwise, the file needs to be opened with parallel access.
  if (FESystem::total_processors > 1)
    {
//      hid_t ierr = H5Pset_fapl_mpio(access_property_list_id, 
//                                    FESystem::COMM_WORLD,
//                                    FESystem::MPI_INFO);
//      AssertThrow (ierr >= 0, ExcInternalError());
		AssertThrow (false, ExcInternalError());
    }
  
  
  // now, using the creation and access property list, create a file
  this->hdf5_file_id = H5Fcreate(this->database_name.c_str(),
                                 H5F_ACC_TRUNC,
                                 creation_property_list_id,
                                 access_property_list_id);
  
  // make sure that the file was created without problems.
  AssertThrow (this->hdf5_file_id > 0, ExcInternalError());
  
  
  // now destroy the property lists
  int ierr = 0;
  ierr = H5Pclose(creation_property_list_id);
  AssertThrow (ierr >= 0, ExcInternalError());
  
  ierr = H5Pclose(access_property_list_id);
  AssertThrow (ierr >= 0, ExcInternalError());	
}



FESystemDatabase::GlobalDataStorage::~GlobalDataStorage()
{
  // close the HDF5 file 
  int ierr = 0;
  ierr = H5Fclose(this->hdf5_file_id);
  AssertThrow (ierr >= 0, ExcInternalError());
}




bool 
FESystemDatabase::GlobalDataStorage::checkIfQtyExists
(const FESystemDatabase::DataInfoBase& data_info)
{
  std::string file_name = data_info.getPath();
  
  unsigned int count = this->stored_quantity.count(file_name);
  
  if (count == 0)
    return false;
  else 
    return true;
}





void FESystemDatabase::GlobalDataStorage::storeMatrixToHDF5(const std::string& file_name, 
                                                            const unsigned int m,
                                                            const unsigned int n,
                                                            const double *mat_data,
                                                            const bool if_overwrite)
{	
	// write this matrix to the hdf5 database file
  
  hid_t property_id = 0, dataspace_id=0, dataset_id=0, transfer_property_id =0,
    ierr = 0;
  
  // create dataset creation property list
  property_id = H5Pcreate(H5P_DATASET_CREATE);
  AssertThrow (property_id > 0, ExcInternalError());
  
  // set layout to contiguous data
  ierr = H5Pset_layout(property_id,H5D_CONTIGUOUS);
  AssertThrow(ierr >= 0, ExcInternalError());
  
  // fill value is left at 0.0, which is default
  // this also lets HDF5 handle the fill time at default
  // alloc time is also left to default  
  
  // create dataspace
  hsize_t dims[2];
  dims[0] = m;
  dims[1] = n;
  
  dataspace_id = H5Screate_simple(2, dims, NULL);
  AssertThrow (dataspace_id > 0, ExcInternalError());
  
  // create dataset
  if (!if_overwrite)
    {
    dataset_id = H5Dcreate1(this->hdf5_file_id, file_name.c_str(), H5T_NATIVE_DOUBLE, 
                           dataspace_id, property_id);
    AssertThrow (dataset_id > 0, ExcInternalError());
    }
  else
    {
    dataset_id  = H5Dopen1(this->hdf5_file_id, file_name.c_str());
    AssertThrow (dataset_id > 0, ExcInternalError());
    }
  
  // create transfer property list
  transfer_property_id = H5Pcreate(H5P_DATASET_XFER);
  AssertThrow (transfer_property_id > 0, ExcInternalError());  
  
  // write dataset
  ierr = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                  transfer_property_id, mat_data);
  AssertThrow (ierr >= 0, ExcInternalError());
  
  // destroy the objects created
  ierr = H5Pclose(property_id);   AssertThrow (ierr >= 0, ExcInternalError());
  ierr = H5Pclose(transfer_property_id);   AssertThrow (ierr >= 0, ExcInternalError());
  ierr = H5Sclose(dataspace_id);   AssertThrow (ierr >= 0, ExcInternalError());
  ierr = H5Dclose(dataset_id);   AssertThrow (ierr >= 0, ExcInternalError());
	
}





void FESystemDatabase::GlobalDataStorage::fillMatrixFromHDF5(const std::string& file_name,
                                           unsigned int &m,
                                           unsigned int &n,
                                           double **mat_data)
{
//	std::string file_name = this->createName(load_case_ID, quantity);
	
  hid_t ndims = 0, dataspace_id=0, dataset_id=0, transfer_property_id =0,
    ierr = 0;
    
  // open dataset
  dataset_id = H5Dopen1(this->hdf5_file_id, file_name.c_str());
  AssertThrow (dataset_id > 0, ExcInternalError());
  
  // get the dataspace for this dataset
  dataspace_id = H5Dget_space(dataset_id);
  AssertThrow (dataspace_id > 0, ExcInternalError());

  // get the dimesnsion of the stored matrix, and make sure that the dimensions agree 
  // with the matrix
  ndims = H5Sget_simple_extent_ndims(dataspace_id);
  AssertThrow (ndims == 2, ExcInternalError());
  
  hsize_t dims[2];

  ierr = H5Sget_simple_extent_dims(dataspace_id, dims,NULL);
  AssertThrow (ierr >= 0, ExcInternalError());

  m = dims[0];
  n = dims[1];
  *mat_data = new double[m*n];
  
  // create transfer property list
  transfer_property_id = H5Pcreate(H5P_DATASET_XFER);
  AssertThrow (transfer_property_id > 0, ExcInternalError());
    
  // read dataset
  ierr = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                  transfer_property_id, *mat_data);
  AssertThrow (ierr >= 0, ExcInternalError());
  
  // destroy the objects created
  ierr = H5Pclose(transfer_property_id);   AssertThrow (ierr >= 0, ExcInternalError());
  ierr = H5Sclose(dataspace_id);   AssertThrow (ierr >= 0, ExcInternalError());
  ierr = H5Dclose(dataset_id);   AssertThrow (ierr >= 0, ExcInternalError());
}





