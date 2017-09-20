// $Id: GMVOutputProcessor.C,v 1.22.6.2 2008-08-21 00:55:21 manav Exp $

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>


// FESystem includes
#include "OutputProcessor/GMVOutputProcessor.h"
#include "FESystem/FESystemController.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "Solvers/FESystemSolverBase.h"
#include "DesignData/DesignParameter.h"
#include "DesignData/DesignDatabase.h"
#include "PostProcess/ElemPostProcessQty.h"
#include "FESystem/AnalysisCase.h"
#include "Mesh/FEMeshData.h"
#include "Mesh/MeshList.h"
#include "PostProcess/PostProcessQtyDatabase.h"
#include "Database/GlobalDataStorage.h"
#include "Discipline/ThermalDisciplineInfo.h"
#include "Discipline/StructuralDisciplineInfo.h"
#include "OutputProcessor/OutputInfo.h"


// libmesh includes
#include "mesh/mesh.h"
#include "geom/elem.h"
#include "geom/node.h"
#include "numerics/numeric_vector.h"


GMVOutputProcessor::GMVOutputProcessor(FESystem::FESystemController& fesys_controller,
                                       const OutputInfo& info):
OutputProcessor(fesys_controller, info, GMV_OUTPUT_PROCESSOR::num()),
mesh(NULL),
mesh_data(NULL),
dummy_elem_set(NULL)
{
	
}




GMVOutputProcessor::~GMVOutputProcessor()
{
	
}



void GMVOutputProcessor::writeData()
{
//   abort();
//   // depending on the discipline, write the output
//   unsigned int discipline = 0;
// //  this->fesystem_controller.analysis_case->getAnalysisDisciplineEnumID();
  
//   unsigned int n_disciplines = 1;
  
//   if (discipline == Discipline::THERMOELASTICITY_DISCIPLINE::num())
//     {
//     n_disciplines = 2;
//     }
  
//   unsigned int mesh_ID =0;
  
//   for (unsigned int i=0; i<n_disciplines; i++)
//     {
//     if (n_disciplines > 1)
//       {
//       switch (i)
//         {
//         case 0:
//           // write the thermal mesh
//           mesh_ID = this->fesystem_controller.analysis_case->getAnalysisDisciplineInfo
//           (Discipline::THERMAL_DISCIPLINE::num()).getMeshID();
//           break;
          
//         case 1:
//           // write the structural mesh
//           mesh_ID = this->fesystem_controller.analysis_case->getAnalysisDisciplineInfo
//           (Discipline::STRUCTURAL_DISCIPLINE::num()).getMeshID();
//           break;
          
//         default:
//           abort();
//           break;
//         }
//       }
    
        
//     // get the mesh and meshdata from mesh list
//     this->mesh = 
//       this->fesystem_controller.mesh_list->getMeshFromID(mesh_ID);
//     this->mesh_data = 
//       this->fesystem_controller.mesh_list->getMeshDataFromID(mesh_ID);
//     this->dummy_elem_set = 
//       this->fesystem_controller.mesh_list->getDummyElemSetFromID(mesh_ID);
    
//     // get the name of the file to be written to
//     std::string file_name = Discipline::AnalysisDisciplineEnum::enumName(discipline);
        
//     std::string name = this->output_info.getOutputFileName();
    
//     file_name += name;
    
    
//     // open the file, and call the functions to write the data to this stream
//     std::fstream output_file;
//     output_file.open(file_name.c_str(), std::fstream::out);
    
//     this->writeHeader(output_file);
//     this->writeNodes(output_file);
//     this->writeElems(output_file);
//     this->writeDisplayNodeIDs(output_file);
//     this->writeDisplayElemIDs(output_file);
//     this->writeVariables(output_file, discipline);
//     this->writeFooter(output_file);
//     }
}



void GMVOutputProcessor::writeHeader(std::ostream& output)
{
//   // input is ascii formatted
//   output << "gmvinput ascii" << std::endl;
}




void GMVOutputProcessor::writeNodes(std::ostream& output)
{
//   // write the xyz locations of all the nodes
//   output << "nodev " <<  this->mesh->n_nodes() << std::endl;
  
//   output.setf(std::ios::showpoint);
	
//   // iterate over the nodes and write the data
//   for (unsigned int node_it=0; node_it < this->mesh->n_nodes(); node_it++)
//     {
//     const Node* node = this->mesh_data->getNodeFromInternalID(node_it);
		
//     for (unsigned int i=0; i<3; i++)
//       output << std::setprecision(8) << (*node)(i) << "  ";
//     output << std::endl;
//     }
}




void GMVOutputProcessor::writeElems(std::ostream& output)
{
//   output << "cells " << (this->mesh->n_elem() - 
//                          this->dummy_elem_set->size()) << std::endl;
  
//   // get the iterator for the elems
//   MeshBase::const_element_iterator elem_it = this->mesh->elements_begin();
//   MeshBase::const_element_iterator elem_end = this->mesh->elements_end();
  
//   std::set<Elem*>::const_iterator dummy_elem_end =
//   this->dummy_elem_set->end();
  
  
	
//   for (; elem_it != elem_end; elem_it++)
//     {
//     Elem* elem = *elem_it;
    
//     if (this->dummy_elem_set->find(elem) != dummy_elem_end)
//       continue;
    
//     switch (elem->type())
//       {
//       case EDGE2:
//         {
//           output << "line 2 " << std::endl;
//         }
//         break;
        
//       case TRI3:
//         {
//           output << "tri 3 " << std::endl;
//         }	
//         break;
        
//       case QUAD4:
//         {
//           output << "quad 4 " << std::endl;
//         }	
//         break;
        
//       case HEX8:
//         {
//           output << "hex 8 " << std::endl;
//         }	
//         break;
				
        
//       case PRISM6:
//         {
//           output << "pprism6 8 " << std::endl;
//         }	
//         break;
        
//       default:
//         abort();
//         break;
//       }
		
//     for (unsigned int i=0; i<elem->n_nodes(); i++)
//       {
//       Node* node = elem->get_node(i);
//       output << this->mesh_data->getInternalIDFromNode(node) << "  " ;
//       }
		
//     output << std::endl;
//     }
}





void GMVOutputProcessor::writeDisplayNodeIDs(std::ostream& output)
{
//   // write the foreign id of all the nodes
//   output << "nodeids " <<  std::endl;
	
//   // iterate over the nodes and write the data
//   for (unsigned int node_it=1; node_it <= this->mesh->n_nodes(); node_it++)
//     {
//     output << this->mesh_data->getNodeForeignIDFromInternalID(node_it) << "  ";
//     }
	
//   output << std::endl;
}




void GMVOutputProcessor::writeDisplayElemIDs(std::ostream& output)
{
//   output << "cellids " << std::endl;
	
//   // get the iterator for the elems
//   MeshBase::const_element_iterator elem_it = this->mesh->elements_begin();
//   MeshBase::const_element_iterator elem_end = this->mesh->elements_end();
  
//   std::set<Elem*>::const_iterator dummy_elem_end =
//   this->dummy_elem_set->end();
  
  
	
//   for (; elem_it != elem_end; elem_it++)
//     {
//     Elem* elem = *elem_it;
    
//     if (this->dummy_elem_set->find(elem) != dummy_elem_end)
//       continue;
		
//     output << this->mesh_data->getForeignIDFromElem(elem) << "  ";
//     }	
	
//   output << std::endl;
}







void GMVOutputProcessor::writeVariables(std::ostream& output,
                                        const unsigned int discipline)
{	
//   output << "variable" << std::endl;
	
//   output.setf(std::ios::showpoint);
	
//   Discipline::AnalysisDisciplineBase* analysis_discipline =
//     this->fesystem_controller.getAnalysisDiscipline(discipline);
  
//   // this is the vector in which the data will be obtained and then written from
//   std::auto_ptr<NumericVector<double> > vector
//     (NumericVector<double>::build().release()); 
	
// //  // iterate over the number of load cases
// //  // get the load case vector from the analysis load case
//   const std::vector<unsigned int> load_cases ;//= 
// //    this->fesystem_controller.analysis_case->getLoadCaseIDs();
	
//   std::vector<unsigned int>::const_iterator load_case_it = load_cases.begin();
//   std::vector<unsigned int>::const_iterator load_case_end = load_cases.end();	
  
//   std::string discipline_name, solution_name;
  
//   // get the name of the file to be written to
//   discipline_name = analysis_discipline->getDisciplineEnumName();

//   solution_name = discipline_name;
//   solution_name += "Solution";
  
//   // iterate over each load case
//   for (; load_case_it != load_case_end; load_case_it++)
//     {
//     std::ostringstream num_to_str;
//     num_to_str << "_lc_" <<*load_case_it;
    
    
//     this->fesystem_controller.global_data_storage->fillVector(*load_case_it,
//                                                               solution_name,
//                                                               *vector);
    
    
//     // iterate over the number of variables at nodes
//     unsigned int n_vars = 
//       analysis_discipline->getNVars(), dof=0;
		
		
//     for (unsigned int var_it =0; var_it < n_vars; var_it++)
//       {
//       std::string name;
//       name = analysis_discipline->getVariableName(var_it);
//       name += num_to_str.str();
      
//       output << name << "  1 " << std::endl;
			
//       // iterate over all the nodes
//       for (unsigned int node_it=1; node_it <= this->mesh->n_nodes(); node_it++)
//         {
// 	      dof = this->mesh_data->getNodeFromInternalID(node_it)->dof_number(0,var_it,0);
// 	      output << std::setprecision(8) << (*vector)(dof) << "  ";
//         }
			
//       output << std::endl;
//       }
    
//     // next, iterate over the DVs, and follow the same procedure as obove
//     // first the property DVs
//     std::auto_ptr<std::vector<DesignData::DesignParameter*> > dv_vector =
//       this->fesystem_controller.design_database->getParameters();
//     std::vector<DesignData::DesignParameter*>::const_iterator dv_it = dv_vector->begin();
//     std::vector<DesignData::DesignParameter*>::const_iterator dv_end = dv_vector->end();
				
//     for (; dv_it != dv_end; dv_it++)
//       {		
//       std::string quantity; 
			
//       std::ostringstream dv_ID_to_str;
//       dv_ID_to_str << (*dv_it)->getID();
      
//       quantity = discipline_name;
//       quantity += "dSolution_dDV_"; 
//       quantity += dv_ID_to_str.str(); // add the DV ID to the variable name
			
//       vector->zero();
			
//       this->fesystem_controller.global_data_storage->fillVector(*load_case_it,
//                                                                 quantity,
//                                                                 *vector);
			
//       // next, write the variable to the output
//       for (unsigned int var_it =0; var_it < n_vars; var_it++)
//         {
// 	      std::string name;
// 	      name += analysis_discipline->getVariableName(var_it);
// 	      name += "_sens_DV";
// 	      name += dv_ID_to_str.str();
// 	      name += num_to_str.str();
				
// 	      output << name << "  1 " << std::endl;
				
// 	      // iterate over all the nodes
// 	      for (unsigned int node_it=1; node_it <= this->mesh->n_nodes(); node_it++)
//           {
//           dof = this->mesh_data->getNodeFromInternalID(node_it)->dof_number(0,var_it,0);
//           output << std::setprecision(8) << (*vector)(dof) << "  ";
//           }
				
// 	      output << std::endl;
//         }
			
//       }				
		
//     }
  
  
//   // next, write the element values
//   // if discipline is structural, thermal gradients and flux will be written
//   // else, if structural analysis is performed, strain and stress will be written 
//   // iterate over all elements
//   // get the post process qty, and write the variables
//   PostProcessQtyDatabase& database = 
//     *(this->fesystem_controller.post_process_qty_database.get());
  
//   std::auto_ptr<std::vector<DesignData::DesignParameter*> > dv_vector = 
//     this->fesystem_controller.design_database->getParameters();
//   unsigned int elem_ID = 0, load_case_ID = 0;
  
//   std::vector<unsigned int>::const_iterator lc_it, lc_begin, lc_end;
//   lc_begin = load_cases.begin();
//   lc_end = load_cases.end();
  
//   std::vector<DesignData::DesignParameter*>::const_iterator dv_it, dv_begin, dv_end;
//   dv_begin = dv_vector->begin();
//   dv_end = dv_vector->end();
  
//   switch (discipline)
//     {
//     case THERMAL_DISCIPLINE_ENUM_ID:
//       {
//         // for now, do nothing. This will be reimplemented though
//       }
//       break;
      
//     case STRUCTURAL_DISCIPLINE_ENUM_ID:
//       {
//         // there will be 6 string streams, 3 for strains and 3 for stresses
        
        
//         lc_it = lc_begin;
//         for ( ; lc_it != lc_end; lc_it++)
//           {
//           load_case_ID = *lc_it;
          
//           std::ostringstream lc_num;
//           lc_num << load_case_ID;
          
//           std::ostringstream 
//             stream_vector_0, stream_vector_1, stream_vector_2, stream_vector_3, 
//             stream_vector_4, stream_vector_5, stream_vector_6, stream_vector_7;;
            
//             // write the names of the quantities to the streams
//             stream_vector_0 << "Strain-XX_lc_" << lc_num.str() << "   " << 0 << std::endl; 
//             stream_vector_1 << "Strain-YY_lc_" << lc_num.str() << "   " << 0 << std::endl; 
//             stream_vector_2 << "Strain-XY_lc_" << lc_num.str() << "   " << 0 << std::endl; 
            
//             stream_vector_3 << "Stress-XX_lc_" << lc_num.str() << "   " << 0 << std::endl; 
//             stream_vector_4 << "Stress-YY_lc_" << lc_num.str() << "   " << 0 << std::endl; 
//             stream_vector_5 << "Stress-XY_lc_" << lc_num.str() << "   " << 0 << std::endl; 
            
            
//             stream_vector_6 << "Strain-Energy-Density_lc_" << lc_num.str()  << "  "  << 0 << std::endl;
//             stream_vector_7 << "VonMises_lc_" << lc_num.str()  << "  "  << 0 << std::endl;
            
            
//             // get the iterator for the elems
//             MeshBase::const_element_iterator elem_it = this->mesh->elements_begin();
//             MeshBase::const_element_iterator elem_end = this->mesh->elements_end();
            
//             std::set<Elem*>::const_iterator dummy_elem_end =
//               this->dummy_elem_set->end();
            
//             for (; elem_it != elem_end; elem_it++)
//               {
//               Elem* elem = *elem_it;
              
//               if (this->dummy_elem_set->find(elem) != dummy_elem_end)
//                 continue;
              
              
//               elem_ID = this->mesh_data->getForeignIDFromElem(elem);
              
//               ElemPostProcessQty& qty =  
//                 dynamic_cast<ElemPostProcessQty&>(database.getElementPostProcessQty
//                                                   (discipline,
//                                                    elem_ID));
              
//               // get the strain and stress tensor
//               TensorValue<double>& strain = qty.getStrainTensor(load_case_ID);
//               TensorValue<double>& stress = qty.getStressTensor(load_case_ID);
              
//               // write the values to the streams
//               stream_vector_0 << strain(0,0) << "  ";
//               stream_vector_1 << strain(1,1) << "  ";
//               stream_vector_2 << strain(0,1) << "  ";
              
//               stream_vector_3 << stress(0,0) << "  ";
//               stream_vector_4 << stress(1,1) << "  ";
//               stream_vector_5 << stress(0,1) << "  ";
              
//               stream_vector_6 << qty.getStrainEnergyDensity(*lc_it) << "  ";
//               stream_vector_7 << qty.getVonMisesStress(*lc_it) << "  " ;
//               }	
            
//             // now write the variables to the output
//             output << stream_vector_0.str() << std::endl;
//             output << stream_vector_1.str() << std::endl;
//             output << stream_vector_2.str() << std::endl;
//             output << stream_vector_3.str() << std::endl;
//             output << stream_vector_4.str() << std::endl;
//             output << stream_vector_5.str() << std::endl;
//             output << stream_vector_6.str() <<  std::endl;
//             output << stream_vector_7.str() << std::endl;
//           }
        
//         // for now, the DVs are not handles. However, this will e changed
//       }
//       break;
      
//     default:
//       abort();
//       break;
//     }
  
  
  
	
//   output << "endvars" << std::endl;
}





void GMVOutputProcessor::writeFooter(std::ostream& output)
{
//   output << "endgmv" << std::endl;
}

