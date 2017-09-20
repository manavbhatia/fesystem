// $Id: plate_MITC4.C,v 1.11.6.2 2008-02-25 04:32:54 manav Exp $

// C++ includes

// FESystem includes
#include "StructuralElems/plate_MITC4.h"
#include "Discipline/AnalysisDisciplineBase.h"
#include "AnalysisDriver/AnalysisDriver.h"
#include "PostProcess/ElemPostProcessQty.h"
#include "DesignData/ShapeParameter.h"
#include "DesignData/PropertyParameter.h"
#include "Properties/Isotropic2DElemDataCard.h"
#include "FESystem/AnalysisCase.h"


// libMesh includes
#include "geom/face_quad4.h"



FESystemElem::PlateMITC4::PlateMITC4(Discipline::AnalysisDisciplineBase& discipline):
FESystemElem::StructuralElem(2, FESystemElem::STRUCTURAL_PLATE_MITC4_QUAD4::num(), discipline)
{
  
}





FESystemElem::PlateMITC4::~PlateMITC4()
{
	
}





void
FESystemElem::PlateMITC4::getFETypes(std::vector<FEType>& fetypes)
{
  fetypes.clear();
  
  fetypes.push_back(FEType());
  FEType& fe = fetypes.back();
  fe.order = FIRST;
  fe.family = LAGRANGE;
}





void
FESystemElem::PlateMITC4::getQuadratureRules
(std::map<FEFamily, std::pair<QuadratureType, Order> >& quadratures)
{
  quadratures.clear();
  
  bool insert_success = 
    quadratures.insert(std::map<FEFamily, std::pair<QuadratureType, Order> >::
                       value_type(LAGRANGE, std::make_pair(QGAUSS, FIFTH))).second;
  
  Assert(insert_success, ExcInternalError());
}




void
FESystemElem::PlateMITC4::calculate_M(DenseMatrix<double>* matrix, 
                                      const unsigned int design_point_enum_ID,
                                      bool sensitivity_calculation)
{
  // params not used here
  (void) matrix;
  (void) design_point_enum_ID;
  (void) sensitivity_calculation;
  
  abort();
}



void
FESystemElem::PlateMITC4::calculate_K_G(DenseMatrix<double>* matrix, 
                                        const unsigned int design_point_enum_ID,
                                        bool sensitivity_calculation)
{
  // params not used here
  (void) matrix;
  (void) design_point_enum_ID;
  (void) sensitivity_calculation;

  abort();
}




void
FESystemElem::PlateMITC4::calculate_K(DenseMatrix<double>* matrix, 
                                      const unsigned int design_point_enum_ID,
                                      bool sensitivity_calculation)
{
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  const double k = 5.0/6.0;
	
	static DenseMatrix<double> factor(3,3), factor_bend(3,3);
  factor.zero(); factor_bend.zero();
  
  switch (sensitivity_calculation && (this->sensitivity_parameter == DesignData::PROPERTY_PARAMETER::num()))
    {
    case true:
      {
        this->elem_property_card->getFactorSensitivityForGlobalParameter
        (factor, STIFFNESS_A_MATRIX_FACTOR::num(), this->sensitivity_parameter_ID);
        this->elem_property_card->getFactorSensitivityForGlobalParameter
          (factor_bend, STIFFNESS_D_MATRIX_FACTOR::num(), this->sensitivity_parameter_ID);
      }
      break;
      
    default:
      {
        this->elem_property_card->getFactor(factor, STIFFNESS_A_MATRIX_FACTOR::num());
        this->elem_property_card->getFactor(factor_bend, STIFFNESS_D_MATRIX_FACTOR::num());
      }
    }
  
  static DenseMatrix<double> shear_stiff(6*n_nodes, 6*n_nodes), Nx_Nx(n_nodes,n_nodes),
    Ny_Ny(n_nodes,n_nodes), Nx_Ny(n_nodes,n_nodes);
  
  shear_stiff.zero(); Nx_Nx.zero(); Ny_Ny.zero(); Nx_Ny.zero();
	
  switch (sensitivity_calculation && (this->sensitivity_parameter == DesignData::SHAPE_PARAMETER::num()))
    {
    case true:
      {
        this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_X_N_X_FACTOR::num(), &Nx_Nx,
                                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_Y_N_Y_FACTOR::num(), &Ny_Ny,
                                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_X_N_Y_FACTOR::num(), &Nx_Ny,
                                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        this->getPlateMITC4Qty(FESystemElem::PLATE_MITC4_SHEAR_STIFFNESS_FACTOR::num(), &shear_stiff,
                               design_point_enum_ID, FESystemElem::ELEM_VOLUME::num(), true);
      }
      break;
      
    default:
      {
        this->getFESystemElemQty(FESystemElem::N_X_N_X_FACTOR::num(), &Nx_Nx,
                                 design_point_enum_ID,
                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        this->getFESystemElemQty(FESystemElem::N_Y_N_Y_FACTOR::num(), &Ny_Ny,
                                 design_point_enum_ID,
                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        this->getFESystemElemQty(FESystemElem::N_X_N_Y_FACTOR::num(), &Nx_Ny,
                                 design_point_enum_ID,
                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        this->getPlateMITC4Qty(FESystemElem::PLATE_MITC4_SHEAR_STIFFNESS_FACTOR::num(),
                               &shear_stiff, design_point_enum_ID,  
                               FESystemElem::ELEM_VOLUME::num(), false);
      }
    }
  
	
  // calculate the stiffness matrix
  for (unsigned int i=0; i<n_nodes; i++)
    for (unsigned int j=0; j<n_nodes; j++)
      {
      (*matrix)(i,j) = (factor(0,0)* Nx_Nx(i,j) + factor(2,2)* Ny_Ny(i,j));
      (*matrix)(n_nodes+i,j) = (factor(0,1) * Nx_Ny(j,i) + factor(2,2) * Nx_Ny(i,j));
      (*matrix)(i,n_nodes+j) = (factor(0,1) * Nx_Ny(i,j) + factor(2,2) * Nx_Ny(j,i));
      (*matrix)(n_nodes+i,n_nodes+j) = (factor(0,0)* Ny_Ny(i,j) + factor(2,2) * Nx_Nx(i,j));
      
      (*matrix)(3*n_nodes+i,3*n_nodes+j) += (factor_bend(0,0)*Ny_Ny(i,j) + 
                                             factor_bend(2,2)* Nx_Nx(i,j));
      (*matrix)(4*n_nodes+i,3*n_nodes+j) += -(factor_bend(0,1)*Nx_Ny(i,j) +
                                              factor_bend(2,2)* Nx_Ny(j,i));
      (*matrix)(3*n_nodes+i,4*n_nodes+j) += -(factor_bend(0,1)*Nx_Ny(j,i) + 
                                              factor_bend(2,2)* Nx_Ny(i,j));
      (*matrix)(4*n_nodes+i,4*n_nodes+j) += (factor_bend(0,0)*Nx_Nx(i,j) +
                                             factor_bend(2,2)* Ny_Ny(i,j));
      }
      
      // now add the effects of the shear strains
      shear_stiff.scale(factor(2,2));
  shear_stiff.scale(k);
  matrix->add(1.0, shear_stiff);
		
  //  std::cout << *(shear_stiff.get()) << std::endl;
  
  // finally, place a unity at the diagonal for the theta_z rows, to remove singularity 
  // of the element
  // this will be changed in future, since the value of this fictitious stiffness is dependent on the 
  // value of the other diagonal terms 
  for (unsigned int i=0; i<4; i++)
    (*matrix)(5*n_nodes+i,5*n_nodes+i) = 1.0;
}





void 
FESystemElem::PlateMITC4::calculateShearStiffnessFactor(DenseMatrix<double>* matrix,
                                                        const unsigned int design_point_enum_ID)
{
  assert (matrix != NULL);
  
  switch (matrix->m() +  matrix->n() - 48)
    {
    case 0:
      {
        // keep going
      }
      break;
      
    default:
      matrix->resize(24,24);
      break;
    }
  
  matrix->zero();
  
  // initialize the element if it is not already initialized
  switch (!this->local_elem_is_initialized_for_DV_map[design_point_enum_ID])
      this->initialize_element(design_point_enum_ID);
  
  //   unsigned int num_nodes = this->base_elem->n_nodes();
  
  //   DenseMatrix<double> Nx_Nx(num_nodes,num_nodes),
  //     Ny_Ny(num_nodes,num_nodes),
  //     Nx_N(num_nodes,num_nodes), Ny_N(num_nodes,num_nodes),
  //     N_N(num_nodes,num_nodes);
  
  //   this->getFactor(&Nx_Nx, StructuralElem::N_X_N_X_FACTOR);
  //   this->getFactor(&Ny_Ny, StructuralElem::N_Y_N_Y_FACTOR);
  //   this->getFactor(&Nx_N, StructuralElem::N_X_N_FACTOR);
  //   this->getFactor(&Ny_N, StructuralElem::N_Y_N_FACTOR);
  //   this->getFactor(&N_N, StructuralElem::N_N_FACTOR);
  
  
  //   for (unsigned int i=0; i<num_nodes; i++)
  //     for (unsigned int j=0; j<num_nodes; j++)
  //       {
  // 	(*matrix)(2*num_nodes+i,2*num_nodes+j) += (Nx_Nx(i,j) + Ny_Ny(i,j));
  // 	(*matrix)(2*num_nodes+i,3*num_nodes+j) += -Ny_N(i,j);
  // 	(*matrix)(2*num_nodes+i,4*num_nodes+j) += Nx_N(i,j);
  
  
  // 	(*matrix)(3*num_nodes+i,2*num_nodes+j) += -Ny_N(j,i);
  // 	(*matrix)(4*num_nodes+i,2*num_nodes+j) += Nx_N(j,i);
	
  // 	(*matrix)(3*num_nodes+i,3*num_nodes+j) += N_N(i,j);
  // 	(*matrix)(4*num_nodes+i,4*num_nodes+j) += N_N(i,j);
  //       }
  
  
  //   return;
	
  double x12, x43, x14, x23, y12, y43, y14, y23; // data for tensorial shear strain integration
  double Ax, Bx, Cx, Ay, By, Cy;
  //  double sin_alpha, cos_alpha, sin_beta, cos_beta;
	
  // from the local element, get the nodes
  // the computaion of the shear factor requires an element 
  // where the (x,y) = (0,0) is at node 2, and 
  // (xi, eta) = (1,1) is at node 0. For the calculation in 
  // FESystemElement, these points are node 0 and 2 respectively. 
  // Also, for this shear factor, required x-axis is from node 2 to node 3,
  // and for the local axis, it is from node 0 to node 1. Hence, 
  // if we simply renumber the nodes as 
  // local elem        bathe's elem
  //    0                   2    
  //    1                   3
  //    2                   0 
  //    3                   1
  // we can use the same axis system for calculation.
  
  static Elem * local_elem;
  local_elem = this->local_elem_for_DV_map[design_point_enum_ID];
  
  const Point& new_point0 = local_elem->point(2);
  const Point& new_point1 = local_elem->point(3);
  const Point& new_point2 = local_elem->point(0);
  const Point& new_point3 = local_elem->point(1);
  
  // calculate data for tensorial strain inegration
  x12 = new_point0(0) - new_point1(0);
  x14 = new_point0(0) - new_point3(0);
  x23 = new_point1(0) - new_point2(0);
  x43 = new_point3(0) - new_point2(0);
	
  y12 = new_point0(1) - new_point1(1);
  y14 = new_point0(1) - new_point3(1);
  y23 = new_point1(1) - new_point2(1);
  y43 = new_point3(1) - new_point2(1);
	
  Ax = new_point0(0) - new_point1(0) - new_point2(0) + new_point3(0);
  Bx = new_point0(0) - new_point1(0) + new_point2(0) - new_point3(0);
  Cx = new_point0(0) + new_point1(0) - new_point2(0) - new_point3(0); 
	
  Ay = new_point0(1) - new_point1(1) - new_point2(1) + new_point3(1); 
  By = new_point0(1) - new_point1(1) + new_point2(1) - new_point3(1); 
  Cy = new_point0(1) + new_point1(1) - new_point2(1) - new_point3(1);
	
  //   // calculate the locations of A,B,C,D, and use that to calculate the 
  //   // cosine and sine for the angles alpha and beta. 
  //   // alpha is the angle betweem r and x
  //   // beta is the angle between s and x
  //   Point A_loc, B_loc, C_loc, D_loc, CA_vec, BD_vec, x_vec;
	
  //   A_loc.assign(new_point0+new_point1); A_loc*0.5;
  //   B_loc.assign(new_point1+new_point2); B_loc*0.5;
  //   C_loc.assign(new_point2+new_point3); C_loc*0.5;
  //   D_loc.assign(new_point3+new_point0); D_loc*0.5;
  
  //   x_vec.assign(new_point3-new_point2); 
  //   CA_vec.assign(A_loc - C_loc);
  //   BD_vec.assign(D_loc - B_loc);
	
  //   cos_beta = (x_vec * CA_vec)/ (x_vec.size() * CA_vec.size());
  //   cos_alpha = (x_vec * BD_vec)/ (x_vec.size() * BD_vec.size());
  //   // following has been done to avoind getting NANs in the code. 
  //   // the problem is that cos_alpha or cos_beta, sometimes ends up being
  //   // equal to 1, and 1.0 - cos^2 gives a very small negative number, which 
  //   // is practically zero. Hence, if this happend, they will be practically 
  //   // set to zero.
  //   // in both cases, if the test does not pass, the code will abort, sicne the 
  //   // cos of an angle cannot be > 1.
  //   double sin_alpha_square = 1.0 - cos_alpha*cos_alpha,
  //     sin_beta_square = 1.0 - cos_beta*cos_beta;
  
  //   if (sin_alpha_square < 0.0)
  //     {
  //       if (fabs(sin_alpha_square) < 1.0e-12)
  // 	sin_alpha_square = 0.0;
  //       else
  // 	abort();
  //     }
  
  //   if (sin_beta_square < 0.0)
  //     {
  //       if (fabs(sin_beta_square) < 1.0e-12)
  // 	sin_beta_square = 0.0;
  //       else
  // 	abort();
  //     }
  
  //   sin_alpha = sqrt(sin_alpha_square);
  //   sin_beta = sqrt(sin_beta_square);
  
  
	
  const std::vector<Real>& JxW = this->fe_base_map_for_DV[design_point_enum_ID][LAGRANGE]->get_JxW();
  const unsigned int n_q_points = 
    this->fe_base_map_for_DV[design_point_enum_ID][LAGRANGE]->n_quadrature_points();
  unsigned int n_nodes = this->getNNodes();
  
  static DenseVector<double> gamma_rz(12), gamma_sz(12), gamma_xz(12), gamma_yz(12);
  gamma_rz.zero(); gamma_sz.zero(); gamma_xz.zero(); gamma_yz.zero();
  
  double xi, eta, ABx_factor, ABy_factor, CBx_factor, CBy_factor, det_J;
  matrix->zero();
  
  for (unsigned int qp=0; qp < n_q_points; qp++)
    {
    gamma_rz.zero(); gamma_sz.zero(); gamma_xz.zero(); gamma_yz.zero();
    
    xi = (this->qbase_map_for_DV[design_point_enum_ID][LAGRANGE]->get_points())[qp](0);
    eta = (this->qbase_map_for_DV[design_point_enum_ID][LAGRANGE]->get_points())[qp](1);
    //            det_J = (JxW[qp])/((this->_qbase_2d->get_weights())[qp]);
    det_J = (1.0 / 16.0) * ((Cy + xi * By)*(Ax + eta * Bx) - 
                            (Cx + xi * Bx)*(Ay + eta * By) );
    
    // numerators of the dot product of 
    ABy_factor = (Ay + eta * By) / 16.0 / det_J;
    ABx_factor = (Ax + eta * Bx) / 16.0 / det_J;
    CBy_factor = (Cy + xi * By) / 16.0 / det_J;
    CBx_factor = (Cx + xi * Bx) / 16.0 / det_J;
    
    gamma_rz(0) =  (1+eta);
    gamma_rz(1) = -(1+eta);
    gamma_rz(2) = -(1-eta);
    gamma_rz(3) =  (1-eta);
		
    gamma_rz(4) = -.5 * (1+eta) * y12;
    gamma_rz(5) = -.5 * (1+eta) * y12;
    gamma_rz(6) = -.5 * (1-eta) * y43;
    gamma_rz(7) = -.5 * (1-eta) * y43;
		
    gamma_rz(8) =  .5 * (1+eta) * x12;
    gamma_rz(9) =  .5 * (1+eta) * x12;
    gamma_rz(10) = .5 * (1-eta) * x43;
    gamma_rz(11) = .5 * (1-eta) * x43;
		
    gamma_sz(0) =  (1+xi);
    gamma_sz(1) =  (1-xi);
    gamma_sz(2) = -(1-xi);
    gamma_sz(3) = -(1+xi);
		
    gamma_sz(4) = -.5 * (1+xi) * y14;
    gamma_sz(5) = -.5 * (1-xi) * y23;
    gamma_sz(6) = -.5 * (1-xi) * y23;
    gamma_sz(7) = -.5 * (1+xi) * y14;
		
    gamma_sz(8) =  .5 * (1+xi) * x14;
    gamma_sz(9) =  .5 * (1-xi) * x23;
    gamma_sz(10) = .5 * (1-xi) * x23;
    gamma_sz(11) = .5 * (1+xi) * x14;
    
    // now initialize the values of gamma_xz, and gamma_yz
    for (unsigned int i=0; i < 12; i++)
      {
      gamma_xz(i) = gamma_rz(i) * CBy_factor - gamma_sz(i) * ABy_factor;
      gamma_yz(i) = -gamma_rz(i) * CBx_factor + gamma_sz(i) * ABx_factor;
      }
    
    //       gamma_rz.scale(1.0/16.0);
    //       gamma_sz.scale(1.0/16.0);
    
    //       std::cout << "Quad point: xi-> " << xi << "  eta-> " << eta << std::endl;
    //       std::cout << " Det J = " << det_J << std::endl; 
    //       std::cout << " gamma_r_z : \n" << gamma_rz << std::endl;
    //       std::cout << " gamma_s_z : \n" << gamma_sz << std::endl;
    //       std::cout << " gamma_x_z : \n" << gamma_xz << std::endl;
    //       std::cout << " gamma_y_z : \n" << gamma_yz << std::endl;
    
    //       std::cout << "sin beta : " << sin_beta << "  " 
    // 		<< (CBy_factor/ sqrt(pow(CBy_factor,2) + pow(CBx_factor,2))) 
    // 		<< std::endl;
    //       std::cout << "sin alpha : " << sin_alpha << "  " 
    // 		<< (ABy_factor/ sqrt(pow(ABy_factor,2) + pow(ABx_factor,2))) 
    // 		<< std::endl;
    //       std::cout << "cos beta : " << cos_beta << "  " 
    // 		<< (CBx_factor/ sqrt(pow(CBy_factor,2) + pow(CBx_factor,2))) 
    // 		<< std::endl;
    //       std::cout << "cos_alpha : " << cos_alpha << "  " 
    // 		<< (ABx_factor/ sqrt(pow(ABy_factor,2) + pow(ABx_factor,2))) 
    // 		<< std::endl;
		
    // Now add the shear contribution to the stiffness matrix.
    for (unsigned int i=0; i<4; i++)
      for (unsigned int j=0; j<4; j++)
        {
        (*matrix)(2*n_nodes+i,2*n_nodes+j) += (JxW[qp])*(gamma_xz(i)*gamma_xz(j) +
                                                         gamma_yz(i)*gamma_yz(j));
        (*matrix)(2*n_nodes+i,3*n_nodes+j) += (JxW[qp])*(gamma_xz(i)*gamma_xz(n_nodes+j) + 
                                                         gamma_yz(i)*gamma_yz(n_nodes+j));
        (*matrix)(2*n_nodes+i,4*n_nodes+j) += (JxW[qp])*(gamma_xz(i)*gamma_xz(2*n_nodes+j) + 
                                                         gamma_yz(i)*gamma_yz(2*n_nodes+j));
				
        (*matrix)(3*n_nodes+i,2*n_nodes+j) += (JxW[qp])*(gamma_xz(n_nodes+i)*gamma_xz(j) + 
                                                         gamma_yz(n_nodes+i)*gamma_yz(j));
        (*matrix)(3*n_nodes+i,3*n_nodes+j) += (JxW[qp])*(gamma_xz(n_nodes+i)*gamma_xz(n_nodes+j) + 
                                                         gamma_yz(n_nodes+i)*gamma_yz(n_nodes+j));
        (*matrix)(3*n_nodes+i,4*n_nodes+j) += (JxW[qp])*(gamma_xz(n_nodes+i)*gamma_xz(2*n_nodes+j) + 
                                                         gamma_yz(n_nodes+i)*gamma_yz(2*n_nodes+j));
				
        (*matrix)(4*n_nodes+i,2*n_nodes+j) += (JxW[qp])*(gamma_xz(2*n_nodes+i)*gamma_xz(j) + 
                                                         gamma_yz(2*n_nodes+i)*gamma_yz(j));
        (*matrix)(4*n_nodes+i,3*n_nodes+j) += (JxW[qp])*(gamma_xz(2*n_nodes+i)*gamma_xz(n_nodes+j) + 
                                                         gamma_yz(2*n_nodes+i)*gamma_yz(n_nodes+j));
        (*matrix)(4*n_nodes+i,4*n_nodes+j) += (JxW[qp])*(gamma_xz(2*n_nodes+i)*gamma_xz(2*n_nodes+j) + 
                                                         gamma_yz(2*n_nodes+i)*gamma_yz(2*n_nodes+j));
        }
    }
}




// bool FESystemElem::PlateMITC4::calculate_F_Pressure(DenseVector<double>* vector, 
// 				      bool sensitivity_calculation)
// {
//   //	const std::vector<std::vector<Real> >& phi = *(this->_phi);	
//   //	
//   //	std::vector<Real>& JxW = *(this->_JxW);
//   //	unsigned int n_nodes = phi.size();
//   //	
//   //	std::vector<Load>::const_iterator load_it = load->begin();
//   //	std::vector<Load>::const_iterator load_end = load->end();
//   //	
//   //	for (; load_it != load_end; load_it++)
//   //	{
//   //		double load_value = load_it->value();
//   //		
//   //		for (unsigned int qp=0; qp<this->_qbase_2d->n_points(); qp++)
//   //			for (unsigned int i=0; i<n_nodes; i++)
//   //			{	
//   //				(*vector)(2*n_nodes+i) += (JxW[qp]) * phi[i][qp] * load_value;
//   //			}	
//   //	}			
//   //		
// }



void
FESystemElem::PlateMITC4::calculate_F_T(DenseVector<double>* vector, 
                                        const unsigned int design_point_enum_ID,
                                        bool sensitivity_calculation)
{
  static double ref_temp = this->analysis_discipline.getFESystemController().
  analysis_case->getRealParameter("REFERENCE_TEMP");
  
  static unsigned int n_nodes;
  n_nodes = this->getNNodes();
  
  // initialize the data structures for the load
  static DenseVector<double> nodal_temp(n_nodes), nodal_temp_sens(n_nodes);

  switch (nodal_temp.size() - n_nodes)
    {
    case 0:
      // keep going
      break;
      
    default:
      nodal_temp.resize(n_nodes);
      nodal_temp_sens.resize(n_nodes);
    }
  
  nodal_temp.zero(); nodal_temp_sens.zero();
  
  // get the loads for this element
  this->extractNodalTemperatureVectorFromLoads(nodal_temp, false);
  
  static double factor = 0.0, factor_sens = 0.0;
  factor = 0.0; factor_sens = 0.0;
  this->elem_property_card->getFactor(factor, THERMAL_EXPANSION_FACTOR::num());  
  
  switch (sensitivity_calculation)
    {
    case true:
      {
        // get the nodal temperature vector sensitivity
        this->extractNodalTemperatureVectorFromLoads(nodal_temp_sens, true);
        
        switch (this->sensitivity_parameter)
          {
          case PROPERTY_PARAMETER_ENUM_ID:
            this->elem_property_card->getFactorSensitivityForGlobalParameter
            (factor_sens, THERMAL_EXPANSION_FACTOR::num(), this->sensitivity_parameter_ID);
            break;
            
          default:
            {
              // keep going
            }
            break;
          }
      }
      break;
      
    default:
      {
        // keep going
      }
    }
  
  static DenseMatrix<double> nx_n_matrix(n_nodes,n_nodes), 
    ny_n_matrix(n_nodes,n_nodes), 
    nx_n_matrix_sens(n_nodes,n_nodes), 
    ny_n_matrix_sens(n_nodes,n_nodes);
  
  nx_n_matrix.zero(); ny_n_matrix.zero(); 
  nx_n_matrix_sens.zero(); ny_n_matrix_sens.zero(); 
  
  // now, scale the matrix with the right factors to take into 
  // account the sensitivity if needed
  switch (sensitivity_calculation)
    {
    case true:
      {
        switch (this->sensitivity_parameter)
          {
          case PROPERTY_PARAMETER_ENUM_ID:
            {
              this->getFESystemElemQty(FESystemElem::N_X_N_FACTOR::num(), &nx_n_matrix,
                                       design_point_enum_ID,
                                       FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
              this->getFESystemElemQty(FESystemElem::N_Y_N_FACTOR::num(), &ny_n_matrix,
                                       design_point_enum_ID,
                                       FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
              
              double sum1 = 0.0, sum2 = 0.0;
              for (unsigned int i=0; i<n_nodes; i++)
                {
                sum1 = 0.0; sum2 = 0.0;
                for (unsigned int j=0; j<n_nodes; j++)
                  {
                  sum1 += nx_n_matrix(i,j) * (factor_sens* (nodal_temp(j)-ref_temp) +
                                              factor * nodal_temp_sens(j));
                  sum2 += ny_n_matrix(i,j) * (factor_sens* (nodal_temp(j)-ref_temp) + 
                                              factor * nodal_temp_sens(j));
                  }
                (*vector)(i) = sum1;
                (*vector)(n_nodes + i) = sum2;
                }
            }
            break;
            
          case SHAPE_PARAMETER_ENUM_ID:
            {
              this->getFESystemElemQty(FESystemElem::N_X_N_FACTOR::num(), &nx_n_matrix,
                                       design_point_enum_ID,
                                       FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
              this->getFESystemElemQty(FESystemElem::N_Y_N_FACTOR::num(), &ny_n_matrix,
                                       design_point_enum_ID,
                                       FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
              this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_X_N_FACTOR::num(),
                                                       &nx_n_matrix_sens,
                                                       FESystemElem::ELEM_VOLUME::num(),
                                                       LAGRANGE);
              this->getFESystemElemQtyShapeSensitivity(FESystemElem::N_Y_N_FACTOR::num(), 
                                                       &ny_n_matrix_sens,
                                                       FESystemElem::ELEM_VOLUME::num(),
                                                       LAGRANGE);
              
              double sum1 = 0.0, sum2 = 0.0;
              for (unsigned int i=0; i<n_nodes; i++)
                {
                sum1 = 0.0; sum2 = 0.0;
                for (unsigned int j=0; j<n_nodes; j++)
                  {
                  sum1 += factor* ( nx_n_matrix_sens(i,j) * (nodal_temp(j)-ref_temp) +
                                    nx_n_matrix(i,j) * nodal_temp_sens(j));
                  sum2 += factor* ( ny_n_matrix_sens(i,j) * (nodal_temp(j)-ref_temp) +
                                    ny_n_matrix(i,j) * nodal_temp_sens(j));
                  }
                (*vector)(i) = sum1;
                (*vector)(n_nodes + i) = sum2;
                }
            }
            break;
            
          default:
            Assert(false, ExcInternalError());
          }
      }
      break;
      
    case false:
    default:
      {
        this->getFESystemElemQty(FESystemElem::N_X_N_FACTOR::num(), &nx_n_matrix, 
                                 design_point_enum_ID,
                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        this->getFESystemElemQty(FESystemElem::N_Y_N_FACTOR::num(), &ny_n_matrix, 
                                 design_point_enum_ID,
                                 FESystemElem::ELEM_VOLUME::num(), LAGRANGE);
        
        double sum1 = 0.0, sum2 = 0.0;
        for (unsigned int i=0; i<n_nodes; i++)
          {
          sum1 = 0.0; sum2 = 0.0;
          for (unsigned int j=0; j<n_nodes; j++)
            {
            sum1 += nx_n_matrix(i,j) * (factor* (nodal_temp(j)-ref_temp));
            sum2 += ny_n_matrix(i,j) * (factor* (nodal_temp(j)-ref_temp));
            }
          (*vector)(i) = sum1;
          (*vector)(n_nodes + i) = sum2;
          }
      }
      break;
    }
}



std::auto_ptr<ElemPostProcessQty> 
FESystemElem::PlateMITC4::getElementPostProcessQty
(std::vector<unsigned int> load_cases, 
 std::vector<DesignData::DesignParameter*> dv_vector)
{
  (void) load_cases;
  (void) dv_vector;

//  std::auto_ptr<ElemPostProcessQty> return_qty(new ElemPostProcessQty(this->elem_ID));
//  
//  TensorValue<double> tensor;
//  
//  // create an fe that will be used for calculationi of the element strains and
//  // stresses. The vector fe_sens is for shape sensitivity
////  std::auto_ptr<FEBase> fe(new FE<2,LAGRANGE>(FEType())),
////    fe_sens (new FE<2,LAGRANGE>(FEType()));
//    
//  
//  unsigned int n_nodes = this->base_elem->n_nodes();
//  
//  static DenseMatrix<double> B_strain(3, 6*n_nodes) , B_strain_sens(3, 6*n_nodes);
//  B_strain.zero(); B_strain_sens.zero();
//  
//  // create rest of the necessary data structures
//  static DenseVector<double>  strain(3), strain_sens(3),
//    stress(3), stress_sens(3), strain_plus_thermal(3), 
//    strain_plus_thermal_sens(3), strain_scratch(3);
//  
//  strain.zero(); strain_sens.zero();
//  stress.zero(); stress_sens.zero(); strain_plus_thermal.zero(); 
//  strain_plus_thermal_sens.zero(); strain_scratch.zero();
//  
//  static DenseMatrix<double> material_factor(3,3), material_factor_sens(3,3);
//  material_factor.zero(); 
//  material_factor_sens.zero();
//  
//  double temp_strain_material_factor = 0.0, temp_strain_material_factor_sens = 0.0,
//    ref_temp = 0.0, thickness = 0.0, thickness_sens = 0.0;
//    
//  // point at which this element will be initialized
//  // for now, this is only the centroid of the element
//  std::vector<Point> point_vec;
//  point_vec.push_back(Point(0.0,0.0,0.0));
//  
//  // then, reinit it at the base elem, at the element centroid, where the 
//  // element strains and stresses will be evaluated
//  this->initialize_element(&point_vec);
//  
//  // get the necessary data out of this element 
//  const std::vector< std::vector<Real> >& dphi_dx = this->fe_base_map[LAGRANGE]->get_dphidx();
//  const std::vector< std::vector<Real> >& dphi_dy = this->fe_base_map[LAGRANGE]->get_dphidy();
//  const std::vector<std::vector<Real> >& phi = this->fe_base_map[LAGRANGE]->get_phi();
//  
//  // calculate B_strain matrix. Since the calculation is being performed at a single point, 
//  // the vector of shape factor derivatives should have size = 1
//  assert (dphi_dx[0].size() == 1);
//  for (unsigned int i=0; i < n_nodes; i++)
//    {
//    B_strain(0,i) = dphi_dx[i][0];  // epsilon_xx
//    B_strain(0,4*n_nodes+i) = 0.5 * thickness * dphi_dx[i][0]; // bending contribution
//    B_strain(1,n_nodes+i) = dphi_dy[i][0];  // epsilon_yy
//    B_strain(1,3*n_nodes+i) = - 0.5 * thickness * dphi_dy[i][0]; // bending contribution
//    B_strain(2,i) = dphi_dy[i][0];  // epsilon  xy 
//    B_strain(2,3*n_nodes+i) = - 0.5 * thickness * dphi_dx[i][0]; // bending contribution
//    B_strain(2,n_nodes+i) = dphi_dx[i][0];  // epsilon xy
//    B_strain(2,4*n_nodes+i) = 0.5 * thickness * dphi_dy[i][0]; // bending contribution
//    }
//  
//  
//  // iterate over all the load cases, to calculate the strains and stresses
//  std::vector<unsigned int>::const_iterator load_case_it = load_cases.begin();
//  std::vector<unsigned int>::const_iterator load_case_end = load_cases.end();	
//	
//  std::vector<DesignData::DesignParameter*>::const_iterator dv_it, dv_begin, dv_end;
//  dv_begin = dv_vector.begin();
//  dv_end = dv_vector.end();
//  
//  
//  static DenseMatrix<double> transform_mat, transform_mat_sens;
//  transform_mat.zero(); transform_mat_sens.zero();
//  
//  this->getStructuralT_matrix(&transform_mat, false);
//  
//  
//  // also, create dof_vectors
//  static DenseVector<double>  dof_values(6*n_nodes), dof_value_sens(6*n_nodes),
//    local_dof(6*n_nodes), local_dof_sens(6*n_nodes), scratch_vec(6*n_nodes),
//    nodal_temp(n_nodes),nodal_temp_sens(n_nodes),
//    dphi_dx_sens(n_nodes), dphi_dy_sens(n_nodes);
//  
//  dof_values.zero(); dof_value_sens.zero();
//  local_dof.zero(); local_dof_sens.zero(); 
//  nodal_temp.zero(); nodal_temp_sens.zero();
//  dphi_dx_sens.zero(); dphi_dy_sens.zero(); 
//  scratch_vec.zero();
//
//  
//  // initialize the material factors
//  this->elem_property_card->getPropertyValueFromMaterialCard(ALPHA_EXPANSION::num(),
//                                                             temp_strain_material_factor);
//  this->elem_property_card->getPropertyValueFromMaterialCard(TEMP_REF::num(),
//                                                             ref_temp);
//  this->elem_property_card->getPropertyValueFromMaterialCard(THICKNESS_2D_ELEM::num(),
//                                                             thickness);
//  
//  this->elem_property_card->getFactor(material_factor, STRESS_STRAIN_FACTOR::num());
//  
//  double temp = 0.0, temp_sens = 0.0;
//  
//  // iterate over each load case, and ask solver to solve for it
//  for (; load_case_it != load_case_end; load_case_it++)
//    {
//    
//    // get the DOF vector for this elem
//    this->analysis_discipline.getElemDofValues(this->base_elem, 
//                                               dof_values,
//                                               *load_case_it);
//    
//    // transform these dofs to the local coordinate system
//    local_dof.zero();
//    transform_mat.right_multiply_vector(dof_values, local_dof);
//    
//    // get the loads for this element
//    this->getNodalTemperatureVector(nodal_temp, *load_case_it);
//    
//    // calculate the factor and the strain and stress
//    B_strain.right_multiply_vector(local_dof, strain);
//    
//    // add the thermal strain component and then calculate the stress
//    for (unsigned int i=0; i<phi.size(); i++)
//      temp += phi[i][0] * nodal_temp(i);
//    
//    strain_plus_thermal.zero();
//    strain_plus_thermal(0) = strain(0) - temp_strain_material_factor * (temp - ref_temp);
//    strain_plus_thermal(1) = strain(1) - temp_strain_material_factor * (temp - ref_temp);
//    strain_plus_thermal(2) = strain(2);
//    
//    // now multiply to calculate the stresses
//    material_factor.right_multiply_vector(strain_plus_thermal, stress);
//    
//    
//    // add the two tensors
//    tensor.zero();
//    tensor(0,0) = strain(0);  // xx
//    tensor(1,1) = strain(1);  // yy
//    tensor(0,1) = strain(2);  // xy
//    tensor(1,0) = strain(2);  // yx
//    return_qty->addStrainTensor(tensor, *load_case_it);
//    
//    tensor.zero();
//    tensor(0,0) = stress(0);  // xx
//    tensor(1,1) = stress(1);  // yy
//    tensor(0,1) = stress(2);  // xy
//    tensor(1,0) = stress(2);  // yx
//    return_qty->addStressTensor(tensor, *load_case_it);
//    
//    tensor.zero();
//    tensor(0,0) = strain_plus_thermal(0);  // xx
//    tensor(1,1) = strain_plus_thermal(1);  // yy
//    tensor(0,1) = strain_plus_thermal(2);  // xy
//    tensor(1,0) = strain_plus_thermal(2);  // yx
//    return_qty->addMechanicalStrainTensor(tensor, *load_case_it);
//    
//    unsigned int dv_ID = 0;
//
//    // now calculate the sensitivities
//    dv_it = dv_begin;
//    for (; dv_it != dv_end; dv_it++)
//      {
//      dv_ID = (*dv_it)->getID();
//      this->clearSensitivityInitialization();
//      
//      // get the dof value sensitivity for this case
//      this->analysis_discipline.getElemDofValues(this->base_elem,
//                                                 dof_value_sens,
//                                                 *load_case_it,
//                                                 true,
//                                                 dv_ID);
//      
//      // get the load sensitivity vectors
//      this->getNodalTemperatureVector(nodal_temp_sens, *load_case_it, true, dv_ID);
//      
//      material_factor_sens.zero();
//      temp_strain_material_factor_sens = 0.0;
//      
//      switch ((*dv_it)->getParameterTypeEnumID())
//        {
//        case PROPERTY_PARAMETER_ENUM_ID:
//          {
//            this->reinitForPropertySensitivity(dv_ID);
//            
//            // for a property DV, the strain operator sensitivity will be zero
//            B_strain_sens.zero();
//            
//            // if the property ID is the same as this elements property, then set the 
//            // value of the property sensitivity. Else the value be zero
//            if ( this->elem_property_card->checkGlobalParameterDependence(dv_ID))
//              {
//              this->elem_property_card->getFactorSensitivityForGlobalParameter
//              (material_factor_sens, STRESS_STRAIN_FACTOR::num(), dv_ID);
//              this->elem_property_card->
//                getPropertyValueDerivativeForGlobalParameterFromMaterialCard
//                (ALPHA_EXPANSION::num(), dv_ID, temp_strain_material_factor_sens);
//              this->elem_property_card->
//                getPropertyValueDerivativeForGlobalParameterFromMaterialCard
//                (THICKNESS_2D_ELEM::num(), dv_ID, thickness_sens);
//              }
//            else 
//              {
//              material_factor_sens.zero();
//              
//              temp_strain_material_factor_sens = 0.0;
//              }
//            
//            // also, calculate the sensitivity of the dof vector. The transformation matrix sensitivity 
//            // will be zero for this case, since shape sensitivity will be zero
//            transform_mat.right_multiply_vector(dof_value_sens, local_dof_sens);
//          }
//          break;
//          
//        case SHAPE_PARAMETER_ENUM_ID:
//          {
//            Elem* pert_elem = NULL;
//            pert_elem = this->analysis_discipline.getPerturbedElemForShapeParameter(this->elem_ID,
//                                                                                    *dv_it);
//            this->reinitForShapeSensitivity(dv_ID, pert_elem, 
//                                            (*dv_it)->getPerturbationStepSize());
//            
//            // set the factor_sens to zero
//            material_factor_sens.zero();
//            temp_strain_material_factor_sens = 0.0;
//            
//            // get the transformation matrix sensitivity and calculate the local_dof_sens
//            this->getStructuralT_matrix(&transform_mat_sens, true);
//            
//            local_dof_sens.zero(); scratch_vec.zero();
//            transform_mat_sens.right_multiply_vector(dof_values, local_dof_sens);
//            
//            scratch_vec.zero();
//            transform_mat.right_multiply_vector(dof_value_sens, scratch_vec);
//            
//            local_dof_sens.add(1.0, scratch_vec);
//            
//            // now init the fe_sens to the perturbed DV and calculate the shape sensitivity of the
//            // B matrix
//            // then, reinit it at the base elem, at the element centroid, where the 
//            // element strains and stresses will be evaluated
//            this->attachNewElem(this->perturbed_elem);
//            this->initialize_element(&point_vec);
//
//            // get the necessary data out of this element 
//            const std::vector< std::vector<Real> >& dphi_dx_perturbed = 
//              this->fe_base_map[LAGRANGE]->get_dphidx();
//            const std::vector< std::vector<Real> >& dphi_dy_perturbed = 
//              this->fe_base_map[LAGRANGE]->get_dphidy();
//            
//            dphi_dx_sens.zero();
//            dphi_dy_sens.zero();
//            
//            for (unsigned int i=0; i<n_nodes; i++)
//              {
//              dphi_dx_sens(i) = dphi_dx_perturbed[i][0] - dphi_dx[i][0];
//              dphi_dy_sens(i) = dphi_dy_perturbed[i][0] - dphi_dy[i][0];
//              }
//            
//            dphi_dx_sens.scale(1.0/this->perturbation);
//            dphi_dy_sens.scale(1.0/this->perturbation);
//            
//            // calculate B_strain matrix. Since the calculation is being performed at a single point, 
//            // the vector of shape factor derivatives should have size = 1
//            
//            for (unsigned int i=0; i < n_nodes; i++)
//              {
//              B_strain_sens(0,i) = dphi_dx_sens(i);  // epsilon_xx
//              B_strain_sens(0,4*n_nodes+i) = 0.5 * (thickness * dphi_dx_sens(i) + 
//                                                    thickness_sens * dphi_dx[i][0]); // bending contribution
//              B_strain_sens(1,n_nodes+i) = dphi_dy_sens(i);  // epsilon_yy
//              B_strain_sens(1,3*n_nodes+i) = - 0.5 * (thickness * dphi_dy_sens(i) + 
//                                                      thickness_sens * dphi_dy[i][0]); // bending contribution
//              B_strain_sens(2,i) = dphi_dy_sens(i);  // epsilon  xy 
//              B_strain_sens(2,3*n_nodes+i) = -0.5 * (thickness * dphi_dx_sens(i) + 
//                                                     thickness_sens * dphi_dx[i][0]); // bending contribution
//              B_strain_sens(2,n_nodes+i) = dphi_dx_sens(i);  // epsilon xy
//              B_strain_sens(2,4*n_nodes+i) = 0.5 * (thickness * dphi_dy_sens(i) + 
//                                                    thickness_sens * dphi_dy[i][0]); // bending contribution
//              }
//            
//          }
//          break;
//          
//        default:
//          abort();
//          break;
//        }
//      
//      // calculate the strain and stress
//      // calculate the factor and the strain and stress
//      strain_scratch.zero();
//      B_strain.right_multiply_vector(local_dof_sens, strain_sens);
//      B_strain_sens.right_multiply_vector(local_dof, strain_scratch);
//      strain_sens.add(1.0, strain_scratch);
//      
//      // add the thermal strain component and then calculate the stress
//      temp_sens = 0.0;
//      for (unsigned int i=0; i<phi.size(); i++)
//        temp_sens += phi[i][0] * nodal_temp_sens(i);
//      
//      // add the temperture sensitivity effects to the strain sensitivity
//      strain_plus_thermal_sens.zero();
//      strain_plus_thermal_sens(0) = strain_sens(0) - temp_strain_material_factor * temp_sens -
//        temp_strain_material_factor_sens * (temp - ref_temp);
//      strain_plus_thermal_sens(1) = strain_sens(1) - temp_strain_material_factor * temp_sens -
//        temp_strain_material_factor_sens * (temp - ref_temp);
//      strain_plus_thermal_sens(2) = strain_sens(2);
//      
//      // now multiply to calculate the stresses
//      strain_scratch.zero();
//      material_factor_sens.right_multiply_vector(strain_plus_thermal, stress_sens);
//      material_factor.right_multiply_vector(strain_plus_thermal_sens, strain_scratch);
//      stress_sens.add(1.0, strain_scratch);
//      
//      // add the two tensors
//      tensor.zero();
//      tensor(0,0) = strain_sens(0);  // xx
//      tensor(1,1) = strain_sens(1);  // yy
//      tensor(0,1) = strain_sens(2);  // xy
//      tensor(1,0) = strain_sens(2);  // yx
//      return_qty->addStrainTensor(tensor, *load_case_it,  dv_ID);
//      
//      tensor.zero();
//      tensor(0,0) = stress_sens(0);  // xx
//      tensor(1,1) = stress_sens(1);  // yy
//      tensor(0,1) = stress_sens(2);  // xy
//      tensor(1,0) = stress_sens(2);  // yx
//      return_qty->addStressTensor(tensor, *load_case_it, dv_ID);
//      
//      tensor.zero();
//      tensor(0,0) = strain_plus_thermal_sens(0);  // xx
//      tensor(1,1) = strain_plus_thermal_sens(1);  // yy
//      tensor(0,1) = strain_plus_thermal_sens(2);  // xy
//      tensor(1,0) = strain_plus_thermal_sens(2);  // yx
//      return_qty->addMechanicalStrainTensor(tensor, *load_case_it,  dv_ID);
//      
//      }
//    }  
//  
//  return return_qty;
}




void 
FESystemElem::PlateMITC4::calculateMITC4ElemQty(DenseMatrix<double>* qty, 
                                                const unsigned int name, 
                                                const unsigned int design_point_enum_ID,
                                                const unsigned int domain)
{
  // unused parameter
  (void) domain;

  assert (qty != NULL);
	
  switch(name)
    {
    case PLATE_MITC4_SHEAR_STIFFNESS_FACTOR_ENUM_ID:
      this->calculateShearStiffnessFactor(qty, design_point_enum_ID);
      break;
      
    default:
      abort();
      break;
    }
}




