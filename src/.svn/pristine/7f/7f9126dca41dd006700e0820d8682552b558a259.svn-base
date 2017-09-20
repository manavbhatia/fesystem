// $Id: PlateDKQ.C,v 1.2 2006-11-13 10:13:42 manav Exp $

// C++ includes

// FESystem includes
#include "PlateDKQ.h"
#include "ElemPostProcessQty.h"

// libMesh includes
#include "face_quad4.h"





PlateDKQ::PlateDKQ(const unsigned int elem_ID,
		       const Elem* elem, 
		       AnalysisDriver* analysis_driver):
  StructuralElem(elem_ID,elem,analysis_driver)
{
  assert(elem->type() == QUAD4);
}





PlateDKQ::~PlateDKQ()
{
	
}






bool PlateDKQ::calculate_M(DenseMatrix<double>* matrix, 
			     bool sensitivity_calculation)
{
  abort();
}







bool PlateDKQ::calculate_K(DenseMatrix<double>* matrix, 
			     bool sensitivity_calculation)
{
  Property::PropertyName prpty_name = Property::INVALID_PROPERTY_NAME;
	
  if (sensitivity_calculation == true)
    {
      switch (this->DV->type())
	{
	case PROPERTY_DV:
	  prpty_name = (dynamic_cast<PropertyDesignVariable*>(this->DV))->propertyName();
	  break;
				
	case SHAPE_DV:
	  break;
				
	default:
	  abort();
	  break;
	}
    }
	
	
  PropertyCard& elem_property = *(this->_elem_property_card);
	
  double E = elem_property[Property::E_11];
  double h = elem_property[Property::THICKNESS];
  double nu = elem_property[Property::NU];
	
  const double k = 5.0/6.0;
	
  double Eh_nu2, Eh_nu, Ehnu_nu2, Eh3_12nu2, Eh3_12nu, Eh3nu_12nu2, Ehk_2nu;
	
  switch (prpty_name)
    {
    case Property::INVALID_PROPERTY_NAME:
      Eh_nu2 = E*h/(1.0-nu*nu);
      Eh_nu = E*h/(1.0+nu);
      Ehnu_nu2 = E*nu*h/(1.0-nu*nu);
      Eh3_12nu2 = E*h*h*h/12.0/(1.0-nu*nu);
      Eh3_12nu = E*h*h*h/12.0/(1.0+nu);
      Eh3nu_12nu2 = E*nu*h*h*h/12.0/(1.0-nu*nu);
      Ehk_2nu = E*h*k/2.0/(1.0+nu);
      break;
			
    case Property::E_11:
      Eh_nu2 = h/(1.0-nu*nu);
      Eh_nu = h/(1.0+nu);
      Ehnu_nu2 = nu*h/(1.0-nu*nu);
      Eh3_12nu2 = h*h*h/12.0/(1.0-nu*nu);
      Eh3_12nu = h*h*h/12.0/(1.0+nu);
      Eh3nu_12nu2 = nu*h*h*h/12.0/(1.0-nu*nu);
      Ehk_2nu = h*k/2.0/(1.0+nu);
      break;
			
    case Property::THICKNESS:
      Eh_nu2 = E/(1.0-nu*nu);
      Eh_nu = E/(1.0+nu);
      Ehnu_nu2 = E*nu/(1.0-nu*nu);
      Eh3_12nu2 = E*h*h/4.0/(1.0-nu*nu);
      Eh3_12nu = E*h*h/4.0/(1.0+nu);
      Eh3nu_12nu2 = E*nu*h*h/4.0/(1.0-nu*nu);
      Ehk_2nu = E*k/2.0/(1.0+nu);
      break;
			
    case Property::NU:
      Eh_nu2 = 2.0*E*h*nu/(1.0-nu*nu)/(1.0-nu*nu);
      Eh_nu = -E*h/(1.0+nu)/(1.0+nu);
      Ehnu_nu2 = E*h*(1.0+nu*nu)/(1.0-nu*nu)/(1.0-nu*nu);
      Eh3_12nu = E*h*h*h*nu/6.0/(1.0-nu*nu)/(1.0-nu*nu);
      Eh3_12nu2 = -E*h*h*h/12.0/(1.0+nu)/(1.0+nu);
      Eh3nu_12nu2 = E*h*h*h*(1.0+nu*nu)/12.0/(1.0-nu*nu)*(1.0-nu*nu);
      Ehk_2nu = -E*h*k/2.0/(1.0+nu)/(1.0+nu);
      break;
			
    default:
      Eh_nu2 = 0.0;
      Eh_nu = 0.0;
      Ehnu_nu2 = 0.0;
      Eh3_12nu = 0.0;
      Eh3_12nu2 = 0.0;
      Eh3nu_12nu2 = 0.0;
      Ehk_2nu = 0.0;
      break;
    }
	
	
  unsigned int n_nodes = this->base_elem->n_nodes();
	
  std::auto_ptr<DenseMatrix<double> > shear_stiff(new DenseMatrix<double>());
  this->resizeQty(shear_stiff.get());
	
  DenseMatrix<double> Nx_Nx(n_nodes,n_nodes),
    Ny_Ny(n_nodes,n_nodes),
    Nx_Ny(n_nodes,n_nodes);
	
	
  if (sensitivity_calculation == false || ( sensitivity_calculation == true && this->DV->type() == PROPERTY_DV))
    {
      this->getFactor(&Nx_Nx, StructuralElem::N_X_N_X_FACTOR);
      this->getFactor(&Ny_Ny, StructuralElem::N_Y_N_Y_FACTOR);
      this->getFactor(&Nx_Ny, StructuralElem::N_X_N_Y_FACTOR);
      this->getFactor(shear_stiff.get(), StructuralElem::PLATE_QUAD4_SHEAR_STIFF_FACTOR);
    }
  else if (sensitivity_calculation == true && this->DV->type() == SHAPE_DV)
    {
      this->getFactorShapeSensitivity(&Nx_Nx, StructuralElem::N_X_N_X_FACTOR);
      this->getFactorShapeSensitivity(&Ny_Ny, StructuralElem::N_Y_N_Y_FACTOR);
      this->getFactorShapeSensitivity(&Nx_Ny, StructuralElem::N_X_N_Y_FACTOR);
      this->getFactorShapeSensitivity(shear_stiff.get(), StructuralElem::PLATE_QUAD4_SHEAR_STIFF_FACTOR);
    }

	
  // calculate the stiffness matrix
  for (unsigned int i=0; i<n_nodes; i++)
    for (unsigned int j=0; j<n_nodes; j++)
      {
	(*matrix)(i,j) += (Eh_nu2*Nx_Nx(i,j) + .5*Eh_nu* Ny_Ny(i,j));
	(*matrix)(n_nodes+i,j) += (Ehnu_nu2 * Nx_Ny(j,i) + .5*Eh_nu * Nx_Ny(i,j));
	(*matrix)(i,n_nodes+j) += (Ehnu_nu2 * Nx_Ny(i,j) + .5*Eh_nu * Nx_Ny(j,i));
	(*matrix)(n_nodes+i,n_nodes+j) += (Eh_nu2*Ny_Ny(i,j) + .5*Eh_nu* Nx_Nx(i,j));
	
	(*matrix)(3*n_nodes+i,3*n_nodes+j) += (Eh3_12nu2*Ny_Ny(i,j) + .5*Eh3_12nu* Nx_Nx(i,j));
	(*matrix)(4*n_nodes+i,3*n_nodes+j) += -(Eh3nu_12nu2*Nx_Ny(i,j) + .5*Eh3_12nu* Nx_Ny(j,i));
	(*matrix)(3*n_nodes+i,4*n_nodes+j) += -(Eh3nu_12nu2*Nx_Ny(j,i) + .5*Eh3_12nu* Nx_Ny(i,j));
	(*matrix)(4*n_nodes+i,4*n_nodes+j) += (Eh3_12nu2*Nx_Nx(i,j) + .5*Eh3_12nu* Ny_Ny(i,j));
      }
	
  // now add the effects of the shear strains
  shear_stiff->scale(Ehk_2nu);
  matrix->add(1.0, *(shear_stiff.get()));
		
  //  std::cout << *(shear_stiff.get()) << std::endl;

  // finally, place a unity at the diagonal for the theta_z rows, to remove singularity 
  // of the element
  // this will be changed in future, since the value of this fictitious stiffness is dependent on the 
  // value of the other diagonal terms 
  for (unsigned int i=0; i<4; i++)
    (*matrix)(5*n_nodes+i,5*n_nodes+i) = 1.0;

  return true;
}





void PlateDKQ::calculateShearStiffnessFactor(DenseMatrix<double>* matrix)
{
  // initialize the element if it is not already initialized
  if (this->local_elem_is_initialized == false)
    this->initialize_element();

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
  double Ax, Bx, Cx, Ay, By, Cy, sin_alpha, cos_alpha, sin_beta, cos_beta;
	
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

  const Point& new_point0 = this->_local_elem->point(2);
  const Point& new_point1 = this->_local_elem->point(3);
  const Point& new_point2 = this->_local_elem->point(0);
  const Point& new_point3 = this->_local_elem->point(1);
  
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

  
	
  std::vector<Real>& JxW = *(this->_JxW);
  unsigned int n_nodes = this->base_elem->n_nodes();

  DenseVector<double> gamma_rz(12), gamma_sz(12), gamma_xz(12), gamma_yz(12);
  double xi, eta, ABx_factor, ABy_factor, CBx_factor, CBy_factor, det_J;
  matrix->zero();

  for (unsigned int qp=0; qp<this->_qbase_2d->n_points(); qp++)
    {
      gamma_rz.zero(); gamma_sz.zero(); gamma_xz.zero(); gamma_yz.zero();
      
      xi = (this->_qbase_2d->get_points())[qp](0);
      eta = (this->_qbase_2d->get_points())[qp](1);
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




// bool PlateDKQ::calculate_F_Pressure(DenseVector<double>* vector, 
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



bool PlateDKQ::calculate_F_T(DenseVector<double>* vector, 
			       bool sensitivity_calculation)
{
  PropertyCard& elem_property = *(this->_elem_property_card);
  Property::PropertyName prpty_name = Property::INVALID_PROPERTY_NAME;
  	
  unsigned int load_case_num = this->analysis_driver->getCurrentLoadCase();
  unsigned int n_nodes = this->base_elem->n_nodes();

  // initialize the data structures for the load
  DenseVector<double> nodal_temp(n_nodes), nodal_temp_sens(n_nodes);

  // get the loads for this element
  this->getNodalTemperatureVector(nodal_temp, load_case_num);
  
  // if sensitivity is being asked for, then get the iterators for the
  // sensitivity load vector
  if (sensitivity_calculation == true)
    {
      switch (this->DV->type())
	{
	case PROPERTY_DV:
	  prpty_name = (dynamic_cast<PropertyDesignVariable*>(this->DV))->propertyName();
	  break;
  				
	case SHAPE_DV:
	  break;
  				
	default:
	  abort();
	  break;
	}
  		
      // get the load sensitivity vectors
      this->getNodalTemperatureVector(nodal_temp_sens, load_case_num, true, this->DV->ID());
    }

  	
  double E = elem_property[Property::E_11];
  double h = elem_property[Property::THICKNESS];
  double nu = elem_property[Property::NU];
  double alpha = elem_property[Property::ALPHA_EXPANSION];
  double ref_temp = this->analysis_discipline.getFESystemController().
    analysis_case->getRealParameter("REFERENCE_TEMP");
  	
  double factor = 0.0, factor_sens = 0.0;
  	
  factor = E*h*alpha/(1.0-nu);
 
  switch (prpty_name)
    {
    case Property::E_11:
      factor_sens = h*alpha/(1.0-nu);
      break;
  			
    case Property::THICKNESS:
      factor_sens = E*alpha/(1.0-nu);
      break;
  			
    case Property::NU:
      factor_sens = E*h*alpha/pow((1.0-nu),2);
      break;
  			
    case Property::ALPHA_EXPANSION:
      factor_sens = E*h/(1.0-nu);
      break;
  
    default:
      factor_sens = 0.0;
      break;
    }
  
  // get the factors necessary for calculation of the load vector
  
  
  if (sensitivity_calculation == true)
    {
      // the reference temperature is assumed to be constant
      switch(this->DV->type())
	{
	case PROPERTY_DV:
	  {
	    DenseMatrix<double> nx_n_matrix(n_nodes,n_nodes), 
	      ny_n_matrix(n_nodes,n_nodes);

	    this->getFactor(&nx_n_matrix, StructuralElem::N_X_N_FACTOR);
	    this->getFactor(&ny_n_matrix, StructuralElem::N_Y_N_FACTOR);

	    double sum1 = 0.0, sum2 = 0.0;
	    for (unsigned int i=0; i<n_nodes; i++)
	      {
		sum1 = 0.0; sum2 = 0.0;
		for (unsigned int j=0; j<n_nodes; j++)
		  {
		    sum1 += nx_n_matrix(i,j) * (factor* (nodal_temp(j)-ref_temp) + factor_sens * nodal_temp_sens(j));
		    sum2 += ny_n_matrix(i,j) * (factor* (nodal_temp(j)-ref_temp) + factor_sens * nodal_temp_sens(j));
		  }
		(*vector)(i) = sum1;
		(*vector)(n_nodes + i) = sum2;
	      }
	  }
	  break;

	case SHAPE_DV:
	  {	  
	    DenseMatrix<double> nx_n_matrix(n_nodes,n_nodes), 
	      ny_n_matrix(n_nodes,n_nodes), 
	      nx_n_matrix_sens(n_nodes,n_nodes), 
	      ny_n_matrix_sens(n_nodes,n_nodes);

	    this->getFactor(&nx_n_matrix, StructuralElem::N_X_N_FACTOR);
	    this->getFactor(&ny_n_matrix, StructuralElem::N_Y_N_FACTOR);
	    this->getFactorShapeSensitivity(&nx_n_matrix_sens, StructuralElem::N_X_N_FACTOR);
	    this->getFactorShapeSensitivity(&ny_n_matrix_sens, StructuralElem::N_Y_N_FACTOR);
	    
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
	  abort();
	  break;
	}
    }
  else 
    {
      DenseMatrix<double> nx_n_matrix(n_nodes,n_nodes), 
	ny_n_matrix(n_nodes,n_nodes);

      this->getFactor(&nx_n_matrix, StructuralElem::N_X_N_FACTOR);
      this->getFactor(&ny_n_matrix, StructuralElem::N_Y_N_FACTOR);

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
  return true;
}


std::auto_ptr<ElemPostProcessQty> 
PlateDKQ::getElementPostProcessQty(std::vector<unsigned int> load_cases, 
				     std::vector<DesignData::DesignParameter*> dv_vector)
{
  std::auto_ptr<ElemPostProcessQty> return_qty(new ElemPostProcessQty(this->_elem_ID));

  TensorValue<double> tensor;

  // create an fe that will be used for calculationi of the element strains and
  // stresses. The vector fe_sens is for shape sensitivity
  std::auto_ptr<FEBase> fe(new FE<2,LAGRANGE>(FEType())),
    fe_sens (new FE<2,LAGRANGE>(FEType()));

  double temp_strain_material_factor = 0.0, temp_strain_material_factor_sens = 0.0,
    ref_temp = 0.0, E_11 =0.0, nu = 0.0, thickness = 0.0;
 

  // initialize the material factors
  PropertyCard& property_card = *(this->_elem_property_card);
  temp_strain_material_factor = property_card[Property::ALPHA_EXPANSION];
  ref_temp = property_card[Property::TEMP_REF];
  E_11 = property_card[Property::E_11];
  nu = property_card[Property::NU];
  thickness = property_card[Property::THICKNESS];


  unsigned int n_nodes = this->base_elem->n_nodes();

  DenseMatrix<double> B_strain(3, 6*n_nodes) , B_strain_sens(3, 6*n_nodes);

  // create rest of the necessary data structures
  DenseVector<double>  strain(3), strain_sens(3),
    stress(3), stress_sens(3), strain_plus_thermal(3), 
    strain_plus_thermal_sens(3), strain_scratch(3);
  
  DenseMatrix<double> material_factor(3,3), material_factor_sens(3,3);


  
  // point at which this element will be initialized
  // for now, this is only the centroid of the element
  std::vector<Point> point_vec;
  point_vec.push_back(Point(0.0,0.0,0.0));

  // then, reinit it at the base elem, at the element centroid, where the 
  // element strains and stresses will be evaluated
  fe->reinit(this->base_elem, &point_vec);

  // get the necessary data out of this element 
  const std::vector< std::vector<Real> >& dphi_dx = fe->get_dphidx();
  const std::vector< std::vector<Real> >& dphi_dy = fe->get_dphidy();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  // calculate B_strain matrix. Since the calculation is being performed at a single point, 
  // the vector of shape factor derivatives should have size = 1
  assert (dphi_dx[0].size() == 1);
  for (unsigned int i=0; i < n_nodes; i++)
    {
      B_strain(0,i) = dphi_dx[i][0];  // epsilon_xx
      B_strain(0,4*n_nodes+i) = 0.5 * thickness * dphi_dx[i][0]; // bending contribution
      B_strain(1,n_nodes+i) = dphi_dy[i][0];  // epsilon_yy
      B_strain(1,3*n_nodes+i) = - 0.5 * thickness * dphi_dy[i][0]; // bending contribution
      B_strain(2,i) = dphi_dy[i][0];  // epsilon  xy 
      B_strain(2,3*n_nodes+i) = - 0.5 * thickness * dphi_dx[i][0]; // bending contribution
      B_strain(2,n_nodes+i) = dphi_dx[i][0];  // epsilon xy
      B_strain(2,4*n_nodes+i) = 0.5 * thickness * dphi_dy[i][0]; // bending contribution
    }
  

  // iterate over all the load cases, to calculate the strains and stresses
  std::vector<unsigned int>::const_iterator load_case_it = load_cases.begin();
  std::vector<unsigned int>::const_iterator load_case_end = load_cases.end();	
	
  std::vector<DesignData::DesignParameter*>::const_iterator dv_it, dv_begin, dv_end;
  dv_begin = dv_vector.begin();
  dv_end = dv_vector.end();
    	
  
  unsigned int property_card_ID = 0;
  Property::PropertyName prpty_name;

  std::auto_ptr<DenseMatrix<double> > transform_mat(new DenseMatrix<double>), 
    transform_mat_sens(new DenseMatrix<double>);
  this->resizeQty(transform_mat.get());
  this->resizeQty(transform_mat_sens.get());

  this->getFactor(transform_mat.get(), StructuralElem::TRANSFORM_MATRIX);


  // also, create dof_vectors
  DenseVector<double>  dof_values(6*n_nodes), dof_value_sens(6*n_nodes),
    local_dof(6*n_nodes), local_dof_sens(6*n_nodes), scratch_vec(6*n_nodes),
    nodal_temp(n_nodes),nodal_temp_sens(n_nodes),
    dphi_dx_sens(n_nodes), dphi_dy_sens(n_nodes);
  

  // set the material factor matrix
  material_factor(0,0) = E_11 / (1.0 - nu * nu);
  material_factor(0,1) = E_11 * nu / (1.0 - nu * nu);
  material_factor(1,0) = E_11 * nu / (1.0 - nu * nu);
  material_factor(1,1) = E_11 / (1.0 - nu * nu);
  material_factor(2,2) = E_11 * 0.5 / (1.0 + nu);
  
  double temp = 0.0, temp_sens = 0.0;

  // iterate over each load case, and ask solver to solve for it
  for (; load_case_it != load_case_end; load_case_it++)
    {
      
      // get the DOF vector for this elem
      dof_values = *(this->getElemDofValuesForLoadCase(*load_case_it).get());

      // transform these dofs to the local coordinate system
      local_dof.zero();
      transform_mat->right_multiply_vector(dof_values, local_dof);

      // get the loads for this element
      this->getNodalTemperatureVector(nodal_temp, *load_case_it);

      // calculate the factor and the strain and stress
      B_strain.right_multiply_vector(local_dof, strain);

      // add the thermal strain component and then calculate the stress
      for (unsigned int i=0; i<phi.size(); i++)
	temp += phi[i][0] * nodal_temp(i);

      strain_plus_thermal.zero();
      strain_plus_thermal(0) = strain(0) - temp_strain_material_factor * (temp - ref_temp);
      strain_plus_thermal(1) = strain(1) - temp_strain_material_factor * (temp - ref_temp);
      strain_plus_thermal(2) = strain(2);

      // now multiply to calculate the stresses
      material_factor.right_multiply_vector(strain_plus_thermal, stress);


      // add the two tensors
      tensor.zero();
      tensor(0,0) = strain(0);  // xx
      tensor(1,1) = strain(1);  // yy
      tensor(0,1) = strain(2);  // xy
      tensor(1,0) = strain(2);  // yx
      return_qty->addStrainTensor(tensor, *load_case_it);
      
      tensor.zero();
      tensor(0,0) = stress(0);  // xx
      tensor(1,1) = stress(1);  // yy
      tensor(0,1) = stress(2);  // xy
      tensor(1,0) = stress(2);  // yx
      return_qty->addStressTensor(tensor, *load_case_it);

      tensor.zero();
      tensor(0,0) = strain_plus_thermal(0);  // xx
      tensor(1,1) = strain_plus_thermal(1);  // yy
      tensor(0,1) = strain_plus_thermal(2);  // xy
      tensor(1,0) = strain_plus_thermal(2);  // yx
      return_qty->addMechanicalStrainTensor(tensor, *load_case_it);

      // now calculate the sensitivities
      dv_it = dv_begin;
      for (; dv_it != dv_end; dv_it++)
	{
	  // get the dof value sensitivity for this case
	  dof_value_sens = *(this->getElemDofValueSensitivityForLoadCaseAndDV(*load_case_it, (*dv_it)->ID()).get());
	  
	  // get the load sensitivity vectors
	  this->getNodalTemperatureVector(nodal_temp_sens, *load_case_it, true, (*dv_it)->ID());

	  material_factor_sens.zero();
	  temp_strain_material_factor_sens = 0.0;
	  
	  // set the DV for this element 
	  this->DV = *dv_it;

	  switch ((*dv_it)->type())
	    {
	    case PROPERTY_DV:
	      {
		// for a property DV, the strain operator sensitivity will be zero
		B_strain_sens.zero();

		PropertyDesignVariable* prpty_DV = dynamic_cast<PropertyDesignVariable*>(*dv_it);

		// get the property ID and name
		property_card_ID = prpty_DV->propertyCardID();
		prpty_name = prpty_DV->propertyName();
		
		// if the property ID is the same as this elements property, then set the 
		// value of the property sensitivity. Else the value be zero
		if ( property_card_ID == this->_elem_property_card->id())
		  {
		    switch (prpty_name)
		      {
		      case Property::E_11:
			{
			  material_factor_sens(0,0) = 1.0 / (1.0- nu * nu);
			  material_factor_sens(0,1) = nu / (1.0- nu * nu);
			  material_factor_sens(1,0) =  nu / (1.0- nu * nu);
			  material_factor_sens(1,1) = 1.0 / (1.0- nu * nu);
			  material_factor_sens(2,2) =  0.5 / (1 + nu);

			  temp_strain_material_factor_sens = 0.0;
			}
			break;
								

		      case Property::NU:
			{
			  material_factor_sens(0,0) = 2.0 * nu * E_11 / pow((1.0- nu * nu),2);
			  material_factor_sens(0,1) = E_11 * (1.0+ nu*nu) / pow((1.0 - nu * nu),2);
			  material_factor_sens(1,0) = E_11 * (1.0+ nu*nu) / pow((1.0 - nu * nu),2);
			  material_factor_sens(1,1) = 2.0 * nu * E_11 / pow((1.0- nu * nu),2);
			  material_factor_sens(2,2) = - E_11 * 0.5 / pow((1 + nu),2);

			  temp_strain_material_factor_sens = 0.0;
			}
			break;


		      case Property::ALPHA_EXPANSION:
			{
			  material_factor_sens.zero();

			  temp_strain_material_factor_sens = 1.0;
			}
			break;

		      case Property::THICKNESS:
			{
			  material_factor_sens.zero();

			  temp_strain_material_factor_sens = 0.0;
			  
			  // set the value of the 
			  for (unsigned int i=0; i < n_nodes; i++)
			    {
			      B_strain_sens(0,4*n_nodes+i) = 0.5 * dphi_dx[i][0]; // bending contribution
			      B_strain_sens(1,3*n_nodes+i) = - 0.5 * dphi_dy[i][0]; // bending contribution
			      B_strain_sens(2,3*n_nodes+i) = -0.5 * dphi_dx[i][0]; // bending contribution
			      B_strain_sens(2,4*n_nodes+i) = 0.5 * dphi_dy[i][0]; // bending contribution
			    }

			}
			break;


		      default:
			{
			  material_factor_sens.zero();

			  temp_strain_material_factor_sens = 0.0;
			}
			break;
		      }
		  }
		else 
		  {
		    material_factor_sens.zero();

		    temp_strain_material_factor_sens = 0.0;
		  }

		// also, calculate the sensitivity of the dof vector. The transformation matrix sensitivity 
		// will be zero for this case, since shape sensitivity will be zero
		transform_mat->right_multiply_vector(dof_value_sens, local_dof_sens);
	      }
	      break;
				
	    case SHAPE_DV:
	      {
		// calculate the length for the perturbed element, and find the sensitivity
		// calculate the length of the element
		// first init the perturbed elem
		this->initPerturbedElemForShapeDV();

		// set the factor_sens to zero
		material_factor_sens.zero();

		temp_strain_material_factor_sens = 0.0;

		// get the transformation matrix sensitivity and calculate the local_dof_sens
		this->getFactorShapeSensitivity(transform_mat_sens.get(), StructuralElem::TRANSFORM_MATRIX);
		
		local_dof_sens.zero(); scratch_vec.zero();
		transform_mat_sens->right_multiply_vector(dof_values, local_dof_sens);
		
		scratch_vec.zero();
		transform_mat->right_multiply_vector(dof_value_sens, scratch_vec);
		
		local_dof_sens.add(1.0, scratch_vec);

		
		ShapeDesignVariable* shape_DV = dynamic_cast<ShapeDesignVariable*>(this->DV);
		double perturbation = shape_DV->perturbation();
	  				
		// now init the fe_sens to the perturbed DV and calculate the shape sensitivity of the
		// B matrix
		// then, reinit it at the base elem, at the element centroid, where the 
		// element strains and stresses will be evaluated
		fe_sens->reinit(this->perturbed_elem, &point_vec);

		// get the necessary data out of this element 
		const std::vector< std::vector<Real> >& dphi_dx_perturbed = fe_sens->get_dphidx();
		const std::vector< std::vector<Real> >& dphi_dy_perturbed = fe_sens->get_dphidy();
		
		dphi_dx_sens.zero();
		dphi_dy_sens.zero();
		
		for (unsigned int i=0; i<n_nodes; i++)
		  {
		    dphi_dx_sens(i) = dphi_dx_perturbed[i][0] - dphi_dx[i][0];
		    dphi_dy_sens(i) = dphi_dy_perturbed[i][0] - dphi_dy[i][0];
		  }

		dphi_dx_sens.scale(1.0/perturbation);
		dphi_dy_sens.scale(1.0/perturbation);

		// calculate B_strain matrix. Since the calculation is being performed at a single point, 
		// the vector of shape factor derivatives should have size = 1

		for (unsigned int i=0; i < n_nodes; i++)
		  {
		    B_strain_sens(0,i) = dphi_dx_sens(i);  // epsilon_xx
		    B_strain_sens(0,4*n_nodes+i) = 0.5 * thickness * dphi_dx_sens(i); // bending contribution
		    B_strain_sens(1,n_nodes+i) = dphi_dy_sens(i);  // epsilon_yy
		    B_strain_sens(1,3*n_nodes+i) = - 0.5 * thickness * dphi_dy_sens(i); // bending contribution
		    B_strain_sens(2,i) = dphi_dy_sens(i);  // epsilon  xy 
		    B_strain_sens(2,3*n_nodes+i) = -0.5 * thickness * dphi_dx_sens(i); // bending contribution
		    B_strain_sens(2,n_nodes+i) = dphi_dx_sens(i);  // epsilon xy
		    B_strain_sens(2,4*n_nodes+i) = 0.5 * thickness * dphi_dy_sens(i); // bending contribution
		  }
		
	      }
	      break;
	      
	    default:
	      abort();
	      break;
	    }
	  
	  // calculate the strain and stress
	  // calculate the factor and the strain and stress
	  strain_scratch.zero();
	  B_strain.right_multiply_vector(local_dof_sens, strain_sens);
	  B_strain_sens.right_multiply_vector(local_dof, strain_scratch);
	  strain_sens.add(1.0, strain_scratch);

	  // add the thermal strain component and then calculate the stress
	  temp_sens = 0.0;
	  for (unsigned int i=0; i<phi.size(); i++)
	    temp_sens += phi[i][0] * nodal_temp_sens(i);
	  
	  // add the temperture sensitivity effects to the strain sensitivity
	  strain_plus_thermal_sens.zero();
	  strain_plus_thermal_sens(0) = strain_sens(0) - temp_strain_material_factor * temp_sens -
	    temp_strain_material_factor_sens * (temp - ref_temp);
	  strain_plus_thermal_sens(1) = strain_sens(1) - temp_strain_material_factor * temp_sens -
	    temp_strain_material_factor_sens * (temp - ref_temp);
	  strain_plus_thermal_sens(2) = strain_sens(2);

	  // now multiply to calculate the stresses
	  strain_scratch.zero();
	  material_factor_sens.right_multiply_vector(strain_plus_thermal, stress_sens);
	  material_factor.right_multiply_vector(strain_plus_thermal_sens, strain_scratch);
	  stress_sens.add(1.0, strain_scratch);

	  // add the two tensors
	  tensor.zero();
	  tensor(0,0) = strain_sens(0);  // xx
	  tensor(1,1) = strain_sens(1);  // yy
	  tensor(0,1) = strain_sens(2);  // xy
	  tensor(1,0) = strain_sens(2);  // yx
	  return_qty->addStrainTensor(tensor, *load_case_it,  (*dv_it)->ID());
      
	  tensor.zero();
	  tensor(0,0) = stress_sens(0);  // xx
	  tensor(1,1) = stress_sens(1);  // yy
	  tensor(0,1) = stress_sens(2);  // xy
	  tensor(1,0) = stress_sens(2);  // yx
	  return_qty->addStressTensor(tensor, *load_case_it, (*dv_it)->ID());

	  tensor.zero();
	  tensor(0,0) = strain_plus_thermal_sens(0);  // xx
	  tensor(1,1) = strain_plus_thermal_sens(1);  // yy
	  tensor(0,1) = strain_plus_thermal_sens(2);  // xy
	  tensor(1,0) = strain_plus_thermal_sens(2);  // yx
	  return_qty->addMechanicalStrainTensor(tensor, *load_case_it,  (*dv_it)->ID());

	}
    }  

  return return_qty;
}

