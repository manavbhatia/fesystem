// $Id: RadiationElement.C,v 1.10.4.1 2007-03-04 03:37:12 manav Exp $



// C++ includes


// FESystem includes
#include "RadiationElement.h"
#include "RadiationCavity.h"
#include "FESystem/FESystemExceptions.h"
#include "Numerics/PetscSeqVector.h"

// libMesh includes
#include "node.h"
#include "elem.h"
#include "face_quad4.h"
#include "cell_hex8.h"
#include "fe.h"
#include "quadrature_gauss.h"

RadiationElement::RadiationElement(RadiationCavityFiniteElem& rad_cav_elem,
				   const Elem* fe_elem):
  fe_rad_elem(rad_cav_elem),
  fe_elem(fe_elem),
  elem_ID(0),
  x_division_number(0),
  y_division_number(0),
  geometric_rad_elem(NULL),
  local_elem(NULL),
  area(0.0),
  initialized(false)
{

}


RadiationElement::~RadiationElement()
{
  // iterate over the nodes and delete them
  std::vector<Node*>::iterator it, end;
  it = this->global_coordinate_nodes.begin();
  end = this->global_coordinate_nodes.end();

  for (; it != end; it++)
    delete *it;
  
  it = this->local_coordinate_nodes.begin();
  end = this->local_coordinate_nodes.end();
  for (; it != end; it++)
    delete *it;

  it = this->local_nodes.begin();
  end = this->local_nodes.end();
  for (; it != end; it++)
    delete *it;

}


void RadiationElement::init()
{
 
  // get the finite element data from the cavity
  unsigned int  n_x_divs = 0, n_y_divs = 0;
      
  n_x_divs = this->fe_rad_elem.n_x_divs;
  n_y_divs = this->fe_rad_elem.n_y_divs;
      
  // calculate the local coordinates of the local elem
  double dx=2.0/n_x_divs, dy=2.0/n_y_divs;
  static Point point;
  static Elem* fe_elem_ptr = NULL;
  // first point 
      
  std::vector<Point> point_vec;
     
  switch (this->fe_elem->type())
    {
    case TRI3:
      {
        // for now, no divisions are supported for triangle elements
        Assert (n_x_divs == 1, ExcInternalError());
        Assert (n_y_divs == 1, ExcInternalError());
        
        /// first point
        point.zero();
        point(0) = 0.0;
        point(1) = 0.0;
        point_vec.push_back(point);
        this->local_coordinate_nodes.push_back(new Node(point,0));
        
        // second point
        point.zero();
        point(0) = 1.0;
        point(1) = 0.0;
        point_vec.push_back(point);
        this->local_coordinate_nodes.push_back(new Node(point,1));
        
        // third point
        point.zero();
        point(0) = 0.0;
        point(1) = 1.0;
        point_vec.push_back(point);
        this->local_coordinate_nodes.push_back(new Node(point,2));
        
        fe_elem_ptr = const_cast<Elem*>(this->fe_elem);
      }
    break;
      
    case QUAD4:
      {
        /// first point
        point.zero();
        point(0) = -1.0 + dx * this->x_division_number;
        point(1) = -1.0 + dy * this->y_division_number;
        point_vec.push_back(point);
        this->local_coordinate_nodes.push_back(new Node(point,0));
        
        // second point
        point.zero();
        point(0) = -1.0 + dx * (this->x_division_number + 1);
        point(1) = -1.0 + dy * (this->y_division_number);
        point_vec.push_back(point);
        this->local_coordinate_nodes.push_back(new Node(point,1));
        
        // third point
        point.zero();
        point(0) = -1.0 + dx * (this->x_division_number + 1);
        point(1) = -1.0 + dy * (this->y_division_number + 1);
        point_vec.push_back(point);
        this->local_coordinate_nodes.push_back(new Node(point,2));
        
        // fourth point
        point.zero();
        point(0) = -1.0 + dx * (this->x_division_number);
        point(1) = -1.0 + dy * (this->y_division_number + 1);
        point_vec.push_back(point);
        this->local_coordinate_nodes.push_back(new Node(point,3));
       
        fe_elem_ptr = const_cast<Elem*>(this->fe_elem);
      }
      break;
      
    default:
      Assert(false, ExcInternalError());
    }


      
  // now, init a finite element with the FE geometric elem, 
  // and init it as the points calculated above
  std::auto_ptr<FEBase> fe(new FE<2,LAGRANGE>(FEType()));

//  else if (this->fe_elem->type() == HEX8)
//    {
//      elem_for_init.reset(this->fe_elem->build_side(this->fe_rad_elem.side_num).release());
//      fe_elem_ptr = elem_for_init.get();
//    }

  fe->reinit(fe_elem_ptr, &point_vec);

  // get the phi map
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
      
  // now iterate over the points and calculate the nodal values, 
  // create the nodes and add them to the nodal vector
  for (unsigned int i=0; i< this->fe_elem->n_nodes(); i++)
    {
      point.zero();
      for (unsigned int j=0; j < this->fe_elem->n_nodes(); j++)
        {
        const Node& elem_node = *(this->fe_elem->get_node(j));
        for (unsigned int k = 0; k < 3; k++)
          point(k) += phi[j][i] * elem_node(k);
        }
	  
      std::auto_ptr<Node> node(new Node(point, i));
      this->global_coordinate_nodes.push_back(node.release());
    }
      
  // now create an elem
  this->geometric_rad_elem.reset(Elem::build(this->fe_elem->type()).release());
  for (unsigned int i=0; i < this->fe_elem->n_nodes(); i++)
    this->geometric_rad_elem->set_node(i) = 
      this->global_coordinate_nodes[i];

  // now calculate the area and surface normal of the element 
  // get the location of the four nodes of this element
  const Point& point0 = this->geometric_rad_elem->point(0);
  const Point& point1 = this->geometric_rad_elem->point(1);
  const Point& point2 = this->geometric_rad_elem->point(2);
		
  //calculate the position vectors from point 0
  static Point vector02, vector01, tmp_vector1;
  vector02.assign(point2 - point0);
  vector01.assign(point1 - point0);
		
  // calculate the cross product and scale it by the side (for +ve or -ve)
  tmp_vector1.assign(vector01.cross(vector02));
  
  this->area = this->geometric_rad_elem->volume();
  
  this->normal.assign(tmp_vector1.unit());

  // now, based on the sign of the side, scale the normal to indicate a +ve or -ve
  // element normal
  if (fe_rad_elem.side_num < 0)
  this->normal *= -1.0;
  

  this->initLocalRadiationElement();
  
  this->initialized = true;

}




Elem* RadiationElement::getRadiationGeometricElem()
{
  // if the elem has not been initialized, init it and 
  // then return the elem
  if (!this->initialized)
    this->init();
    
  return this->geometric_rad_elem.get();
}



void RadiationElement::printInfo(std::ostream& output)
{
  // if the element is not initialized, initialize it and then print 
  // the info
  if (!this->initialized)
    this->init();

  output << "Radiation Element ID: " << this->elem_ID << std::endl
	 << "        FE Elem: " <<  this->fe_rad_elem.elem_ID 
	 << "  with Face: " << this->fe_rad_elem.side_num << std::endl
	 << "        Division number in FE Elem: x = " << this->x_division_number 
	 << "  y = " << this->y_division_number << std::endl
	 << "        Area = " << this->area << std::endl
	 << "        Surface Normal = " << this->normal << std::endl
	 << "        Node Locations: " << std::endl;
  
  // now print node locations
  std::vector<Node*>::const_iterator it, end;
  it = this->global_coordinate_nodes.begin();
  end = this->global_coordinate_nodes.end();

  for ( ; it != end; it++)
    output << (*(*it));
    output << std::endl;

}



void RadiationElement::initLocalRadiationElement()
{
  // this method will initialize the Elem for this radiation element in its 
  // own local coordinates. This coordinate is different from the 
  // global coordinate or the local coordinate system in the FE elem
  //get the four nodes of this elem
  switch (this->fe_elem->type())
    {
    case TRI3:
      {
        //get the three nodes of this elem
        const Point& point0 = this->geometric_rad_elem->point(0);
        const Point& point1 = this->geometric_rad_elem->point(1);
        const Point& point2 = this->geometric_rad_elem->point(2);
        
        //calculate the position vectors from point 0
        static Point vector01, vector02, tmp_vector1, tmp_vector2, 
          new_unit_vec_x, new_unit_vec_y;
        
        vector01.assign(point1 - point0);
        vector02.assign(point2 - point0);
        
        new_unit_vec_x = vector01.unit();
        
        tmp_vector1.assign(new_unit_vec_x.cross(vector02));
        tmp_vector2 = tmp_vector1.unit();
        
        new_unit_vec_y = tmp_vector2.cross(new_unit_vec_x);
        
        //create a new element and set its two nodes at the
        //origin and the elem length
        static Point new_point0, new_point1, new_point2;
        new_point0.zero(); new_point1.zero(); new_point2.zero();
        new_point1(0)  = vector01.size();
        new_point2(0) = vector02 * new_unit_vec_x;
        new_point2(1) = vector02 * new_unit_vec_y;
        
        this->local_nodes.push_back(new Node(new_point0, 0));
        this->local_nodes.push_back(new Node(new_point1, 1));
        this->local_nodes.push_back(new Node(new_point2, 2));
      }			
      break;
			
    case QUAD4:
      {
        //get the four nodes of this elem
        const Point& point0 = this->geometric_rad_elem->point(0);
        const Point& point1 = this->geometric_rad_elem->point(1);
        const Point& point2 = this->geometric_rad_elem->point(2);
        const Point& point3 = this->geometric_rad_elem->point(3);
        
        //calculate the position vectors from point 0
        static Point vector01, vector02, vector03, tmp_vector1, tmp_vector2,
          new_unit_vec_x, new_unit_vec_y;
        vector01.assign(point1 - point0);
        vector02.assign(point2 - point0);
        vector03.assign(point3 - point0);
        
        new_unit_vec_x = vector01.unit();
        
        tmp_vector1.assign(new_unit_vec_x.cross(vector03));
        tmp_vector2 = tmp_vector1.unit();
        
        new_unit_vec_y = tmp_vector2.cross(new_unit_vec_x);
        
        //create a new element and set its two nodes at the
        //origin and the elem length
        static Point new_point0, new_point1, new_point2, new_point3;
        new_point0.zero(); new_point1.zero(); new_point2.zero(); new_point3.zero();
        
        new_point1(0) = vector01.size();
        
        new_point2(0) = vector02 * new_unit_vec_x; 
        new_point2(1) = vector02 * new_unit_vec_y; 
        
        new_point3(0) = vector03 * new_unit_vec_x; 
        new_point3(1) = vector03 * new_unit_vec_y;
        
        this->local_nodes.push_back(new Node(new_point0, 0));
        this->local_nodes.push_back(new Node(new_point1, 1));
        this->local_nodes.push_back(new Node(new_point2, 2));
        this->local_nodes.push_back(new Node(new_point3, 3));
      }
      break;
			
    default:
      // cases for TRI6, QUAD8, QUAD9 are not handled
      abort();
      break;
    }
  				
		
  // now create the local element and set its nodes
  this->local_elem.reset(Elem::build(this->fe_elem->type()).release());
  for (unsigned int i=0; i<this->fe_elem->n_nodes(); i++)
    this->local_elem->set_node(i) = this->local_nodes[i];
}




void RadiationElement::getRadElemG_Vector(FESystemNumerics::PetscSeqVector<double>& vector)
{
  // for now, this method will handle only the 2-D elems
  assert (this->fe_elem->dim() == 2);
  
  unsigned int n_nodes = this->fe_elem->n_nodes();

  // resize the vector if necessary
  if (vector.size() != n_nodes)
    vector.resize(n_nodes);
    
  // now process and insert the data
  // init a gauss quadrature and calculate the local coordinates for this elem
  std::auto_ptr<QBase> quadrature(new QGauss(2,FIFTH));
  std::auto_ptr<FEBase> fe_rad(new FE<2, LAGRANGE>(FEType()));

  // to calculate the JxW, an element in its local coordinate needs to be 
  // used, since in general, an element can exist in a 3-D space, and this will 
  // create problems for 3-D elems
  // 
  
  fe_rad->attach_quadrature_rule(quadrature.get());
  fe_rad->reinit(local_elem.get());

  const std::vector<std::vector<Real> >& phi_rad = fe_rad->get_phi();
  const std::vector<Real>& JxW_rad = fe_rad->get_JxW();

  std::vector<Point> points;

  for (unsigned int qp = 0; qp < phi_rad[0].size(); qp++)
    {
      Point pt;
      for (unsigned int i=0; i < phi_rad.size(); i++)
	{
	  Node& node = *(this->local_coordinate_nodes[i]);
	  pt(0) += phi_rad[i][qp] * node(0);
	  pt(1) += phi_rad[i][qp] * node(1);
	}
      points.push_back(pt);
    }

  // the point localtions calculated above are the nodal values of the
  // quad points in the local coordinates of the FE to which this 
  // element is attached
  std::auto_ptr<FEBase> global_fe(new FE<2,LAGRANGE>(FEType()));
  global_fe->reinit(this->fe_elem, &points);

  const std::vector<std::vector<Real> >& phi_global = global_fe->get_phi();
  
  // now, calculate the G_vector
  vector.zero();
  for (unsigned int qp=0; qp < JxW_rad.size(); qp++)
    for (unsigned int i=0; i < phi_global.size(); i++)
      vector.add(i, phi_global[i][qp] * JxW_rad[qp]);
}
