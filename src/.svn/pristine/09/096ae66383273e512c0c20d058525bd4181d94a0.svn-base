// $Id: RadiationElementPair.C,v 1.14.6.2 2008-04-06 04:06:28 manav Exp $

// C++ includes


// FESystem includes
#include "Radiation/RadiationElement.h"
#include "Radiation/RadiationElementPair.h"
#include "FESystem/FESystemController.h"
#include "Radiation/ContourIntegration.h"
#include "Radiation/MitalasIntegration.h"
 
// Shape factor includes
extern "C" {
#include "Radiation/ShapeFactors/ff.h"
#include "Radiation/ShapeFactors/claussen.h"
}

// libmesh includes
#include "geom/face_quad4.h"
#include "geom/point.h"
#include "fe/fe.h"
#include "quadrature/quadrature_gauss.h"
#include "quadrature/quadrature_grid.h"

#ifdef WITH_MATHEMATICA
double MathematicaFormFactor(RadiationElement* elem1, RadiationElement* elem2, MLINK mathlink);
static void read_and_print_expression( MLINK lp);
static void print_packet_type(int type);
#endif


RadiationElementPair::RadiationElementPair():
initialized(false),
rad_elem1(NULL),
rad_elem2(NULL)
#ifdef WITH_MATHEMATICA
,mathematica_link(NULL)
#endif
{
}





RadiationElementPair::~RadiationElementPair()
{

}



std::auto_ptr<Elem>
RadiationElementPair::createLocalElem(Elem* elem, 
                                      std::vector<Node>& nodes)
{
		
  //get the four nodes of this elem
  const Point& point0 = elem->point(0);
  const Point& point1 = elem->point(1);
  const Point& point2 = elem->point(2);
  const Point& point3 = elem->point(3);
		
  //calculate the position vectors from point 0
  Point vector01, vector02, vector03, tmp_vector1, tmp_vector2;
  vector01.assign(point1 - point0);
  vector02.assign(point2 - point0);
  vector03.assign(point3 - point0);
		
  Point new_unit_vec_x = vector01.unit();
		
  tmp_vector1.assign(new_unit_vec_x.cross(vector03));
  tmp_vector2 = tmp_vector1.unit();
		
  Point new_unit_vec_y = tmp_vector2.cross(new_unit_vec_x);
		
  //create a new element and set its two nodes at the
  //origin and the elem length
  Point new_point0(0. ,0. ,0. );
  Point new_point1(vector01.size() ,0. ,0.);
  Point new_point2(vector02 * new_unit_vec_x, 
		   vector02 * new_unit_vec_y, 
		   0.);
  Point new_point3(vector03 * new_unit_vec_x, 
		   vector03 * new_unit_vec_y, 
		   0.);
		
  nodes.push_back(Node(new_point0, 0));
  nodes.push_back(Node(new_point1, 1));
  nodes.push_back(Node(new_point2, 2));
  nodes.push_back(Node(new_point3, 3));
		
  // now create the local element and set its nodes
  std::auto_ptr<Elem> local_elem(new Quad4);
  local_elem->set_node(0) = &(nodes[0]);
  local_elem->set_node(1) = &(nodes[1]);
  local_elem->set_node(2) = &(nodes[2]);
  local_elem->set_node(3) = &(nodes[3]);
  
  return local_elem;
}



std::pair<double, double>
RadiationElementPair::getShapeFactor(ShapeFactorMethod method)
{
  switch (method)
    {
    case DOUBLE_AREA_INT:
      {
	return this->doubleAreaIntShapeFactor();
      }
      break;
      
      
    case ANALYTIC:
      {
	return this->analyticShapeFactor();
      }
      break;
      

    case CONTOUR_INT:
      {
        return this->contourIntShapeFactor();
      }
      break;
      
    case MITALAS_CONTOUR_INT:
      {
        return this->mitalasContourIntShapeFactor();
      }
      break;

    default:
      abort();
    }
}



std::pair<double, double> 
RadiationElementPair::analyticShapeFactor()
{
  // get the geometric element for this element pair, and call the librarry function 
  // to get the shape factors
  Elem* elem1 = this->rad_elem1->getRadiationGeometricElem();
  Elem* elem2 = this->rad_elem2->getRadiationGeometricElem();

  double factor_12 = 0.0, factor_21 = 0.0;
  const unsigned int n_nodes1 = elem1->n_nodes();
  const unsigned int n_nodes2 = elem2->n_nodes();
  
  double * node1_coords[3];
  double * node2_coords[3];
  // create the matrices for passing to the FormFactor routine
//  double *x_coords1 = new double[n_nodes1], *y_coords1 = new double[n_nodes1],
//  *z_coords1 = new double[n_nodes1];
//  double *x_coords2 = new double[n_nodes2], *y_coords2 = new double[n_nodes2],
//  *z_coords2 = new double[n_nodes2];

//  for (unsigned int i=0; i < n_nodes1; i++)
//    std::cout << elem1->point(i);
//  
//  std::cout << " \n \n";
//  
//  for (unsigned int i=0; i < n_nodes2; i++)
//    std::cout << elem2->point(i);

  for (unsigned int j=0; j < 3; j++)
    {
      node1_coords[j] = new double[n_nodes1];
      node2_coords[j] = new double[n_nodes2];
      for (unsigned int i=0; i < n_nodes1; i++)
        {
          if (this->rad_elem1->FeElementFaceID() < 0)
            node1_coords[j][n_nodes1-1-i] = (elem1->point(i))(j);
          else
            node1_coords[j][i] = (elem1->point(i))(j);
        }

      for (unsigned int i=0; i < n_nodes2; i++)
        {
          if (this->rad_elem2->FeElementFaceID() < 0)
            node1_coords[j][n_nodes2-1-i] = (elem2->point(i))(j);
          else
            node1_coords[j][i] = (elem2->point(i))(j);
        }
    }
  
//  for (unsigned int i=0; i < n_nodes1; i++)
//    for (unsigned int j=0; j < 3; j++)
//      {
//      if (this->rad_elem1->FeElementFaceID() < 0)
//        elem1_points[((n_nodes1-1)-i)*3+j] = (elem1->point(i))(j);
//      else
//        elem1_points[i*3+j] = (elem1->point(i))(j);
//      }  
	  
//  for (unsigned int i=0; i < n_nodes2; i++)
//    for (unsigned int j=0; j < 3; j++)
//      {
//      if (this->rad_elem2->FeElementFaceID() < 0)
//        elem2_points[((n_nodes2-1)-i)*3+j] = (elem2->point(i))(j);
//      else
//        elem2_points[i*3+j] = (elem2->point(i))(j);
//      }

  Assert(false, ExcInternalError());
  
//      FormFactor( node1_coords, n_nodes1, node2_coords, n_nodes2, 
//                  &factor_12, &factor_21);

#ifdef WITH_MATHEMATICA
  // get the factor from mathematica
  double factor = MathematicaFormFactor(this->rad_elem1, this->rad_elem2, this->mathematica_link);
  //  factor_12 = factor / this->rad_elem1->getArea();
//  factor_21 = factor / this->rad_elem2->getArea();
  std::cout << factor_12 <<  "  " << factor / this->rad_elem1->getArea() 
    << "  " << factor_21 << "  " << factor / this->rad_elem2->getArea() << std::endl;
  //  std::cout << factor_12 <<  "  " << factor << std::endl; 
#endif
  
  if (factor_12 < 0)
    factor_12 = 0.0;
  if (factor_21 < 0)
    factor_21 = 0.0;


  // if the factor is 1.0, then certain checks need to be made, since the
  // analytic method is unable to catch the error of surfaces lying on top 
  // of each other, but with opposite orientation, in which case, the factor
  // is 0.
  if (fabs(factor_12 -1.0) <= 1.0e-10)
    {
      const Point& normal1 = this->rad_elem1->getSurfaceNormal(), 
	normal2 = this->rad_elem2->getSurfaceNormal();
      
      if (normal1 * normal2 < 0.0)
	{
	  factor_12 = 0.0;
	  
	  // also, check the other shape factor for the pair.
	  if (fabs(factor_21 -1.0) <= 1.0e-10)
		factor_21 = 0.0;
	}
    }

  for (unsigned int j=0; j < 3; j++)
    {
      delete[] node1_coords[j];
      delete[] node2_coords[j];
    }      

  return std::pair<double, double>(factor_12, factor_21);
}





std::pair<double, double> 
RadiationElementPair::doubleAreaIntShapeFactor()
{
  bool initialized = false;

  FEType fetype(FIRST, LAGRANGE);
//  std::auto_ptr<Elem> local_elem1(NULL), local_elem2(NULL);
//  std::vector<Node> local_elem1_nodes, local_elem2_nodes;
  std::auto_ptr<FEBase > 
    fe_elem1(FEBase::build(2,fetype).release()), 
    fe_elem2(FEBase::build(2,fetype).release()); 
  std::auto_ptr<QBase> 
    qbase_elem1(QBase::build(QGAUSS, 2, ELEVENTH).release()), 
    qbase_elem2(QBase::build(QGAUSS, 2, ELEVENTH).release());
  
  // get the elems from the rad elems
  Elem* elem1 = this->rad_elem1->getRadiationGeometricElem();
  Elem* elem2 = this->rad_elem2->getRadiationGeometricElem();

  Point& normal1 = this->rad_elem1->getSurfaceNormal();
  Point& normal2 = this->rad_elem2->getSurfaceNormal();

  // since the cavity is assumed to be convex, if the normals are parallel, 
  // and if the dot product of the normals is not negative, then they lie in
  // the same plane, and their shape factor should be zero.
  if (fabs(normal1 * normal2 - 1.0) <= FESystemNumbers::Epsilon)
    return std::pair<double, double>(0.0, 0.0);
  
  // now, initialize the elements
	if (!initialized)
    {
    fe_elem1->attach_quadrature_rule(qbase_elem1.get());
    fe_elem2->attach_quadrature_rule(qbase_elem2.get());
    initialized = true;
    }

  // this has been commented out, since the local element is now being 
  // obtained from the radiation element
//  // now create local elements
//  local_elem1.reset(this->createLocalElem(elem1,
//					  local_elem1_nodes).release());
//  
//  local_elem2.reset(this->createLocalElem(elem2,
//					  local_elem2_nodes).release());
  
  fe_elem1->reinit(this->rad_elem1->getLocalRadiationGeometricElem());
  fe_elem2->reinit(this->rad_elem2->getLocalRadiationGeometricElem());

  double factor = 0.0, area1 = 0.0, area2 = 0.0;
  double length = 0.0, cos_theta1 = 0.0, cos_theta2 = 0.0;
  Point point1, point2, point12;

  area1 = this->rad_elem1->getArea();
  area2 = this->rad_elem2->getArea();

  const std::vector<Real>& JxW1 = fe_elem1->get_JxW();
  const std::vector<std::vector<Real> >& phi1 = fe_elem1->get_phi();
  const std::vector<Real>& JxW2 = fe_elem2->get_JxW();
  const std::vector<std::vector<Real> >& phi2 = fe_elem2->get_phi();
  
  
  // get the location of the nodes of the two elements
  // this is available from the local_elem_node pointers
  
  // iterate  over the quadrature points of both the elements
  for (unsigned int qp1=0; qp1<qbase_elem1->n_points(); qp1++)
    {
      point1.zero();
		
    for (unsigned int node_it=0; node_it<phi1.size(); node_it++)
      for (unsigned int coord_it=0; coord_it<3; coord_it++)
        point1(coord_it) += 
	      phi1[node_it][qp1] * (elem1->point(node_it))(coord_it);
				
    for (unsigned int qp2=0; qp2<qbase_elem2->n_points(); qp2++)
      {			
      point2.zero();
      point12.zero(); length = 0.0; cos_theta1 = 0.0; cos_theta2= 0.0;
      // get location of the quadrature points 1 and 2 using the nodal location
      // and shape function values at the quadrature points
      for (unsigned int node_it=0; node_it<phi2.size(); node_it++)
        for (unsigned int coord_it=0; coord_it<3; coord_it++)
          point2(coord_it) += 
          phi2[node_it][qp2] * (elem2->point(node_it))(coord_it);
          
          
      // calculate the vector from point 1 to 2
      point12.assign(point2);
      point12 -= point1;
      length = point12.size();
      
      // calculate the value of cos(theta) at point 1 and 2
      if (length >  1.0e-12)
        {
        cos_theta1 = fabs((normal1) * point12) / length;
        cos_theta2 = fabs((normal2) * point12) / length;
        
        factor += (cos_theta1 * cos_theta2)/(FESystemNumbers::Pi* length * length) * 
          JxW1[qp1] * JxW2[qp2];
        }
      }
    }
    
    if (factor < 0.0)
      factor = 0.0;
    
    return std::pair<double, double>(factor/area1, factor/area2);
}






std::pair<double, double> 
RadiationElementPair::contourIntShapeFactor()
{
  static double factor;
  static unsigned int n_sides1, n_sides2, j1, j2;
  factor = 0.0;
    
  Point& normal1 = this->rad_elem1->getSurfaceNormal();
  Point& normal2 = this->rad_elem2->getSurfaceNormal();
  
  // since the cavity is assumed to be convex, if the normals are parallel, 
  // and if the dot product of the normals is not negative, then they lie in
  // the same plane, and their shape factor should be zero.
  if (fabs(normal1 * normal2 - 1.0) <= FESystemNumbers::Epsilon)
    return std::pair<double, double>(0.0, 0.0);

  // get the elems from the rad elems
  const Elem* elem1 = this->rad_elem1->getRadiationGeometricElem();
  const Elem* elem2 = this->rad_elem2->getRadiationGeometricElem();

  n_sides1 = elem1->n_sides();
  n_sides2 = elem2->n_sides();

  for (unsigned int i=0; i < n_sides1; i++)
    {
    // get the nodes for this side
    if (i == (n_sides1-1))
      j1 = 0;
    else 
      j1 = i + 1;

    const Point& p1 = elem1->point(i);
    const Point& q1 = elem1->point(j1);

    for (unsigned int j=0; j < elem2->n_sides(); j++)
      {
      if (j == (n_sides2-1))
        j2 = 0;
      else
        j2 = j + 1;

      const Point& p2 = elem2->point(j);
      const Point& q2 = elem2->point(j2);
      
      factor += EdgeFactor(p1, q1, p2, q2);
      }
    }
  
  factor /= (2.0 * FESystemNumbers::Pi);
  
  double area1 = this->rad_elem1->getArea(),
    area2 = this->rad_elem2->getArea();
  
  if (factor < 0.0)
    factor = 0.0;
  
  return std::pair<double, double>(factor/area1, factor/area2);
}






std::pair<double, double> 
RadiationElementPair::mitalasContourIntShapeFactor()
{
  static double factor;
  static unsigned int n_sides1, n_sides2, j1, j2;
  factor = 0.0;

//  double val1, val2;
  
  Point& normal1 = this->rad_elem1->getSurfaceNormal();
  Point& normal2 = this->rad_elem2->getSurfaceNormal();
  
  // since the cavity is assumed to be convex, if the normals are parallel, 
  // and if the dot product of the normals is not negative, then they lie in
  // the same plane, and their shape factor should be zero.
  if (fabs(normal1 * normal2 - 1.0) <= FESystemNumbers::Epsilon)
    return std::pair<double, double>(0.0, 0.0);
  
  // get the elems from the rad elems
  const Elem* elem1 = this->rad_elem1->getRadiationGeometricElem();
  const Elem* elem2 = this->rad_elem2->getRadiationGeometricElem();
  
  n_sides1 = elem1->n_sides();
  n_sides2 = elem2->n_sides();
  
//  for (unsigned int i=0; i < elem1->n_nodes(); i++)
//    std::cout << elem1->point(i);
//  std::cout << std::endl;
//  for (unsigned int i=0; i < elem2->n_nodes(); i++)
//    std::cout << elem2->point(i);
//  std::cout << std::endl;
  
  for (unsigned int i=0; i < n_sides1; i++)
    {
    // get the nodes for this side
    if (i == (n_sides1-1))
      j1 = 0;
    else 
      j1 = i + 1;
    
    const Point& p1 = elem1->point(i);
    const Point& q1 = elem1->point(j1);
    
    for (unsigned int j=0; j < n_sides2; j++)
      {
      if (j == (n_sides2-1))
        j2 = 0;
      else
        j2 = j + 1;
      
      const Point& p2 = elem2->point(j);
      const Point& q2 = elem2->point(j2);
      
//      val1 = MitalasEdgeFactor(p1, q1, p2, q2);
//      val2 = NumIntEdgeFactor(p1, q1, p2, q2);
//      if (fabs(val1 - val2) > 1.0e-5)
//        {
//        std::cout << p1 << q1;
//        std::cout << p2 << q2;
//        std::cout << "milatas : " << val1 << "  numeric: " <<  val2<< std::endl; 
//        }
  
      factor += MitalasEdgeFactor(p1, q1, p2, q2);
      }
    }

  if (this->rad_elem1->FeElementFaceID() < 0)
    factor *= -1.0;
  if (this->rad_elem2->FeElementFaceID() < 0)
    factor *= -1.0;
  
  factor /= (2.0 * FESystemNumbers::Pi);

  
  double area1 = this->rad_elem1->getArea(),
    area2 = this->rad_elem2->getArea();
  
//  std::cout << this->rad_elem1->FeElementID() << "," << 
//    this->rad_elem2->FeElementID() << " : " << 
//    factor/area1 << "," << factor/area2 << std::endl;

  //  std::cout << "factor : " << factor << std::endl;
//  std::cout << std::endl;
  
  if (factor < 0.0)
    factor = 0.0;
  
  return std::pair<double, double>(factor/area1, factor/area2); 
}






#ifdef WITH_MATHEMATICA
void
RadiationElementPair::setMathLink(MLINK mathlink)
{
  this->mathematica_link = mathlink;
}





double 
MathematicaFormFactor(RadiationElement* rad_elem1, RadiationElement* rad_elem2, MLINK mathlink)
{
  assert(rad_elem1 != NULL);
  assert(rad_elem2 != NULL);
  
  int pkt, n_nodes1, n_nodes2; long int lenp;
  double re, im;
  Elem* elem1 = rad_elem1->getRadiationGeometricElem();
  Elem* elem2 = rad_elem2->getRadiationGeometricElem();
  n_nodes1 = elem1->n_nodes();
  n_nodes2 = elem2->n_nodes();
  
    
  MLPutFunction( mathlink, "EvaluatePacket", 1);
  
  // build the function call for form factor
  MLPutFunction( mathlink, "FormFactor", 2);
  
  // now give coordinates of the first element
  MLPutFunction( mathlink, "List", n_nodes1);
  unsigned int lower_index, upper_index, increment;
  
  if (rad_elem1->FeElementFaceID() < 0)
    {
    lower_index = n_nodes1 -1;
    upper_index = -1;
    increment = -1;
    }
  else
    {
    lower_index = 0;
    upper_index = n_nodes1;
    increment = 1;
    }
  

  for (int i=lower_index; i != upper_index; i+=increment)
    {
    // list of nodal coordinates for each node
    MLPutFunction( mathlink, "List", 3);
    for (unsigned int j=0; j < 3; j++)
      MLPutReal(mathlink, elem1->point(i)(j));
    }

    // put the list for second elemen
    MLPutFunction( mathlink, "List", n_nodes2);

  if (rad_elem2->FeElementFaceID() < 0)
    {
    lower_index = n_nodes2 -1;
    upper_index = -1;
    increment = -1;
    }
  else
    {
    lower_index = 0;
    upper_index = n_nodes2;
    increment = 1;
    }

  for (int i=lower_index; i != upper_index; i+=increment)
      {
      // list of nodal coordinates for each node
      MLPutFunction( mathlink, "List", 3);
      for (unsigned int j=0; j < 3; j++)
        MLPutReal(mathlink, elem2->point(i)(j));
      }
  
//  MLPutFunction(mathlink, "ToExpression",1L);
//		MLPutString(mathlink, "FormFactor[{{6.24627, 4.5849, 0}, {6.44958, 4.91239 , 0}, {6.44958, 4.91239, 0.539143}},{{1.83183, 1.96496, 0.827657}, {2.13714,2.29245, 0}, {1.83183, 1.83396, -0.84048}}]");

      MLEndPacket(mathlink);
    MLFlush(mathlink);
          
    while( (pkt = MLNextPacket(mathlink), pkt) && pkt != RETURNPKT) {
//      print_packet_type(pkt);
//      read_and_print_expression( mathlink);
      MLNewPacket( mathlink);
      if( MLError( FESystem::mathlink_link))
        fprintf( stderr, "Error detected by MathLink: %s.\n",
                 MLErrorMessage(FESystem::mathlink_link));
    }
    
//    print_packet_type(pkt);
//    read_and_print_expression( mathlink);
    switch( MLGetNext( mathlink)) 
      {
      case MLTKREAL:
        {
          MLGetReal( mathlink, &re);
          im = 0.0;
//          printf( "%g \n", re);
        }
        break;
      case MLTKFUNC:
        {
          if (MLCheckFunction( mathlink, "Complex", &lenp)
              &&  lenp == 2
              &&  MLGetReal( mathlink, &re)
              &&  MLGetReal( mathlink, &im)
              )
            {
//            printf( "%g +i %g\n", re, im);
            }
          else
            {
            if( MLError( FESystem::mathlink_link))
              fprintf( stderr, "Error detected by MathLink: %s.\n",
                       MLErrorMessage(FESystem::mathlink_link));
            
            }  
        }
        break;
      default:
        {
          printf("did not find appropriate return from MathLink");
          abort();
        }
      }
    
    
    return re;
}

static void print_packet_type(int type)
{
  switch(type)
    {
    case RETURNPKT:
      printf("RETURNPKT\n");
      break;
    case RETURNTEXTPKT:
      printf("RETURNTEXTPKT\n");
      break;
    case INPUTNAMEPKT:
      printf("INPUTNAMEPKT\n");
      break;
    case OUTPUTNAMEPKT:
      printf("OUTPUTNAMEPKT\n");
      break;
    case TEXTPKT:
      printf("TEXTPKT\n");
      break;
    case MESSAGEPKT:
      printf("MESSAGEPKT\n");
      break;
    case DISPLAYPKT:
      printf("DISPLAYPKT\n");
      break;
    case DISPLAYENDPKT:
      printf("DISPLAYENDPKT\n");
      break;
    case INPUTPKT:
      printf("INPUTPKT\n");
      break;
    case CALLPKT:
      printf("CALLPKT\n");
      break;
    default:
      {
        printf("could not find anything appropriate\n");
        MLErrorMessage(FESystem::mathlink_link);
      }
    }  
}

static void read_and_print_expression( MLINK lp)
{
	kcharp_ct  s;
	int n;
	long i, len;
	double r;
	static int indent;
  
	switch( MLGetNext( lp)) {
    case MLTKSYM:
      {
        printf("symbol\n");
        MLGetSymbol( lp, &s);
        printf( "%s ", s);
        MLDisownSymbol( lp, s);
      }
      break;
    case MLTKSTR:
      {
        printf("string\n");
        MLGetString( lp, &s);
        printf( "\"%s\" ", s);
        MLDisownString( lp, s);
      }
      break;
    case MLTKINT:
      {
        printf("integer\n");
        MLGetInteger( lp, &n);
        printf( "%d ", n);
      }
      break;
    case MLTKREAL:
      {
        printf("real\n");
        MLGetReal( lp, &r);
        printf( "%g ", r);
      }
      break;
    case MLTKFUNC:
      {
        printf("function\n");
        indent += 3;
        printf( "\n %*.*s", indent, indent, "");
        if( MLGetArgCount( lp, &len) == 0){
          MLErrorMessage(FESystem::mathlink_link);
        }else{
          read_and_print_expression( lp);
          printf( "[");
          for( i = 1; i <= len; ++i){
            read_and_print_expression( lp);
            if( i != len) printf( ", ");
          }
          printf( "]");
        }
        indent -= 3;
      }
      break;
    case MLTKERROR:
    default:
      {
        printf("error\n");
        MLErrorMessage(FESystem::mathlink_link);
      }
	}
}
#endif

