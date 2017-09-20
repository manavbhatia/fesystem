#include<iostream>
#include<iomanip>
#include<fstream>
#include <cassert>
#include <math.h>

using std::cout;
using std::cin;
using std::endl;
using std::setw;
using std::ios;
int main()
{
  const int max_grid_pt = 25000;
  double grid_cord[max_grid_pt][3] = {{0.}}; //maximum grid points is 500
  int n_grid_x, n_grid_y, n_grid_total; // number of grid points 
  double x_side_len,y_side_len; // plate side lengths
  double dx, dy, Nx=0.0, Ny=0.0, Nxy=0.0, skew_angle; 
  double x_init = 0.,y_init= 0., temperature = 0.0;
  int bc_1=0, bc_2=0, bc_3=0, bc_4=0, use_geometric_stiffness=0; // number of grid points 
  
  cout << "\nEnter the side length along x- axis: ";
  cin >> x_side_len;
  
  cout << "\nEnter the side length along y- axis: ";
  cin >> y_side_len;

  cout << "\nSkew angle: ";
  cin >> skew_angle;
  
  cout << "\nGrid points along x- axis :";
  cin >> n_grid_x; 
  
  cout << "\nGrid points along y- axis :";
  cin >> n_grid_y;
  
  cout << "\n Temperature :";
  cin >> temperature;
  
  cout << "\n Use geometric stiffness (0/1) :";
  cin >> use_geometric_stiffness;
  
  // following statements calculate the data
  dx = x_side_len / double(n_grid_x -1);
  dy = y_side_len / double(n_grid_y -1);
  n_grid_total = n_grid_x * n_grid_y;
  assert(n_grid_total < max_grid_pt);
  
  int grid_number = 0;
  for (int i = 0; i < n_grid_y; i++ )
    {
      y_init = dy * double(i);
      for (int j = 0; j < n_grid_x; j++)
        { x_init = dx * double(j) + y_init  * sin(skew_angle * 3.14159 / 180.0);
	        grid_cord[grid_number][1] = x_init;
		grid_cord[grid_number][2] = y_init * cos(skew_angle * 3.14159 / 180.0);  	
		grid_number = grid_number + 1; 
        } 
    }
  
  
  unsigned int n_conds = 10;
  double density = 0.040, start_mach = 1.5, mach_incr = .1, speed_sound = 400.0;

  // write the fesystem input file
  {
    
    std::fstream output_file;
    output_file.open("fesys_input.txt", std::fstream::out);
    
    output_file << "FUNCTION_DATABASE BEGIN " << std::endl
    <<" N_FUNCTIONS 0 " << std::endl
    <<"FUNCTION_DATABASE END " << std::endl
    <<"PROPERTY_DATABASE BEGIN " << std::endl
    <<"MATERIAL_PROPERTY_CARDS BEGIN " << std::endl
    <<"N_MATERIAL_PROPERTY_CARDS 1 " << std::endl
    <<"ISOTROPIC_MATERIAL_PROPERTY_CARD BEGIN " << std::endl
    <<"ID 1 " << std::endl
    <<"N_PROPERTIES 4 " << std::endl
    <<"DENSITY 2000.0 " << std::endl
    <<"YOUNGS_MODULUS 27.0e+9 " << std::endl
    <<"POISSONS_RATIO 0.33 " << std::endl
    <<"ALPHA_EXPANSION 13.0e-6 " << std::endl
    <<"N_LOCAL_PARAMETERS 0 " << std::endl
    <<"N_GLOBAL_PARAMETERS 0 " << std::endl
    <<"N_FUNCTION_IDS 0 " << std::endl
    <<"ISOTROPIC_MATERIAL_PROPERTY_CARD END " << std::endl
    <<"MATERIAL_PROPERTY_CARDS END " << std::endl
    <<"ELEM_DATA_CARDS BEGIN " << std::endl
    <<"N_ELEM_DATA_CARDS 1 " << std::endl
    <<"ISOTROPIC_2D_ELEM_DATA_CARD BEGIN " << std::endl
    <<"ID 1 " << std::endl
    <<"MATERIAL_PROPERTY_CARD_ID 1 " << std::endl
    <<"N_PROPERTIES 1 " << std::endl
    <<"THICKNESS 0.0075 " << std::endl
    <<"N_LOCAL_PARAMETERS 0 " << std::endl
    <<"N_GLOBAL_PARAMETERS 0 " << std::endl
    <<"N_FUNCTION_IDS 0 " << std::endl
    <<"ISOTROPIC_2D_ELEM_DATA_CARD END " << std::endl
    <<"ELEM_DATA_CARDS END " << std::endl
    <<"PROPERTY_DATABASE END " << std::endl
    <<"MESH_LIST BEGIN " << std::endl
    <<"N_MESH 1 " << std::endl
    <<"MESH BEGIN " << std::endl
    <<"MESH_ID 1 " << std::endl
    << "DIMENSION 3" << std::endl
    <<"NODES BEGIN" << std::endl
    <<"N_NODES " << n_grid_total << std::endl;
    
    for ( unsigned int i=0; i < n_grid_total; i++)
      {
        output_file << i+1	<< "   " 
        << std::setprecision(7) << grid_cord[i][1]
        << "   " << std::setprecision(7) << grid_cord[i][2] 
        << "   0.0" << std::endl; 
      }
    
    output_file << "NODES END" << std::endl;
    int n_elem = (n_grid_y-1)*(n_grid_x-1) * 2;
    output_file << "ELEMENTS BEGIN" <<  std::endl
    << "N_ELEMENTS " << n_elem << std::endl;
    
    int n1 = 1;
    int elem_i = 1;
    std::string elem_name;
    if (use_geometric_stiffness)
      elem_name = "STRUCTURAL_TRI3_VON_KARMAN ";
    else
      elem_name = "STRUCTURAL_PLATE_DKT ";
    
    for (int i = 0; i < n_grid_y-1; i++ )
      {
        for (int j = 0; j < n_grid_x-1; j++)
          {
            output_file << elem_name <<setw(8)<<elem_i++<<
            setw(8)<<"1"<<setw(8)<<n1<<setw(8)<<n1 + 1<<
            setw(8)<<n1+1+n_grid_x<<setw(8)<<"\n"; 
            
            output_file << elem_name <<setw(8)<<elem_i++<<
            setw(8)<<"1"<<setw(8)<<n1 <<
            setw(8)<<n1+1+n_grid_x<<setw(8)<<n1+n_grid_x<<"\n"; 
            n1 = n1 + 1;
          }
        n1 = n1 + 1 ;
      }	
    
    output_file << "ELEMENTS END " << std::endl
    <<"MESH END " << std::endl
    <<"MESH_LIST END " << std::endl
    <<"RADIATION_CAVITY_LIST BEGIN " << std::endl
    <<"N_RADIATION_CAVITIES 0 " << std::endl
    <<"RADIATION_CAVITY_LIST END " << std::endl
    <<"ELEM_SET_LIST BEGIN " << std::endl
    <<"N_ELEM_SETS 0 " << std::endl
    <<"ELEM_SET_LIST END " << std::endl
    <<"LOAD_DATABASE BEGIN " << std::endl
    <<"LOAD_SET_DATA BEGIN " << std::endl
    <<"N_LOAD_SETS 5 " << std::endl;
    
    // piston theory surface
    output_file << "SURFACE_LOAD_SET BEGIN" << std::endl
    <<"NAME PISTON_SURFACES " << std::endl
    <<"ID 4" << std::endl
    <<"KIND PISTON_THEORY_SURFACE" << std::endl
    <<"N_LOADS " << n_elem << std::endl;
    {
      for (unsigned int i=1; i<= n_elem; i++)
        output_file << i << "  " << i << "  3  0.0" << std::endl;
    }
    output_file << "SURFACE_LOAD_SET END" << std::endl;
    
    // point load
    unsigned int n_nodal_force = 2*(n_grid_y + n_grid_x) ; 
    output_file << "NODAL_LOAD_SET BEGIN" << std::endl
    <<"NAME NODAL_FORCE " << std::endl
    <<"ID 2" << std::endl
    <<"KIND NODAL_FORCE" << std::endl
    <<"N_DOFS 6" << std::endl
    <<"N_LOADS " << n_nodal_force << std::endl;
    {
      unsigned int load_id=1;
      
      // write the Ny and Nxy loads
      double val1 = Ny * dx;
      double val2 = Nxy * dx;
      for (unsigned int i=0; i < n_grid_x; i++)
        {
          output_file << load_id << "  " << (i + 1) << "  " 
          << -val2 << "  " << -val1 <<  "  0.0 0.0 0.0 0.0" << std::endl;
          load_id++;
          output_file << load_id << "  " << ((n_grid_y-1)*n_grid_x + i+1)<< "  " 
          << val2 << "  " << val1 <<  "  0.0 0.0 0.0 0.0" << std::endl;
          load_id++;
        }
      
      // write the Nx and Nxy loads
      val1 = Nx * dy;
      val2 = Nxy * dy;
      for (unsigned int i=0; i < n_grid_y; i++)
        {
          output_file << load_id << "  " << ((i+1)* n_grid_x)<< "  " 
          << val1 << "  " << val2 << "  0.0 0.0 0.0 0.0" << std::endl;
          load_id++;
          output_file << load_id << "  " << (i* n_grid_x + 1)<< "  " 
          << -val1 << "  " << -val2 << "  0.0 0.0 0.0 0.0" << std::endl;
          load_id++;
        }
    }
    output_file <<"NODAL_LOAD_SET END" << std::endl;
    
    
    // displacement boundary conditions	
    output_file << "BOUNDARY_CONDITION_LOAD_SET BEGIN " << std::endl;
    output_file << "NAME DISPLACEMENT_BOUNDARY_CONDITIONS" << std::endl;
    output_file << "ID 5" << std::endl;
    output_file << "KIND DISPLACEMENT_BOUNDARY_CONDITION" << std::endl;
    output_file << "N_LOADS " << (n_grid_total+(n_grid_x+(n_grid_y-1))*6) << std::endl;
    unsigned int bc_id=1 ;
    for (unsigned int i=1; i<=n_grid_y; i++)
      for (unsigned int j=1; j<=n_grid_x; j++)
        {
          unsigned int n_node = (i-1)*n_grid_x + j ;
          output_file <<bc_id << "  " << n_node << "  6  0.0" << std::endl;
          bc_id++;
        }
    
    for (unsigned int j=1; j<=n_grid_x; j++)
      {
        output_file << bc_id << "  " << j << "  1  0.0" << std::endl; 
        bc_id++;
        output_file << bc_id << "  " << j << "  2  0.0" << std::endl; 
        bc_id++;
        output_file << bc_id << "  " << j << "  3  0.0" << std::endl; 
        bc_id++;
        
        output_file << bc_id << "  " <<  (j+n_grid_x*(n_grid_y-1)) << "  1  0.0" << std::endl; 
        bc_id++;
        output_file << bc_id << "  " <<  (j+n_grid_x*(n_grid_y-1)) << "  2  0.0" << std::endl; 
        bc_id++;
        output_file << bc_id << "  " <<  (j+n_grid_x*(n_grid_y-1)) << "  3  0.0" << std::endl; 	
        bc_id++;
      }
    
    for (unsigned int i=2; i<=(n_grid_y); i++)
      {
        output_file  << bc_id << "  " << ((i-1)*n_grid_x+1) << "  1  0.0" << std::endl; 
        bc_id++;
        output_file  << bc_id << "  " << ((i-1)*n_grid_x+1) << "  2  0.0" << std::endl; 
        bc_id++;
        output_file  << bc_id << "  " << ((i-1)*n_grid_x+1) << "  3  0.0" << std::endl; 
        bc_id++;
        
        output_file  << bc_id << "  " << ((i-1)*n_grid_x) << "  1  0.0" << std::endl; 
        bc_id++;
        output_file  << bc_id << "  " << ((i-1)*n_grid_x) << "  2  0.0" << std::endl; 
        bc_id++;
        output_file  << bc_id << "  " << ((i-1)*n_grid_x) << "  3  0.0" << std::endl; 	
        bc_id++;
      }
    
    
    output_file << "BOUNDARY_CONDITION_LOAD_SET END " << std::endl;
    
    
    // plate_transverse_load
    output_file << "SURFACE_LOAD_SET BEGIN"  << std::endl;
    output_file << "NAME SURFACE_PRESSURE_LOAD" << std::endl;
    output_file << "ID 3" << std::endl;
    output_file << "KIND SURFACE_PRESSURE" << std::endl;
    output_file << "N_LOADS " << (n_grid_x -1) << std::endl;
    elem_i = 1;
    for (int i = 0; i < (n_grid_x -1); i++ )
      {
        output_file << i+1 << "  " <<(i*2)+1<< "  0   " << 
        setw(8)<<"100.0"  << std::endl;
      }	
    output_file << "SURFACE_LOAD_SET END"  << std::endl;

    // plate_temperature load
    output_file << "NODAL_LOAD_SET BEGIN" << std::endl
    <<"NAME NODAL_TEMPERATURE " << std::endl
    <<"ID 6" << std::endl
    <<"KIND NODAL_TEMPERATURE" << std::endl
    <<"N_DOFS 1" << std::endl
    <<"N_LOADS " << n_grid_x*n_grid_y << std::endl;
    for ( int i = 0; i < n_grid_x*n_grid_y; i++ )
      output_file << i+1 << "  " << i+1 << "  " << setw(8)<< temperature << std::endl;
    
    output_file << "NODAL_LOAD_SET END" << std::endl    
    <<"LOAD_SET_DATA END " << std::endl
    <<"LOAD_CASE_DATA BEGIN " << std::endl
    <<"N_LOAD_CASES 1 " << std::endl
    <<"STATIC_LOAD_CASE BEGIN " << std::endl
    <<"ID 1 " << std::endl
    <<"N_DESIGN_VARIABLES 0 " << std::endl
    <<"N_LOAD_SETS 3 " << std::endl
    <<"DISPLACEMENT_BOUNDARY_CONDITION 1.0 5 " << std::endl
    << "PISTON_THEORY_SURFACE 1.0 4 " << std::endl
    << "NODAL_TEMPERATURE 1.0 6 " << std::endl; 
    
    output_file <<"STATIC_LOAD_CASE END " << std::endl
    <<"LOAD_CASE_DATA END " << std::endl
    <<"LOAD_DATABASE END " << std::endl
    <<"ANALYSIS_CASE BEGIN " << std::endl
    <<"ANALYSIS_TITLE TEST_CASE_FOR_NEW_FESYSTEM " << std::endl
    <<"N_ANALYSIS_DISCIPLINE_INFO 2 " << std::endl
    <<"STRUCTURAL_DISCIPLINE_INFO BEGIN " << std::endl
    <<"ID 1 " << std::endl
    <<"MESH_ID 1 " << std::endl
    <<"ANALYSIS_TYPE LINEAR_ANALYSIS " << std::endl
    <<"N_LOCAL_PARAMETERS 0 " << std::endl
    <<"STRUCTURAL_DISCIPLINE_INFO END " << std::endl
    <<"PISTON_THEORY_INFO BEGIN  " << std::endl
    <<"ID 1  " << std::endl
    <<"MESH_ID 1  " << std::endl
    <<"ORDER 1 " << std::endl
    <<"PISTON_THEORY_INFO END " << std::endl 
    <<"N_SOLUTIONS 1  " << std::endl
    <<"AEROELASTICITY_SOLUTION BEGIN " << std::endl 
    <<"ID 1  " << std::endl
    <<"AERODYNAMIC_DISCIPLINE PISTON_THEORY " << std::endl
    <<"SOLUTION_METHOD ROGER_APPROXIMATION_AEROELASTICITY_DRIVER " << std::endl
    <<"COUPLED_THERMAL_STRUCTURAL FALSE  " << std::endl
    <<"GEOMETRIC_EFFECTS   "; 
    if (use_geometric_stiffness == 0)
      output_file << "FALSE" << std::endl;
    else
      output_file << "TRUE" << std::endl;
    
    output_file <<"N_ANALYSIS_LOAD_CASES 1  " << std::endl
    <<"1  " << std::endl;
    output_file <<"N_FLIGHT_CONDITIONS " << n_conds << std::endl;
    for (unsigned int i=0; i<n_conds; i++)
      output_file <<"MACH_NUMBER " << (start_mach)/* + i * mach_incr)*/ << std::endl 
      <<"DENSITY " << density << std::endl
      <<"DYNAMIC_PRESSURE " << (.5*density*pow(speed_sound*(start_mach+i*mach_incr),2)) << std::endl
      <<"FLUID_FLOW_VECTOR 1.0 0.0 0.0 " << std::endl;
    output_file <<"N_SOLVER_INFO_ID 2  " << std::endl
    <<"EIGEN_SOLVER 2  " << std::endl
    <<"LINEAR_SOLVER 1  " << std::endl
    <<"AEROELASTICITY_SOLUTION END  " << std::endl
    <<"SOLVER_INFO BEGIN " << std::endl
    <<"N_SOLVER_INFO 2 " << std::endl
    <<"LINEAR_SOLVER_INFO BEGIN " << std::endl
    <<"ID 1 " << std::endl
    << "SOLVER_PACKAGE PETSC_LINEAR_SOLVER " << std::endl
    <<"KSP_TYPE KSP_PRE_ONLY " << std::endl
    <<"PC_TYPE PC_LU " << std::endl
    <<"LINEAR_SOLVER_INFO END " << std::endl
    <<"EIGEN_SOLVER_INFO BEGIN " << std::endl
    <<"ID 2 " << std::endl
    <<"EIGEN_SOLVER_TYPE EPS_ARNOLDI_EIGEN_SOLVER " << std::endl
    << "EIGEN_PROBLEM_TYPE  GENERALIZED_NON_HERMITIAN_EIGENPROBLEM" << std::endl
    << "SWAP_A_AND_B_MATRICES TRUE" << std::endl
    << "LINEAR_SOLVER_INFO_ID 1" << std::endl
    << "N_EIGEN_PAIRS   30" << std::endl
    << "CALCULATE_EIGEN_VECTORS TRUE" << std::endl
    << "EIGEN_SHIFT_TYPE SHIFT_AND_INVERT" << std::endl
    << "EIGEN_SHIFT_VALUE -0.5E-1" << std::endl
    << "EIGEN_SPECTRUM_END_TO_CALCULATE  LARGEST_IMAGINARY" << std::endl
    << "TOLERANCE  1.0e-6" << std::endl
    << "MAX_ITERATIONS 100" << std::endl
    <<"EIGEN_SOLVER_INFO END " << std::endl
    <<"SOLVER_INFO END " << std::endl
    << "INTERPOLATION_CASES BEGIN" << std::endl
    << "N_INTERPOLATION_CASES 0" << std::endl
    << "INTERPOLATION_CASES END" << std::endl
    <<"INTEGER_PARAMETERS BEGIN " << std::endl
    <<"N_INTEGER_PARAMETERS 0 " << std::endl
    <<"INTEGER_PARAMETERS END " << std::endl
    <<"REAL_PARAMETERS BEGIN " << std::endl
    <<"N_REAL_PARAMETERS 1 " << std::endl
    <<"REFERENCE_TEMP 30.00 " << std::endl
    <<"REAL_PARAMETERS END " << std::endl
    <<"STRING_PARAMETERS BEGIN " << std::endl
    <<"N_STRING_PARAMETERS 0 " << std::endl
    <<"STRING_PARAMETERS END " << std::endl
    <<"OUTPUT_INFO BEGIN " << std::endl
    <<"OUTPUT_FORMAT GMSH_OUTPUT_PROCESSOR " << std::endl
    <<"OUTPUT_FILE_NAME output_file " << std::endl
    <<"OUTPUT_INFO END " << std::endl
    <<"ANALYSIS_CASE END " << std::endl
    <<"DESIGN_DATABASE BEGIN " << std::endl
    <<"DESIGN_PARAMETERS BEGIN " << std::endl
    <<"N_PARAMETERS 0 " << std::endl
    <<"DESIGN_PARAMETERS END " << std::endl
    <<"DESIGN_DATABASE END " << std::endl
    <<"END_LOAD_SET" << std::endl;	
  }
  
  
  // write the nastran input file
  {
    
    std::fstream output_file;
    output_file.open("nastran_input.bdf", std::fstream::out);
    
    output_file << "SOL 145 " << std::endl
    <<"CEND" << std::endl
    << "$POST TOFILE 12 DISPLACE\n"
    << "ECHO=NONE\n"
    << "SPC=2\n"
    << "AEROF=NONE\n"
    << "APRESSURE=NONE\n"
    << "DISPLACEMENT=ALL\n"
    << "FMETHOD=1\n"
    << "METHOD=1\n"
    << "CMETHOD=2\n"
    << "SVECTOR=ALL\n"
    << "BEGIN BULK\n" 
    << "PARAM,POST,0\n"
    << "PARAM,LMODES,50\n"
    << "MKAERO2,2.0,.01,2.0,.03,2.0,0.06,2.0,.09\n"
    << "MKAERO2,2.0,0.12,2.0,0.15\n"
    << "AERO,,,5.0,0.04\n"
    << "$AEROS,,,10.0,50.0,50.0\n"
    << "FLUTTER,1,PK,1,2,3,S,20,1.0e-6\n"
    << "$ density ratio\n"
    << "FLFACT,1,1.0\n"
    << "$ Mach numbers for analysis\n"
    << "FLFACT,2,2.0\n"
    << "$ velocities\n"
    << "FLFACT,3,-400.0,-420.0,-440.0,-460.0,-480.0,-500.0,-520.0,\n"
    << ",-540.0,-560.0,-580.0,-600.0,-620.0,-640.0,-660.0,680.0\n"
    << "EIGRL,1,0.0,,50\n"
    << "EIGC,2,CLAN,,,,,20\n"
    << "$\n" 
    << "MAT1,1,27.0e+9,,0.33,2000.0" << std::endl
    << "PSHELL,1,1,0.0075,1" << std::endl
    << "PAERO5,1,1,4,\n,0.0" << std::endl
    << "AEFACT,4,2.0,0.0997,2.5,0.0997,3.0,0.0997,3.5\n,0.0997" << std::endl
    << "AEFACT,5,0.0,0.0,0.0,0.0,0.0,0.0"<<std::endl;
    
    for ( unsigned int i=0; i < n_grid_total; i++)
      {
        output_file << std::showpoint << "GRID," << i+1	<< ",," 
        << std::setprecision(7) << grid_cord[i][1]
        << "," << std::setprecision(7) << grid_cord[i][2] 
        << ",0.0" << std::endl; 
      }
    
    int n1 = 1;
    int elem_i = 1;
    for (int i = 0; i < n_grid_y-1; i++ )
      {
        for (int j = 0; j < n_grid_x-1; j++)
          {
            output_file << "SET1,"<<2000000+elem_i<<","<<n1<<","<<n1+1<<","<<n1+n_grid_x<<","
            <<n1+1+n_grid_x<<std::endl;
            output_file << "CTRIA3,"<<elem_i++<<",1,"<<n1<<","<<n1 + 1<<","<<n1+1+n_grid_x<<"\n"; 
            output_file << "CTRIA3,"<<elem_i++<<",1,"<<n1<<","<<n1+1+n_grid_x<<","<<n1+n_grid_x<<"\n"; 
            n1 = n1 + 1;
          }
        n1 = n1 + 1 ;
      }	
    
    // piston theory surface
    n1 = 1;
    elem_i = 1;
    for (int i = 0; i < n_grid_y-1; i++ )
      for (int j = 0; j < n_grid_x-1; j++)
        {
          output_file << std::showpoint << "CAERO5," <<100000+elem_i << ",1,,1,,0,5,," << std::endl << ","
          << dx*j <<","<<dy*i<<",0.0,"<<dx<<","<<dx*j<<","<<dy*(i+1)<<",0.0,"<<dx<<std::endl;
          output_file << "SPLINE1," <<1000000+elem_i<<","<<100000+elem_i<<","<<100000+elem_i<<","
          <<100000+elem_i<<","<<2000000+elem_i<< ",0.0,IPS,BOTH" << std::endl;
          elem_i++; elem_i++;
        }

    // displacement boundary conditions	
    for (unsigned int i=1; i<=n_grid_y; i++)
      for (unsigned int j=1; j<=n_grid_x; j++)
        {
          unsigned int n_node = (i-1)*n_grid_x + j ;
          output_file <<"SPC,2,"<< n_node << ",6,0.0" << std::endl;
        }
    
    for (unsigned int j=1; j<=n_grid_x; j++)
      {
        output_file << "SPC,2,"<< j <<",123,0.0" << std::endl; 
        output_file << "SPC,2,"<<(j+n_grid_x*(n_grid_y-1)) << ",123,0.0" << std::endl; 
      }
    
    for (unsigned int i=2; i<=(n_grid_y); i++)
      {
        output_file  << "SPC,2," << ((i-1)*n_grid_x+1) << ",123,0.0" << std::endl; 
        output_file  << "SPC,2," << ((i-1)*n_grid_x) << ",123,0.0" << std::endl; 
      }
    
    output_file<<"ENDDATA" << std::endl;
  }
  
  return 0;
}
