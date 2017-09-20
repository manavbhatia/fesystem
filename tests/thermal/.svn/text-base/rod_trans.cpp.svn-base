#include<iostream>
#include<iomanip>
#include<fstream>

using std::cout;
using std::cin;
using std::endl;
using std::setw;
using std::ios;
int main()
{
 const int max_grid_pt = 2500;
 int grid_no[max_grid_pt]={0}; // maximum grid points is 500
 float grid_cord[max_grid_pt][3] = {{0.}}; //maximum grid points is 500
 int n_grid_x, n_grid_y, n_grid_total; // number of grid points 
 float x_side_len,y_side_len; // plate side lengths
 float dx,dy; 
 float x_init = 0.,y_init= 0.;

 cout << "\nEnter the side length along x- axis: ";
 cin >> x_side_len;

 cout << "\nEnter the side length along y- axis: ";
 cin >> y_side_len;

 cout << "\nGrid points along x- axis :";
 cin >> n_grid_x; 
  
 cout << "\nGrid points along y- axis :";
 cin >> n_grid_y;

// following statements calculate the data
 dx = x_side_len / float(n_grid_x -1);
 dy = y_side_len / float(n_grid_y -1);
 n_grid_total = n_grid_x * n_grid_y;

 int grid_number = 0;
 for (int i = 0; i < n_grid_y; i++ )
	{
	y_init = dy * float(i);
	for (int j = 0; j < n_grid_x; j++)
		{ x_init = dx * float(j);
	        grid_cord[grid_number][1] = x_init;
		grid_cord[grid_number][2] = y_init;  	
                grid_number = grid_number + 1; 
		} 
	}

 std::fstream output_file;
output_file.open("fesys_input.txt", std::fstream::out);

 output_file << "FUNCTION_DATABASE BEGIN " << std::endl
	     <<" N_FUNCTIONS 1 " << std::endl
   << "MULTILINEAR_FUNCTION BEGIN" << std::endl
   << "ID 1 " << std::endl
   << "N_PAIRS 2" << std::endl
   << "0.0 1.0" << std::endl
   << "1.0 1.0" << std::endl
  << "MULTILINEAR_FUNCTION END" << std::endl
<<"FUNCTION_DATABASE END " << std::endl
<<"PROPERTY_DATABASE BEGIN " << std::endl
<<"MATERIAL_PROPERTY_CARDS BEGIN " << std::endl
<<"N_MATERIAL_PROPERTY_CARDS 1 " << std::endl
<<"ISOTROPIC_MATERIAL_PROPERTY_CARD BEGIN " << std::endl
<<"ID 1 " << std::endl
<<"N_PROPERTIES 6 " << std::endl
<<"YOUNGS_MODULUS 27.0e+9 " << std::endl
<<"POISSONS_RATIO 0.33 " << std::endl
<<"ALPHA_EXPANSION 13.0e-6 " << std::endl
   << "THERMAL_CONDUCTIVITY 2.0" << std::endl
   << "SPECIFIC_HEAT 2.0" << std::endl
   << "DENSITY 15.0" << std::endl
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
<<"THICKNESS 0.01 " << std::endl
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
for (int i = 0; i < n_grid_y-1; i++ )
{
	for (int j = 0; j < n_grid_x-1; j++)
	{
		output_file << "THERMAL_CONDUCTION_TRI3 "	<<setw(8)<<elem_i++<<
		setw(8)<<"1"<<setw(8)<<n1<<setw(8)<<n1 + 1<<
		  setw(8)<<n1+1+n_grid_x<<setw(8)<<"\n"; 

		output_file << "THERMAL_CONDUCTION_TRI3 "	<<setw(8)<<elem_i++<<
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
<<"N_LOAD_SETS 4 " << std::endl
<<"NODAL_LOAD_SET BEGIN " << std::endl
<<"NAME NODAL_TEMP " << std::endl
<<"ID 1 " << std::endl
<<"KIND NODAL_TEMPERATURE " << std::endl
<<"N_DOFS 1 " << std::endl
<<"N_LOADS 0 " << std::endl
<<"NODAL_LOAD_SET END" << std::endl ;

// point load
output_file << "NODAL_LOAD_SET BEGIN" << std::endl
<<"NAME NODAL_FORCE " << std::endl
<<"ID 2" << std::endl
<<"KIND NODAL_FORCE" << std::endl
<<"N_DOFS 6" << std::endl
<<"N_LOADS 0" << std::endl
<<"NODAL_LOAD_SET END" << std::endl;


// plate_transverse_load
output_file << "SURFACE_LOAD_SET BEGIN"  << std::endl;
output_file << "NAME SURFACE_HEAT_LOAD" << std::endl;
output_file << "ID 3" << std::endl;
output_file << "KIND SURFACE_HEAT_LOAD" << std::endl;
output_file << "N_LOADS " << (n_grid_x -1) << std::endl;
elem_i = 1;
unsigned int load_id=1;
for (int i = 0; i < (n_grid_x -1); i++ )
{
  output_file << i+1 << "  " <<(i*2)+1<< "  0   " << 
  setw(8)<<"1.0"  << std::endl;
}	

output_file << "SURFACE_LOAD_SET END" << std::endl;	

// displacement boundary conditions	
output_file << "BOUNDARY_CONDITION_LOAD_SET BEGIN " << std::endl;
output_file << "NAME DISPLACEMENT_BOUNDARY_CONDITIONS" << std::endl;
output_file << "ID 4" << std::endl;
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


output_file << "BOUNDARY_CONDITION_LOAD_SET END " << std::endl
<<"LOAD_SET_DATA END " << std::endl
<<"LOAD_CASE_DATA BEGIN " << std::endl
<<"N_LOAD_CASES 1 " << std::endl
<<"TRANSIENT_LOAD_CASE BEGIN " << std::endl
<<"ID 1 " << std::endl
<<"N_TRANSIENT_LOAD_FACTOR_DEPENDENT_LOAD_SETS 1 " << std::endl
<<"N_TIME_INSTANTS 0 " << std::endl
<<"N_DESIGN_VARIABLES 0 " << std::endl
<<"SURFACE_HEAT_LOAD 1 3 " << std::endl
<<"TRANSIENT_LOAD_CASE END " << std::endl
<<"LOAD_CASE_DATA END " << std::endl
<<"LOAD_DATABASE END " << std::endl
<<"ANALYSIS_CASE BEGIN " << std::endl
<<"ANALYSIS_TITLE TEST_CASE_FOR_NEW_FESYSTEM " << std::endl
<<"N_ANALYSIS_DISCIPLINE_INFO 1 " << std::endl
<<"THERMAL_DISCIPLINE_INFO BEGIN " << std::endl
<<"ID 1 " << std::endl
<<"MESH_ID 1 " << std::endl
<<"ANALYSIS_TYPE LINEAR_ANALYSIS " << std::endl
<<"N_RADIATION_CAVITIES 0 " << std::endl
<<"N_LOCAL_PARAMETERS 0 " << std::endl
<<"THERMAL_DISCIPLINE_INFO END " << std::endl
<<"N_SOLUTIONS 1 " << std::endl
<<"TRANSIENT_THERMAL_SOLUTION BEGIN " << std::endl
<<"ID 1 " << std::endl
<<"N_ANALYSIS_LOAD_CASES 1 " << std::endl
<<"1 " << std::endl
<<"N_SOLVER_INFO_ID 2 " << std::endl
<<"LINEAR_TRANSIENT_SOLVER 3 " << std::endl
<<"NONLINEAR_TRANSIENT_SOLVER 4 " << std::endl
<<"TRANSIENT_THERMAL_SOLUTION END " << std::endl
<<"SOLVER_INFO BEGIN " << std::endl
<<"N_SOLVER_INFO 5 " << std::endl
<<"LINEAR_SOLVER_INFO BEGIN " << std::endl
<<"ID 1 " << std::endl
<<"KSP_TYPE KSP_PRE_ONLY " << std::endl
<<"PC_TYPE PC_LU " << std::endl
<<"LINEAR_SOLVER_INFO END " << std::endl
<<"EIGEN_SOLVER_INFO BEGIN " << std::endl
<<"ID 2 " << std::endl
<<"EIGEN_SOLVER_TYPE ARPACK_EIGEN_SOLVER " << std::endl
<< "EIGEN_PROBLEM_TYPE  GENERALIZED_HERMITIAN_EIGENPROBLEM" << std::endl
<< "LINEAR_SOLVER_INFO_ID 1" << std::endl
<< "N_EIGEN_PAIRS   10" << std::endl
<< "CALCULATE_EIGEN_VECTORS TRUE" << std::endl
<< "EIGEN_SHIFT_TYPE SHIFT_AND_INVERT" << std::endl
<< "EIGEN_SHIFT_VALUE 1.0e-4" << std::endl
<< "EIGEN_SPECTRUM_END_TO_CALCULATE  LARGEST_MAGNITUDE" << std::endl
<< "TOLERANCE  1.0e-6" << std::endl
<< "MAX_ITERATIONS 100" << std::endl
<<"EIGEN_SOLVER_INFO END " << std::endl
<<"NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO BEGIN " << std::endl
<<"ID 3 " << std::endl
<<"LINEAR_SOLVER_INFO_ID 1 " << std::endl
<< "MAX_ITERATIONS 100" << std::endl
<< "INITIAL_TIME 0.0" << std::endl
<< "TOTAL_DURATION 20.0" << std::endl
<< "IF_ADAPTIVE_SIZE FALSE" << std::endl
<< "INITIAL_STEP_SIZE 1.0" << std::endl
<< "GLOBAL_TOLERANCE 1.0e-6" << std::endl
<< "SYSTEM_ORDER 1" << std::endl
<< "CONSTANTS .5" << std::endl
<<"NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO END " << std::endl
<<"NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO BEGIN " << std::endl
<<"ID 4 " << std::endl
<<"NONLINEAR_SOLVER_INFO_ID 5 " << std::endl
<< "MAX_ITERATIONS 100" << std::endl
<< "INITIAL_TIME 0.0" << std::endl
<< "TOTAL_DURATION 20.0" << std::endl
<< "IF_ADAPTIVE_SIZE FALSE" << std::endl
<< "INITIAL_STEP_SIZE 1.0" << std::endl
<< "GLOBAL_TOLERANCE 1.0e-6" << std::endl
<< "SYSTEM_ORDER 1" << std::endl
<< "CONSTANTS .5" << std::endl
<<"NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO END " << std::endl
<< "NONLINEAR_SOLVER_INFO BEGIN" << std::endl
<< "ID 5" << std::endl
<< "SOLVER_TYPE FESYSTEM_NEWTON_NONLINEAR_SOLVER" << std::endl
<< "LINEAR_SOLVER_INFO_ID 1" << std::endl
<< "MAXIMUM_ITERATIONS 20" << std::endl
<< "CONVERGENCE_TOLERANCE 1.0e-6" << std::endl
<< "NONLINEAR_SOLVER_INFO END" << std::endl
<<"SOLVER_INFO END " << std::endl
<<"INTEGER_PARAMETERS BEGIN " << std::endl
<<"N_INTEGER_PARAMETERS 0 " << std::endl
<<"INTEGER_PARAMETERS END " << std::endl
<<"REAL_PARAMETERS BEGIN " << std::endl
<<"N_REAL_PARAMETERS 0 " << std::endl
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






return 0;
}
