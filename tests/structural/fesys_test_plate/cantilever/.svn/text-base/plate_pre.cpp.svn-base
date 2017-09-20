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

output_file << "BEGIN_PROPERTY_LIST 1" << std::endl;
output_file << "BEGIN_PROPERTY_CARD" << std::endl;
output_file << "ID 100" << std::endl;
output_file << "THICKNESS 0.005" << std::endl;
output_file << "E_11 2.0e9" << std::endl;
output_file << "NU 0.33" << std::endl;
output_file << "ALPHA_EXPANSION 1.0e-3" << std::endl;
output_file << "TEMP_REF 0.0" << std::endl;
output_file << "END_PROPERTY_CARD" << std::endl;
output_file << "END_PROPERTY_LIST" << std::endl;
output_file << "BEGIN_MESHLIST 1" << std::endl;
output_file << "BEGIN_MESH 100" << std::endl;
output_file << "BEGIN_NODES " << n_grid_total << std::endl;

for ( unsigned int i=0; i < n_grid_total; i++)
{
	output_file << i+1	<< "   " 
		    << std::setprecision(7) << grid_cord[i][1]
		    << "   " << std::setprecision(7) << grid_cord[i][2] 
		    << "   0.0" << std::endl; 
}

output_file << "END_NODES" << std::endl;
int n_elem = (n_grid_y-1)*(n_grid_x-1);
output_file << "BEGIN_ELEMENTS " <<  n_elem << std::endl;

int n1 = 1;
int elem_i = 1;
for (int i = 0; i < n_grid_y-1; i++ )
{
	for (int j = 0; j < n_grid_x-1; j++)
	{
		output_file << "PLATE_QUAD4 "	<<setw(8)<<elem_i++<<
		setw(8)<<"100"<<setw(8)<<n1<<setw(8)<<n1 + 1<<
		setw(8)<<n1+1+n_grid_x<<setw(8)<<n1+n_grid_x<<"\n"; 
		n1 = n1 + 1;
	}
	n1 = n1 + 1 ;
}	

output_file << "END_ELEMENTS" << std::endl;
output_file << "END_MESH" << std::endl;
output_file << "END_MESHLIST" << std::endl;
output_file << "BEGIN_RADIATION_CAVITY_LIST 0" << std::endl;
output_file << "END_RADIATION_CAVITY_LIST" << std::endl 
<< "BEGIN_ELEMENT_SET_LIST 0"<< std::endl
<< "END_ELEMENT_SET_LIST" << std::endl;
output_file << "BEGIN_LOAD_SET_DATA 4" << std::endl;


// temperature load
output_file << "BEGIN_LOAD_SET " << n_grid_total << std::endl;
output_file << "LOAD_SET_NAME NODAL_TEMPERATURE_LOAD" << std::endl;
output_file << "LOAD_SET_ID 101" << std::endl;
output_file << "LOAD_SET_KIND POINT_LOAD" << std::endl;
output_file << "LOAD_SET_N_DOFS 1" << std::endl;

for ( unsigned int i=0; i < n_grid_total; i++)
  output_file << i+1 << "  "<<  i+1	<< "  0.0 " << std::endl; 

output_file << "END_LOAD_SET" << std::endl;


// point load
output_file << "BEGIN_LOAD_SET 0" << std::endl;
output_file << "LOAD_SET_NAME NODAL_LOADS" << std::endl;
output_file << "LOAD_SET_ID 102" << std::endl;
output_file << "LOAD_SET_KIND POINT_LOAD" << std::endl;
output_file << "LOAD_SET_N_DOFS 6" << std::endl;
output_file << "END_LOAD_SET" << std::endl;



// plate_transverse_load
output_file << "BEGIN_LOAD_SET " << n_elem << std::endl;
output_file << "LOAD_SET_NAME SURFACE_PRESSURE_LOAD" << std::endl;
output_file << "LOAD_SET_ID 103" << std::endl;
output_file << "LOAD_SET_KIND SURFACE_LOAD" << std::endl;
output_file << "LOAD_SET_N_DOFS 1" << std::endl;
elem_i = 1;
unsigned int load_id=1;
for (int i = 0; i < n_grid_y-1; i++ )
{
	for (int j = 0; j < n_grid_x-1; j++)
	{
	  output_file << load_id << "  " <<elem_i++<< "  4   " << 
		setw(8)<<"100.0"  << std::endl;
	  load_id++;
	}
}	

output_file << "END_LOAD_SET" << std::endl;	

// displacement boundary conditions	
output_file << "BEGIN_LOAD_SET " << (n_grid_y*6) << std::endl;
output_file << "LOAD_SET_NAME DISPLACEMENT_BOUNDARY_CONDITIONS" << std::endl;
output_file << "LOAD_SET_ID 104" << std::endl;
output_file << "LOAD_SET_KIND BOUNDARY_CONDITION" << std::endl;
output_file << "LOAD_SET_N_DOFS 1" << std::endl;
/*for (unsigned int i=1; i<=n_grid_y; i++)
	for (unsigned int j=1; j<=n_grid_x; j++)
	{
		unsigned int n_node = (i-1)*n_grid_x + j ;
		output_file << n_node << "  6  0.0" << std::endl; 
	}
*/
load_id=1;
for (unsigned int j=1; j<=n_grid_y; j++)
{
  output_file << load_id << "   " <<(n_grid_x*(j-1) + 1) << "  1  0.0" << std::endl; 
  load_id++;
  output_file << load_id << "   " <<(n_grid_x*(j-1) + 1) << "  2  0.0" << std::endl; 
  load_id++;
  output_file << load_id << "   " <<(n_grid_x*(j-1) + 1) << "  3  0.0" << std::endl; 
  load_id++;
  output_file << load_id << "   " <<(n_grid_x*(j-1) + 1) << "  4  0.0" << std::endl; 
  load_id++;
  output_file << load_id << "   " <<(n_grid_x*(j-1) + 1) << "  5  0.0" << std::endl; 
  load_id++;
  output_file << load_id << "   " <<(n_grid_x*(j-1) + 1) << "  6  0.0" << std::endl; 
  load_id++;
}


output_file << "END_LOAD_SET" << std::endl;	



output_file << "END_LOAD_SET_DATA" << std::endl;



output_file << "BEGIN_LOAD_CASE_DATA 1" << std::endl;

output_file << "BEGIN_LOAD_CASE" << std::endl;
output_file << "LOAD_CASE_ID 1" << std::endl;
output_file << "N_TIME_INSTANTS 1" << std::endl;
output_file << "N_LOAD_SETS 4" << std::endl;
output_file << "TIME_INSTANT 0.0" << std::endl;
output_file << "STRUCTURAL_BOUNDARY_CONDITION 104" << std::endl;
output_file << "NODAL_LOAD 102" << std::endl;
output_file << "NODAL_TEMPERATURE_LOAD 101" << std::endl;
output_file << "SURFACE_PRESSURE_LOAD 103" << std::endl;
output_file << "END_LOAD_CASE" << std::endl;

output_file << "END_LOAD_CASE_DATA" << std::endl;


output_file << "BEGIN_ANALYSIS_CASE " << std::endl;
output_file << 	"ANALYSIS_TITLE test_case" << std::endl;
output_file << 	"ANALYSIS_DISCIPLINE STRUCTURAL" << std::endl;
output_file << 	"MESH_ID 100 " << std::endl;
output_file << 	"ANALYSIS_TYPE LINEAR" << std::endl;
output_file << 	"INTEGER_PARAMETERS 1" << std::endl;
output_file << 	"NONLINEAR_SOLVER_MAXIMUM_ITERATIONS 20" << std::endl;
output_file << 	"REAL_PARAMETERS 1" << std::endl;
output_file << 	"NONLINEAR_SOLVER_TOLERANCE 1.0e-3" << std::endl;
output_file << 	"ANALYSIS_LOAD_CASES 1" << std::endl;
output_file << 	"1 " << std::endl;
output_file << 	"BEGIN_DESIGN_VARIABLES" << std::endl;
output_file << 	"PROPERTY_DESIGN_VARIABLES 0" << std::endl;
output_file << 	"SHAPE_DESIGN_VARIABLES 0" << std::endl;
output_file << 	"END_DESIGN_VARIABLES" << std::endl;
output_file << 	"BEGIN_OUTPUT_REQUEST" << std::endl;
output_file << 	"DATA_VISUALIZER_INPUT_FILE " << std::endl;
output_file << "FORMAT TECPLOT" << std::endl;
output_file << "NAME output_plots.dat" << std::endl;
output_file << 	"END_OUTPUT_REQUEST" << std::endl;
output_file << 	"END_ANALYSIS_CASE" << std::endl;



return 0;
}
