Include "Cylinder.geo";
Include "EllipticShell.geo";


Geometry.AutoCoherence = 0;

// set the number of stages and then all the arrays of DVs will 
// have this dimension
n_stages = 2;



// 
// ********initialize the number of div variables **************
//
n_circumference_div = 0;
n_circum_elems_per_div = 0;

n_divs_along_length_nose = 0;
n_elems_along_length_per_div_nose = 0;

n_divs_along_length_cylinder1_upper_end[n_stages] = 0;
n_elems_along_length_per_div_cylinder1_upper_end[n_stages] = 0;

n_divs_along_length_cylinder1[n_stages] = 0;
n_elems_along_length_per_div_cylinder1[n_stages] = 0;

n_divs_along_length_cylinder1_lower_end[n_stages] = 0 ;
n_elems_along_length_per_div_cylinder1_lower_end[n_stages] = 0;

n_divs_along_length_intertank_skirt[n_stages] = 0;
n_elems_along_length_per_div_intertank_skirt[n_stages] = 0;

n_divs_along_length_cylinder2_upper_end[n_stages] = 0;
n_elems_along_length_per_div_cylinder2_upper_end[n_stages] = 0;

n_divs_along_length_cylinder2[n_stages] = 0;
n_elems_along_length_per_div_cylinder2[n_stages] = 0;

n_divs_along_length_cylinder2_lower_end[n_stages] = 0;
n_elems_along_length_per_div_cylinder2_lower_end[n_stages] = 0;

n_divs_along_length_interstage_skirt = 0;
n_elems_along_length_per_div_interstage_skirt = 0;



//
// set the variables for number of divisions 
//
n_circumference_div = 8;
n_circum_elems_per_div = 4;

// 
// nose
// 
n_divs_along_length_nose = 2;
n_elems_along_length_per_div_nose = 4;


// 
// stage 1
// 

stage_number = 0;

// 
// cylinder end1
// 
n_divs_along_length_cylinder1_upper_end[stage_number] = 2;
n_elems_along_length_per_div_cylinder1_upper_end[stage_number] = 4;


// 
// cylinder
// 
n_divs_along_length_cylinder1[stage_number] = 2;
n_elems_along_length_per_div_cylinder1[stage_number] = 4;


// 
// cylinder end 2
// 
n_divs_along_length_cylinder1_lower_end[stage_number] = 2;
n_elems_along_length_per_div_cylinder1_lower_end[stage_number] = 4;


// 
// skirt
// 
n_divs_along_length_intertank_skirt[stage_number] = 1;
n_elems_along_length_per_div_intertank_skirt[stage_number] = 4;


// 
// cylinder end1
// 
n_divs_along_length_cylinder2_upper_end[stage_number] = 2;
n_elems_along_length_per_div_cylinder2_upper_end[stage_number] = 4;


// 
// cylinder
// 
n_divs_along_length_cylinder2[stage_number] = 2;
n_elems_along_length_per_div_cylinder2[stage_number] = 4;


// 
// cylinder end 2
// 
n_divs_along_length_cylinder2_lower_end[stage_number] = 2;
n_elems_along_length_per_div_cylinder2_lower_end[stage_number] = 4;


// 
// skirt
// 
n_divs_along_length_interstage_skirt = 1;
n_elems_along_length_per_div_interstage_skirt  = 4;



// 
// stage 2
// 

stage_number = 1;

// 
// cylinder end1
// 
n_divs_along_length_cylinder1_upper_end[stage_number] = 2;
n_elems_along_length_per_div_cylinder1_upper_end[stage_number] = 4;


// 
// cylinder
// 
n_divs_along_length_cylinder1[stage_number] = 2;
n_elems_along_length_per_div_cylinder1[stage_number] = 4;


// 
// cylinder end 2
// 
n_divs_along_length_cylinder1_lower_end[stage_number] = 2;
n_elems_along_length_per_div_cylinder1_lower_end[stage_number] = 4;


// 
// skirt
// 
n_divs_along_length_intertank_skirt[stage_number] = 1;
n_elems_along_length_per_div_intertank_skirt[stage_number] = 4;


// 
// cylinder end1
// 
n_divs_along_length_cylinder2_upper_end[stage_number] = 2;
n_elems_along_length_per_div_cylinder2_upper_end[stage_number] = 4;


// 
// cylinder
// 
n_divs_along_length_cylinder2[stage_number] = 2;
n_elems_along_length_per_div_cylinder2[stage_number] = 4;


// 
// cylinder end 2
// 
n_divs_along_length_cylinder2_lower_end[stage_number] = 2;
n_elems_along_length_per_div_cylinder2_lower_end[stage_number] = 4;






//
// ********* initialize all the variables ********
//
// stage independent variables
// nose
nose_height_DV = 0;


// stage dependent variables 
stage_radius_DV[n_stages-1] = 0;
stage_cylinder1_upper_shell_height_DV[n_stages-1] = 0;
stage_cylinder1_height_DV[n_stages-1] = 0;
stage_cylinder1_lower_shell_height_DV[n_stages-1] = 0;

stage_skirt_height_DV[n_stages-1] = 0;

stage_cylinder2_upper_shell_height_DV[n_stages-1] = 0;
stage_cylinder2_height_DV[n_stages-1] = 0;
stage_cylinder2_lower_shell_height_DV[n_stages-1] = 0;




//
// *********** now set the values of the variables *******
//
nose_height_DV = 1.5;

// 
// stage 1 
//
stage_number = 0;

stage_radius_DV[stage_number] = 0.6;
stage_cylinder1_upper_shell_height_DV[stage_number] = 0.25;
stage_cylinder1_height_DV[stage_number] = 1.5;
stage_cylinder1_lower_shell_height_DV[stage_number] = 0.25;

stage_skirt_height_DV[stage_number] = 0.6;

stage_cylinder2_upper_shell_height_DV[stage_number] = 0.25;
stage_cylinder2_height_DV[stage_number] = 0.75;
stage_cylinder2_lower_shell_height_DV[stage_number] = 0.25;


interstage_skirt_height_DV = 0.6;

// stage 2
stage_number = 1;

stage_radius_DV[stage_number] = 0.75;
stage_cylinder1_upper_shell_height_DV[stage_number] = 0.25;
stage_cylinder1_height_DV[stage_number] = 2.0;
stage_cylinder1_lower_shell_height_DV[stage_number] = 0.25;

stage_skirt_height_DV[stage_number] = .6;

stage_cylinder2_upper_shell_height_DV[stage_number] = 0.25;
stage_cylinder2_height_DV[stage_number] = 1.25;
stage_cylinder2_lower_shell_height_DV[stage_number] = 0.25;




// 
// ********* now create the points **********
// 
// there will be 8 control points for each stage
//


//
// nose
//
Point(1) = {0.0,0.0,0.0,0.1};


//
// stage 1
//
stage_number = 0;
stage_starting_control_point_ID = 2;
stage_starting_x_coord = 0.0 + nose_height_DV;


// point for cylinder 1 upper end
next_point_ID = stage_starting_control_point_ID;
next_x_coord = stage_starting_x_coord;
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for cylinder 1 upper shell end
next_point_ID += 1;
next_x_coord -= stage_cylinder1_upper_shell_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for cylinder 1 lower end
next_point_ID += 1;
next_x_coord = stage_starting_x_coord + stage_cylinder1_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for cylinder 1 lower shell end
next_point_ID += 1;
next_x_coord += stage_cylinder1_lower_shell_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for inter-cylinder skirt lower end, and cylinder 2 upper end
next_point_ID += 1;
next_x_coord = stage_starting_x_coord + stage_cylinder1_height_DV[stage_number] + stage_skirt_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for cylinder 2 upper shell end
next_point_ID += 1;
next_x_coord -= stage_cylinder2_upper_shell_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for cylinder 2 lower end
next_point_ID += 1;
next_x_coord = stage_starting_x_coord + stage_cylinder1_height_DV[stage_number] + stage_skirt_height_DV[stage_number] + stage_cylinder2_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for cylinder 2 lower shell end
next_point_ID += 1;
next_x_coord += stage_cylinder2_lower_shell_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};



// 
//  repeat the same process for the second stage
//  
//
stage_number = 1;
stage_starting_control_point_ID = next_point_ID + 1;
stage_starting_x_coord = nose_height_DV + stage_cylinder1_height_DV[stage_number-1] + stage_skirt_height_DV[stage_number-1] +  stage_cylinder2_height_DV[stage_number-1] + interstage_skirt_height_DV;

// point for cylinder 1 upper end
next_point_ID = stage_starting_control_point_ID;
next_x_coord = stage_starting_x_coord;
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for cylinder 1 upper shell end
next_point_ID += 1;
next_x_coord -= stage_cylinder1_upper_shell_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for cylinder 1 lower end
next_point_ID += 1;
next_x_coord = stage_starting_x_coord + stage_cylinder1_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for cylinder 1 lower shell end
next_point_ID += 1;
next_x_coord += stage_cylinder1_lower_shell_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for inter-cylinder skirt lower end, and cylinder 2 upper end
next_point_ID += 1;
next_x_coord = stage_starting_x_coord + stage_cylinder1_height_DV[stage_number] + stage_skirt_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for cylinder 2 upper shell end
next_point_ID += 1;
next_x_coord -= stage_cylinder2_upper_shell_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for cylinder 2 lower end
next_point_ID += 1;
next_x_coord = stage_starting_x_coord + stage_cylinder1_height_DV[stage_number] + stage_skirt_height_DV[stage_number] + stage_cylinder2_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};


// point for cylinder 2 lower shell end
next_point_ID += 1;
next_x_coord += stage_cylinder2_lower_shell_height_DV[stage_number];
Point(next_point_ID) = {next_x_coord, 0.0, 0.0, 0.1};





// now create the geometries



//n_divs_along_length_cylinder1[n_stages] = 0;
//n_elems_along_length_per_div_cylinder1[n_stages] = 0;

//n_divs_along_length_cylinder1_lower_end[n_stages] = 0 ;
//n_elems_along_length_per_div_cylinder1_lower_end[n_stages] = 0;

//n_divs_along_length_intertank_skirt[n_stages] = 0;
//n_elems_along_length_per_div_intertank_skirt[n_stages] = 0;

//n_divs_along_length_cylinder2_upper_end[n_stages] = 0;
//n_elems_along_length_per_div_cylinder2_upper_end[n_stages] = 0;

//n_divs_along_length_cylinder2[n_stages] = 0;
//n_elems_along_length_per_div_cylinder2[n_stages] = 0;

//n_divs_along_length_cylinder2_lower_end[n_stages] = 0;
//n_elems_along_length_per_div_cylinder2_lower_end[n_stages] = 0;

//n_divs_along_length_interstage_skirt = 0;
//n_elems_along_length_per_div_interstage_skirt = 0;


//
// nose
//

Printf("FESystem_Message: GEOMETRY_BEGIN NOSE");

stage_number = 0;
func_ellshell_end1_pointID = 1;
func_ellshell_end2_pointID = 2;
func_ellshell_end2_radius = stage_radius_DV[stage_number];
func_ellshell_n_divs_along_circumference = n_circumference_div;
func_ellshell_n_angular_divs_along_length = n_divs_along_length_nose;
func_ellshell_elems_along_circum_per_div = n_circum_elems_per_div;
func_ellshell_elems_along_length_per_div = n_elems_along_length_per_div_nose;
func_ellshell_local_y_axis[2] = 0;

func_ellshell_local_y_axis[0] = 0;
func_ellshell_local_y_axis[1] = 1;
func_ellshell_local_y_axis[2] = 0;

Call ELLIPTICSHELL;

Printf("FESystem_Message: GEOMETRY_END NOSE");

//
// stage 1
//
stage_number = 0;

//
// clyinder1 end 1
// 

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_1_CYLINDER_1_END_1");

func_ellshell_end1_pointID = 3;
func_ellshell_end2_pointID = 2;
func_ellshell_end2_radius = stage_radius_DV[stage_number];
func_ellshell_n_divs_along_circumference = n_circumference_div;
func_ellshell_n_angular_divs_along_length = n_divs_along_length_cylinder1_upper_end[stage_number];
func_ellshell_elems_along_circum_per_div = n_circum_elems_per_div;
func_ellshell_elems_along_length_per_div = n_elems_along_length_per_div_cylinder1_upper_end[stage_number];
func_ellshell_local_y_axis[2] = 0;

func_ellshell_local_y_axis[0] = 0;
func_ellshell_local_y_axis[1] = 1;
func_ellshell_local_y_axis[2] = 0;

Call ELLIPTICSHELL;

Printf("FESystem_Message: GEOMETRY_END STAGE_1_CYLINDER_1_END_1");

//
// clyinder1
//

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_1_CYLINDER_1");

func_cylinder_end1_pointID = 2;
func_cylinder_end2_pointID = 4;
func_cylinder_end1_radius = stage_radius_DV[stage_number];
func_cylinder_end2_radius = stage_radius_DV[stage_number];
func_cylinder_n_divs_along_circumference = n_circumference_div;
func_cylinder_n_divs_along_length = n_divs_along_length_cylinder1[stage_number];
func_cylinder_elems_along_circum_per_div = n_circum_elems_per_div ;
func_cylinder_elems_along_length_per_div = n_elems_along_length_per_div_cylinder1[stage_number];
func_cylinder_local_y_axis[2] = 0;

func_cylinder_local_y_axis[0] = 0;
func_cylinder_local_y_axis[1] = 1;
func_cylinder_local_y_axis[2] = 0;

Call CYLINDER;

Printf("FESystem_Message: GEOMETRY_END STAGE_1_CYLINDER_1");

// 
// clyinder1 end 2
// 

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_1_CYLINDER_1_END_2");

func_ellshell_end1_pointID = 5;
func_ellshell_end2_pointID = 4;
func_ellshell_end2_radius = stage_radius_DV[stage_number];
func_ellshell_n_divs_along_circumference = n_circumference_div;
func_ellshell_n_angular_divs_along_length = n_divs_along_length_cylinder1_lower_end[stage_number];
func_ellshell_elems_along_circum_per_div = n_circum_elems_per_div;
func_ellshell_elems_along_length_per_div = n_elems_along_length_per_div_cylinder1_lower_end[stage_number];
func_ellshell_local_y_axis[2] = 0;

func_ellshell_local_y_axis[0] = 0;
func_ellshell_local_y_axis[1] = 1;
func_ellshell_local_y_axis[2] = 0;

Call ELLIPTICSHELL;

Printf("FESystem_Message: GEOMETRY_END STAGE_1_CYLINDER_1_END_2");


// 
// skirt
// 

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_1_INTERCYLINDER_SKIRT");

func_cylinder_end1_pointID = 4;
func_cylinder_end2_pointID = 6;
func_cylinder_end1_radius = stage_radius_DV[stage_number];
func_cylinder_end2_radius = stage_radius_DV[stage_number];
func_cylinder_n_divs_along_circumference = n_circumference_div;
func_cylinder_n_divs_along_length = n_divs_along_length_intertank_skirt[stage_number];
func_cylinder_elems_along_circum_per_div = n_circum_elems_per_div ;
func_cylinder_elems_along_length_per_div = n_elems_along_length_per_div_intertank_skirt[stage_number];
func_cylinder_local_y_axis[2] = 0;

func_cylinder_local_y_axis[0] = 0;
func_cylinder_local_y_axis[1] = 1;
func_cylinder_local_y_axis[2] = 0;

Call CYLINDER;

Printf("FESystem_Message: GEOMETRY_END STAGE_1_INTERCYLINDER_SKIRT");

//
// clyinder2 end 1
// 

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_1_CYLINDER_2_END_1");

func_ellshell_end1_pointID = 7;
func_ellshell_end2_pointID = 6;
func_ellshell_end2_radius = stage_radius_DV[stage_number];
func_ellshell_n_divs_along_circumference = n_circumference_div;
func_ellshell_n_angular_divs_along_length = n_divs_along_length_cylinder2_upper_end[stage_number];
func_ellshell_elems_along_circum_per_div = n_circum_elems_per_div;
func_ellshell_elems_along_length_per_div = n_elems_along_length_per_div_cylinder2_upper_end[stage_number];
func_ellshell_local_y_axis[2] = 0;

func_ellshell_local_y_axis[0] = 0;
func_ellshell_local_y_axis[1] = 1;
func_ellshell_local_y_axis[2] = 0;

Call ELLIPTICSHELL;

Printf("FESystem_Message: GEOMETRY_END STAGE_1_CYLINDER_2_END_1");

//
// clyinder2
//

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_1_CYLINDER_2");

func_cylinder_end1_pointID = 6;
func_cylinder_end2_pointID = 8;
func_cylinder_end1_radius = stage_radius_DV[stage_number];
func_cylinder_end2_radius = stage_radius_DV[stage_number];
func_cylinder_n_divs_along_circumference = n_circumference_div;
func_cylinder_n_divs_along_length = n_divs_along_length_cylinder2[stage_number];
func_cylinder_elems_along_circum_per_div = n_circum_elems_per_div ;
func_cylinder_elems_along_length_per_div = n_elems_along_length_per_div_cylinder2[stage_number];
func_cylinder_local_y_axis[2] = 0;

func_cylinder_local_y_axis[0] = 0;
func_cylinder_local_y_axis[1] = 1;
func_cylinder_local_y_axis[2] = 0;

Call CYLINDER;

Printf("FESystem_Message: GEOMETRY_END STAGE_1_CYLINDER_2");

// 
// clyinder2 end 2
// 

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_1_CYLINDER_2_END_2");

func_ellshell_end1_pointID = 9;
func_ellshell_end2_pointID = 8;
func_ellshell_end2_radius = stage_radius_DV[stage_number];
func_ellshell_n_divs_along_circumference = n_circumference_div;
func_ellshell_n_angular_divs_along_length = n_divs_along_length_cylinder2_lower_end[stage_number];
func_ellshell_elems_along_circum_per_div = n_circum_elems_per_div;
func_ellshell_elems_along_length_per_div = n_elems_along_length_per_div_cylinder2_lower_end[stage_number];
func_ellshell_local_y_axis[2] = 0;

func_ellshell_local_y_axis[0] = 0;
func_ellshell_local_y_axis[1] = 1;
func_ellshell_local_y_axis[2] = 0;

Call ELLIPTICSHELL;

Printf("FESystem_Message: GEOMETRY_END STAGE_1_CYLINDER_2_END_2");


// 
// skirt
// 

Printf("FESystem_Message: GEOMETRY_BEGIN INTERSTAGE_SKIRT");

func_cylinder_end1_pointID = 8;
func_cylinder_end2_pointID = 10;
func_cylinder_end1_radius = stage_radius_DV[stage_number];
func_cylinder_end2_radius = stage_radius_DV[stage_number+1];
func_cylinder_n_divs_along_circumference = n_circumference_div;
func_cylinder_n_divs_along_length = n_divs_along_length_interstage_skirt;
func_cylinder_elems_along_circum_per_div = n_circum_elems_per_div ;
func_cylinder_elems_along_length_per_div = n_elems_along_length_per_div_interstage_skirt;
func_cylinder_local_y_axis[2] = 0;

func_cylinder_local_y_axis[0] = 0;
func_cylinder_local_y_axis[1] = 1;
func_cylinder_local_y_axis[2] = 0;

Call CYLINDER;

Printf("FESystem_Message: GEOMETRY_END INTERSTAGE_SKIRT");


//
// stage 2
// 

stage_number = 1;

//
// clyinder1 end 1
// 

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_2_CYLINDER_1_END_1");

func_ellshell_end1_pointID = 11;
func_ellshell_end2_pointID = 10;
func_ellshell_end2_radius = stage_radius_DV[stage_number];
func_ellshell_n_divs_along_circumference = n_circumference_div;
func_ellshell_n_angular_divs_along_length = n_divs_along_length_cylinder1_upper_end[stage_number];
func_ellshell_elems_along_circum_per_div = n_circum_elems_per_div;
func_ellshell_elems_along_length_per_div = n_elems_along_length_per_div_cylinder1_upper_end[stage_number];
func_ellshell_local_y_axis[2] = 0;

func_ellshell_local_y_axis[0] = 0;
func_ellshell_local_y_axis[1] = 1;
func_ellshell_local_y_axis[2] = 0;

Call ELLIPTICSHELL;

Printf("FESystem_Message: GEOMETRY_END STAGE_2_CYLINDER_1_END_1");

//
// clyinder1
//

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_2_CYLINDER_1");

func_cylinder_end1_pointID = 10;
func_cylinder_end2_pointID = 12;
func_cylinder_end1_radius = stage_radius_DV[stage_number];
func_cylinder_end2_radius = stage_radius_DV[stage_number];
func_cylinder_n_divs_along_circumference = n_circumference_div;
func_cylinder_n_divs_along_length = n_divs_along_length_cylinder1[stage_number];
func_cylinder_elems_along_circum_per_div = n_circum_elems_per_div ;
func_cylinder_elems_along_length_per_div = n_elems_along_length_per_div_cylinder1[stage_number];
func_cylinder_local_y_axis[2] = 0;

func_cylinder_local_y_axis[0] = 0;
func_cylinder_local_y_axis[1] = 1;
func_cylinder_local_y_axis[2] = 0;

Call CYLINDER;

Printf("FESystem_Message: GEOMETRY_END STAGE_2_CYLINDER_1");

// 
// clyinder1 end 2
// 

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_2_CYLINDER_1_END_2");

func_ellshell_end1_pointID = 13;
func_ellshell_end2_pointID = 12;
func_ellshell_end2_radius = stage_radius_DV[stage_number];
func_ellshell_n_divs_along_circumference = n_circumference_div;
func_ellshell_n_angular_divs_along_length = n_divs_along_length_cylinder1_lower_end[stage_number];
func_ellshell_elems_along_circum_per_div = n_circum_elems_per_div;
func_ellshell_elems_along_length_per_div = n_elems_along_length_per_div_cylinder1_lower_end[stage_number];
func_ellshell_local_y_axis[2] = 0;

func_ellshell_local_y_axis[0] = 0;
func_ellshell_local_y_axis[1] = 1;
func_ellshell_local_y_axis[2] = 0;

Call ELLIPTICSHELL;

Printf("FESystem_Message: GEOMETRY_END STAGE_2_CYLINDER_1_END_2");

// 
// skirt
// 

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_2_INTERCYLINDER_SKIRT");

func_cylinder_end1_pointID = 12;
func_cylinder_end2_pointID = 14;
func_cylinder_end1_radius = stage_radius_DV[stage_number];
func_cylinder_end2_radius = stage_radius_DV[stage_number];
func_cylinder_n_divs_along_circumference = n_circumference_div;
func_cylinder_n_divs_along_length = n_divs_along_length_intertank_skirt[stage_number];
func_cylinder_elems_along_circum_per_div = n_circum_elems_per_div ;
func_cylinder_elems_along_length_per_div = n_elems_along_length_per_div_intertank_skirt[stage_number];
func_cylinder_local_y_axis[2] = 0;

func_cylinder_local_y_axis[0] = 0;
func_cylinder_local_y_axis[1] = 1;
func_cylinder_local_y_axis[2] = 0;

Call CYLINDER;

Printf("FESystem_Message: GEOMETRY_END STAGE_2_INTERCYLINDER_SKIRT");


//
// clyinder2 end 1
// 

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_2_CYLINDER_2_END_1");

func_ellshell_end1_pointID = 15;
func_ellshell_end2_pointID = 14;
func_ellshell_end2_radius = stage_radius_DV[stage_number];
func_ellshell_n_divs_along_circumference = n_circumference_div;
func_ellshell_n_angular_divs_along_length = n_divs_along_length_cylinder2_upper_end[stage_number];
func_ellshell_elems_along_circum_per_div = n_circum_elems_per_div;
func_ellshell_elems_along_length_per_div = n_elems_along_length_per_div_cylinder2_upper_end[stage_number];
func_ellshell_local_y_axis[2] = 0;

func_ellshell_local_y_axis[0] = 0;
func_ellshell_local_y_axis[1] = 1;
func_ellshell_local_y_axis[2] = 0;

Call ELLIPTICSHELL;

Printf("FESystem_Message: GEOMETRY_END STAGE_2_CYLINDER_2_END_1");

//
// clyinder2
//

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_2_CYLINDER_2");

func_cylinder_end1_pointID = 14;
func_cylinder_end2_pointID = 16;
func_cylinder_end1_radius = stage_radius_DV[stage_number];
func_cylinder_end2_radius = stage_radius_DV[stage_number];
func_cylinder_n_divs_along_circumference = n_circumference_div;
func_cylinder_n_divs_along_length = n_divs_along_length_cylinder2[stage_number];
func_cylinder_elems_along_circum_per_div = n_circum_elems_per_div ;
func_cylinder_elems_along_length_per_div = n_elems_along_length_per_div_cylinder2[stage_number];
func_cylinder_local_y_axis[2] = 0;

func_cylinder_local_y_axis[0] = 0;
func_cylinder_local_y_axis[1] = 1;
func_cylinder_local_y_axis[2] = 0;

Call CYLINDER;

Printf("FESystem_Message: GEOMETRY_END STAGE_2_CYLINDER_2");

// 
// clyinder2 end 2
// 

Printf("FESystem_Message: GEOMETRY_BEGIN STAGE_2_CYLINDER_2_END_2");

func_ellshell_end1_pointID = 17;
func_ellshell_end2_pointID = 16;
func_ellshell_end2_radius = stage_radius_DV[stage_number];
func_ellshell_n_divs_along_circumference = n_circumference_div;
func_ellshell_n_angular_divs_along_length = n_divs_along_length_cylinder2_lower_end[stage_number];
func_ellshell_elems_along_circum_per_div = n_circum_elems_per_div;
func_ellshell_elems_along_length_per_div = n_elems_along_length_per_div_cylinder2_lower_end[stage_number];
func_ellshell_local_y_axis[2] = 0;

func_ellshell_local_y_axis[0] = 0;
func_ellshell_local_y_axis[1] = 1;
func_ellshell_local_y_axis[2] = 0;

Call ELLIPTICSHELL;

Printf("FESystem_Message: GEOMETRY_END STAGE_2_CYLINDER_2_END_2");

