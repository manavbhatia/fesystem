Function THICK_CYLINDER

// This function builds up on the  CYLINDER function in the sense that it uses the 
// cylinder function to create the cylinderical points and then creates the subsurfaces 
// to make the solid elements in the thickness of the cylinder

//
// Following are the input parameters that the function needs

// func_thick_cylinder_end1_pointID 
// func_thick_cylinder_end2_pointID
// func_thick_cylinder_end1_radius
// func_thick_cylinder_end2_radius
// func_thick_cylinder_n_divs_along_circumference
// func_thick_cylinder_n_divs_along_length
// func_thick_cylinder_elems_along_circum_per_div
// func_thick_cylinder_elems_along_length_per_div
// func_thick_cylinder_local_y_axis
//
// following are the return variables
// func_thick_cylinder_n_points
// func_thick_cylinder_n_axial_lines
// func_thick_cylinder_n_circum_lines
// func_thick_cylinder_n_surfaces
// func_thick_cylinder_point_ID_vector
// func_thick_cylinder_axial_line_ID_vector
// func_thick_cylinder_circum_line_ID_vector
// func_thick_cylinder_surface_ID_vector
// 
//
//


func_cylinder_end1_pointID = func_thick_cylinder_end1_pointID;
func_cylinder_end2_pointID = func_thick_cylinder_end2_pointID;
func_cylinder_end1_radius = func_thick_cylinder_end1_radius;
func_cylinder_end2_radius = ;
func_cylinder_n_divs_along_circumference = n_circumference_div;
func_cylinder_n_divs_along_length = n_divs_along_length_cylinder1[stage_number];
func_cylinder_elems_along_circum_per_div = n_circum_elems_per_div ;
func_cylinder_elems_along_length_per_div = n_elems_along_length_per_div_cylinder1[stage_number];
func_cylinder_local_y_axis[2] = 0;

func_cylinder_local_y_axis[0] = 0;
func_cylinder_local_y_axis[1] = 1;
func_cylinder_local_y_axis[2] = 0;


func_cylinder_n_points =  local_n_points;
func_cylinder_n_axial_lines = local_n_axial_lines;
func_cylinder_n_circum_lines = local_n_circum_lines;
func_cylinder_n_surfaces = local_n_surfaces;

func_cylinder_point_ID_vector[func_cylinder_n_points] = 0;
func_cylinder_circum_line_ID_vector[func_cylinder_n_circum_lines] = 0;
func_cylinder_axial_line_ID_vector[func_cylinder_n_axial_lines] = 0;
func_cylinder_surface_ID_vector[func_cylinder_n_surfaces] = 0;

// now copy the variables
For local_incr In {0: func_cylinder_n_points-1}
	func_cylinder_point_ID_vector[local_incr] = local_point_ID_vector[local_incr];
EndFor

For local_incr In {0: func_cylinder_n_axial_lines-1}
	func_cylinder_axial_line_ID_vector[local_incr] = local_axial_line_ID_vector[local_incr]; 
EndFor

For local_incr In {0: func_cylinder_n_circum_lines-1}
	func_cylinder_circum_line_ID_vector[local_incr] = local_circum_line_ID_vector[local_incr]; 
EndFor

For local_incr In {0: func_cylinder_n_surfaces-1}
	func_cylinder_surface_ID_vector[local_incr] = local_surf_ID_vector[local_incr]; 
EndFor





RETURN
