Function CYLINDER


// This is a function to create a cylinder and its surface in the gmsh
// The function allows for different radii at the two ends of the cylinder
// 
// For the cylinder, the local axis are such that z-axis forms the longitudinal axis, and
// the local y-axis is given as an input. This is used to create the first point for the circumferencial 
// points
// 
//
// Following are the input parameters that the function needs

// func_cylinder_end1_pointID 
// func_cylinder_end2_pointID
// func_cylinder_end1_radius
// func_cylinder_end2_radius
// func_cylinder_n_divs_along_circumference
// func_cylinder_n_divs_along_length
// func_cylinder_elems_along_circum_per_div
// func_cylinder_elems_along_length_per_div
// func_cylinder_local_y_axis
//
// following are the return variables
// func_cylinder_n_points
// func_cylinder_n_axial_lines
// func_cylinder_n_circum_lines
// func_cylinder_n_surfaces
// func_cylinder_point_ID_vector
// func_cylinder_axial_line_ID_vector
// func_cylinder_circum_line_ID_vector
// func_cylinder_surface_ID_vector
// 

// get the coordinates of the first and the last point
local_end1_coord[] = Point{func_cylinder_end1_pointID};
local_end2_coord[] = Point{func_cylinder_end2_pointID};

// vector of point IDs
local_n_points = (func_cylinder_n_divs_along_length + 1) * (func_cylinder_n_divs_along_circumference +1);
local_point_ID_vector[local_n_points] = 0;

// vector of line IDs
local_n_circum_lines = func_cylinder_n_divs_along_circumference * (func_cylinder_n_divs_along_length + 1);
local_circum_line_ID_vector[local_n_circum_lines] = 0;

local_n_axial_lines = func_cylinder_n_divs_along_circumference * func_cylinder_n_divs_along_length ;
local_axial_line_ID_vector[local_n_axial_lines] = 0;

// vector of surface IDs
local_n_surfaces = func_cylinder_n_divs_along_circumference * func_cylinder_n_divs_along_length;
local_surf_ID_vector[local_n_surfaces] = 0;


// this will keep track of the element in the point ID vector to write to
local_point_coord_tmp[2] = 0;

// calculate the position vector from end 1 to 2
local_vector_end1_to_end2[] = Point{func_cylinder_end2_pointID};
local_point_coord_tmp[] = Point{func_cylinder_end1_pointID};
local_vector_end1_to_end2[0] -= local_point_coord_tmp[0];
local_vector_end1_to_end2[1] -= local_point_coord_tmp[1];
local_vector_end1_to_end2[2] -= local_point_coord_tmp[2];


// calculate the increment in rotation for each circumferencial point
local_rotate_incr = 2*Pi/ func_cylinder_n_divs_along_circumference;
// claculate the translation unit vector for the circumferencial point
local_circum_translate[2] = 0;
local_circum_translate_unit[2] = 0;
local_tmp_value = func_cylinder_local_y_axis[0] * func_cylinder_local_y_axis[0] ;
local_tmp_value += func_cylinder_local_y_axis[1] * func_cylinder_local_y_axis[1] ;
local_tmp_value += func_cylinder_local_y_axis[2] * func_cylinder_local_y_axis[2] ;
local_tmp_value = Sqrt(local_tmp_value);
local_circum_translate_unit[0] = func_cylinder_local_y_axis[0] / local_tmp_value;
local_circum_translate_unit[1] = func_cylinder_local_y_axis[1] / local_tmp_value;
local_circum_translate_unit[2] = func_cylinder_local_y_axis[2] / local_tmp_value;

local_incr = 0;
local_circum_line_incr = 0;
local_axial_line_incr = 0;
// loop over each station and create all the necessary points
For local_axial_station In {0 : func_cylinder_n_divs_along_length}

	local_fraction = local_axial_station/ func_cylinder_n_divs_along_length;
	// calculate the local_radius
	local_station_radius = func_cylinder_end1_radius + (func_cylinder_end2_radius - func_cylinder_end1_radius) * local_fraction;

	// create the axial point if necessary and record its point ID
	local_create_axial_point = 1;
	If (local_axial_station == 0 )
		local_point_ID_vector[local_incr] = func_cylinder_end1_pointID;
		local_create_axial_point = 0;
	EndIf

	If (local_axial_station == func_cylinder_n_divs_along_length)
		local_point_ID_vector[local_incr] = func_cylinder_end2_pointID;
		local_create_axial_point = 0;
	EndIf

	If (local_create_axial_point == 1)
		local_point_ID_vector[local_incr] = newp;
		local_point_coord_tmp[0] =  local_vector_end1_to_end2[0] * local_fraction;
		local_point_coord_tmp[1] =  local_vector_end1_to_end2[1] * local_fraction;
		local_point_coord_tmp[2] =  local_vector_end1_to_end2[2] * local_fraction;
		local_pointID[] = Translate {local_point_coord_tmp[0],local_point_coord_tmp[1],local_point_coord_tmp[2]} { Duplicata{ Point{func_cylinder_end1_pointID};} } ;
		local_point_ID_vector[local_incr] = local_pointID[0];
	EndIf

	local_incr += 1;
	
	local_axial_point_ID = local_point_ID_vector[local_incr-1];
	
	For local_circum_station In {0 : func_cylinder_n_divs_along_circumference-1}
		
		
		// get the point ID of the axial point for this station
		local_rotation_base[] = Point{local_axial_point_ID};
		local_previous_point = local_point_ID_vector[local_incr - 1];
		
		// if this is the first circumferencial point, then translate the axial point to get 
		// the first circumferencial point
		If (local_circum_station == 0)
			local_point_coord_tmp[0] = local_circum_translate_unit[0] * local_station_radius;
			local_point_coord_tmp[1] = local_circum_translate_unit[1] * local_station_radius;
			local_point_coord_tmp[2] = local_circum_translate_unit[2] * local_station_radius;
			local_pointID[] = Translate {local_point_coord_tmp[0],local_point_coord_tmp[1],local_point_coord_tmp[2]} { Duplicata{ Point{local_previous_point};} } ;
			local_circum_first_point = local_pointID[0];

		EndIf
		
		// otherwise, rotate the previous circumferencial point
		If (local_circum_station != 0)
			
			local_pointID[] = Rotate {{local_vector_end1_to_end2[0],local_vector_end1_to_end2[1],local_vector_end1_to_end2[2]}, {local_rotation_base[0],local_rotation_base[1],local_rotation_base[2]}, local_rotate_incr} { Duplicata{ Point{local_previous_point};} } ;

			
			// create an arc connecting the previous point to this point
			local_lineID = newl;
			local_circum_line_ID_vector[local_circum_line_incr] = local_lineID;
			Circle(local_lineID) = {local_previous_point, local_axial_point_ID, local_pointID[0]};
			local_circum_line_incr += 1;
			Transfinite Line {local_lineID} = func_cylinder_elems_along_circum_per_div + 1;

			
		EndIf
		

		// finally, if the last point is being created, then a line connecting the last and the first point need to be created
		If (local_circum_station == (func_cylinder_n_divs_along_circumference -1))
			
			local_lineID = newl;
			local_circum_line_ID_vector[local_circum_line_incr] = local_lineID;
			Circle(local_lineID) = {local_pointID[0], local_axial_point_ID, local_circum_first_point};
			local_circum_line_incr += 1;
			Transfinite Line {local_lineID} = func_cylinder_elems_along_circum_per_div + 1;
						
		EndIf
		
		
		// If this is not the first axial station, create the axial lines 
		If (local_axial_station != 0)
		
			local_lineID = newl;
			local_point_ID_previous_station = local_point_ID_vector[(local_incr - func_cylinder_n_divs_along_circumference - 1)];
			Line(local_lineID) = {local_point_ID_previous_station, local_pointID[0]};
			local_axial_line_ID_vector[local_axial_line_incr] = local_lineID;
			Transfinite Line {local_lineID} = func_cylinder_elems_along_length_per_div + 1;
			local_axial_line_incr += 1;
		
		EndIf
		
		// insert the point ID
		local_point_ID_vector[local_incr] = local_pointID[0];
		// increment the index
		local_incr += 1;
				
	EndFor
	
EndFor


// now iterate over the stations and create the line loops and surfaces

local_surf_incr = 0;
For local_axial_station In {0 : func_cylinder_n_divs_along_length-1}
	For local_circum_station In {0 : func_cylinder_n_divs_along_circumference-1}
	
		// get the IDs of the 4 lines that will make up this surface
		local_line_index = func_cylinder_n_divs_along_circumference  * local_axial_station + local_circum_station;
		local_lineID = newl;
		
		If (local_circum_station != (func_cylinder_n_divs_along_circumference-1))
			Line Loop(local_lineID) = { local_circum_line_ID_vector[local_line_index], local_axial_line_ID_vector[local_line_index + 1], -local_circum_line_ID_vector[local_line_index + func_cylinder_n_divs_along_circumference], -local_axial_line_ID_vector[local_line_index]};
		EndIf
		
		If (local_circum_station == (func_cylinder_n_divs_along_circumference-1))
			Line Loop(local_lineID) = { local_circum_line_ID_vector[local_line_index], local_axial_line_ID_vector[local_line_index - func_cylinder_n_divs_along_circumference+1], -local_circum_line_ID_vector[local_line_index + func_cylinder_n_divs_along_circumference], -local_axial_line_ID_vector[local_line_index]};		
		EndIf

		local_surfaceID = news;
		Ruled Surface(local_surfaceID) = {local_lineID};
		
		// find the point IDs for the four corners of the surface and create a transfinite surface
		local_point_ID_index = local_axial_station * (func_cylinder_n_divs_along_circumference + 1) + 1 + local_circum_station;
		local_pointID_station1 = local_point_ID_vector[local_point_ID_index];
		local_pointID_station2 = local_point_ID_vector[local_point_ID_index + func_cylinder_n_divs_along_circumference + 1];
		
		If (local_circum_station != (func_cylinder_n_divs_along_circumference-1))
			Transfinite Surface{local_surfaceID} = {local_pointID_station1, (local_pointID_station1+1), (local_pointID_station2+1) ,local_pointID_station2};
		EndIf
		
		If (local_circum_station == (func_cylinder_n_divs_along_circumference-1))
			Transfinite Surface{local_surfaceID} = {local_pointID_station1, (local_pointID_station1-func_cylinder_n_divs_along_circumference +1), (local_pointID_station2-func_cylinder_n_divs_along_circumference+1) ,local_pointID_station2};
		EndIf
		
		Recombine Surface {local_surfaceID};
		
		// insert the surface ID;
		local_surf_ID_vector[local_surf_incr] = local_surfaceID;
		local_surf_incr += 1;
		
	EndFor
EndFor


// now initialize the return variables and copy the data into them
//local_point_ID_vector[local_n_points] = 0;
//local_circum_line_ID_vector[local_n_circum_lines] = 0;
//local_axial_line_ID_vector[local_n_axial_lines] = 0;
//local_surf_ID_vector[local_n_surfaces] = 0;

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



Return
