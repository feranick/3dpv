#include "Cell_3.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "geometry.h"

//int Cell_3::cell_counter=0;                 // global class variable to hold # of cells instantiated
//int Cell_3::cell_shape=0;                   // common base shape to all instantiated cells
//double Cell_3::cell_radius=1;                  // common to cells: distance from position to vertices

Cell_3::Cell_3() { // default constructor
//	cell_counter++;
	
	// from header it is difficult to tell quickly, but
	// should be sure x,y,z are member variables that can be used this way
//	position.x=0;
//	position.y=0;
//	position.z=0;
	
//	rot_x=0;
//	rot_y=0;
//	rot_z=0;
//
//	rsize=2.0;
	
//	numShadowingEdges=0;
}

//Cell_3::Cell_3(Position * new_pos, double new_rotx, double new_roty, double new_rotz) { // constructor
////	cell_counter++;
//
//	position.x = new_pos->x;
//	position.y = new_pos->y;
//	position.z = new_pos->z;
//
//	rot_x = new_rotx;
//	rot_y = new_roty;
//	rot_z = new_rotz;
//
//	rsize=2.0;
//
//	numShadowingEdges=0;
//}


Cell_3::~Cell_3()   // default destructor
{  
//	--cell_counter;
}

//void Cell_3::getNormal(double pn[3]) {
//	double p1[3], p2[3], p3[3];
//	p1[0] = vertices[0]->x; p1[1] = vertices[0]->y; p1[2] = vertices[0]->z;
//	p2[0] = vertices[1]->x; p2[1] = vertices[1]->y; p2[2] = vertices[1]->z;
//	p3[0] = vertices[2]->x; p3[1] = vertices[2]->y; p3[2] = vertices[2]->z;
//
//	plane_exp_normal_3d(p1, p2, p3, pn);
//}
