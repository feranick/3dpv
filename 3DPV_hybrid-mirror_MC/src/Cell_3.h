#ifndef CELL_3_H_
#define CELL_3_H_

// probably don't need to inlcude here

struct Position
{
	double x;
	double y;
	double z;
};

struct Point
{
	double pt[3];
};

struct Vector
{
	double vec[3];
};

struct Refl_vector
{
	Vector ref_dir;     // direction of reflected ray
	double ref_power;   // power carried here (was calculated as releflection coefficient times solar constant or whatever
	                    // power was multiplied by reflection coefficient
};

struct Plane
{
	// define by an 3 size array of pointers to CCW-defined points
	// and its norm
	Point *points[3];
	Vector Normal;
};

struct Edge
{
	int vertexIndex[2];
	int triangleIndex[2];
	bool is_shadower_edge;         // true if one of its triangles faces sun and other doesn't
	bool extrusion_intersects[2]; /*= {false, false}*/; // true when extruded rays are intersecting the current triangle
												  // note there is unfortunate redundancy here.... since each
												  // ray will be shared by two objects of type Edge :[
	Plane *shadowPlane;
};

struct Triangle
{
	int index[3];                  // indices of the vertices for the cell
	double vertList[3*3];                     // list of the actual 3 vertices in world space (3 coordiantes each)
	//Vector_3 Normal;                         // the (yet "un-normalized") orthogonal vector to surface
	double Normal[3];
	double SunScalarProduct;                 // scalar product of this face's normal with a sun ray directed to cell
	double Area;                             // scalar corresponding to total surface area of this triangle
											 // note: area should be the same for ALL the analagous componets
											 // of every cell IF they are identical in shape and size
	bool is_facing_sun; /*= false*/;              // switched from dot product signs
	bool is_positive;			// if this triangle is the postiveone
	int numVerticesShadowed; /*= 0*/;             // how many of triangles vertices are shadowed by ANY shadow (0,1,2,3)
	bool is_shadowedIndex[3]; /*= {false, false, false}*/; // which of the three vertices are currently shadowed
	bool is_deactivated;
	
	//double vert2sun_dists[3]; // each timestep, calculate the distance from each vertex to sun
	//int least2sun_vert;  // which of these vertices are closest to sun

	int MainAxis; // 0=x, 1=y, 2=z
};

// CELL COMPONENTS CLASS
class Cell_3 {
public:
	Cell_3();                                // default constructor
	~Cell_3();                               // destructor
	
	//int setupCell(); // no longer a member function
	
	// buildedges will no longer be a member function
	//int BuildEdges(int triangleCount, const Triangle *triangleArray, Edge **edgeArray);

	int numVertices;                         // # of points actually in the cell
	int numTriangles;

	Position *vertices[8];                   // make room for up to 8 vertices
	Triangle *triangleList[12];               // array of pointers to each of the cell's triangle instances


//	static int cell_counter;                 // global class variable to hold # of cells instantiated
//	static int cell_shape;                   // common base shape to all instantiated cells
//	static double cell_radius;                  // common to cells: distance from position to vertices
												 // thus, all vertices lie on the same spherical shell (For now)

//	Position position;                        // 3 component cell position (x,y,z)
//	int numEdges;                           // how many total edges
//	int numShadowingEdges;                  // how many of its edges have extrusions
	// numShadowingEdges is also how many shadowplanes it has
//	Edge /* **edgeList*/ *edgeList[12];                        // pointer to location where pointer to array of edges will be returned
//	Plane *shadowPlaneList[12];           // array of pointers to planes containing extruding rays (size 12 is bigger than needed)
//	double rot_x;                            // theta (rotation angle about x axis)
//	double rot_y;                            // phi (rotation angle about y axis)
//	double rot_z;                             // mainly for testing
//	double rsize;                            // RADIUS, individual
	
//	void getNormal(double pn[3]);
};

//// EXTRA JUST FOR RAYTRACING/////////////////////
struct RAYTRI
{
	double org[3];
	double end[3];
	double dir[3];
	double v0[3],v1[3],v2[3];

	struct PLANE
	{
		double x, y, z, d;
		enum MAIN_AXIS { X, Y, Z };
		MAIN_AXIS type;
	};
	PLANE plane;
};
////////////////////////////////////////////////


#endif /*CELL_3_H_*/
