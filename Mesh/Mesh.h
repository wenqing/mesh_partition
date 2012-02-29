#ifndef Mesh_INC
#define Mesh_INC

#include<string>
#include<vector>
#include<iostream>

//------------------------------------------------------ 
//   Topology declartion of geometrical element.
//   WW. 06.2005
//------------------------------------------------------ 

namespace Mesh_Group
{
class Node;
class Edge;
class Elem;


/*!
   \class Mesh
*/
class Mesh
{
   public:
     Mesh(bool quad = false);
     ~Mesh();
   
     void ConstructSubDomain_by_Elements(const std::string fname,  const int num_parts);
     void ConstructSubDomain_by_Nodes(const std::string fname,  const int num_parts, const bool is_quad);
    
     void ReadGrid(std::istream& is = std::cin);
     void ReadGridGeoSys(std::istream& is = std::cin);

	 void Write2METIS(std::ostream& os); 
    void ConnectedNodes(bool quadratic); 
    void ConnectedElements2Node(bool quadratic=false);

    void ConstructGrid();
    void GenerateHighOrderNodes();
   private:
      // The following can be members of grid class
      long NodesNumber_Linear;
      long NodesNumber_Quadratic;
    bool useQuadratic;
    bool axisymmetry;
   
    // Coordinate indicator
	// 1:  X component only
	// 12: Y component only
	// 13: Z component only
	// 2:  X, Y component
	// 23:  X, Z component
	// 3:  X, Y, Z component
    int coordinate_system; 
    int max_ele_dim; 

    // All nodes
    std::vector<Node*> node_vector;
    // All edges
    std::vector<Edge*> edge_vector;
    // All surface feces
    std::vector<Elem*> face_vector;
    // All elements 
    std::vector<Elem*> elem_vector;

    long msh_no_line;
    long msh_no_quad;
    long msh_no_hexs;
    long msh_no_tris;
    long msh_no_tets;
    long msh_no_pris;
    int msh_max_dim;   
};

}

#endif