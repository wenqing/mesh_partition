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
     Mesh() {};
     ~Mesh();
   
     void ConstructGrid( const bool quadratic=false);
     void ConstructDomain(char *fname,  const int num_parts);
    
     void ReadGrid(std::istream& is = std::cin);
     void ReadGridGeoSys(std::istream& is = std::cin);

	 void Write2METIS(std::ostream& os); 
   private:
      // The following can be members of grid class
      long NodesNumber_Linear;
      long NodesNumber_Quadratic;

      // All nodes
      std::vector<Node*> NodesVector;
      // All edges
      std::vector<Edge*> EdgeVector;
      // All surface feces
      std::vector<Elem*> SurfaceFaces;
      //  All elements 
      std::vector<Elem*> ElementsVector;

};

}

#endif