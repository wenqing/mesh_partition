#ifndef Mesh_INC
#define Mesh_INC

#include<string>
#include<vector>
#include<iostream>

//------------------------------------------------------
//   Topology declartion of geometrical element.
//   WW. 06.2005
//------------------------------------------------------

#include "Elem.h"

namespace Mesh_Group
{
class Node;

struct MeshPartConfig
{
   std::string fname;  /// file name
   std::string fpath;  /// file path
   std::string mat_fname; /// file name of material data (permeability and porosity)
   int num_parts;  /// number of parittions
   bool osdom;     /// flag to output subdomain mesh in vtk
   bool out_renum_gsmsh;  /// flag to output the node index renumbered ogs mesh
   bool is_vtk_out; ;  /// flag to output the node index renumbered vtk mesh
   bool binary_output;  /// output binary file of the partitioned mesh
   bool out_cct;  /// output of exended partitions for the flux corrected transport
};
/*!
   \class Mesh
*/
class Mesh
{
   public:
      Mesh(bool quad = false);
      ~Mesh();


      void ReadGrid(std::istream& is, const bool high_order);
      void ReadGridGeoSys(std::istream& is, const bool high_order);

      void Write2METIS(std::ostream& os);

      void WriteVTK_Nodes(std::ostream& os);

      void WriteVTK_Head(std::ostream& os, const size_t number_of_nodes);

      void ConstructSubDomain_by_Elements(const std::string fname,  const int num_parts, const bool osdom);
      void ConstructSubDomain_by_Nodes(const MeshPartConfig mpc);

      void ConnectedNodes(bool quadratic);
      void ConnectedElements2Node(bool quadratic=false);

      void ConstructGrid();
      void GenerateHighOrderNodes();

      void setOrder(const bool is_quad)
      {
         useQuadratic = is_quad;
      }
   private:
      // The following can be members of grid class
      MyInt NodesNumber_Linear;
      MyInt NodesNumber_Quadratic;
      bool useQuadratic;
      bool axisymmetry;

      bool ouput_vtk_part_info;
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

      // All surface feces
#ifdef BUILD_MESH_FACE
      std::vector<Elem*> face_vector;
#endif
      // All elements
      std::vector<Elem*> elem_vector;

      MyInt msh_no_line;
      MyInt msh_no_quad;
      MyInt msh_no_hexs;
      MyInt msh_no_tris;
      MyInt msh_no_tets;
      MyInt msh_no_pris;
      MyInt msh_no_pyra;
      int msh_max_dim;

      void WriteVTK_Nodes(std::ostream& os, std::vector<Node*>& nod_vec, const size_t start, const size_t end);

      void WriteVTK_Elements_of_Subdomain(std::ostream& os, std::vector<Elem*>& ele_vec,
                                          const int sbd_index);
      void writeSubdomainElementsVTK(std::ostream& os, const std::vector<Elem*>& sdom_elems,
                                     const std::vector<Elem*>& sdom_elems_ghost, const int sbd_index);

      void fillNodeVector4BinaryOuput(const std::vector<Node*> &sdom_nodes,
                                      std::vector<Node_Str> &sdom_nodes4bin,
                                      const size_t start, const size_t end,  size_t &counter);

      void writeSubDomainNodes(std::ostream& os, const std::vector<Node*>& sdom_nodes, const size_t start, const size_t end);

      friend class Elem;
};

}

#endif
