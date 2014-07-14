#ifndef node_INC
#define node_INC

#include<vector>

#include "Grain.h"

//------------------------------------------------------
//   Topology declartion of geometrical element.
//   WW. 06.2005
//   WW. 02.2012
//------------------------------------------------------
namespace Mesh_Group
{
class Edge;
class Elem;
class Mesh;

typedef struct
{
   MyInt id;
   MyInt id_l;  // My node ID for linear elements
   double x;
   double y;
   double z;
} Node_Str;

//2.  Node declaration
class Node:public Grain
{
   public:
      Node(const MyInt Index):Grain(Index), global_index(Index), global_index_l( -1 )
      {
         Coordinate = new double[3];
         index_org = index;
      }
      Node(const MyInt Index, const double x,
           const double y, const double z=0.0);
      ~Node()
      {
         delete [] Coordinate;
         Coordinate = NULL;
         ElementsRelated.resize(0);
      }

      // Operator
      void operator = (const Node& n);
      bool operator == (const Node & n);

      // Change members;
      // By component
      void setX(const double argX)
      {
         Coordinate[0] = argX;
      }
      void setY(const double argY)
      {
         Coordinate[1] = argY;
      }
      void setZ(const double argZ)
      {
         Coordinate[2] = argZ;
      }
      void SetCoordinates(const double* argCoord);

      // Access to members
      double X() const
      {
         return Coordinate[0];
      }
      double Y() const
      {
         return Coordinate[1];
      }
      double Z() const
      {
         return Coordinate[2];
      }
      double *getCoordinates() const
      {
         return  Coordinate;
      }

      void setLocalIndex(const MyInt l_index)
      {
         local_index = l_index;
      }

      MyInt getLocalIndex() const
      {
         return local_index;
      }

      MyInt getGlobalIndex4LinearElement() const
      {
         return global_index_l;
      }

      void setGlobalIndex4LinearElement(const MyInt node_id)
      {
         global_index_l = node_id;
      }

      // Output
      void Write(std::ostream& os = std::cout) const;
      void WriteBIN(std::ostream& os = std::cout) const;
      void WriteCoordinates(std::ostream& os = std::cout) const;

   private:
      double *Coordinate;
      MyInt local_index; // For domain decomposition
      MyInt global_index;
      MyInt global_index_l; // global index for nodes that make up linear elements
      size_t index_org; // For quad elements in ddc
      std::vector<size_t>  ElementsRelated;
      std::vector<size_t>  NodesRelated;
      friend class Mesh_Group::Edge;
      friend class Mesh_Group::Elem;
      friend class Mesh_Group::Mesh;

};

}

#endif
