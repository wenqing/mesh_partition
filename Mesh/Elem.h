#ifndef Elem_INC
#define Elem_INC

#include<vector>

#include "Grain.h"
#include "Node.h"

//------------------------------------------------------
//   Topology declartion of geometrical element.
//   WW. 06.2005
//   WW. 02.2012
//------------------------------------------------------

namespace Mesh_Group
{
class Mesh;
enum ElemType {line, quadri, hex, tri, tet, prism, pyramid};

//3.  Element declaration
class Elem:public Grain
{
   public:
      explicit Elem(const int Index) : Grain(Index),
         quadratic(false), sub_dom(0), no_faces_on_surface(0),
         PatchIndex(0), ele_Type(line), locnodes_index(NULL), nodes(NULL), neighbors(NULL)
      {}

      ~Elem();

      int getNodesNumber() const
      {
         switch(ele_Type)
         {
            case line:
               return 2;
            case quadri:
               return 4;
            case hex:
               return 8;
            case tri:
               return 3;
            case tet:
               return 4;
            case prism:
               return 6;
            case pyramid:
               return 5;
            default:
               break;
         }
         return 0;
      }

      int getNodesNumberHQ() const
      {
         switch(ele_Type)
         {
            case line:
               return 3;
            case quadri:
               return 9;
            case hex:
               return 20;
            case tri:
               return 6;
            case tet:
               return 10;
            case prism:
               return 15;
            case pyramid:
               return 14;
            default:
               break;
         }
         return 0;
      }

      int getNodesNumber(bool quad) const
      {
         if(quad) return getNodesNumberHQ();
         else return getNodesNumber();
      }

      int getEdgesNumber() const
      {
         switch(ele_Type)
         {
            case line:
               return 0;
            case quadri:
               return 4;
            case hex:
               return 12;
            case tri:
               return 3;
            case tet:
               return 6;
            case prism:
               return 9;
            case pyramid:
               return 8;
            default:
               break;
         }
         return 0;
      }

      int getFacesNumber() const
      {
         switch(ele_Type)
         {
            case line:
               return 2;
            case quadri:
               return 4;
            case hex:
               return 6;
            case tri:
               return 3;
            case tet:
               return 4;
            case prism:
               return 5;
            case pyramid:
               return 5;
            default:
               break;
         }
         return 0;
      }

      int Dim() const
      {
         switch(ele_Type)
         {
            case line:
               return 1;
            case quadri:
               return 2;
            case hex:
               return 3;
            case tri:
               return 2;
            case tet:
               return 3;
            case prism:
               return 3;
            case pyramid:
               return 3;
            default:
               break;
         }
         return 0;
      }

      int getElementType() const
      {
         return ele_Type;
      }

      int getPatchIndex() const
      {
         return PatchIndex;
      }

      std::string getName() const;

      // Nodes
      void getNodeIndeces(MyInt *node_index, const bool quad = false) const
      {
         for (int i=0; i< getNodesNumber(quad); i++)
            node_index[i]= nodes[i]->index;
      }

      MyInt getNodeIndex(const int loc_lndex) const
      {
         return nodes[loc_lndex]->index;
      }

      Node **getNodes()
      {
         return nodes;
      }

      Node* getNode(const int i)
      {
         return nodes[i];
      }

      void MarkingNodes(bool maker);

      //
      void setLocalNodeIndex(const int li, const MyInt n_lindex);

      //
      MyInt getLocalNodeIndex(const int li) const;

      // Neighbors
      void setNeighbors(Elem **ele_neighbors)
      {
         for (int i=0; i<getFacesNumber(); i++)
            neighbors[i] = ele_neighbors[i];
      }

      void setNeighbor(const int LocalIndex, Elem* ele_neighbor)
      {
         neighbors[LocalIndex] = ele_neighbor;
      }

      Elem** getNeighbors()
      {
         return neighbors;
      }

      Elem* getNeighbor(const int local_id) const
      {
         return neighbors[local_id];
      }

      //Domain partition
      MyInt getDomNodeIndex(const int loc_index)
      {
         return locnodes_index[loc_index];
      }

      void setDomNodeIndex(const int loc_index, const MyInt dom_nindex)
      {
         locnodes_index[loc_index] = dom_nindex;
      }

      // For partitioning by element
      void AllocateLocalIndexVector()
      {
         if(!locnodes_index)
            locnodes_index = new MyInt[getNodesNumber()];
      }

      void setDomainIndex(const int dom)
      {
         sub_dom = dom;
      }

      int getDomainIndex() const
      {
         return sub_dom;
      }

      // Local indicies
      void getLocalIndices_EdgeNodes(const int Edge, int *EdgeNodes);

      int getElementFaceNodes(const int Face, int *FacesNode);

      int getFaceType();

      void setOrder(const bool order)
      {
         quadratic = order;
      }

      // Output
      void Read(std::istream& is, Mesh_Group::Mesh *mesh, int fileType, const bool high_order);

      void WriteIndex(std::ostream& os = std::cout) const;
      void WriteGmsh(std::ostream& os, const int sdom_idx = 0) const;
      void WriteGSmsh(std::ostream& os, bool quad = false) const;

      void WriteSubDOM(std::ostream& os, bool quad = false) const;
      int getDataArray4BinaryOut(MyInt *ivar, bool quad = false) const;

      void WriteVTK_Type(std::ostream& os, bool isquad) const;
      void Write_index(std::ostream& os = std::cout) const;
      void WriteAll(std::ostream& os = std::cout) const;
      void WriteNeighbors(std::ostream& os = std::cout) const;

   private:
      // High order
      bool quadratic;

      int sub_dom;
      int no_faces_on_surface;
      //

      int PatchIndex;
      // Element type
      // 1 Line, 2 Quad, 3 Hex, 4 Tri, 5 Tet, 6 Pris
      ElemType ele_Type;

      MyInt *locnodes_index;
      Node **nodes;
      Elem **neighbors;

      int nnodes_ngl; //> number of non ghost nodes for linear element
      std::vector<int>  non_ghost_nodes;

      // Private methods
      int getElementFaces1D(int *FaceNode);
      int getElementFacesTri(const int Face, int *FaceNode);
      int getElementFacesQuad(const int Face, int *FaceNode);
      int getElementFacesHex(const int Face, int *FaceNode);
      int getElementFacesTet(const int Face, int *FaceNode);
      int getElementFacesPri(const int Face, int *FaceNode);
      int getElementFacesPyramid(const int Face, int *FaceNode);
      friend class Mesh;
};

} //end namespace

#endif
