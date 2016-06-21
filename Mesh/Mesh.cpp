#include "Mesh.h"

#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <limits>
#include <algorithm>

#include <cstdio> // for binary output

#include "Node.h"

//#define BUILD_MESH_EDGE

//------------------------------------------------------
//   Topology definition of geometrical element.
//    WW. 10.01.2005
//------------------------------------------------------
namespace Mesh_Group
{
//using namespace std;

Mesh::Mesh(bool quad)
{
   useQuadratic = quad;
   coordinate_system = 1;
   axisymmetry=false;
   max_ele_dim = 0;

   ouput_vtk_part_info = true;
}

Mesh::~Mesh()
{
   MyInt i;
   // Nodes
   for(i=0; i<(MyInt)node_vector.size(); i++)
      delete node_vector[i];
   node_vector.clear();

#ifdef BUILD_MESH_FACE
   // Surface faces
   for(i=0; i<(MyInt)face_vector.size(); i++)
      delete face_vector[i];
   face_vector.clear();
#endif

   // Element
   for(i=0; i<(MyInt)elem_vector.size(); i++)
      delete elem_vector[i];
   elem_vector.clear();
}


// Construct grid
//
/**************************************************************************
ConnectedNodes
**************************************************************************/
void Mesh::ConnectedNodes(bool quadratic)
{
   size_t i, j, k, n;
   Node* m_nod = NULL;
   Elem* m_ele = NULL;
   bool exist = false;
   //----------------------------------------------------------------------
   for(i=0; i<node_vector.size(); i++)
   {
      m_nod = node_vector[i];
      for(j=0; j<m_nod->ElementsRelated.size(); j++)
      {
         m_ele = elem_vector[m_nod->ElementsRelated[j]];
         for(int l=0; l<m_ele->getNodesNumber(quadratic); l++)
         {
            exist = false;
            const size_t nidx = m_ele->nodes[l]->index;
            for(k=0; k<m_nod->NodesRelated.size(); k++)
            {
               if(m_nod->NodesRelated[k] == nidx)
               {
                  exist = true;
                  break;
               }
            }
            if(!exist)
               m_nod->NodesRelated.push_back(nidx);
         }
      }
   }
   for(i=0; i<node_vector.size(); i++)
   {
      m_nod = node_vector[i];
      j = m_nod->NodesRelated.size();
      for(k=0; k<j; k++)
      {
         for(size_t ll=k; ll<j; ll++)
         {
            if(m_nod->NodesRelated[ll]<m_nod->NodesRelated[k])
            {
               n = m_nod->NodesRelated[k];
               m_nod->NodesRelated[k] = m_nod->NodesRelated[ll];
               m_nod->NodesRelated[ll] = n;
            }
         }
      }

   }
}

/**************************************************************************
ConnectedElements2Node
**************************************************************************/
void Mesh::ConnectedElements2Node(bool quadratic)
{
   size_t j, e, ni;
   Elem* thisElem0=NULL;
   Node * node = NULL;
   bool done = false;
   // set neighbors of node
   for(e=0; e<node_vector.size(); e++)
      node_vector[e]->ElementsRelated.clear();
   for(e=0; e<elem_vector.size(); e++)
   {
      thisElem0 = elem_vector[e];
      if(!thisElem0->getStatus()) continue;      // Not marked for use
      for(int i=0; i<thisElem0->getNodesNumber(quadratic); i++)
      {
         done = false;
         ni = thisElem0->getNodeIndex(i);
         node = node_vector[ni];
         for(j=0; j<node->ElementsRelated.size(); j++)
         {
            if(e==node->ElementsRelated[j])
            {
               done = true;
               break;
            }
         }
         if(!done)
            node->ElementsRelated.push_back(e);
      }
   }
}
void Mesh::ConstructGrid()
{
   int counter;
   int i, j, k, ii, jj, m0, m, n0, n;
   int nnodes0;
   size_t e,  e_size;
   bool done;
   double x_sum,y_sum,z_sum;

   int faceIndex_loc0[10];
   int faceIndex_loc[10];
   MyInt node_index_glb[20];
   MyInt node_index_glb0[20];

   Node **e_nodes0;
   Elem **Neighbors;
   Elem **Neighbors0;

   Elem* thisElem0=NULL;
   Elem* thisElem=NULL;

   clock_t start, finish;
   start = clock();

   //Elem->nodes not initialized

   e_size = (MyInt)elem_vector.size();
   NodesNumber_Linear= (MyInt)node_vector.size();

   //----------------------------------------------------------------------
   // set neighbors of node
   ConnectedElements2Node();
   //----------------------------------------------------------------------

   //----------------------------------------------------------------------
   // Compute neighbors and edges
   for(e=0; e<e_size; e++)
   {
      thisElem0 = elem_vector[e];
      nnodes0 = thisElem0->getNodesNumber(); // Number of nodes for linear element
      thisElem0->getNodeIndeces(node_index_glb0);

      e_nodes0 = thisElem0->getNodes();
      for(i=0; i<nnodes0; i++) // Nodes
         e_nodes0[i] = node_vector[node_index_glb0[i]];

      Neighbors0 = thisElem0->getNeighbors();
      m0 = thisElem0->getFacesNumber();
      // neighbors
      for(i=0; i<m0; i++) // Faces
      {
         if(Neighbors0[i])
            continue;
         n0 = thisElem0->getElementFaceNodes(i, faceIndex_loc0);

         done = false;
         for(k=0; k<n0; k++)
         {
            const size_t e_size_l = e_nodes0[faceIndex_loc0[k]]->ElementsRelated.size();
            for(size_t ei=0; ei<e_size_l; ei++)
            {
               const size_t ee = e_nodes0[faceIndex_loc0[k]]->ElementsRelated[ei];
               if(ee==e)
                  continue;

               thisElem = elem_vector[ee];
               thisElem->getNodeIndeces(node_index_glb);
               Neighbors = thisElem->getNeighbors();
               m = thisElem->getFacesNumber();

               for(ii=0; ii<m; ii++) // Faces
               {
                  n = thisElem->getElementFaceNodes(ii, faceIndex_loc);
                  if(n0!=n)
                     continue;

                  if(Neighbors0[ii])
                     continue;

                  counter = 0;
                  for(j=0; j<n0; j++)
                  {
                     for(jj=0; jj<n; jj++)
                     {
                        if(node_index_glb0[faceIndex_loc0[j]]
                              ==node_index_glb[faceIndex_loc[jj]])
                        {
                           counter++;
                           break;
                        }
                     }
                  }
                  if(counter==n)
                  {
                     Neighbors0[i] = thisElem;
                     Neighbors[ii] = thisElem0;
                     done = true;
                     break;
                  }
               }

               if(done)
                  break;
            }

            if(done)
               break;
         }
      }

      // set nodes
      thisElem0->setOrder(false);
      //
   }// Over elements

   // set faces on surfaces and others
   msh_no_line=0;  // Should be members of mesh
   msh_no_quad=0;
   msh_no_hexs=0;
   msh_no_tris=0;
   msh_no_tets=0;
   msh_no_pris=0;
   msh_no_pyra=0;
   for(e=0; e<e_size; e++)
   {
      thisElem0 = elem_vector[e];
      switch(thisElem0->getElementType())
      {
         case line:
            msh_no_line++;
            break;
         case quadri:
            msh_no_quad++;
            break;
         case hex:
            msh_no_hexs++;
            break;
         case tri:
            msh_no_tris++;
            break;
         case tet:
            msh_no_tets++;
            break;
         case prism:
            msh_no_pris++;
            break;
         case pyramid:
            msh_no_pyra++;
            break;
      }
      // Compute volume meanwhile
      //thisElem0->ComputeVolume();

      if(thisElem0->getElementType() == line)
         continue; // line element

      thisElem0->getNodeIndeces(node_index_glb0);
      Neighbors0 = thisElem0->getNeighbors();
      m0 = thisElem0->getFacesNumber();

#ifdef BUILD_MESH_FACE
      // Check face on surface
      for(i=0; i<m0; i++) // Faces
      {
         if(Neighbors0[i])
            continue;
         Elem* newFace = new Elem((MyInt)face_vector.size(), thisElem0, i);
//          thisElem0->boundary_type='B';
         thisElem0->no_faces_on_surface++;
         face_vector.push_back(newFace);
         Neighbors0[i] = newFace;
      }
#endif
      thisElem0->setNeighbors(Neighbors0);

   }
   NodesNumber_Quadratic= (MyInt)node_vector.size();
   if((msh_no_hexs+msh_no_tets+msh_no_pris+msh_no_pyra)>0) max_ele_dim=3;
   else if((msh_no_quad+msh_no_tris)>0) max_ele_dim=2;
   else max_ele_dim=1;
   //----------------------------------------------------------------------
   // Node information
   // 1. Default node index <---> eqs index relationship
   // 2. Coordiate system flag
   x_sum=0.0;
   y_sum=0.0;
   z_sum=0.0;
   for(e=0; e<node_vector.size(); e++)
   {
      x_sum += fabs(node_vector[e]->X());
      y_sum += fabs(node_vector[e]->Y());
      z_sum += fabs(node_vector[e]->Z());
   }
   if(x_sum>0.0&&y_sum<DBL_MIN&&z_sum<DBL_MIN)
      coordinate_system = 10;
   else if(y_sum>0.0&&x_sum<DBL_MIN&&z_sum<DBL_MIN)
      coordinate_system = 11;
   else if(z_sum>0.0&&x_sum<DBL_MIN&&y_sum<DBL_MIN)
      coordinate_system = 12;
   else if(x_sum>0.0&&y_sum>0.0&&z_sum<DBL_MIN)
      coordinate_system = 21;
   else if(x_sum>0.0&&z_sum>0.0&&y_sum<DBL_MIN)
      coordinate_system = 22;
   else if(x_sum>0.0&&y_sum>0.0&&z_sum>0.0)
      coordinate_system = 32;
   /*  // 23.05.2008. WW. Futher test is needed
   // 1D in 2D
   if(msh_no_line>0)   //
   {
     if(x_sum>0.0&&y_sum>0.0&&z_sum<DBL_MIN)
        coordinate_system = 22;
     if(x_sum>0.0&&z_sum>0.0&&y_sum<DBL_MIN)
        coordinate_system = 22;
   }
   */
   //----------------------------------------------------------------------

   // For sparse matrix
   ConnectedNodes(false);
   //
   finish = clock();
   std::cout<<"\nCPU time elapsed in constructing topology of grids: "
       <<(double)(finish - start) / CLOCKS_PER_SEC<<"s"<<std::endl<<std::endl;

}

/**************************************************************************
Programing:
07/2007 WW Implementation
**************************************************************************/
void Mesh::GenerateHighOrderNodes()
{
   int i, k, ii;
   int nnodes0, nedges0, nedges;
   size_t e, ei,  e_size,  e_size_l;
   int edgeIndex_loc0[3];
   int edgeIndex_loc1[3];
   bool done;
   double x0,y0,z0;

   clock_t start, finish;
   start = clock();

   //
   Node *aNode=NULL;
   Node **e_nodes0;
   Elem *thisElem0=NULL;
   Elem *thisElem=NULL;

   //----------------------------------------------------------------------
   // Loop over elements
   e_size = elem_vector.size();

   for(e=0; e<e_size; e++)
   {
      elem_vector[e]->Marking(false);
   }

   for(e=0; e<e_size; e++)
   {
      thisElem0 = elem_vector[e];
      nnodes0 = thisElem0->getNodesNumber(); // Number of nodes for linear element
//      thisElem0->GetNodeIndeces(node_index_glb0);

      e_nodes0 = thisElem0->getNodes();
      for(i=0; i<nnodes0; i++) // Nodes
         e_nodes0[i] = thisElem0->getNode(i);
      // --------------------------------

      // Edges
      nedges0 = thisElem0->getEdgesNumber();
      // Check if there is any neighbor that has new middle points
      for(i=0; i<nedges0; i++)
      {
         thisElem0->getLocalIndices_EdgeNodes(i, edgeIndex_loc0);
         const size_t ena0 = thisElem0->getNodeIndex(edgeIndex_loc0[0]);
         const size_t ena1 = thisElem0->getNodeIndex(edgeIndex_loc0[1]);
         // Check neighbors
         done = false;
         for(k=0; k<2; k++)
         {
            e_size_l = (MyInt)e_nodes0[edgeIndex_loc0[k]]->ElementsRelated.size();
            for(ei=0; ei<e_size_l; ei++)
            {
               const size_t ee = e_nodes0[edgeIndex_loc0[k]]->ElementsRelated[ei];
               if(ee==e)
                  continue;
               thisElem = elem_vector[ee];

               // If this element has already been proccessed
               if( thisElem->getStatus() )
               {
                  nedges = thisElem->getEdgesNumber();
                  // Edges of neighbors
                  for(ii=0; ii<nedges; ii++)
                  {
                     thisElem->getLocalIndices_EdgeNodes(ii, edgeIndex_loc1);

                     const size_t enb0 = thisElem->getNodeIndex(edgeIndex_loc1[0]);
                     const size_t enb1 = thisElem->getNodeIndex(edgeIndex_loc1[1]);

                     if(   (( ena0 == enb0 ) && (ena1 == enb1))
                           || (( ena0 == enb1 ) && (ena1 == enb0))
                       )
                     {
                        aNode = thisElem->getNode(edgeIndex_loc1[2]);
                        e_nodes0[edgeIndex_loc0[2]] = aNode;
                        done = true;
                        break;
                     }

                  } //  for(ii=0; ii<nedges; ii++)

               }

               if(done) break;
            } // for(ei=0; ei<e_size_l; ei++)
            if(done) break;
         }//for(k=0;k<2;k++)
         if(!done)
         {
            const Node *na = thisElem0->getNode(edgeIndex_loc0[0]);
            const Node *nb = thisElem0->getNode(edgeIndex_loc0[1]);
            aNode = new Node((MyInt)node_vector.size());
            aNode->setX(0.5*(na->X() + nb->X()));
            aNode->setY(0.5*(na->Y() + nb->Y()));
            aNode->setZ(0.5*(na->Z() + nb->Z()));

            e_nodes0[edgeIndex_loc0[2]] = aNode;
            node_vector.push_back(aNode);
         }
      } //  for(i=0; i<nedges0; i++)

      // No neighors or no neighbor has new middle point
      //
      if(thisElem0->getElementType()==quadri) // Quadrilateral
      {
         x0=y0=z0=0.0;
         for(i=0; i<4; i++) // Nodes
         {
            x0 += e_nodes0[i]->X();
            y0 += e_nodes0[i]->Y();
            z0 += e_nodes0[i]->Z();
         }
         x0 /= 4.;
         y0 /= 4.;
         z0 /= 4.;

         aNode = new Node((MyInt)node_vector.size());
         aNode->setX(x0);
         aNode->setY(y0);
         aNode->setZ(z0);

         e_nodes0[8] = aNode;
         node_vector.push_back(aNode);
      }
      // Set edges and nodes
      thisElem0->setOrder(true);
      //
      thisElem0->Marking(true);
   }// Over elements

   //
   NodesNumber_Quadratic= (MyInt)node_vector.size();
   for(e=0; e<e_size; e++)
   {
      thisElem0 = elem_vector[e];
      for(i=thisElem0->getNodesNumber(); i<thisElem0->getNodesNumberHQ(); i++)
      {
         done = false;
         aNode = thisElem0->getNode(i);
         for(k=0; k<(int)aNode->ElementsRelated.size(); k++)
         {
            if(e==aNode->ElementsRelated[k])
            {
               done = true;
               break;
            }
         }
         if(!done)
            aNode->ElementsRelated.push_back(e);
      }
   }

   // For sparse matrix
   ConnectedNodes(true);
   //ConnectedElements2Node(true);
   //

   finish = clock();
   std::cout<<"\n\tCPU time elapsed in generating high oder elements: "
       <<(double)(finish - start) / CLOCKS_PER_SEC<<"s"<<std::endl;

}

void Mesh::ConstructSubDomain_by_Elements(const std::string fname, const int num_parts, const bool osdom)
{
   std::string str;
   std::string stro;
   int dom;
   int max_dom;
   int k,kk;
   MyInt i,j;
   //  int ntags = 3;

   std::fstream gmsh_out;

   std::string deli = " ";
   //

   std::string s_nparts;
   std::stringstream ss;
   ss << num_parts;
   ss >> s_nparts;
   ss.clear();

   str = fname + ".mesh.epart." + s_nparts;
   stro = fname + "." + s_nparts +"ddc";

   //namef = ".mesh.epart."; //+str_buf;
   std::ifstream part_in;
   std::fstream part_out;
   part_out.open(stro.c_str(), std::ios::out | std::ios::trunc );
   // Output for gmsh

   if(osdom)
   {
      stro = fname + "_gmsh.msh";
      gmsh_out.open(stro.c_str(), std::ios::out );
      //gmsh_out<<"$NOD" << std::endl;
      gmsh_out<<"$MeshFormat\n2 0 8\n$EndMeshFormat\n$Nodes"<<std::endl;
      gmsh_out<<node_vector.size() << std::endl;
      Node *node;
      for(i=0; i<(MyInt)node_vector.size(); i++)
      {
         gmsh_out<<i+1<<" ";
         node = node_vector[i];
         gmsh_out<<node->X()<<" ";
         gmsh_out<<node->Y()<<" ";
         gmsh_out<<node->Z() << std::endl;
      }
      //gmsh_out<<"$ENDNOD" << std::endl;
      //gmsh_out<<"$ELM" << std::endl;
      gmsh_out<<"$EndNodes\n$Elements" << std::endl;
      gmsh_out<<(MyInt)elem_vector.size() << std::endl;
   }

   //
   part_in.open(str.c_str());
   if(!part_in.is_open())
   {
      std::cerr<<("Error: cannot open .epart file . It may not exist !");
      abort();
   }

   max_dom=0;
   Elem *ele = NULL;

   for(i=0; i<(MyInt)elem_vector.size(); i++)
   {
      part_in>>dom>>std::ws;
      ele = elem_vector[i];
      ele->setDomainIndex(dom);

      if(dom>max_dom) max_dom = dom;

      if(osdom)
      {
         ele->WriteGmsh(gmsh_out, dom+1);
      }

   }
   max_dom++;
   part_in.close();
   std::remove(str.c_str()); // Remove the metis file.

   if(osdom)
   {
      gmsh_out<<"$EndElements"<<std::endl;
      gmsh_out.close();
   }
   //

   //Output ddc file
   // MyInt *nod_dom = new MyInt[max_dom];
   MyInt *ele_dom = new MyInt[max_dom];
   for(k=0; k<max_dom; k++)
   {
      ele_dom[k]=0;
      //nod_dom[k]=0;
      for(j=0; j<(MyInt)elem_vector.size(); j++)
      {
         if(elem_vector[j]->getDomainIndex()==k)
            ele_dom[k] += 1;
      }
   }

   bool done = false;
   MyInt n_index=0;
   std::vector<std::size_t> nodes_dom;
   std::vector<Elem*> eles_dom;

   //
   for(k=0; k<max_dom; k++)
   {
      part_out<<"#DOMAIN "<<k<<std::endl;
      part_out<<"$ELEMENTS "<<ele_dom[k]<<std::endl;
      nodes_dom.clear();
      eles_dom.clear();
      for(j=0; j<(MyInt)elem_vector.size(); j++)
      {
         ele = elem_vector[j];
         for(kk=0; kk<ele->getNodesNumber(); kk++)
         {
            ele->setLocalNodeIndex(kk, -1);
            ele->AllocateLocalIndexVector();
            ele->setDomNodeIndex(kk, -1);
         }
      }
      for(j=0; j<(MyInt)elem_vector.size(); j++)
      {
         ele = elem_vector[j];

         if(ele->getDomainIndex()==k)
         {
            for(kk=0; kk<ele->getNodesNumber(); kk++)
            {
               done = false;
               n_index = ele->getLocalNodeIndex(kk);
               if(n_index>-1)
               {
                  ele->setDomNodeIndex(kk, n_index);
                  done = true;
               }
               if(!done)
               {
                  ele->setDomNodeIndex(kk, (MyInt)nodes_dom.size()); //For test output
                  ele->setLocalNodeIndex(kk, (MyInt)nodes_dom.size());
                  nodes_dom.push_back(ele->getNodeIndex(kk));
               }
            }
            part_out<<ele->getIndex() << std::endl;
            eles_dom.push_back(ele); //TEST OUT
         }
      }
      part_out<<"$NODES_INNER "<<(MyInt)nodes_dom.size() << std::endl;
      for(j=0; j<(MyInt)nodes_dom.size(); j++)
         part_out<<nodes_dom[j] << std::endl;

      if(osdom)
      {
         std::string i_nparts;
         ss << k;
         ss >> i_nparts;
         ss.clear();

         std::string name_f = fname+"_"+i_nparts+"_of_"+s_nparts+"subdomains.msh";
         std::fstream test_out;
         test_out.open(name_f.c_str(), std::ios::out|std::ios::trunc );

         Node *nod = 0;
         //GMSH test_out<<"$NOD" << std::endl;
         //GMSH test_out<<(MyInt)nodes_dom.size() << std::endl;
         test_out<<"#0#0#0#1#0.0#0#################################################################" << std::endl;
         test_out<<"0 "<<(MyInt)nodes_dom.size()<<" "<<(MyInt)eles_dom.size() << std::endl;
         for(j=0; j<(MyInt)nodes_dom.size(); j++)
         {
            nod = node_vector[nodes_dom[j]];
            //GMSH  test_out<<j+1<<"  "
            test_out<<j<<deli
                    << nod->X()<<deli<< nod->Y()<<deli<< nod->Z()  << std::endl;
         }
         //GMSH test_out<<"$ENDNOD" << std::endl;
         //GMSH test_out<<"$ELE" << std::endl;
         //GMSG test_out<<(MyInt)eles_dom.size() << std::endl;
         for(j=0; j<(MyInt)eles_dom.size(); j++)
         {
            ele = eles_dom[j];

            //GMSH  ele->WriteGmsh(test_out, k+1);

            test_out<<j<<deli<<ele->getPatchIndex()<<deli<<ele->getName()<<deli;
            for(kk=0; kk<ele->getNodesNumber(); kk++)
               test_out<< ele->getDomNodeIndex(kk)<<deli;
            test_out << std::endl;
         }
         test_out.clear();
         test_out.close();
      }

   }
   part_out<<"#STOP " << std::endl;
   part_out.clear();
   part_out.close();

   //
   delete ele_dom;
   //delete nod_dom;
}

// NW
struct ConnEdge
{
   Node* first;
   Node* second;
   ConnEdge(Node* i, Node* j) : first(i), second(j) {};
};

inline bool operator<(const ConnEdge& lhs, const ConnEdge& rhs)
{
   if (lhs.first->getIndex() != rhs.first->getIndex()) return lhs.first->getIndex() < rhs.first->getIndex();
   else return lhs.second->getIndex() < rhs.second->getIndex();
}

/**
* Converts a number (double, float, int, ...) into a string
* \param d The number to be converted
* \return The number as string
*/
template<typename T> std::string number2str(T d)
{
   std::stringstream out;
   out << d;
   return out.str();
}

/*!
\brief void Mesh::ConstructSubDomain_by_Nodes

Partition a mesh ny nodes

02.2012 WW
*/

void Mesh::ConstructSubDomain_by_Nodes(const MeshPartConfig mpc)
{
   // Material data partitioning
   int num_data = 0;
   std::vector<std::string> m_headers;
   std::vector<std::size_t> m_header_marker_per_data;
   std::vector<std::string> m_datanames;
   //std::vector<MyInt> m_ele_idx;
   std::vector<double> m_ele_val;

   if(mpc.mat_fname.size() !=0)
   {
      std::string line_buffer;
      std::string mat_fname_abs = mpc.fpath + mpc.mat_fname;

      std::ifstream is_mat(mat_fname_abs.c_str());
      if(!is_mat.good())
      {
         std::cout<<"Material data file "<<mat_fname_abs<<" does not exist"<<std::endl;
         exit(1);
      }
      is_mat >> num_data;
      m_datanames.resize(num_data);
      m_header_marker_per_data.resize(num_data +1 );
      m_header_marker_per_data[0] = 0;
      for(int k=0; k<num_data; k++)
      {
         std::string data_name;
         is_mat >> m_datanames[k];
         m_datanames[k] = mpc.fpath + m_datanames[k];
      }
      is_mat.close();

      // Read each data file
      for(int k=0; k<num_data; k++)
      {
         is_mat.open(m_datanames[k].c_str());
         if(!is_mat.good())
         {
            std::cout<<"Material data file "<<m_datanames[k]<<" does not exist"<<std::endl;
            exit(1);
         }
         while(!is_mat.eof())
         {
            getline(is_mat, line_buffer);
            if(line_buffer.find("$DATA")!=std::string::npos)
            {
               m_headers.push_back(line_buffer);
               const size_t ne = elem_vector.size();
               for(size_t ie = 0; ie<ne; ie++)
               {
                  MyInt index;
                  double m_val;
                  is_mat >> index >> m_val;
                  //m_ele_idx.push_back(index);
                  m_ele_val.push_back(m_val);
               }
            }
            else if(line_buffer.find("#STOP")!=std::string::npos)
            {
               break;
            }
            else if(line_buffer.size()>0)
            {
               m_headers.push_back(line_buffer);
            }
         }
         m_header_marker_per_data[k+1] = m_headers.size();
         is_mat.clear();
         is_mat.close();
      }
   }

   std::string f_iparts;
   std::string o_part_msh;
   std::string fname = mpc.fname;

   // Number of integer variables of subdomain elements
   MyInt nmb_element_idxs;
   MyInt nmb_element_idxs_g;

   std::string deli = " ";

   //
   const MyInt num_headers = 13;
   MyInt head[13];
   MyInt ivar[23];
   // For binary output: all data in on string
   std::vector<MyInt> ele_info;
   // node rank offset
   head[10] = 0;
   // element rank offset
   head[11] = 0;
   // ghost element rank offset
   head[12] = 0;

   const int num_parts = mpc.num_parts;
   // Read a METIS node partitioned file
   //--------------------------------------------------------------------------
   const std::string s_nparts = number2str(num_parts);
   f_iparts = fname + ".mesh.npart." + s_nparts;
   //o_part_msh = fname + "." + s_nparts +"mesh";

   std::ifstream npart_in(f_iparts.c_str());
   if(!npart_in.is_open())
   {
      std::cerr<<("Error: cannot open .npart file . It may not exist !");
      exit(1);
   }

   const MyInt nn = static_cast<MyInt>(node_vector.size());

   std::vector<bool> sdom_marked(nn);
   std::vector<MyInt> dom_idx(nn);

   for(MyInt i=0; i<nn; i++)
   {
      npart_in >> dom_idx[i] >>std::ws;
      sdom_marked[i] = false;
   }
   npart_in.close();
   //remove(f_iparts.c_str());

   //--------------------------------------------------------------------------
   // Output a partitioned mesh
   //--------------------------------------------------------------------------
#define OUTPUT_TO_SINGLE_FILE
#ifdef OUTPUT_TO_SINGLE_FILE
   std::string name_f;
   std::fstream os_subd;
   std::fstream os_subd_head;
   std::fstream os_subd_node;
   FILE *of_bin_cfg = 0; // ostream os_bin_cfg
   FILE *of_bin_nod = 0; // ostream os_bin_nod
   FILE *of_bin_ele = 0; // ostream os_bin_ele
   FILE *of_bin_ele_g = 0; // ostream os_bin_ele_g
   const bool binary_output = mpc.binary_output;

   if(!binary_output)
   {
      name_f = fname+"_partitioned_cfg"+ s_nparts + ".msh";
      os_subd_head.open(name_f.c_str(), std::ios::out|std::ios::trunc );
      name_f = "Subdomain mesh "
               "(Nodes;  Nodes_linear; Elements; Ghost elements; Nodes of Linear elements; Nodes of quadratic elements) "
               "Nodes of whole linear elements; Nodes of whole quadratic elements; "
               "Total integer variables of elements;Total integer variables of ghost elements  ";
      os_subd_head << name_f << std::endl;
      os_subd_head<<num_parts << std::endl;

      name_f = fname+"_partitioned_elems_"+ s_nparts + ".msh";
      os_subd.open(name_f.c_str(), std::ios::out|std::ios::trunc );

      name_f = fname+"_partitioned_nodes_"+ s_nparts + ".msh";
      os_subd_node.open(name_f.c_str(), std::ios::out|std::ios::trunc );

      std::setw(14);
      os_subd.precision(14);
      //os_subd.setf(ios::fixed, ios::scientific);
      os_subd.setf(std::ios::scientific);
   }
   else
   {
      name_f = fname+"_partitioned_msh_cfg"+ s_nparts + ".bin";
      of_bin_cfg = fopen(name_f.c_str(), "wb");

      name_f = fname+"_partitioned_msh_ele"+ s_nparts + ".bin";
      of_bin_ele = fopen(name_f.c_str(), "wb");

      name_f = fname+"_partitioned_msh_ele_g"+ s_nparts + ".bin";
      of_bin_ele_g = fopen(name_f.c_str(), "wb");

      name_f = fname+"_partitioned_msh_nod"+ s_nparts + ".bin";
      of_bin_nod = fopen(name_f.c_str(), "wb");

   }

#endif

   // for linear
   std::vector<std::size_t> sdom_start_node(num_parts);
   std::vector<std::size_t> sdom_end_act_node(num_parts);
   std::vector<std::size_t> sdom_end_node(num_parts);

   // quadratic
   std::vector<std::size_t> sdom_start_node_hq(num_parts);
   std::vector<std::size_t> sdom_end_act_node_hq(num_parts);
   std::vector<std::size_t> sdom_end_node_hq(num_parts);

   const MyInt ne_total = static_cast<MyInt>(elem_vector.size());

   std::vector<Node*> sbd_nodes;        // nodes of a partition for linear element
   std::vector<Node*> sbd_nodes_hq;     // nodes of a partition for quadratic element
   std::vector<std::vector<std::set<ConnEdge> > > vec_neighbors(num_parts); //NW

   // Should be removed for openmp 
   for (std::size_t i=0; i<elem_vector.size(); i++)
   {
      elem_vector[i]->Marking(false);
   }

   MyInt node_id_offset = 0;
   MyInt node_id_offset_h = NodesNumber_Linear;
   for(int idom=0; idom<num_parts; idom++)
   {

      sdom_start_node[idom] = sbd_nodes.size();
      sdom_start_node_hq[idom] = sbd_nodes_hq.size();

      if (mpc.out_cct)
      {
         vec_neighbors[idom].resize(num_parts);
      }
      std::cout << "Process partition: " << idom << std::endl;

      nmb_element_idxs = 0;
      nmb_element_idxs_g = 0;

      // Release all nodes of the whole mesh
      for(MyInt j=0; j<nn; j++)
         node_vector[j]->Marking(false);

      MyInt num_nodes_active_l = 0; // acitve nodes for linear elements
      for(MyInt j=0; j<NodesNumber_Linear; j++)
      {
         if(dom_idx[j] == idom && (!sdom_marked[j]))
         {
            sbd_nodes.push_back(node_vector[j]);
            sdom_marked[j] = true;  // avoid other subdomain use this node
            node_vector[j]->Marking(true);
            num_nodes_active_l++;
         }
      }
      sdom_end_act_node[idom] = sbd_nodes.size();
      // Find the elements in this subdomain.
      std::vector<Elem*> in_subdom_elements;
      std::vector<Elem*> ghost_subdom_elements;

      MyInt num_nodes_active_h = 0; // acitve nodes for quadratic elements
      if (useQuadratic)
      {
         for(MyInt j=NodesNumber_Linear; j< nn; j++)
         {
            if(dom_idx[j] == idom && (!sdom_marked[j]))
            {
               sbd_nodes_hq.push_back(node_vector[j]);
               sdom_marked[j] = true;  // avoid other subdomain use this node
               node_vector[j]->Marking(true);
               num_nodes_active_h++;
            }
         }
      }
      sdom_end_act_node_hq[idom] = sbd_nodes_hq.size();

      for (std::size_t i=0; i<elem_vector.size(); i++)
      {
         Elem* elem = elem_vector[i];

         // Should be removed for openmp 
         if ( elem->getStatus() )
            continue;

         elem->non_ghost_nodes.clear();
         int nn_ngl = 0;
         std::vector<int> g_nodes; // ghost nodes in ghost elements
         for(int kk=0; kk<elem->getNodesNumber(useQuadratic); kk++)
         {
            if (elem->getNode(kk)->getStatus())
            {
               elem->non_ghost_nodes.push_back(kk);
               if (kk < elem->getNodesNumber())
                  nn_ngl++;
            }
            else
               g_nodes.push_back(kk);
         }

         if (elem->non_ghost_nodes.size() == 0)
            continue;

         if (elem->non_ghost_nodes.size()
               == static_cast<std::size_t>(elem->getNodesNumber(useQuadratic)) )
         {
            in_subdom_elements.push_back(elem);
            elem->Marking(true);    // Should be removed for openmp 
            elem->non_ghost_nodes.clear();
         }
         else
         {
            elem->nnodes_ngl = nn_ngl;
            ghost_subdom_elements.push_back(elem);
         }

         // overlapping edges
         if (mpc.out_cct)
         {
            for (std::size_t ig=0; ig<g_nodes.size(); ig++)
            {
               MyInt ig_id = elem->getNode(g_nodes[ig])->global_index; //index; //;
               MyInt ig_dom = dom_idx[ig_id];
               for (size_t ing=0; ing<elem->non_ghost_nodes.size(); ing++)
               {
                  //MyInt ing_id = a_elem->getNode(ng_nodes[ing])->index; //global_index;
                  vec_neighbors[idom][ig_dom].insert(ConnEdge(elem->getNode(elem->non_ghost_nodes[ing]),
                                                     elem->getNode(g_nodes[ig]))); // inner node - ghost node
                  //vec_neighbors[ig_dom].insert(ConnEdge(ing_id, ig_id)); // inner node - ghost node
               }
            }
         }
      }

      //-----------------------------------------------
      // Add inactive nodes in ghost elements
      const MyInt ne_g = static_cast<MyInt>(ghost_subdom_elements.size());
      for(MyInt j=0; j<ne_g; j++)
      {
         Elem* a_elem = ghost_subdom_elements[j];
         for(int k=0; k<a_elem->getNodesNumber(useQuadratic); k++)
            a_elem->nodes[k]->Marking(false);

         // Existing nodes
         const int ngh_nodes = static_cast<int>(a_elem->non_ghost_nodes.size());
         for(int k=0; k<ngh_nodes; k++)
            a_elem->nodes[a_elem->non_ghost_nodes[k]]->Marking(true);
      }
      //
      for(MyInt j=0; j<ne_g; j++)
      {
         Elem*  a_elem = ghost_subdom_elements[j];
         for(int k=0; k<a_elem->getNodesNumber(); k++)
         {
            Node* a_node = a_elem->nodes[k];
            if(a_node->getStatus())
               continue;
            a_node->Marking(true);
            sbd_nodes.push_back(a_node);
         }
         //
         for(int k=a_elem->getNodesNumber(); k<a_elem->getNodesNumber(useQuadratic); k++)
         {
            Node* a_node = a_elem->nodes[k];
            if(a_node->getStatus())
               continue;
            a_node->Marking(true);
            sbd_nodes_hq.push_back(a_node);
         }
      }

      // Ending mark
      sdom_end_node[idom] = sbd_nodes.size();
      sdom_end_node_hq[idom] = sbd_nodes_hq.size();

      // Renumber
      size_t start_id = sdom_start_node[idom];
      size_t end_id = sdom_end_act_node[idom];
      MyInt local_id = 0;
      for(size_t in = start_id; in < end_id; in++)
      {
         Node *node = sbd_nodes[in];
         node->index = node_id_offset;
         node->local_index = local_id;
         local_id++;
         node_id_offset++;
      }

      //Inactive nodes, local node IDs
      start_id = sdom_end_act_node[idom];
      end_id = sdom_end_node[idom];
      for(size_t in = start_id; in < end_id; in++)
      {
         sbd_nodes[in]->local_index = local_id;
         local_id++;
      }

      // Active nodes, higher order
      start_id = sdom_start_node_hq[idom];
      end_id = sdom_end_act_node_hq[idom];
      for(size_t in = start_id; in < end_id; in++)
      {
         Node *node = sbd_nodes_hq[in];
         node->index = node_id_offset_h;
         node->local_index = local_id;
         local_id++;
         node_id_offset_h++;
      }

      //Inactive nodes, local node IDs,  higher order
      start_id = sdom_end_act_node_hq[idom];
      end_id = sdom_end_node_hq[idom];
      for(size_t in = start_id; in < end_id; in++)
      {
         sbd_nodes_hq[in]->local_index = local_id;
         local_id++;
      }

      // Count the total integer variables of this subdomain
      const MyInt nei = static_cast<MyInt>(in_subdom_elements.size());
      nmb_element_idxs =  3*nei;
      for(MyInt j=0; j<nei; j++)
      {
         nmb_element_idxs += in_subdom_elements[j]->getNodesNumber(useQuadratic);
      }
      const MyInt neg = static_cast<MyInt>( ghost_subdom_elements.size());
      //  mat index, element type, number of element, number of ghost nodes, number of ghost nodes of high order elements
      nmb_element_idxs_g = 5*neg;
      for(MyInt j=0; j<neg; j++)
      {
         nmb_element_idxs_g += ghost_subdom_elements[j]->getNodesNumber(useQuadratic);
         nmb_element_idxs_g += static_cast<MyInt>(ghost_subdom_elements[j]->non_ghost_nodes.size());
      }

      std::string dom_str = number2str(idom);
#ifdef OUTPUT_TO_DIFF_FILES
      // Make output of this subdomain for simulation
      //string name_f = fname+"_"+dom_str+"_of_"+s_nparts+"_subdomains.msh";
      std::string name_f = fname+"_"+dom_str+".msh";
      fstream os_subd(name_f.c_str(), ios::out|ios::trunc );
      //os_subd<<"#FEM_MSH\n   $PCS_TYPE\n    NULL" << std::endl;
      //os_subd<<" $NODES\n"<<size_sbd_nodes << std::endl;


      name_f = "Subdomain mesh "
               "(Nodes;  Nodes_linear; Elements; Ghost elements; Nodes of Linear elements; Nodes of quadratic elements) "
               "Nodes of Linear whole elements; Nodes of whole quadratic elements; "
               "Total integer variables of elements;Total integer variables of ghost elements  ";
      os_subd<<name_f << std::endl;
#endif
      // Number of active elements
      const MyInt nei_size = static_cast<MyInt>(in_subdom_elements.size());
      const MyInt offset_e = nei_size + nmb_element_idxs;
      const MyInt offset_e_g = ne_g + nmb_element_idxs_g;

      const MyInt sdom_nodes_l = static_cast<MyInt>(sdom_end_node[idom] - sdom_start_node[idom]);
      const MyInt sdom_nodes_h = static_cast<MyInt>(sdom_end_node_hq[idom] - sdom_start_node_hq[idom]);

      if(!binary_output)
      {
         os_subd_head<<sdom_nodes_l + sdom_nodes_h <<deli<<sdom_nodes_l<<deli<<nei_size
                     <<deli<<ne_g<<deli<<num_nodes_active_l<<deli<<num_nodes_active_l + num_nodes_active_h
                     <<deli<<NodesNumber_Linear<<deli<<NodesNumber_Quadratic
                     <<deli<<nmb_element_idxs<<deli<<nmb_element_idxs_g << std::endl;
      }
      else
      {
         head[0] = sdom_nodes_l + sdom_nodes_h;
         head[1] = sdom_nodes_l;
         head[2] = nei_size;
         head[3] = ne_g;
         head[4] = num_nodes_active_l; // active nodes
         head[5] = num_nodes_active_l + num_nodes_active_h; // active nodes
         head[6] = NodesNumber_Linear;
         head[7] = NodesNumber_Quadratic;
         head[8] = nmb_element_idxs;
         head[9] = nmb_element_idxs_g;

         //os_bin_cfg.write( reinterpret_cast <const char*> (head), num_headers*sizeof(MyInt));
         fwrite(head, 1, num_headers*sizeof(MyInt), of_bin_cfg);

         head[10] += head[0] * sizeof(Node_Str);
         head[11] += offset_e * sizeof(MyInt);
         head[12] += offset_e_g * sizeof(MyInt);

      }

      //os_subd<<"Elements" << std::endl;
      MyInt counter = 0;
      if(binary_output)
      {
         ele_info.clear();
         ele_info.resize(offset_e);
      }
      counter = nei_size;

      for(MyInt j=0; j<nei_size; j++)
      {
         if(binary_output)
         {
            const int nvar = in_subdom_elements[j]->getDataArray4BinaryOut(ivar, useQuadratic);

            ele_info[j] = counter;
            for(int m=0; m<nvar; m++)
            {
               ele_info[counter] = ivar[m];
               counter++;
            }

            //os_bin_ele.write( reinterpret_cast <const char*> (ivar), nvar*sizeof(MyInt));
            //fwrite(ivar, 1, nvar*sizeof(MyInt), of_bin_ele);
         }
         else
            in_subdom_elements[j]->WriteSubDOM(os_subd, useQuadratic);
      }
      // in one MyInt string buffer
      if(binary_output)
      {
         fwrite(&ele_info[0], 1, (offset_e)*sizeof(MyInt), of_bin_ele);
         ele_info.clear();

      }

      //os_subd<<"Ghost elements" << std::endl;
      counter = 0;
      if(binary_output)
      {
         ele_info.clear();
         ele_info.resize(offset_e_g);
      }
      counter = ne_g;

      for(MyInt j=0; j<ne_g; j++)
      {
         Elem* a_elem = ghost_subdom_elements[j];
         const int ngh_nodes = static_cast<int>(a_elem->non_ghost_nodes.size());

         if(binary_output)
         {

            ele_info[j] = counter;

            const int nvar = a_elem->getDataArray4BinaryOut(ivar, useQuadratic);
            //os_bin_ele_g.write( reinterpret_cast <const char*> (ivar), nvar*sizeof(MyInt));
            //fwrite(ivar, 1, nvar*sizeof(MyInt), of_bin_ele);

            for(int m=0; m<nvar; m++)
            {
               ele_info[counter] = ivar[m];
               counter++;
            }

            ivar[0] = a_elem->nnodes_ngl;
            ivar[1] = ngh_nodes;
            for(int kk=0; kk<ngh_nodes; kk++)
            {
               ivar[kk+2] = a_elem->non_ghost_nodes[kk];
            }
            // os_bin_ele_g.write( reinterpret_cast <const char*> (ivar), (ngh_nodes+2)*sizeof(MyInt));
            // fwrite(ivar, 1, (ngh_nodes+2)*sizeof(MyInt), of_bin_ele);
            for(int m=0; m<ngh_nodes+2; m++)
            {
               ele_info[counter] = ivar[m];
               counter++;
            }

         }
         else
         {
            a_elem->WriteSubDOM(os_subd, useQuadratic);
            os_subd<<a_elem->nnodes_ngl<<deli<<ngh_nodes<<deli;
            for(int kk=0; kk<ngh_nodes; kk++)
            {
               os_subd<<a_elem->non_ghost_nodes[kk]<<deli;
            }
            os_subd << std::endl;
         }
      }

      os_subd << std::endl;

      // in one long string buffer
      if(binary_output)
      {
         fwrite(&ele_info[0], 1, (offset_e_g)*sizeof(MyInt), of_bin_ele_g);
         ele_info.clear();

      }

      //----------------------------------------------------------------------------------
      /// Material data partitioning
      if( num_data > 0)
      {
         std::ofstream os_mat;
         for(int mm = 0; mm<num_data; mm++)
         {
            std::string mat_ofile_name = m_datanames[mm] + dom_str;
            os_mat.open(mat_ofile_name.c_str(), std::ios::trunc);
            for(size_t mh = m_header_marker_per_data[mm]; mh<m_header_marker_per_data[mm+1]; mh++ )
            {
               os_mat<<m_headers[mh] << std::endl;
            }

            const MyInt e_shift = ne_total*mm;
            for(MyInt j=0; j<nei_size; j++)
            {
               const size_t entry_index = in_subdom_elements[j]->getIndex() + e_shift;
               os_mat<<j<<deli<<m_ele_val[entry_index] << std::endl;
            }
            for(MyInt j=0; j<ne_g; j++)
            {
               const size_t entry_index  = ghost_subdom_elements[j]->getIndex() + e_shift;
               os_mat<<j+nei_size<<deli<<m_ele_val[entry_index] << std::endl;
            }
            os_mat<<"#STOP" << std::endl;
            os_mat.clear();
            os_mat.close();
         }
      }
      //----------------------------------------------------------------------------------

#ifdef OUTPUT_TO_DIFF_FILES
      os_subd.clear();
      os_subd.close();
#endif
      if(mpc.osdom)
      {
         //-----------------------------------------------------------
         /// VTK output
         ouput_vtk_part_info = true;
         // Elements in this subdomain
         f_iparts = fname+"_"+dom_str+"_of_"+s_nparts+"_subdomains.vtk";
         //f_iparts = fname+"_"+str_buf+".vtk";
         std::ofstream os(f_iparts.c_str(), std::ios::out|std::ios::trunc);
         WriteVTK_Head(os, sdom_nodes_l + sdom_nodes_h );

         WriteVTK_Nodes(os, sbd_nodes, sdom_start_node[idom], sdom_end_act_node[idom]);
         WriteVTK_Nodes(os, sbd_nodes, sdom_end_act_node[idom], sdom_end_node[idom]);
         WriteVTK_Nodes(os, sbd_nodes_hq, sdom_start_node_hq[idom], sdom_end_act_node_hq[idom]);
         WriteVTK_Nodes(os, sbd_nodes_hq, sdom_end_act_node_hq[idom], sdom_end_node_hq[idom]);

         writeSubdomainElementsVTK(os, in_subdom_elements, ghost_subdom_elements,  idom+1);

         /// Material data partitioning
         if( num_data > 0)
         {
            for(int mm = 0; mm<num_data; mm++)
            {
               // Partition
               os<<"SCALARS "<< m_headers[m_header_marker_per_data[mm]+4] <<" double 1\nLOOKUP_TABLE default" << std::endl;
               const MyInt e_shift = ne_total*mm;
               for(MyInt i=0; i<nei_size; i++)
               {
                  const size_t entry_index = in_subdom_elements[i]->getIndex() + e_shift;
                  os<<m_ele_val[entry_index] << std::endl;
               }
               for(MyInt i=0; i<ne_g; i++)
               {
                  const size_t entry_index = ghost_subdom_elements[i]->getIndex() + e_shift;
                  os<<m_ele_val[entry_index] << std::endl;
               }
            }
         }

         os.clear();
         os.close();
         //-----------------------------------------------------------

      }

      //sbd_nodes.clear();
      in_subdom_elements.clear();
      ghost_subdom_elements.clear();
   }

   //----------------------------------------------------------------------------------
   if (mpc.out_cct)
   {
      std::string cct_file_name = fname + "_" + s_nparts + ".cct";
      std::fstream cct_file(cct_file_name.c_str(), std::ios::out|std::ios::trunc );
      for (size_t idom=0; idom<vec_neighbors.size(); idom++)
      {
         size_t n_nei = 0;
         for (size_t ii=0; ii<vec_neighbors[idom].size(); ii++)
         {
            if (!vec_neighbors[idom][ii].empty()) n_nei++;
         }
         cct_file << "#COMMUNICATION_TABLE\n";
         cct_file << "  $MYRANK\n";
         cct_file << "    " << idom << "\n";
         cct_file << "  $NNEIGHBORS\n";
         cct_file << "    " << n_nei << "\n";
         for (size_t ii=0; ii<vec_neighbors[idom].size(); ii++)
         {
            std::set<ConnEdge> &edges = vec_neighbors[idom][ii];
            if (edges.empty()) continue;
            cct_file << "  $NEIGHBOR\n";
            cct_file << "    " << ii << "\n";
            cct_file << "    " << edges.size() << "\n";
            for (std::set<ConnEdge>::iterator itr=edges.begin(); itr!=edges.end(); ++itr)
            {
               cct_file << "    " << itr->first->local_index << " " << itr->second->local_index << "\n";
               //              cct_file << "    " << itr->first << " " << itr->second << "\n";
            }
         }
      }
      cct_file<<"#STOP" << std::endl;
      cct_file.close();
   }

   //----------------------------------------------------------------------------------

#ifdef OUTPUT_TO_SINGLE_FILE
   // Write nodes
   for(int idom=0; idom<num_parts; idom++)
   {

      const size_t start_l =  sdom_start_node[idom];
      const size_t end_l_act =  sdom_end_act_node[idom];
      const size_t end_l = sdom_end_node[idom];
      const size_t start_h =  sdom_start_node_hq[idom];
      const size_t end_h_act =  sdom_end_act_node_hq[idom];
      const size_t end_h = sdom_end_node_hq[idom];

      // For that in one long string buffer
      if(binary_output)
      {
         const size_t nnodes = end_l - start_l + end_h - start_h;
         std::vector<Node_Str> nodes_buffer(nnodes);

         size_t counter = 0;
         fillNodeVector4BinaryOuput(sbd_nodes, nodes_buffer, start_l, end_l_act, counter);
         fillNodeVector4BinaryOuput(sbd_nodes, nodes_buffer, end_l_act, end_l, counter);
         //
         fillNodeVector4BinaryOuput(sbd_nodes_hq, nodes_buffer, start_h, end_h_act, counter);
         fillNodeVector4BinaryOuput(sbd_nodes_hq, nodes_buffer, end_h_act, end_h, counter);

         fwrite(&nodes_buffer[0], sizeof(Node_Str), counter, of_bin_nod);
      }
      else
      {
         std::setw(14);
         os_subd_node.precision(14);
         //os_subd_node.setf(ios::fixed, ios::scientific);
         os_subd_node.setf(std::ios::scientific);

         writeSubDomainNodes(os_subd_node, sbd_nodes, start_l, end_l_act);
         writeSubDomainNodes(os_subd_node, sbd_nodes,  end_l_act, end_l);
         //
         writeSubDomainNodes(os_subd_node, sbd_nodes_hq, start_h, end_h_act);
         writeSubDomainNodes(os_subd_node, sbd_nodes_hq,  end_h_act, end_h);

         os_subd_node << std::endl;
      }
   }

   if(binary_output)
   {
      fclose(of_bin_cfg);
      fclose(of_bin_nod);
      fclose(of_bin_ele);
      fclose(of_bin_ele_g);
   }
   else
   {
      os_subd_head.clear();
      os_subd_head.close();
      os_subd.clear();
      os_subd.close();
      os_subd_node.clear();
      os_subd_node.close();
   }

#endif

   std::ofstream os;
   if(mpc.out_renum_gsmsh)
   {
      f_iparts = fname + "_renum_"+ s_nparts +".msh";
      os.open(f_iparts.c_str(), std::ios::out|std::ios::trunc);

      // Output renumbered mesh
      os<<"#FEM_MSH\n   $PCS_TYPE\n    NULL" << std::endl;
      os<<" $NODES\n"<<NodesNumber_Linear << std::endl;
      for(int idom=0; idom<num_parts; idom++)
      {
         const size_t start_l =  sdom_start_node[idom];
         const size_t end_l_act =  sdom_end_act_node[idom];
         for(size_t i=start_l; i<end_l_act; i++)
         {
            const Node* a_node = sbd_nodes[i];
            a_node->Write(os);
         }
      }

      os<<" $ELEMENTS\n"<<elem_vector.size() << std::endl;
      for(size_t e=0; e<elem_vector.size(); e++)
      {
         elem_vector[e]->WriteGSmsh(os);
      }
      os<<"#STOP" << std::endl;
      os.clear();
      os.close();
   }

   // Output VTK of the whole mesh
   if(mpc.is_vtk_out)
   {
      f_iparts = fname + "_renum_"+ s_nparts +".vtk";
      os.open(f_iparts.c_str(), std::ios::out|std::ios::trunc);
      ouput_vtk_part_info = false;

      os<<"# vtk DataFile Version 4.0\nGrid Partition by WW \nASCII\n" << std::endl;
      os<<"DATASET UNSTRUCTURED_GRID" << std::endl;
      os<<"POINTS "<< NodesNumber_Linear<<" double" << std::endl;
      std::setw(14);
      os.precision(14);

      for(int idom=0; idom<num_parts; idom++)
      {
         const size_t start_l =  sdom_start_node[idom];
         const size_t end_l_act =  sdom_end_act_node[idom];

         for(size_t i=start_l; i<end_l_act; i++)
         {
            Node* a_node = sbd_nodes[i];
            a_node->local_index = a_node->index; // Because WriteVTK_Elements_of_Subdomain uses local node IDs.
            os<<a_node->X()<<" "<<a_node->Y()<<" "<<a_node->Z() << std::endl;
         }
      }
      useQuadratic = false;
      WriteVTK_Elements_of_Subdomain(os, elem_vector, 0);

      os<<"POINT_DATA " << NodesNumber_Linear << std::endl;
   }

   os.clear();
   os.close();

   sbd_nodes.clear();
}

void Mesh::readGrid(const std::string& fname, const bool order)
{
   std::ifstream infile(fname.data());
   if(!infile.is_open())
   {
      std::cerr<<("Error: cannot open msh file . It may not exist !");
      abort();
   }

   bool rfiMesh = true;
   std::string line_string;
   getline(infile,line_string); // The first line
   if(line_string.find("#FEM_MSH")!=std::string::npos)
      rfiMesh = false;
   if(line_string.find("GeoSys-MSH")!=std::string::npos)
      rfiMesh = false;
   infile.seekg(0L, std::ios::beg);

   if(rfiMesh)
      ReadGrid(infile, order);
   else
      ReadGridGeoSys(infile, order);

}

// Read grid for test purpose
void Mesh::ReadGrid(std::istream& is, const bool high_order)
{
   MyInt i, ne, nn, counter;
   int ibuff;
   double x,y,z;
   std::string buffer;
//   is.seekg(position);
   // Read description
   is>>buffer>>std::ws;
   // Read numbers of nodes and elements
   is>>ibuff>>nn>>ne>>std::ws;
   if(nn==0||ne==0)
   {
      std::cout<<"Error: number of elements or nodes is zero" << std::endl;
      exit(1);
   }

   // Read Nodes
   counter = 0;
   node_vector.resize(nn);
   for(i=0; i<nn; i++)
   {
      is>>ibuff>>x>>y>>z>>std::ws;
      Node* newNode = new Node(ibuff,x,y,z);
      newNode->Marking(true);
      node_vector[counter] = newNode;
      counter++;
   }
   if(counter!=nn)
   {
      std::cout<<"Error: number nodes do not match" << std::endl;
      exit(1);
   }
   NodesNumber_Linear = nn;
   NodesNumber_Quadratic = nn;

   // Read Elements
   counter = 0;
   elem_vector.resize(ne);
   for(i=0; i<ne; i++)
   {
      Elem* newElem = new Elem(i);
      newElem->Read(is, this, 1, high_order);
      newElem->Marking(true);
      elem_vector[counter] = newElem;
      counter++;
   }
   if(counter!=ne)
   {
      std::cout<<"Error: number elements do not match" << std::endl;
      exit(1);
   }
//   position = is.tellg();
}

void Mesh::ReadGridGeoSys(std::istream& is, const bool high_order)
{
   std::string sub_line;
   std::string line_string;
   bool new_keyword = false;
   std::string hash("#");
   std::string sub_string,sub_string1;
   MyInt i, ibuff;
   MyInt no_elements;
   MyInt no_nodes;
   double x,y,z;
   Node* newNode = NULL;
   Elem* newElem = NULL;
   //========================================================================
   // Keyword loop
   while (!new_keyword)
   {
      //if(!GetLineFromFile(line,fem_file))
      //  break;
      //line_string = line;
      getline(is, line_string);
      if(is.fail())
         break;
      /*
      if(line_string.find(hash)!=string::npos)
      {
           new_keyword = true;
           break;
       }
      */
      //....................................................................
      //....................................................................
      if(line_string.find("$NODES") != std::string::npos)   // subkeyword found
      {
         is  >> no_nodes>>std::ws;

         node_vector.resize(no_nodes);
         for(i=0; i<no_nodes; i++)
         {
            is>>ibuff>>x>>y>>z>>std::ws;
            newNode = new Node(ibuff,x,y,z);
            node_vector[i] = newNode;
         }
         continue;
      }
      //....................................................................
      if(line_string.find("$ELEMENTS") != std::string::npos)   // subkeyword found
      {
         is >> no_elements>>std::ws;

         elem_vector.resize(no_elements);

         for(i=0; i<no_elements; i++)
         {
            newElem = new Elem(i);
            newElem->Read(is, this, 0, high_order);
            newElem->Marking(true);
            elem_vector[i] = newElem;
         }
         continue;
      }
   }
   //========================================================================
}

// 02.2012. WW
void  Mesh::WriteVTK_Nodes(std::ostream& os)
{
   size_t i;
   Node *a_node = NULL;

   os<<"# vtk DataFile Version 4.0\nGrid Partition by WW \nASCII\n" << std::endl;
   os<<"DATASET UNSTRUCTURED_GRID" << std::endl;
   os<<"POINTS "<<node_vector.size()<<" double" << std::endl;
   std::setw(14);
   os.precision(14);
   for(i=0; i<node_vector.size(); i++)
   {
      a_node = node_vector[i];
      os<<a_node->X()<<" "<<a_node->Y()<<" "<<a_node->Z() << std::endl;
   }
}

void Mesh::WriteVTK_Head(std::ostream& os, const size_t number_of_nodes)
{
   os<<"# vtk DataFile Version 4.0\nGrid Partition by WW \nASCII\n" << std::endl;
   os<<"DATASET UNSTRUCTURED_GRID" << std::endl;
   os<<"POINTS "<< number_of_nodes <<" double" << std::endl;

}

// 03.2012. WW
void Mesh::WriteVTK_Nodes(std::ostream& os, std::vector<Node*>& nod_vec, const size_t start, const size_t end)
{
   std::setw(14);
   os.precision(14);
   for(size_t i=start;  i<end; i++)
   {
      const double *x = nod_vec[i]->getCoordinates();
      os << x[0] << " " << x[1] << " " << x[2]  << std::endl;
   }
}

// 02.2012. WW
void  Mesh::WriteVTK_Elements_of_Subdomain(std::ostream& os, std::vector<Elem*>& ele_vec,
      const int sbd_index)
{
   size_t i;
   int k;
   int nne;

   //-----------------------------------------------------------
   //  VTK output
   // Elements in this subdomain
   size_t ne0 = ele_vec.size();
   size_t size = ne0;

   std::string deli = " ";

   Elem *a_elem = NULL;

   for(i=0; i<ne0; i++)
   {
      a_elem = ele_vec[i];
      nne =  a_elem->getNodesNumber(useQuadratic);
      if(useQuadratic&&a_elem->ele_Type == quadri)
         nne -= 1;

      size += nne;
   }
   os<<"\nCELLS "<<ne0<<deli<<size << std::endl;

   // CELLs
   for(i=0; i<ne0; i++)
   {
      a_elem = ele_vec[i];

      nne =  a_elem->getNodesNumber(useQuadratic);
      if(useQuadratic&&a_elem->ele_Type == quadri)
         nne -= 1;

      os<<nne<<deli;

      for(k=0; k<nne; k++)
      {
         os << a_elem->nodes[k]->local_index<< deli;
      }

      os << "\n";
   }
   os << std::endl;

   // CELL types
   os << "CELL_TYPES " << ne0 << std::endl;
   for(i=0; i<ne0; i++)
   {
      a_elem = ele_vec[i];
      a_elem->WriteVTK_Type(os, useQuadratic);
   }
   os << std::endl;

   if(ouput_vtk_part_info)
   {
      // Partition
      os<<"CELL_DATA "<<ne0 << std::endl;
      os<<"SCALARS Partition int 1\nLOOKUP_TABLE default" << std::endl;
      for(i=0; i<ne0; i++)
         os<<sbd_index << std::endl;
   }
}

void Mesh::Write2METIS(std::ostream& os)
{
   os<<(MyInt)elem_vector.size()<<" ";

#ifdef METIS4_0
   int e_type =0;
   switch(elem_vector[0]->getElementType())
   {
      case line:
         std::cout<<"Not for 1D element" << std::endl;
         exit(1);
      case quadri:
         e_type =4;
         break;
      case hex:
         e_type =3;
         break;
      case tri:
         e_type =1;
         break;
      case tet:
         e_type =2;
         break;
      case 6:
         std::cout<<"Not for prismal element" << std::endl;
         abort();
   }
   os<<e_type;
#endif
   os << std::endl;
   for(MyInt i=0; i<(MyInt)elem_vector.size(); i++)
   {
      elem_vector[i]->setOrder(useQuadratic);
      elem_vector[i]->Write_index(os);
   }
}

void Mesh::writeSubdomainElementsVTK(std::ostream& os, const std::vector<Elem*>& sdom_elems,
                                     const std::vector<Elem*>& sdom_elems_ghost, const int sbd_index)
{
   size_t i;
   int k;
   int nne;

   //-----------------------------------------------------------
   //  VTK output
   // Elements in this subdomain
   const size_t ne0 = sdom_elems.size();
   const size_t ne1 = sdom_elems_ghost.size();
   size_t size = ne0 + ne1;

   std::string deli = " ";

   Elem *a_elem = NULL;

   for(i=0; i<ne0; i++)
   {
      a_elem = sdom_elems[i];
      nne =  a_elem->getNodesNumber(useQuadratic);
      if(useQuadratic&&a_elem->ele_Type == quadri)
         nne -= 1;

      size += nne;
   }
   //
   for(i=0; i<ne1; i++)
   {
      a_elem = sdom_elems_ghost[i];
      nne =  a_elem->getNodesNumber(useQuadratic);
      if(useQuadratic&&a_elem->ele_Type == quadri)
         nne -= 1;

      size += nne;
   }
   os<<"\nCELLS "<<ne0 + ne1 << deli << size << std::endl;

   // CELLs
   for(i=0; i<ne0; i++)
   {
      a_elem = sdom_elems[i];

      nne =  a_elem->getNodesNumber(useQuadratic);
      if(useQuadratic&&a_elem->ele_Type == quadri)
         nne -= 1;

      os<<nne<<deli;

      for(k=0; k<nne; k++)
      {
         os << a_elem->nodes[k]->local_index<< deli;
      }

      os << "\n";
   }
   for(i=0; i<ne1; i++)
   {
      a_elem = sdom_elems_ghost[i];

      nne =  a_elem->getNodesNumber(useQuadratic);
      if(useQuadratic&&a_elem->ele_Type == quadri)
         nne -= 1;

      os<<nne<<deli;

      for(k=0; k<nne; k++)
      {
         os << a_elem->nodes[k]->local_index<< deli;
      }

      os << "\n";
   }
   os << std::endl;

   // CELL types
   os << "CELL_TYPES " << ne0  + ne1 << std::endl;
   for(i=0; i<ne0; i++)
   {
      a_elem = sdom_elems[i];
      a_elem->WriteVTK_Type(os, useQuadratic);
   }
   for(i=0; i<ne1; i++)
   {
      a_elem = sdom_elems_ghost[i];
      a_elem->WriteVTK_Type(os, useQuadratic);
   }
   os << std::endl;

   if(ouput_vtk_part_info)
   {
      // Partition
      os<<"CELL_DATA "<<ne0 + ne1  << std::endl;
      os<<"SCALARS Partition int 1\nLOOKUP_TABLE default" << std::endl;
      for(i=0; i<ne0; i++)
         os<<sbd_index << std::endl;
      const int id_gst = 0;
      for(i=0; i<ne1; i++)
         os<<id_gst << std::endl;
   }

}

void Mesh::writeVTK(const std::string& fname)
{
   std::ofstream os(fname.data(), std::ios::trunc);
   if (!os.is_open())
   {
      std::cout << "Could not open " << fname << std::endl;
      abort();
   }

   WriteVTK_Nodes(os);

   std::size_t ne = elem_vector.size();
   std::size_t size = ne;

   for(std::size_t i=0; i<ne; i++)
   {
      const Elem* a_elem = elem_vector[i];
      int nne =  a_elem->getNodesNumber(useQuadratic);
      if(useQuadratic&&a_elem->ele_Type == quadri)
         nne -= 1;
      size += nne;
   }

   os<<"\nCELLS " << ne << " " << size << std::endl;

   // CELLs
   for(std::size_t i=0; i<ne; i++)
   {
      const Elem* a_elem = elem_vector[i];

      int nne =  a_elem->getNodesNumber(useQuadratic);
      if(useQuadratic&&a_elem->ele_Type == quadri)
         nne -= 1;

      os << nne << " ";

      for(int k=0; k<nne; k++)
      {
         os << a_elem->nodes[k]->getIndex()<< " ";
      }

      os << "\n";
   }
   os << std::endl;

   // CELL types
   os << "CELL_TYPES " << ne << std::endl;
   for(std::size_t i=0; i<ne; i++)
   {
      elem_vector[i]->WriteVTK_Type(os, useQuadratic);
   }
   os << std::endl;

}

void Mesh::fillNodeVector4BinaryOuput(const std::vector<Node*> &sdom_nodes,
                                      std::vector<Node_Str> &sdom_nodes4bin,
                                      const size_t start, const size_t end, size_t &counter)
{
   for(size_t i=start; i<end; i++)
   {
      Node *a_node = sdom_nodes[i];

      Node_Str nd;
      nd.id = a_node->index;
      nd.x = a_node->Coordinate[0];
      nd.y = a_node->Coordinate[1];
      nd.z = a_node->Coordinate[2];
      sdom_nodes4bin[counter] = nd;
      counter++;
   }
}

void Mesh::writeSubDomainNodes(std::ostream& os, const std::vector<Node*>& sdom_nodes,
                               const size_t start, const size_t end)
{
   for(size_t i=start; i<end; i++)
   {
      Node *a_node = sdom_nodes[i];
      a_node->Write(os);
   }
}

void Mesh::writeBinary(const std::string& fname)
{
   std::ofstream os(fname.data(), std::ios::binary | std::ios::trunc);
   if (!os.is_open())
   {
      std::cout << "Could not open " << fname << std::endl;
      abort();
   }

   MyInt mesh_header[15];
   mesh_header[0] = NodesNumber_Linear;
   mesh_header[1] = NodesNumber_Quadratic;
   mesh_header[2] = static_cast<MyInt>(elem_vector.size());
   mesh_header[3] = axisymmetry;

   mesh_header[4] = coordinate_system;
   mesh_header[5] = max_ele_dim;

   mesh_header[6] = msh_no_line;
   mesh_header[7] = msh_no_quad;
   mesh_header[8] = msh_no_hexs;
   mesh_header[9] = msh_no_tris;
   mesh_header[10] = msh_no_tets;
   mesh_header[11] = msh_no_pris;
   mesh_header[12] = msh_no_pyra;
   mesh_header[13] = msh_max_dim;
   mesh_header[14] = useQuadratic;

   os.write( (char*)(mesh_header), 15 * sizeof(MyInt));

   for (std::size_t i=0; i<node_vector.size(); i++)
   {
      MyInt node_info[3];
      const Node* node = node_vector[i];
      node_info[0] = node->getIndex();
      node_info[1] = node->NodesRelated.size();
      node_info[2] = node->ElementsRelated.size();
      os.write( (char*)(node_info), 3 * sizeof(MyInt));
      os.write( (char*)(node->getCoordinates()), 3 * sizeof(double) );

      os.write( (char*)(&node->NodesRelated[0]),
                node_info[1] * sizeof(std::size_t) );
      os.write( (char*)(&node->ElementsRelated[0]),
                node_info[2] * sizeof(std::size_t) );
   }

   MyInt ids[100];
   // Write element nodes
   for (std::size_t i=0; i<elem_vector.size(); i++)
   {
      const Elem* elem = elem_vector[i];
      MyInt elem_info[3];
      elem_info[0] = elem->getIndex();
      elem_info[1] = elem->getPatchIndex();
      elem_info[2] = elem->getElementType();
      os.write( (char*)(elem_info), 3 * sizeof(MyInt));

      for (int j=0; j<elem->getNodesNumber(useQuadratic); j++)
      {
         ids[j] = elem->getNodeIndex(j);
      }
      // Write node index
      os.write( (char*)(ids), elem->getNodesNumber(useQuadratic)
                *sizeof(MyInt) );
   }
   // Write element neighbours
   for (std::size_t i=0; i<elem_vector.size(); i++)
   {
      const Elem* elem = elem_vector[i];
      for (int j=0; j<elem->getFacesNumber(); j++)
      {
         if (!elem->getNeighbor(j))
            ids[j] = -1;
         else
            ids[j] = elem->getNeighbor(j)->getIndex();
      }

      // Write neighbour element indexes
      os.write( (char*)(ids), elem->getFacesNumber() * sizeof(MyInt) );
   }
}

void Mesh::readBinary(const std::string& fname)
{
   std::ifstream ins(fname.data(), std::ios::binary);
   if (!ins.is_open())
   {
      std::cout << "Could not open " << fname << std::endl;
      abort();
   }

   MyInt mesh_header[15];
   ins.read( (char*)(mesh_header), 15 * sizeof(MyInt));

   NodesNumber_Linear = mesh_header[0];
   NodesNumber_Quadratic = mesh_header[1];
   MyInt num_elems = mesh_header[2];
   axisymmetry = static_cast<bool>(mesh_header[3]);

   coordinate_system = mesh_header[4];
   max_ele_dim = mesh_header[5];

   msh_no_line = mesh_header[6];
   msh_no_quad = mesh_header[7];
   msh_no_hexs = mesh_header[8];
   msh_no_tris = mesh_header[9];
   msh_no_tets = mesh_header[10];
   msh_no_pris = mesh_header[11];
   msh_no_pyra = mesh_header[12];
   msh_max_dim = mesh_header[13];
   useQuadratic = static_cast<bool>(mesh_header[14]);

   node_vector.resize(NodesNumber_Quadratic);
   double x[3];
   for (MyInt i=0; i<NodesNumber_Quadratic; i++)
   {
      MyInt node_info[3];
      ins.read( (char*)(node_info), 3 * sizeof(MyInt));
      ins.read( (char*)(x), 3 * sizeof(double) );
      Node* node = new Node(node_info[0], x[0], x[1], x[2]);

      node->NodesRelated.resize(node_info[1]);
      node->ElementsRelated.resize(node_info[2]);

      ins.read( (char*)(&node->NodesRelated[0]),
                node_info[1] * sizeof(std::size_t) );
      ins.read( (char*)(&node->ElementsRelated[0]),
                node_info[2] * sizeof(std::size_t) );
      node_vector[i] = node;
   }

   MyInt ids[100];
   // read element nodes
   elem_vector.resize(num_elems);
   for (MyInt i=0; i<num_elems; i++)
   {
      MyInt elem_info[3];
      ins.read( (char*)(elem_info), 3 * sizeof(MyInt));
      Elem* elem = new Elem(elem_info[0]);
      elem->PatchIndex = elem_info[1];
      elem->ele_Type = static_cast<ElemType>(elem_info[2]);

      // read node index
      ins.read( (char*)(ids), elem->getNodesNumber(useQuadratic)
                *sizeof(MyInt) );
      elem->nodes = new Node*[elem->getNodesNumber(useQuadratic)];
      for (int j=0; j<elem->getNodesNumber(useQuadratic); j++)
      {
         elem->nodes[j] = node_vector[ids[j]];
      }
      elem_vector[i] = elem;
   }
   // read element neighbours
   for(std::size_t i=0; i<elem_vector.size(); i++)
   {
      Elem* elem = elem_vector[i];
      // Write neighbour element indexes
      ins.read( (char*)(ids), elem->getFacesNumber() * sizeof(MyInt) );
      elem->neighbors = new Elem*[elem->getFacesNumber()];
      for (int j=0; j<elem->getFacesNumber(); j++)
      {
         if (ids[j] == -1)
            elem->neighbors[j] = NULL;
         else
            elem->neighbors[j] = elem_vector[ids[j]];
      }
   }
}

}//end namespace


