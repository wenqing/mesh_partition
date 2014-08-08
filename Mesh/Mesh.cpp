#include "Mesh.h"

#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <limits>
#include <set>

#include <stdio.h> // for binary output

#include "Node.h"

//#define BUILD_MESH_EDGE

//------------------------------------------------------
//   Topology definition of geometrical element.
//    WW. 10.01.2005
//------------------------------------------------------
namespace Mesh_Group
{
using namespace std;

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
   cout<<"\nCPU time elapsed in constructing topology of grids: "
       <<(double)(finish - start) / CLOCKS_PER_SEC<<"s"<<endl<<endl;

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
   cout<<"\n\tCPU time elapsed in generating high oder elements: "
       <<(double)(finish - start) / CLOCKS_PER_SEC<<"s"<<endl;

}

void Mesh::ConstructSubDomain_by_Elements(const string fname, const int num_parts, const bool osdom)
{
   string str;
   string stro;
   int dom;
   int max_dom;
   int k,kk;
   MyInt i,j;
   //  int ntags = 3;

   fstream gmsh_out;

   string deli = " ";
   //

   string s_nparts;
   stringstream ss;
   ss << num_parts;
   ss >> s_nparts;
   ss.clear();

   str = fname + ".mesh.epart." + s_nparts;
   stro = fname + "." + s_nparts +"ddc";

   //namef = ".mesh.epart."; //+str_buf;
   ifstream part_in;
   fstream part_out;
   part_out.open(stro.c_str(), ios::out | ios::trunc );
   // Output for gmsh

   if(osdom)
   {
      stro = fname + "_gmsh.msh";
      gmsh_out.open(stro.c_str(), ios::out );
      //gmsh_out<<"$NOD"<<endl;
      gmsh_out<<"$MeshFormat\n2 0 8\n$EndMeshFormat\n$Nodes"<<endl;
      gmsh_out<<node_vector.size()<<endl;
      Node *node;
      for(i=0; i<(MyInt)node_vector.size(); i++)
      {
         gmsh_out<<i+1<<" ";
         node = node_vector[i];
         gmsh_out<<node->X()<<" ";
         gmsh_out<<node->Y()<<" ";
         gmsh_out<<node->Z()<<endl;
      }
      //gmsh_out<<"$ENDNOD"<<endl;
      //gmsh_out<<"$ELM"<<endl;
      gmsh_out<<"$EndNodes\n$Elements"<<endl;
      gmsh_out<<(MyInt)elem_vector.size()<<endl;
   }

   //
   part_in.open(str.c_str());
   if(!part_in.is_open())
   {
      cerr<<("Error: cannot open .epart file . It may not exist !");
      abort();
   }

   max_dom=0;
   Elem *ele = NULL;

   for(i=0; i<(MyInt)elem_vector.size(); i++)
   {
      part_in>>dom>>ws;
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
   remove(str.c_str());

   if(osdom)
   {
      gmsh_out<<"$EndElements"<<endl;
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
   vector<size_t> nodes_dom;
   vector<Elem*> eles_dom;

   //
   for(k=0; k<max_dom; k++)
   {
      part_out<<"#DOMAIN "<<k<<endl;
      part_out<<"$ELEMENTS "<<ele_dom[k]<<endl;
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
            part_out<<ele->getIndex()<<endl;
            eles_dom.push_back(ele); //TEST OUT
         }
      }
      part_out<<"$NODES_INNER "<<(MyInt)nodes_dom.size()<<endl;
      for(j=0; j<(MyInt)nodes_dom.size(); j++)
         part_out<<nodes_dom[j]<<endl;

      if(osdom)
      {
         string i_nparts;
         ss << k;
         ss >> i_nparts;
         ss.clear();

         string name_f = fname+"_"+i_nparts+"_of_"+s_nparts+"subdomains.msh";
         fstream test_out;
         test_out.open(name_f.c_str(), ios::out|ios::trunc );

         Node *nod = 0;
         //GMSH test_out<<"$NOD"<<endl;
         //GMSH test_out<<(MyInt)nodes_dom.size()<<endl;
         test_out<<"#0#0#0#1#0.0#0#################################################################"<<endl;
         test_out<<"0 "<<(MyInt)nodes_dom.size()<<" "<<(MyInt)eles_dom.size()<<endl;
         for(j=0; j<(MyInt)nodes_dom.size(); j++)
         {
            nod = node_vector[nodes_dom[j]];
            //GMSH  test_out<<j+1<<"  "
            test_out<<j<<deli
                    << nod->X()<<deli<< nod->Y()<<deli<< nod->Z() <<endl;
         }
         //GMSH test_out<<"$ENDNOD"<<endl;
         //GMSH test_out<<"$ELE"<<endl;
         //GMSG test_out<<(MyInt)eles_dom.size()<<endl;
         for(j=0; j<(MyInt)eles_dom.size(); j++)
         {
            ele = eles_dom[j];

            //GMSH  ele->WriteGmsh(test_out, k+1);

            test_out<<j<<deli<<ele->getPatchIndex()<<deli<<ele->getName()<<deli;
            for(kk=0; kk<ele->getNodesNumber(); kk++)
               test_out<< ele->getDomNodeIndex(kk)<<deli;
            test_out<<endl;
         }
         test_out.clear();
         test_out.close();
      }

   }
   part_out<<"#STOP "<<endl;
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
   vector<string> m_headers;
   vector<size_t> m_header_marker_per_data;
   vector<string> m_datanames;
   //vector<MyInt> m_ele_idx;
   vector<double> m_ele_val;

   if(mpc.mat_fname.size() !=0)
   {
      string line_buffer;
      string mat_fname_abs = mpc.fpath + mpc.mat_fname;

      ifstream is_mat(mat_fname_abs.c_str());
      if(!is_mat.good())
      {
         cout<<"Material data file "<<mat_fname_abs<<" does not exist"<<endl;
         exit(1);
      }
      is_mat >> num_data;
      m_datanames.resize(num_data);
      m_header_marker_per_data.resize(num_data +1 );
      m_header_marker_per_data[0] = 0;
      for(int k=0; k<num_data; k++)
      {
         string data_name;
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
            cout<<"Material data file "<<m_datanames[k]<<" does not exist"<<endl;
            exit(1);
         }
         while(!is_mat.eof())
         {
            getline(is_mat, line_buffer);
            if(line_buffer.find("$DATA")!=string::npos)
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
            else if(line_buffer.find("#STOP")!=string::npos)
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

   string f_iparts;
   string o_part_msh;
   string fname = mpc.fname;

   // Number of integer variables of subdomain elements
   MyInt nmb_element_idxs;
   MyInt nmb_element_idxs_g;

   string deli = " ";

   //
   const MyInt num_headers = 13;
   MyInt head[13];
   MyInt ivar[23];
   // For binary output: all data in on string
   vector<MyInt> ele_info;
   // node rank offset
   head[10] = 0;
   // element rank offset
   head[11] = 0;
   // ghost element rank offset
   head[12] = 0;

   const int num_parts = mpc.num_parts;
   // Read a METIS node partitioned file
   //--------------------------------------------------------------------------
   const string s_nparts = number2str(num_parts);
   f_iparts = fname + ".mesh.npart." + s_nparts;
   //o_part_msh = fname + "." + s_nparts +"mesh";

   ifstream npart_in(f_iparts.c_str());
   if(!npart_in.is_open())
   {
      cerr<<("Error: cannot open .npart file . It may not exist !");
      exit(1);
   }

   const MyInt nn = static_cast<MyInt>(node_vector.size());

   vector<bool> sdom_marked(nn);
   vector<MyInt> dom_idx(NodesNumber_Linear);

   // Re-ordered nodes of the whole mesh for ouput
   MyInt dom;
   for(MyInt i=0; i<NodesNumber_Linear; i++)
   {
      npart_in>>dom>>ws;
      dom_idx[i] = dom;
      sdom_marked[i] = false;
   }
   npart_in.close();
   //remove(f_iparts.c_str());

   //--------------------------------------------------------------------------
   // Output a partitioned mesh
   //--------------------------------------------------------------------------
#define OUTPUT_TO_SINGLE_FILE
#ifdef OUTPUT_TO_SINGLE_FILE
   string name_f;
   fstream os_subd;
   fstream os_subd_head;
   fstream os_subd_node;
   FILE *of_bin_cfg = 0; // ostream os_bin_cfg
   FILE *of_bin_nod = 0; // ostream os_bin_nod
   FILE *of_bin_ele = 0; // ostream os_bin_ele
   FILE *of_bin_ele_g = 0; // ostream os_bin_ele_g
   const bool binary_output = mpc.binary_output;

   if(!binary_output)
   {
      name_f = fname+"_partitioned_cfg"+ s_nparts + ".msh";
      os_subd_head.open(name_f.c_str(), ios::out|ios::trunc );
      name_f = "Subdomain mesh "
               "(Nodes;  Nodes_linear; Elements; Ghost elements; Nodes of Linear elements; Nodes of quadratic elements) "
               "Nodes of Linear whole elements; Nodes of whole quadratic elements; "
               "Total integer variables of elements;Total integer variables of ghost elements  ";
      os_subd_head << name_f << endl;
      os_subd_head<<num_parts<<endl;

      name_f = fname+"_partitioned_elems_"+ s_nparts + ".msh";
      os_subd.open(name_f.c_str(), ios::out|ios::trunc );

      name_f = fname+"_partitioned_nodes_"+ s_nparts + ".msh";
      os_subd_node.open(name_f.c_str(), ios::out|ios::trunc );

      setw(14);
      os_subd.precision(14);
      //os_subd.setf(ios::fixed, ios::scientific);
      os_subd.setf(ios::scientific);
   }
   else
   {
      name_f = fname+"_partitioned_msh_cfg"+ s_nparts + ".bin";
      //os_bin_cfg.open(name_f.c_str(), ios::binary|ios::out|ios::trunc );
      of_bin_cfg = fopen(name_f.c_str(), "wb");

      name_f = fname+"_partitioned_msh_ele"+ s_nparts + ".bin";
      //os_bin_ele.open(name_f.c_str(), ios::binary|ios::out|ios::trunc );
      of_bin_ele = fopen(name_f.c_str(), "wb");

      name_f = fname+"_partitioned_msh_ele_g"+ s_nparts + ".bin";
      //os_bin_ele.open(name_f.c_str(), ios::binary|ios::out|ios::trunc );
      of_bin_ele_g = fopen(name_f.c_str(), "wb");

      name_f = fname+"_partitioned_msh_nod"+ s_nparts + ".bin";
      //os_bin_nod.open(name_f.c_str(), ios::binary|ios::out|ios::trunc );
      of_bin_nod = fopen(name_f.c_str(), "wb");

   }

#endif

   // for linear
   vector<size_t> sdom_start_node(num_parts);
   vector<size_t> sdom_end_act_node(num_parts);
   vector<size_t> sdom_end_node(num_parts);

   // quadratic
   vector<size_t> sdom_start_node_hq(num_parts);
   vector<size_t> sdom_end_act_node_hq(num_parts);
   vector<size_t> sdom_end_node_hq(num_parts);

   const MyInt ne_total = static_cast<MyInt>(elem_vector.size());

   vector<Node*> sbd_nodes;        // nodes of a partition for linear element
   vector<Node*> sbd_nodes_hq;     // nodes of a partition for quadratic element
   std::vector<std::vector<std::set<ConnEdge> > > vec_neighbors(num_parts); //NW

   MyInt node_id_offset = 0;
   MyInt node_id_offset_h = NodesNumber_Linear;
   Node *a_node = NULL;
   Elem *a_elem = NULL;
   for(int idom=0; idom<num_parts; idom++)
   {

      sdom_start_node[idom] = sbd_nodes.size();
      sdom_start_node_hq[idom] = sbd_nodes_hq.size();

      if (mpc.out_cct)
      {
         vec_neighbors[idom].resize(num_parts);
      }
      cout << "Process partition: " << idom << endl;

      nmb_element_idxs = 0;
      nmb_element_idxs_g = 0;

      // Un-making all nodes of the whole mesh
      for(MyInt j=0; j<nn; j++)
         node_vector[j]->Marking(false);

      MyInt num_nodes_active_l = 0; // acitve nodes for linear elements
      MyInt num_nodes_active_h = 0; // acitve nodes for quadratic elements
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

      // Un-making all elements of the whole mesh
      for(MyInt j=0; j<ne_total; j++)
         elem_vector[j]->Marking(false);

      /// Find the elements in this subdomain.
      vector<Elem*> in_subdom_elements;
      vector<Elem*> ghost_subdom_elements;
      for(MyInt j=0; j<num_nodes_active_l; j++)
      {
         a_node = sbd_nodes[ j +  sdom_start_node[idom] ];

         // Search the elements connected to this nodes
         const MyInt ne_rel = static_cast<MyInt>(a_node->ElementsRelated.size());
         for(int k=0; k<ne_rel; k++)
         {
            a_elem = elem_vector[a_node->ElementsRelated[k]];

            // If checked
            if(a_elem->getStatus())
               continue;


            vector<int> ng_nodes; // non ghost nodes in ghost elements
            vector<int> g_nodes; // ghost nodes in ghost elements
            for(int kk=0; kk<a_elem->getNodesNumber(); kk++)
            {
               if(a_elem->getNode(kk)->getStatus())
               {
                  ng_nodes.push_back(kk);
               }
               else
               {
                  g_nodes.push_back(kk);
               }
            }

            // All nodes of this element are inside this subdomain
            if(g_nodes.size() == 0)
            {
               in_subdom_elements.push_back(a_elem);
            }
            else if(g_nodes.size() != static_cast<size_t>(a_elem->getNodesNumber()))
            {
               ghost_subdom_elements.push_back(a_elem);

               const int nn_gl = static_cast<int>(ng_nodes.size());
               a_elem->nnodes_gl = nn_gl;
               a_elem->ghost_nodes.resize(nn_gl);
               for(int kk=0; kk<nn_gl; kk++)
                  a_elem->ghost_nodes[kk] = ng_nodes[kk];

            }

            // overlapping edges
            if (mpc.out_cct)
            {
               for (size_t ig=0; ig<g_nodes.size(); ig++)
               {
                  MyInt ig_id = a_elem->getNode(g_nodes[ig])->global_index; //index; //;
                  MyInt ig_dom = dom_idx[ig_id];
                  for (size_t ing=0; ing<ng_nodes.size(); ing++)
                  {
                     //MyInt ing_id = a_elem->getNode(ng_nodes[ing])->index; //global_index;
                     vec_neighbors[idom][ig_dom].insert(ConnEdge(a_elem->getNode(ng_nodes[ing]), a_elem->getNode(g_nodes[ig]))); // inner node - ghost node
                     //vec_neighbors[ig_dom].insert(ConnEdge(ing_id, ig_id)); // inner node - ghost node
                  }
               }
            }
            a_elem->Marking(true);

            ng_nodes.clear();
            g_nodes.clear();

         }
      }

      // For quadratic element, add additional nodes here
      // Add  non ghost nodes in ghost elements as well
      if(useQuadratic)
      {
         MyInt nei = static_cast<MyInt>(in_subdom_elements.size());
         for(MyInt j=0; j<nei; j++)
         {
            a_elem = in_subdom_elements[j];
            const int nnodes = a_elem->getNodesNumber();
            const int nnodesHQ = a_elem->getNodesNumberHQ();

            for(int k=nnodes; k<nnodesHQ; k++)
               a_elem->nodes[k]->Marking(false);
         }

         MyInt neg = static_cast<MyInt>(ghost_subdom_elements.size());
         for(MyInt j=0; j<neg; j++)
         {
            a_elem = ghost_subdom_elements[j];
            const int nnodesHQ = a_elem->getNodesNumberHQ();

            for(int k=0; k<nnodesHQ; k++)
               a_elem->nodes[k]->Marking(false);
         }

         // Add nodes for quadrtic elements in this subdomain zone
         for(MyInt j=0; j<nei; j++)
         {
            a_elem = in_subdom_elements[j];
            const int nnodes = a_elem->getNodesNumber();
            const int nnodesHQ = a_elem->getNodesNumberHQ();

            for(int k=nnodes; k<nnodesHQ; k++)
            {
               a_node = a_elem->nodes[k];
               const MyInt i = static_cast<MyInt>(a_elem->nodes[k]->index_org);
               if(sdom_marked[i]) // Already in other subdomains
                  continue;

               if(a_node->getStatus()) // Already added
                  continue;

               sdom_marked[i] = true;

               sbd_nodes_hq.push_back(a_node);
               a_node->Marking(true);
               num_nodes_active_h++;
            }
         }

         //-------------------------------------------
         // Check ghost elements
         // Add nodes for quadrtic elements in ghost zone
         for(MyInt j=0; j<neg; j++)
         {
            a_elem = ghost_subdom_elements[j];
            const int nnodes = a_elem->getNodesNumber();
            const int nnodesHQ = a_elem->getNodesNumberHQ();

            for(int k=nnodes; k<nnodesHQ; k++)
            {
               a_node = a_elem->nodes[k];
               // Since a_elem->nodes_index[k] is not touched
               const MyInt i = static_cast<MyInt>(a_elem->nodes[k]->index_org);
               if(sdom_marked[i]) // Already in other subdomains
                  continue;

               if(a_node->getStatus()) // Already added
                  continue;

               sdom_marked[i] = true;

               sbd_nodes_hq.push_back(a_node);
               a_node->Marking(true);
               num_nodes_active_h++;
            }
         }

         // Make non-ghost nodes in ghost elements
         for(MyInt j=0; j<neg; j++)
         {
            a_elem = ghost_subdom_elements[j];
            const int nnodesHQ = a_elem->getNodesNumberHQ();

            for(int k=0; k<nnodesHQ; k++)
            {
               if(a_elem->nodes[k]->getStatus())
               {
                  a_elem->ghost_nodes.push_back(k);
               }
            }

         }
      }

      sdom_end_act_node_hq[idom] = sbd_nodes_hq.size();

      //-----------------------------------------------
      // Add inactive nodes in ghost elements
      const MyInt ne_g = static_cast<MyInt>(ghost_subdom_elements.size());
      for(MyInt j=0; j<ne_g; j++)
      {
         a_elem = ghost_subdom_elements[j];
         for(int k=0; k<a_elem->getNodesNumber(useQuadratic); k++)
            a_elem->nodes[k]->Marking(false);

         // Existing nodes
         const int ngh_nodes = static_cast<int>(a_elem->ghost_nodes.size());
         for(int k=0; k<ngh_nodes; k++)
            a_elem->nodes[a_elem->ghost_nodes[k]]->Marking(true);
      }
      //
      for(MyInt j=0; j<ne_g; j++)
      {
         a_elem = ghost_subdom_elements[j];
         for(int k=0; k<a_elem->getNodesNumber(); k++)
         {
            a_node = a_elem->nodes[k];
            if(a_node->getStatus())
               continue;
            a_node->Marking(true);
            sbd_nodes.push_back(a_node);
         }
         //
         for(int k=a_elem->getNodesNumber(); k<a_elem->getNodesNumber(useQuadratic); k++)
         {
            a_node = a_elem->nodes[k];
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
         nmb_element_idxs_g += static_cast<MyInt>(ghost_subdom_elements[j]->ghost_nodes.size());
      }

      string dom_str = number2str(idom);
#ifdef OUTPUT_TO_DIFF_FILES
      // Make output of this subdomain for simulation
      //string name_f = fname+"_"+dom_str+"_of_"+s_nparts+"_subdomains.msh";
      string name_f = fname+"_"+dom_str+".msh";
      fstream os_subd(name_f.c_str(), ios::out|ios::trunc );
      //os_subd<<"#FEM_MSH\n   $PCS_TYPE\n    NULL"<<endl;
      //os_subd<<" $NODES\n"<<size_sbd_nodes<<endl;


      name_f = "Subdomain mesh "
               "(Nodes;  Nodes_linear; Elements; Ghost elements; Nodes of Linear elements; Nodes of quadratic elements) "
               "Nodes of Linear whole elements; Nodes of whole quadratic elements; "
               "Total integer variables of elements;Total integer variables of ghost elements  ";
      os_subd<<name_f<<endl;
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
                     <<deli<<nmb_element_idxs<<deli<<nmb_element_idxs_g<<endl;
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

      //os_subd<<"Elements"<<endl;
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

      //os_subd<<"Ghost elements"<<endl;
      counter = 0;
      if(binary_output)
      {
         ele_info.clear();
         ele_info.resize(offset_e_g);
      }
      counter = ne_g;

      for(MyInt j=0; j<ne_g; j++)
      {
         a_elem = ghost_subdom_elements[j];
         const int ngh_nodes = static_cast<int>(a_elem->ghost_nodes.size());

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

            ivar[0] = a_elem->nnodes_gl;
            ivar[1] = ngh_nodes;
            for(int kk=0; kk<ngh_nodes; kk++)
            {
               ivar[kk+2] = a_elem->ghost_nodes[kk];
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
            os_subd<<a_elem->nnodes_gl<<deli<<ngh_nodes<<deli;
            for(int kk=0; kk<ngh_nodes; kk++)
            {
               os_subd<<a_elem->ghost_nodes[kk]<<deli;
            }
            os_subd<<endl;
         }
      }

      os_subd << endl;

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
         ofstream os_mat;
         for(int mm = 0; mm<num_data; mm++)
         {
            string mat_ofile_name = m_datanames[mm] + dom_str;
            os_mat.open(mat_ofile_name.c_str(), ios::trunc);
            for(size_t mh = m_header_marker_per_data[mm]; mh<m_header_marker_per_data[mm+1]; mh++ )
            {
               os_mat<<m_headers[mh]<<endl;
            }

            const MyInt e_shift = ne_total*mm;
            for(MyInt j=0; j<nei_size; j++)
            {
               const size_t entry_index = in_subdom_elements[j]->getIndex() + e_shift;
               os_mat<<j<<deli<<m_ele_val[entry_index]<<endl;
            }
            for(MyInt j=0; j<ne_g; j++)
            {
               const size_t entry_index  = ghost_subdom_elements[j]->getIndex() + e_shift;
               os_mat<<j+nei_size<<deli<<m_ele_val[entry_index]<<endl;
            }
            os_mat<<"#STOP"<<endl;
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
         ofstream os(f_iparts.c_str(), ios::out|ios::trunc);
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
               os<<"SCALARS "<< m_headers[m_header_marker_per_data[mm]+4] <<" double 1\nLOOKUP_TABLE default"<<endl;
               const MyInt e_shift = ne_total*mm;
               for(MyInt i=0; i<nei_size; i++)
               {
                  const size_t entry_index = in_subdom_elements[i]->getIndex() + e_shift;
                  os<<m_ele_val[entry_index]<<endl;
               }
               for(MyInt i=0; i<ne_g; i++)
               {
                  const size_t entry_index = ghost_subdom_elements[i]->getIndex() + e_shift;
                  os<<m_ele_val[entry_index]<<endl;
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
      string cct_file_name = fname + "_" + s_nparts + ".cct";
      fstream cct_file(cct_file_name.c_str(), ios::out|ios::trunc );
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
      cct_file<<"#STOP"<<endl;
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
         vector<Node_Str> nodes_buffer(nnodes);

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
         os_subd_node.setf(ios::scientific);

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

   ofstream os;
   if(mpc.out_renum_gsmsh)
   {
      f_iparts = fname + "_renum_"+ s_nparts +".msh";
      os.open(f_iparts.c_str(), ios::out|ios::trunc);

      // Output renumbered mesh
      os<<"#FEM_MSH\n   $PCS_TYPE\n    NULL"<<endl;
      os<<" $NODES\n"<<NodesNumber_Linear<<endl;
      for(int idom=0; idom<num_parts; idom++)
      {
         const size_t start_l =  sdom_start_node[idom];
         const size_t end_l_act =  sdom_end_act_node[idom];
         for(size_t i=start_l; i<end_l_act; i++)
         {
            a_node = sbd_nodes[i];
            a_node->Write(os);
         }
      }

      os<<" $ELEMENTS\n"<<elem_vector.size()<<endl;
      for(size_t e=0; e<elem_vector.size(); e++)
      {
         elem_vector[e]->WriteGSmsh(os);
      }
      os<<"#STOP"<<endl;
      os.clear();
      os.close();
   }

   // Output VTK of the whole mesh
   if(mpc.is_vtk_out)
   {
      f_iparts = fname + "_renum_"+ s_nparts +".vtk";
      os.open(f_iparts.c_str(), ios::out|ios::trunc);
      ouput_vtk_part_info = false;

      os<<"# vtk DataFile Version 4.0\nGrid Partition by WW \nASCII\n"<<endl;
      os<<"DATASET UNSTRUCTURED_GRID"<<endl;
      os<<"POINTS "<< NodesNumber_Linear<<" double"<<endl;
      setw(14);
      os.precision(14);

      for(int idom=0; idom<num_parts; idom++)
      {
         const size_t start_l =  sdom_start_node[idom];
         const size_t end_l_act =  sdom_end_act_node[idom];

         for(size_t i=start_l; i<end_l_act; i++)
         {
            a_node = sbd_nodes[i];
            a_node->local_index = a_node->index; // Because WriteVTK_Elements_of_Subdomain uses local node IDs.
            os<<a_node->X()<<" "<<a_node->Y()<<" "<<a_node->Z()<<endl;
         }
      }
      useQuadratic = false;
      WriteVTK_Elements_of_Subdomain(os, elem_vector, 0);

      os<<"POINT_DATA " << NodesNumber_Linear << endl;
   }

   os.clear();
   os.close();

   sbd_nodes.clear();
}

// 02.2012. WW
void  Mesh::WriteVTK_Nodes(std::ostream& os)
{
   size_t i;
   Node *a_node = NULL;

   os<<"# vtk DataFile Version 4.0\nGrid Partition by WW \nASCII\n"<<endl;
   os<<"DATASET UNSTRUCTURED_GRID"<<endl;
   os<<"POINTS "<<node_vector.size()<<" double"<<endl;
   setw(14);
   os.precision(14);
   for(i=0; i<node_vector.size(); i++)
   {
      a_node = node_vector[i];
      os<<a_node->X()<<" "<<a_node->Y()<<" "<<a_node->Z()<<endl;
   }
}

void Mesh::WriteVTK_Head(std::ostream& os, const size_t number_of_nodes)
{
   os<<"# vtk DataFile Version 4.0\nGrid Partition by WW \nASCII\n"<<endl;
   os<<"DATASET UNSTRUCTURED_GRID"<<endl;
   os<<"POINTS "<< number_of_nodes <<" double"<<endl;

}

// 03.2012. WW
void Mesh::WriteVTK_Nodes(std::ostream& os, std::vector<Node*>& nod_vec, const size_t start, const size_t end)
{
   setw(14);
   os.precision(14);
   for(size_t i=start;  i<end; i++)
   {
      const double *x = nod_vec[i]->getCoordinates();
      os << x[0] << " " << x[1] << " " << x[2] <<endl;
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

   string deli = " ";

   Elem *a_elem = NULL;

   for(i=0; i<ne0; i++)
   {
      a_elem = ele_vec[i];
      nne =  a_elem->getNodesNumber(useQuadratic);
      if(useQuadratic&&a_elem->ele_Type == quadri)
         nne -= 1;

      size += nne;
   }
   os<<"\nCELLS "<<ne0<<deli<<size<<endl;

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
   os << endl;

   // CELL types
   os << "CELL_TYPES " << ne0 << endl;
   for(i=0; i<ne0; i++)
   {
      a_elem = ele_vec[i];
      a_elem->WriteVTK_Type(os, useQuadratic);
   }
   os << endl;

   if(ouput_vtk_part_info)
   {
      // Partition
      os<<"CELL_DATA "<<ne0<<endl;
      os<<"SCALARS Partition int 1\nLOOKUP_TABLE default"<<endl;
      for(i=0; i<ne0; i++)
         os<<sbd_index<<endl;
   }
}

void Mesh::Write2METIS(ostream& os)
{
   os<<(MyInt)elem_vector.size()<<" ";

#ifdef METIS4_0
   int e_type =0;
   switch(elem_vector[0]->getElementType())
   {
      case line:
         cout<<"Not for 1D element"<<endl;
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
         cout<<"Not for prismal element"<<endl;
         abort();
   }
   os<<e_type;
#endif
   os<<endl;
   for(MyInt i=0; i<(MyInt)elem_vector.size(); i++)
      elem_vector[i]->Write_index(os);
}

// Read grid for test purpose
void Mesh::ReadGrid(istream& is, const bool high_order)
{
   MyInt i, ne, nn, counter;
   int ibuff;
   double x,y,z;
   string buffer;
//   is.seekg(position);
   // Read description
   is>>buffer>>ws;
   // Read numbers of nodes and elements
   is>>ibuff>>nn>>ne>>ws;
   if(nn==0||ne==0)
   {
      cout<<"Error: number of elements or nodes is zero"<<endl;
      exit(1);
   }

   // Read Nodes
   counter = 0;
   node_vector.resize(nn);
   for(i=0; i<nn; i++)
   {
      is>>ibuff>>x>>y>>z>>ws;
      Node* newNode = new Node(ibuff,x,y,z);
      newNode->Marking(true);
      node_vector[counter] = newNode;
      counter++;
   }
   if(counter!=nn)
   {
      cout<<"Error: number nodes do not match"<<endl;
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
      cout<<"Error: number elements do not match"<<endl;
      exit(1);
   }
//   position = is.tellg();
}

void Mesh::ReadGridGeoSys(istream& is, const bool high_order)
{
   string sub_line;
   string line_string;
   bool new_keyword = false;
   string hash("#");
   string sub_string,sub_string1;
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
      if(line_string.find("$NODES")!=string::npos)   // subkeyword found
      {
         is  >> no_nodes>>ws;

         node_vector.resize(no_nodes);
         for(i=0; i<no_nodes; i++)
         {
            is>>ibuff>>x>>y>>z>>ws;
            newNode = new Node(ibuff,x,y,z);
            node_vector[i] = newNode;
         }
         continue;
      }
      //....................................................................
      if(line_string.find("$ELEMENTS")!=string::npos)   // subkeyword found
      {
         is >> no_elements>>ws;

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

   string deli = " ";

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
   os<<"\nCELLS "<<ne0 + ne1 << deli << size << endl;

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
   os << endl;

   // CELL types
   os << "CELL_TYPES " << ne0  + ne1<< endl;
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
   os << endl;

   if(ouput_vtk_part_info)
   {
      // Partition
      os<<"CELL_DATA "<<ne0 + ne1 <<endl;
      os<<"SCALARS Partition int 1\nLOOKUP_TABLE default"<<endl;
      for(i=0; i<ne0; i++)
         os<<sbd_index<<endl;
      const int id_gst = 0;
      for(i=0; i<ne1; i++)
         os<<id_gst<<endl;
   }

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
}//end namespace


