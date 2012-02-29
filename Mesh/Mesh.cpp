#include "Mesh.h"

#include <iomanip>
#include <limits>
#include <list>

#include "Node.h"
#include "Edge.h"
#include "Elem.h"


//------------------------------------------------------ 
//   Topology definition of geometrical element.
//    WW. 10.01.2005
//------------------------------------------------------ 
namespace Mesh_Group
{
using namespace std;
using namespace Math_Group;


Mesh::Mesh(bool quad)
{
   useQuadratic = quad;
   coordinate_system = 1;
   axisymmetry=false;
   max_ele_dim = 0; 

}

Mesh::~Mesh()
{
  long i;
  // Nodes
  for(i=0; i<(long)node_vector.size(); i++)
    delete node_vector[i];
  node_vector.clear();
  // Edges
  for(i=0; i<(long)edge_vector.size(); i++)
     delete edge_vector[i];
  edge_vector.clear();
  // Surface faces
  for(i=0; i<(long)face_vector.size(); i++)
     delete face_vector[i];
  face_vector.clear();
  // Element
  for(i=0; i<(long)elem_vector.size(); i++)
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
  int i, j,l, k, n;
  Node* m_nod = NULL;
  Elem* m_ele = NULL;
  bool exist = false;
  //----------------------------------------------------------------------
  for(i=0;i<(long)node_vector.size();i++){
    m_nod = node_vector[i];
    for(j=0;j<(int)m_nod->ElementsRelated.size();j++){
      m_ele = elem_vector[m_nod->ElementsRelated[j]];
      for(l=0;l<m_ele->getNodesNumber(quadratic);l++){
          exist = false;
          for(k=0;k<(int)m_nod->NodesRelated.size();k++) 
          {
             if(m_nod->NodesRelated[k]==m_ele->nodes_index[l])
             {
                 exist = true;
                 break;
             }
          }
          if(!exist) 
             m_nod->NodesRelated.push_back(m_ele->nodes_index[l]);      
      }
    }
  }
  for(i=0;i<(long)node_vector.size();i++){
    m_nod = node_vector[i];
    j = (int)m_nod->NodesRelated.size();
    for(k=0; k<j; k++)
    {            
       for(l=k; l<j; l++)
       {
          if(m_nod->NodesRelated[l]<m_nod->NodesRelated[k])
          {
             n = m_nod->NodesRelated[k];
             m_nod->NodesRelated[k] = m_nod->NodesRelated[l];
             m_nod->NodesRelated[l] = n;
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
   long i, j, e, ni;
   Elem* thisElem0=NULL;
   Node * node = NULL;
   bool done = false;
   // set neighbors of node
   for(e=0; e<(long)node_vector.size(); e++)
     node_vector[e]->ElementsRelated.clear();
   for(e=0; e<(long)elem_vector.size(); e++)
   {
      thisElem0 = elem_vector[e];   
      if(!thisElem0->getStatus()) continue;      // Not marked for use
      for(i=0; i<thisElem0->getNodesNumber(quadratic); i++)
      {
          done = false;
          ni = thisElem0->getNodeIndex(i);
          node = node_vector[ni];
          for(j=0; j<(int)node->ElementsRelated.size(); j++)
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
   int nnodes0, nedges0, nedges;
   long e, ei, ee,  e_size,  e_size_l;
   bool done;
   double x_sum,y_sum,z_sum;

   int edgeIndex_loc0[2];
   int edgeIndex_loc[2];
   int faceIndex_loc0[10];
   int faceIndex_loc[10];
   vec<Node*> e_nodes0(20);
   vec<long> node_index_glb(20);
   vec<long> node_index_glb0(20);
   vec<int> Edge_Orientation(15);
   vec<Edge*> Edges(15);
   vec<Edge*> Edges0(15);
   vec<Elem*> Neighbors(15);
   vec<Elem*> Neighbors0(15);

   vec<Node*> e_edgeNodes0(3);
   vec<Node*> e_edgeNodes(3);
   Elem* thisElem0=NULL;
   Elem* thisElem=NULL;

   clock_t start, finish;
   start = clock();

   //Elem->nodes not initialized

   e_size = (long)elem_vector.size();  
   NodesNumber_Linear= (long)node_vector.size();

   Edge_Orientation = 1;
   //----------------------------------------------------------------------
   // set neighbors of node
   ConnectedElements2Node();
   //----------------------------------------------------------------------

  //----------------------------------------------------------------------
   // Compute neighbors and edges
   for(e=0; e<e_size; e++)
   {
       thisElem0 = elem_vector[e];   
       nnodes0 = thisElem0->nnodes; // Number of nodes for linear element
       thisElem0->getNodeIndeces(node_index_glb0);
       thisElem0->getNeighbors(Neighbors0);
       for(i=0; i<nnodes0; i++) // Nodes
         e_nodes0[i] = node_vector[node_index_glb0[i]];  
       m0 = thisElem0->getFacesNumber();
	   // neighbors
       for(i=0; i<m0; i++) // Faces
       {
          if(Neighbors0[i])
               continue;
          n0 = thisElem0->getElementFaceNodes(i, faceIndex_loc0);
          done = false;  
          for(k=0;k<n0;k++)
          {    
             e_size_l = (long)e_nodes0[faceIndex_loc0[k]]->ElementsRelated.size();         
             for(ei=0; ei<e_size_l; ei++)
             {
                ee = e_nodes0[faceIndex_loc0[k]]->ElementsRelated[ei];   
                if(ee==e) continue;
                thisElem = elem_vector[ee];   
                thisElem->getNodeIndeces(node_index_glb);
                thisElem->getNeighbors(Neighbors);
                m = thisElem->getFacesNumber();

                for(ii=0; ii<m; ii++) // Faces
                {
                   n = thisElem->getElementFaceNodes(ii, faceIndex_loc);
                   if(n0!=n) continue;
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
                       thisElem->setNeighbor(ii, thisElem0);
                       done = true;
                       break;                 
                   }
                }
                if(done) break;
             }
             if(done) break;
          }
       }
       thisElem0->setNeighbors(Neighbors0);						
       
       // --------------------------------
       // Edges
       nedges0 = thisElem0->getEdgesNumber();
       thisElem0->getEdges(Edges0);
       for(i=0; i<nedges0; i++)
       { 
          thisElem0->getLocalIndices_EdgeNodes(i, edgeIndex_loc0);    
          // Check neighbors 
          done = false; 
          for(k=0;k<2;k++)
          {    
             e_size_l = (long)e_nodes0[edgeIndex_loc0[k]]->ElementsRelated.size();         
             for(ei=0; ei<e_size_l; ei++)
             {
                ee = e_nodes0[edgeIndex_loc0[k]]->ElementsRelated[ei];   
                if(ee==e) continue;
                thisElem = elem_vector[ee];                   
                thisElem->getNodeIndeces(node_index_glb);
                nedges = thisElem->getEdgesNumber();
                thisElem->getEdges(Edges);
                // Edges of neighbors
                for(ii=0; ii<nedges; ii++)
                { 
                    thisElem->getLocalIndices_EdgeNodes(ii, edgeIndex_loc);
                    if((  node_index_glb0[edgeIndex_loc0[0]]==node_index_glb[edgeIndex_loc[0]]
                        &&node_index_glb0[edgeIndex_loc0[1]]==node_index_glb[edgeIndex_loc[1]])
				                 ||(  node_index_glb0[edgeIndex_loc0[0]]==node_index_glb[edgeIndex_loc[1]]
                        &&node_index_glb0[edgeIndex_loc0[1]]==node_index_glb[edgeIndex_loc[0]]) )
                     {
                         if(Edges[ii])
                         {
                            Edges0[i] = Edges[ii]; 
                            Edges[ii]->getNodes(e_edgeNodes); 
                            if(  node_index_glb0[edgeIndex_loc0[0]]==e_edgeNodes[1]->getIndex()
                             && node_index_glb0[edgeIndex_loc0[1]]==e_edgeNodes[0]->getIndex())
			                             Edge_Orientation[i] = -1; 
                            done = true;
                            break;
                        }
		            }                   
                } //  for(ii=0; ii<nedges; ii++)
                if(done) break;
             } // for(ei=0; ei<e_size_l; ei++)
             if(done) break;
          }//for(k=0;k<2;k++)
          if(!done) // new edges and new node
          {
              Edges0[i] = new Edge((long)edge_vector.size()); 
              Edges0[i]->setOrder(false); 
              e_edgeNodes0[0] = e_nodes0[edgeIndex_loc0[0]];
              e_edgeNodes0[1] = e_nodes0[edgeIndex_loc0[1]];
              e_edgeNodes0[2] = NULL;
             Edges0[i]->setNodes(e_edgeNodes0); 
             edge_vector.push_back(Edges0[i]);		               
          } // new edges
   	  } //  for(i=0; i<nedges0; i++)
      //
      // set edges and nodes
      thisElem0->setOrder(false);
      thisElem0->setEdges_Orientation(Edge_Orientation); 
      thisElem0->setEdges(Edges0); 
      // Resize is true
      thisElem0->setNodes(e_nodes0, true);						
   }// Over elements
   // set faces on surfaces and others
   msh_no_line=0;  // Should be members of mesh
   msh_no_quad=0;
   msh_no_hexs=0;
   msh_no_tris=0;
   msh_no_tets=0;
   msh_no_pris=0;
   for(e=0; e<e_size; e++)
   {
       thisElem0 = elem_vector[e];   
	   switch(thisElem0->getElementType())
	   {
	      case 1: msh_no_line++; break;
	      case 2: msh_no_quad++; break;
	      case 3: msh_no_hexs++; break;
	      case 4: msh_no_tris++; break;
	      case 5: msh_no_tets++; break;
	      case 6: msh_no_pris++; break;
	   }
       // Compute volume meanwhile
	   //thisElem0->ComputeVolume();
       
	   if(thisElem0->getElementType()==1) continue; // line element  
       thisElem0->getNodeIndeces(node_index_glb0);
       thisElem0->getNeighbors(Neighbors0);
       m0 = thisElem0->getFacesNumber();

       // Check face on surface
       for(i=0; i<m0; i++) // Faces
       {		  
          if(Neighbors0[i])
             continue;
          Elem* newFace = new Elem((long)face_vector.size(), thisElem0, i);
//          thisElem0->boundary_type='B';
		  thisElem0->no_faces_on_surface++;
          face_vector.push_back(newFace);
          Neighbors0[i] = newFace;        
       }
       thisElem0->setNeighbors(Neighbors0);		

   }
   NodesNumber_Quadratic= (long)node_vector.size();
   if((msh_no_hexs+msh_no_tets+msh_no_pris)>0) max_ele_dim=3;
   else if((msh_no_quad+msh_no_tris)>0) max_ele_dim=2;
   else max_ele_dim=1;
   //----------------------------------------------------------------------
   // Node information 
   // 1. Default node index <---> eqs index relationship
   // 2. Coordiate system flag
   x_sum=0.0;
   y_sum=0.0;
   z_sum=0.0;
   for(e=0; e<(long)node_vector.size(); e++)
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
   e_nodes0.resize(0);
   node_index_glb.resize(0);
   node_index_glb0.resize(0);
   Edge_Orientation.resize(0);
   Edges.resize(0);
   Edges0.resize(0);
   Neighbors.resize(0);
   Neighbors0.resize(0);
   e_edgeNodes0.resize(0);
   e_edgeNodes.resize(0);

   finish = clock();
   cout<<"\n\tCPU time elapsed in constructing topology of grids: "
      <<(double)(finish - start) / CLOCKS_PER_SEC<<"s"<<endl;

}

/**************************************************************************
Programing:
07/2007 WW Implementation
**************************************************************************/
void Mesh::GenerateHighOrderNodes()
{
   int i, k, ii;
   int nnodes0, nedges0, nedges;
   long e, ei, ee,  e_size,  e_size_l;
   int edgeIndex_loc0[2];
   bool done;
   double x0,y0,z0;

   clock_t start, finish;
   start = clock();

   //
   Node *aNode=NULL;
   vec<Node*> e_nodes0(20);
   Elem *thisElem0=NULL;
   Elem *thisElem=NULL;
   Edge *thisEdge0=NULL;
   Edge *thisEdge=NULL;
   //----------------------------------------------------------------------
   // Loop over elements
   e_size = (long)elem_vector.size();  
   for(e=0; e<e_size; e++)
   {
      thisElem0 = elem_vector[e];   
      nnodes0 = thisElem0->nnodes; // Number of nodes for linear element
//      thisElem0->GetNodeIndeces(node_index_glb0);
      for(i=0; i<nnodes0; i++) // Nodes
        e_nodes0[i] = thisElem0->getNode(i);  
      // --------------------------------
      // Edges
      nedges0 = thisElem0->getEdgesNumber();
      // Check if there is any neighbor that has new middle points
      for(i=0; i<nedges0; i++)
      { 
          thisEdge0 = thisElem0->getEdge(i);    
          thisElem0->getLocalIndices_EdgeNodes(i, edgeIndex_loc0);    
          // Check neighbors 
          done = false; 
          for(k=0;k<2;k++)
          {    
             e_size_l = (long)e_nodes0[edgeIndex_loc0[k]]->ElementsRelated.size();         
             for(ei=0; ei<e_size_l; ei++)
             {
                ee = e_nodes0[edgeIndex_loc0[k]]->ElementsRelated[ei];   
                if(ee==e) continue;
                thisElem = elem_vector[ee];                   
                nedges = thisElem->getEdgesNumber();
                // Edges of neighbors
                for(ii=0; ii<nedges; ii++)
                { 
                   thisEdge = thisElem->getEdge(ii);
                   if(*thisEdge0==*thisEdge)
                   {
                      aNode = thisEdge->getNode(2); 
                      if(aNode) // The middle point exist
                      {
                         e_nodes0[nnodes0] = aNode;
                         nnodes0++;
                         done = true; 
                         break;                   
                      }  
                   }
                } //  for(ii=0; ii<nedges; ii++)
                if(done) break;
             } // for(ei=0; ei<e_size_l; ei++)
             if(done) break;
          }//for(k=0;k<2;k++)
         if(!done)
         {
            aNode = new Node((long)node_vector.size());
            aNode->setX(0.5*(thisEdge0->getNode(0)->X()+thisEdge0->getNode(1)->X()));                
            aNode->setY(0.5*(thisEdge0->getNode(0)->Y()+thisEdge0->getNode(1)->Y()));                
            aNode->setZ(0.5*(thisEdge0->getNode(0)->Z()+thisEdge0->getNode(1)->Z()));    
            e_nodes0[nnodes0] = aNode;
            thisEdge0->setNode(2, aNode);
            nnodes0++; 
            node_vector.push_back(aNode);
         }
   	  } //  for(i=0; i<nedges0; i++)
      // No neighors or no neighbor has new middle point     
      //
      if(thisElem0->getElementType()==2) // Quadrilateral
      {
         x0=y0=z0=0.0;
         aNode = new Node((long)node_vector.size());
         e_nodes0[nnodes0] = aNode;
         nnodes0 = thisElem0->nnodes;
         for(i=0; i<nnodes0; i++) // Nodes
         {
            x0 += e_nodes0[i]->X();	
            y0 += e_nodes0[i]->Y();	
            z0 += e_nodes0[i]->Z();	
         }         
         x0 /= (double)nnodes0;
         y0 /= (double)nnodes0;
         z0 /= (double)nnodes0;
         aNode->setX(x0);
         aNode->setY(y0);
         aNode->setZ(z0);
         node_vector.push_back(aNode);         
      }     
      // Set edges and nodes
      thisElem0->setOrder(true);
      // Resize is true
      thisElem0->setNodes(e_nodes0, true);		
   }// Over elements
   //
   NodesNumber_Quadratic= (long)node_vector.size();
   for(e=0; e<e_size; e++)
   {
      thisElem0 = elem_vector[e];   
      for(i=thisElem0->nnodes; i<thisElem0->nnodesHQ; i++)
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
   e_nodes0.resize(0);


   finish = clock();
   cout<<"\n\tCPU time elapsed in generating high oder elements: "
      <<(double)(finish - start) / CLOCKS_PER_SEC<<"s"<<endl;

}


void Mesh::ConstructSubDomain_by_Elements(const string fname, const int num_parts)
{
   string str;
   string stro;
   char str_buf[3];
   int dom;
   int max_dom;
   int k,kk;
   long i,j;
   int ntags = 3;
   string deli = " ";
   //
   sprintf(str_buf, "%d",num_parts);


   string s_nparts = str_buf;
   str = fname + ".mesh.epart." + str_buf;
   stro = fname + "." + str_buf +"ddc";

   //namef = ".mesh.epart."; //+str_buf;
   ifstream part_in;
   fstream part_out;
   part_out.open(stro, ios::out );
   // Output for gmsh
 
   stro = fname + "_gmsh.msh";
   fstream gmsh_out;
   gmsh_out.open(stro, ios::out );
   //gmsh_out<<"$NOD"<<endl;
   gmsh_out<<"$MeshFormat\n2 0 8\n$EndMeshFormat\n$Nodes"<<endl;
   gmsh_out<<node_vector.size()<<endl;
   Node *node;
   for(i=0; i<(long)node_vector.size(); i++)
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
   gmsh_out<<(long)elem_vector.size()<<endl; 
   //
   part_in.open(str);
   if(!part_in.is_open())
   {
       cerr<<("Error: cannot open .epart file . It may not exist !");
       abort();
   } 

   max_dom=0;
   int et;
   Elem  *elem;
   for(i=0; i<(long)elem_vector.size(); i++)
   {
      part_in>>dom>>ws;
      elem = elem_vector[i];
	  elem->setDomainIndex(dom);
//      elem_vector[i]->AllocateLocalIndexVector(); 
      if(dom>max_dom) max_dom = dom;
      // GMSH output
      switch(elem->getElementType())
      {
         case 1: et = 1; break;
         case 4: et = 2; break;
         case 2: et = 3; break;
         case 5: et = 4; break;
         case 3: et = 5; break;
         case 6: et = 6; break;
      }
      //gmsh_out<<i+1<<" "<<et<<" "<<dom+1<<" "<<dom+1<<" ";
      gmsh_out<<i+1<<deli<<et<<deli<<ntags<<deli<<elem->getPatchIndex()+1<<deli<<elem->getPatchIndex()+1<<deli<<dom+1<<deli;


      //gmsh_out<<elem->getNodesNumber()<<"  ";
      for(k=0; k<elem->getNodesNumber(); k++)
        gmsh_out<<elem->getNodeIndex(k)+1<<" ";
      gmsh_out<<endl;
   }
   max_dom++;
   //gmsh_out<<"$ENDELM"<<endl; 
   gmsh_out<<"$EndElements"<<endl;
   gmsh_out.close();
   part_in.close();
   //
 
   //Output ddc file
   // long *nod_dom = new long[max_dom];
   long *ele_dom = new long[max_dom];
   for(k=0; k<max_dom; k++)
   {
      ele_dom[k]=0;
      //nod_dom[k]=0;
      for(j=0; j<(long)elem_vector.size(); j++)
      {
	     if(elem_vector[j]->getDomainIndex()==k)
             ele_dom[k] += 1;
      } 
      /*
      for(j=0; j<(long)node_vector.size(); j++)
      {
	     if(node_dom[j]==k)
             nod_dom[k] += 1;
      } 
      */
   }

   Elem *ele=0;
   bool done = false;
   long n_index=0;
   vector<int> nodes_dom;
   vector<int> eles_dom;

   //TEST
   Node *nod = 0;
   //


   //
   for(k=0; k<max_dom; k++)
   {      	   
	  part_out<<"#DOMAIN "<<k<<endl;   
	  part_out<<"$ELEMENTS "<<ele_dom[k]<<endl; 
	  nodes_dom.clear(); 
	  eles_dom.clear();
      for(j=0; j<(long)elem_vector.size(); j++)
	  {
         ele = elem_vector[j];
         for(kk=0; kk<ele->getNodesNumber(); kk++)
		 {
            ele->setLocalNodeIndex(kk, -1);
            ele->AllocateLocalIndexVector();
            ele->setDomNodeIndex(kk, -1);
		 }
      }
      for(j=0; j<(long)elem_vector.size(); j++)
      {
         ele = elem_vector[j]; 
		 //ele->AllocateLocalIndexVector();
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
				   ele->setDomNodeIndex(kk, (long)nodes_dom.size()); //For test output
                   ele->setLocalNodeIndex(kk, (long)nodes_dom.size());
				   nodes_dom.push_back(ele->getNodeIndex(kk));
			   } 
			}
            part_out<<ele->getIndex()<<endl;
			eles_dom.push_back(ele->getIndex()); //TEST OUT
		 }
      }        
	  part_out<<"$NODES_INNER "<<(long)nodes_dom.size()<<endl;
      for(j=0; j<(long)nodes_dom.size(); j++)
          part_out<<nodes_dom[j]<<endl;
      /*
	  part_out<<"$NODES_INNER "<<nod_dom[k]<<endl;
      for(j=0; j<(long)node_vector.size(); j++)
      {
	     if(node_dom[j]==k)
            part_out<<node_vector[j]->GetIndex()<<endl;
      } 
	  */
      //TEST OUT

      sprintf(str_buf, "%d",k);

      string name_f = fname+"_"+str_buf+"_of_"+s_nparts+"subdomains.msh";
      fstream test_out;
      test_out.open(name_f.c_str(), ios::out );

      //GMSH test_out<<"$NOD"<<endl;
 	  //GMSH test_out<<(long)nodes_dom.size()<<endl;
      test_out<<"#0#0#0#1#0.0#0#################################################################"<<endl;
      test_out<<"0 "<<(long)nodes_dom.size()<<" "<<(long)eles_dom.size()<<endl;
      for(j=0; j<(long)nodes_dom.size(); j++)
	  {
          nod = node_vector[nodes_dom[j]];
        //GMSH  test_out<<j+1<<"  "
          test_out<<j<<"  "
                  << nod->X()<<"  "<< nod->Y()<<"  "<< nod->Z() <<endl;
	  }
      //GMSH test_out<<"$ENDNOD"<<endl;
      //GMSH test_out<<"$ELE"<<endl;
	  //GMSG test_out<<(long)eles_dom.size()<<endl;
      for(j=0; j<(long)eles_dom.size(); j++)
	  {
          ele = elem_vector[eles_dom[j]]; 
          /* //GMSH
          test_out<<j+1<<"  ";
          int e_t=1;       
          switch(ele->getElementType())
          {
            case 1: e_t = 1; break;
            case 4: e_t = 2; break;
            case 2: e_t = 3; break;
            case 5: e_t = 4; break;
            case 3: e_t = 5; break;
            case 6: e_t = 6; break;
          }
          test_out<<e_t<<" 1 23 "<<ele->getNodesNumber()<<" ";

          for(kk=0; kk<ele->getNodesNumber(); kk++)
             test_out<< ele->GetDomNodeIndex(kk)+1<<" ";
          */
		  test_out<<j<<"  "<<ele->getPatchIndex()<<" "<<ele->getName()<<" ";
          for(kk=0; kk<ele->getNodesNumber(); kk++)
             test_out<< ele->getDomNodeIndex(kk)<<" ";
          test_out<<endl;
	  }
     //GMSH test_out<<"$ENDELE"<<endl;
      test_out.clear();
	  test_out.close();     
   }
   part_out<<"#STOP "<<endl;  
   part_out.clear();
   part_out.close();

   //
   delete ele_dom;
   //delete nod_dom;
}

/*!
\brief void Mesh::ConstructSubDomain_by_Nodes

Partition a mesh ny nodes

02.2012 WW
*/

void Mesh::ConstructSubDomain_by_Nodes(const string fname, const int num_parts, const bool is_quad)
{

   string f_iparts;
   string o_part_msh;
   char str_buf[3];
   int dom;
   int k,kk;
   long i,j;
   int ntags = 3;
   string deli = " ";

   Node *a_node = NULL;
   Elem *a_elem = NULL;
   //
   sprintf(str_buf, "%d",num_parts);


   string s_nparts = str_buf;
   f_iparts = fname + ".mesh.npart." + str_buf;
   //o_part_msh = fname + "." + str_buf +"mesh";

   vector<long> node_dom;
   ifstream npart_in(f_iparts.c_str());
   if(!npart_in.is_open())
   {
       cerr<<("Error: cannot open .npart file . It may not exist !");
       exit(1);
   } 


   list <vector<long>> sub_dom_nodes(num_parts); 
   list<vector<long>>::iterator itr_subd_list0 = sub_dom_nodes.begin();
   list<vector<long>>::iterator itr_subd_list;

   for(i=0; i<(long)node_vector.size(); i++)
   {
      npart_in>>dom>>ws;

	  itr_subd_list = itr_subd_list0; 
	  std::advance(itr_subd_list, dom);
	  	
      (*itr_subd_list).push_back(i);  
   }
   npart_in.close();

   list<vector<long>>::iterator iter;

   iter=itr_subd_list0;
   
   int counter = 0;
   for(iter=itr_subd_list0; iter!=sub_dom_nodes.end(); iter++)
   {
       vector<long> sbd_nodes = *iter;
       long size_sbd_nodes = (long)sbd_nodes.size();

	   vector<Elem*> in_subdom_elements;
	   vector<Elem*> ghost_subdom_elements;

	   // Un-making all nodes and elements of the whole mesh 
       for(j=0; j<(long)node_vector.size(); j++)
          node_vector[j]->Marking(false);
       for(j=0; j<(long)elem_vector.size(); j++)
          elem_vector[j]->Marking(false);
	   // Only select nodes in this subdomain
	   for(j=0; j<size_sbd_nodes; j++)
	   {
		   node_vector[sbd_nodes[j]]->Marking(true);            
	   }


       /// Find the elements in this subdomain.
	   for(j=0; j<size_sbd_nodes; j++)
	   {
           a_node = node_vector[sbd_nodes[j]]; 
           
           // Search the elements connected to this nodes
		   for(k=0; k<a_node->ElementsRelated.size(); k++)
		   {
               a_elem = elem_vector[a_node->ElementsRelated[k]];

			   // If checked
			   if(a_elem->getStatus())
                  continue;

               
			   vector<int> g_nodes; // Nodes in ghost elements
			   for(kk=0; kk<a_elem->getNodesNumber(); kk++)
			   {                
				   if(!(a_elem->getNode(kk)->getStatus()))   
                     g_nodes.push_back(kk);
			   }
                  
               // All nodes of this element are inside this subdomain
			   if(g_nodes.size() == 0)
			   {
				   in_subdom_elements.push_back(a_elem);
			   }
			   else
			   {
				   ghost_subdom_elements.push_back(a_elem);
				   a_elem->ghost_nodes.resize((int)g_nodes.size());
                 
                   for(kk=0; kk<g_nodes.size(); kk++)
					   a_elem->ghost_nodes[kk] = g_nodes[kk];
			   }


			   a_elem->Marking(true);

		   }
           
	   }

       // Make output of this subdomain for simulation
	   sprintf(str_buf, "%d", counter);
       counter++;

       //string name_f = fname+"_"+str_buf+"_of_"+s_nparts+"_subdomains.msh";    
       string name_f = fname+"_"+str_buf+".msh";    
	   fstream os_subd(name_f.c_str(), ios::out|ios::trunc );

	   //os_subd<<"#FEM_MSH\n   $PCS_TYPE\n    NULL"<<endl;
       //os_subd<<" $NODES\n"<<size_sbd_nodes<<endl;
       //os_subd<<" $NODES\n"<<size_sbd_nodes<<endl;
       os_subd<<"Nodes   Elements  Ghost elements"<<endl; 
       os_subd<<size_sbd_nodes<<deli<<in_subdom_elements.size()
		      <<deli<<ghost_subdom_elements.size()<<endl;

	   //os_subd<<"Nodes"<<endl;
       for(j=0; j<size_sbd_nodes; j++)
		   node_vector[sbd_nodes[j]]->Write(os_subd);

	   //os_subd<<"Elements"<<endl;
	   for(j=0; j<in_subdom_elements.size(); j++)
		   in_subdom_elements[j]->WriteGSmsh(os_subd);
	  
	   //os_subd<<"Ghost elements"<<endl;
	   for(j=0; j<ghost_subdom_elements.size(); j++)
	   {
           a_elem = ghost_subdom_elements[j];
		   a_elem->WriteGSmsh(os_subd);
		   for(kk=0; kk<a_elem->ghost_nodes.Size(); kk++)
		   {
              os_subd<<a_elem->ghost_nodes[kk]<<deli;  
		   }
		   os_subd<<endl;
	   }
	   os_subd.close();

	   //-----------------------------------------------------------
	   /// VTK output
       // Elements in this subdomain  
       f_iparts = fname+"_"+str_buf+"_of_"+s_nparts+"_subdomains.vtk";
       //f_iparts = fname+"_"+str_buf+".vtk";
       ofstream os(f_iparts.c_str(), ios::out|ios::trunc);

	   WriteVTK_Nodes(os);
	   WriteVTK_Elements_of_Subdomain(os, in_subdom_elements, counter);
	   os.close();

	   // Ghost elements in this subdomain  
       f_iparts = fname+"_"+str_buf+"_ghost_of_"+s_nparts+"_subdomains.vtk";
       //f_iparts = fname+"_"+str_buf+"ghost.vtk";
	   os.open(f_iparts.c_str(), ios::out|ios::trunc);
       WriteVTK_Nodes(os);
	   WriteVTK_Elements_of_Subdomain(os, ghost_subdom_elements, 0);
	   os.close();
	   //-----------------------------------------------------------

   }
 
 
}

// 02.2012. WW
void  Mesh::WriteVTK_Nodes(std::ostream& os)
{
  long i;
  Node *a_node = NULL;

  os<<"# vtk DataFile Version 4.0\nGrid Partition by WW \nASCII\n"<<endl;
  os<<"DATASET UNSTRUCTURED_GRID"<<endl;
  os<<"POINTS "<<node_vector.size()<<" double"<<endl;
  setw(14);
  os.precision(14);
  for(i=0; i<(long)node_vector.size(); i++)
  {
     a_node = node_vector[i];
     os<<a_node->X()<<" "<<a_node->Y()<<" "<<a_node->Z()<<endl; 
  }

}
// 02.2012. WW
void  Mesh::WriteVTK_Elements_of_Subdomain(std::ostream& os, std::vector<Elem*>& ele_vec, const int sbd_index) 
{
   long i;
   int k;
   //-----------------------------------------------------------
   //  VTK output
   // Elements in this subdomain
   long ne0 = (long)ele_vec.size();
   long size = ne0;

   string deli = " ";

   Elem *a_elem = NULL;

   for(i=0; i<ne0; i++)
      size += ele_vec[i]->getNodesNumber(false);
   os<<"\nCELLS "<<ne0<<deli<<size<<endl;

  
   // CELLs
   for(i=0;i<ne0;i++)
   {
      a_elem = ele_vec[i];
      os<<a_elem->getNodesNumber(false)<<deli;
      for(k=0; k<a_elem->getNodesNumber(false); k++)
         os << a_elem->nodes_index[k] << deli;
      os << endl;
   }
   os << endl; 

   // CELL types
   os << "CELL_TYPES " << ne0 << endl; 
   for(i=0;i<ne0;i++)
   {
      a_elem = ele_vec[i];
      a_elem->WriteVTK_Type(os);
   }
   os << endl; 
	  
   // Partition
   os<<"CELL_DATA "<<ne0<<endl;
   os<<"SCALARS Partition int 1\nLOOKUP_TABLE default"<<endl;
   for(i=0; i<ne0; i++)
     os<<sbd_index<<endl;

}

void Mesh::Write2METIS(ostream& os)
{
	 os<<(long)elem_vector.size()<<" ";
     int e_type =0;
	 switch(elem_vector[0]->getElementType())    
	 {
        case 1: 
			cout<<"Not for 1D element"<<endl;
			exit(1); 
        case 2: 
			e_type =4; 
			break; 
        case 3: 
			e_type =3; 
			break;  
        case 4: 
			e_type =1; 
			break;  
        case 5: 
			e_type =2; 
			break; 
        case 6: cout<<"Not for prismal element"<<endl; abort(); 
	 } 
     os<<e_type<<endl;
     for(long i=0; i<(long)elem_vector.size(); i++)
        elem_vector[i]->Write_index(os);
}



// Read grid for test purpose
void Mesh::ReadGrid(istream& is)
{
   long i, ne, nn, counter;
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
       abort();
   }
   
   // Read Nodes
   counter = 0;
   for(i=0; i<nn; i++)
   {
      is>>ibuff>>x>>y>>z>>ws;
      Node* newNode = new Node(ibuff,x,y,z);
      node_vector.push_back(newNode);            
      counter++;
   }
   if(counter!=nn)
   {
       cout<<"Error: number nodes do not match"<<endl;
       abort();
   }
   NodesNumber_Linear = nn;
   NodesNumber_Quadratic = nn;
      
   // Read Elements
   counter = 0;
   for(i=0; i<ne; i++)
   { 
      Elem* newElem = new Elem(i);
      newElem->Read(is);
	  newElem->Marking(true);
      elem_vector.push_back(newElem);           
      counter++;
   }     
   if(counter!=ne)
   {
       cout<<"Error: number elements do not match"<<endl;
       abort();
   }

//   position = is.tellg();
}


void Mesh::ReadGridGeoSys(istream& is)
{
  string sub_line;
  string line_string;
  bool new_keyword = false;
  string hash("#");
  string sub_string,sub_string1;
  long i, ibuff;
  long no_elements;
  long no_nodes;
  double x,y,z;
  Node* newNode = NULL;
  Elem* newElem = NULL;
  //========================================================================
  // Keyword loop
  while (!new_keyword) {
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
    if(line_string.find("$NODES")!=string::npos) { // subkeyword found
      is  >> no_nodes>>ws;
      for(i=0;i<no_nodes;i++){
         is>>ibuff>>x>>y>>z>>ws;
         newNode = new Node(ibuff,x,y,z);
         node_vector.push_back(newNode);            
      }
      continue;
    }
    //....................................................................
    if(line_string.find("$ELEMENTS")!=string::npos) { // subkeyword found
      is >> no_elements>>ws;
      for(i=0;i<no_elements;i++){
         newElem = new Elem(i);
         newElem->Read(is, 0);
         newElem->Marking(true);
		 elem_vector.push_back(newElem);
      }
      continue;
    }
  }
  //========================================================================
}

}//end namespace


