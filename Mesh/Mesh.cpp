#include "Mesh.h"

#include <iomanip>
#include <limits>

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


Mesh::~Mesh()
{
			long i;
			for(i=0; i<(long)NodesVector.size(); i++)
     delete NodesVector[i];
			NodesVector.clear();
			for(i=0; i<(long)EdgeVector.size(); i++)
     delete EdgeVector[i];
			EdgeVector.clear();
			for(i=0; i<(long)SurfaceFaces.size(); i++)
     delete SurfaceFaces[i];
			SurfaceFaces.clear();

			for(i=0; i<(long)ElementsVector.size(); i++)
     delete ElementsVector[i];
			ElementsVector.clear();
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
      NodesVector.push_back(newNode);            
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
      ElementsVector.push_back(newElem);           
      counter++;
   }     
   if(counter!=ne)
   {
       cout<<"Error: number elements do not match"<<endl;
       abort();
   }

//   position = is.tellg();
}

// Construct grid
// 
void Mesh::ConstructGrid( const bool quadratic)
{
   int counter;	 
   int i, j, k, ii, jj, m0, m, n0, n;
			int nnodes0, nedges0, nedges;
   long e, ei, ee,  e_size,  e_size_l;
   bool done;
   double x0,y0,z0;

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

   //Elem->nodes not initialized

   e_size = (long)ElementsVector.size();  
   Edge_Orientation = 1;

   // Set neighbors
   for(e=0; e<e_size; e++)
   {
       thisElem0 = ElementsVector[e];   
       thisElem0->getNodeIndeces(node_index_glb0);
       for(i=0; i<thisElem0->getNodesNumber(); i++)
       {
          done = false;
          for(j=0; j<(int)NodesVector[node_index_glb0[i]]
					         ->ElementsBelonged.size(); j++)
          {
            if(e==NodesVector[node_index_glb0[i]]
                            ->ElementsBelonged[j])
              done = true;
              break;
          }
          if(!done)  
          NodesVector[node_index_glb0[i]]->ElementsBelonged.push_back(e);
      }
   }

   // Compute neighbors and edges
   for(e=0; e<e_size; e++)
   {
       thisElem0 = ElementsVector[e];   
       nnodes0 = thisElem0->getNodesNumber(); // Number of nodes for linear element
       thisElem0->getNodeIndeces(node_index_glb0);
       thisElem0->getNeighbors(Neighbors0);
       for(i=0; i<nnodes0; i++) // Nodes
         e_nodes0[i] = NodesVector[node_index_glb0[i]];  
       m0 = thisElem0->getFacesNumber();
	      // neighbors
       for(i=0; i<m0; i++) // Faces
       {
          if(Neighbors0[i])
               continue;
          n0 = thisElem0->GetElementFaceNodes(i, faceIndex_loc0);
          done = false;  
          for(k=0;k<n0;k++)
          {    
             e_size_l = (long)e_nodes0[faceIndex_loc0[k]]->ElementsBelonged.size();         
             for(ei=0; ei<e_size_l; ei++)
             {
                ee = e_nodes0[faceIndex_loc0[k]]->ElementsBelonged[ei];   
                if(ee==e) continue;
                thisElem = ElementsVector[ee];   
                thisElem->getNodeIndeces(node_index_glb);
                thisElem->getNeighbors(Neighbors);
                m = thisElem->getFacesNumber();

                for(ii=0; ii<m; ii++) // Faces
                {
                   n = thisElem->GetElementFaceNodes(ii, faceIndex_loc);
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
       // 
       // Edges
       nedges0 = thisElem0->getEdgesNumber();
       thisElem0->getEdges(Edges0);
       for(i=0; i<nedges0; i++)
       { 
          thisElem0->GetLocalIndices_EdgeNodes(i, edgeIndex_loc0);    
          // Check neighbors 
          done = false; 
          for(k=0;k<2;k++)
          {    
             e_size_l = (long)e_nodes0[edgeIndex_loc0[k]]->ElementsBelonged.size();         
             for(ei=0; ei<e_size_l; ei++)
             {
                ee = e_nodes0[edgeIndex_loc0[k]]->ElementsBelonged[ei];   
                if(ee==e) continue;
                thisElem = ElementsVector[ee];                   
                thisElem->getNodeIndeces(node_index_glb);
                nedges = thisElem->getEdgesNumber();
                thisElem->getEdges(Edges);
                // Edges of neighbors
                for(ii=0; ii<nedges; ii++)
                { 
                    thisElem->GetLocalIndices_EdgeNodes(ii, edgeIndex_loc);
                    if((  node_index_glb0[edgeIndex_loc0[0]]==node_index_glb[edgeIndex_loc[0]]
                        &&node_index_glb0[edgeIndex_loc0[1]]==node_index_glb[edgeIndex_loc[1]])
				                 ||(  node_index_glb0[edgeIndex_loc0[0]]==node_index_glb[edgeIndex_loc[1]]
                        &&node_index_glb0[edgeIndex_loc0[1]]==node_index_glb[edgeIndex_loc[0]]) )
                     {
                         if(Edges[ii])
                         {
                            Edges0[i] = Edges[ii]; 
                            Edges[ii]->getNodes(e_edgeNodes); 
                            if(  node_index_glb0[edgeIndex_loc0[0]]==e_edgeNodes[1]->GetIndex()
                             && node_index_glb0[edgeIndex_loc0[1]]==e_edgeNodes[0]->GetIndex())
			                             Edge_Orientation[i] = -1; 
                            if(quadratic)  // Get middle node
                            {
                               node_index_glb0[nnodes0] = e_edgeNodes[2]->GetIndex();
                               e_nodes0[nnodes0] = e_edgeNodes[2];
                               nnodes0++;
                            }
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
              Edges0[i] = new Edge((long)EdgeVector.size()); 
              Edges0[i]->SetOrder(quadratic); 
              e_edgeNodes0[0] = e_nodes0[edgeIndex_loc0[0]];
              e_edgeNodes0[1] = e_nodes0[edgeIndex_loc0[1]];
              if(quadratic)  // new node: middle point of edges
              {
                  e_edgeNodes0[2] = new Node((long)NodesVector.size());
                  e_edgeNodes0[2]->setX(0.5*(e_edgeNodes0[0]->X()+e_edgeNodes0[1]->X()));                
                  e_edgeNodes0[2]->setY(0.5*(e_edgeNodes0[0]->Y()+e_edgeNodes0[1]->Y()));                
                  e_edgeNodes0[2]->setZ(0.5*(e_edgeNodes0[0]->Z()+e_edgeNodes0[1]->Z()));    
                  NodesVector.push_back(e_edgeNodes0[2]);
                  node_index_glb0[nnodes0] = e_edgeNodes0[2]->GetIndex();
                  e_nodes0[nnodes0] = e_edgeNodes0[2];
                  nnodes0++;
              }
             Edges0[i]->setNodes(e_edgeNodes0); 
             EdgeVector.push_back(Edges0[i]);		               
          } // new edges
   	  } //  for(i=0; i<nedges0; i++)
      //
      if(quadratic&&thisElem0->getElementType()==2) // Quadrilateral
      {
         x0=y0=z0=0.0;
         Node* newNode = new Node((long)NodesVector.size());
         e_nodes0[nnodes0] = newNode;
         nnodes0 = thisElem0->getNodesNumber();
         for(i=0; i<nnodes0; i++) // Nodes
         {
            x0 += e_nodes0[i]->X();	
            y0 += e_nodes0[i]->Y();	
            z0 += e_nodes0[i]->Z();	
         }         
         x0 /= (double)nnodes0;
         y0 /= (double)nnodes0;
         z0 /= (double)nnodes0;
         newNode->setX(x0);
         newNode->setY(y0);
         newNode->setZ(z0);
         NodesVector.push_back(newNode);         
      }     
      // Set edges and nodes
      thisElem0->SetOrder(quadratic);
      thisElem0->setEdges_Orientation(Edge_Orientation); 
      thisElem0->setEdges(Edges0); 
      // Resize is true
      thisElem0->setNodes(e_nodes0, true);						
   }// Over elements

   // Set faces on surfaces
   for(e=0; e<e_size; e++)
   {
       thisElem0 = ElementsVector[e];   
       thisElem0->getNodeIndeces(node_index_glb0);
       thisElem0->getNeighbors(Neighbors0);
       m0 = thisElem0->getFacesNumber();

       // Check face on surface
       for(i=0; i<m0; i++) // Faces
       {		  
          if(Neighbors0[i])
             continue;
          Elem* newFace = new Elem((long)SurfaceFaces.size(), thisElem0, i);
          SurfaceFaces.push_back(newFace);
          Neighbors0[i] = newFace;        
       }
       thisElem0->setNeighbors(Neighbors0);						
   }
   NodesNumber_Quadratic= (long)NodesVector.size(); 
}

// Read grid for test purpose
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
         NodesVector.push_back(newNode);            
      }
      continue;
    }
    //....................................................................
    if(line_string.find("$ELEMENTS")!=string::npos) { // subkeyword found
      is >> no_elements>>ws;
      for(i=0;i<no_elements;i++){
         newElem = new Elem(i);
         newElem->Read(is, 0);
		 ElementsVector.push_back(newElem);
      }
      continue;
    }
  }
  //========================================================================
}


void Mesh::ConstructDomain(char *fname, const int num_parts)
{
   char str[1028];
   char stro[1028];
   char str_buf[3];
   int dom;
   int max_dom;
   int k,kk;
   long i,j;
   int ntags = 3;
   string deli = " ";
   //
   sprintf(str_buf, "%d",num_parts);

   /////////////////vector<int> node_dom;
   strcpy(str,fname);
   strcpy(stro,fname);
   strcat(str,".mesh.epart.");
   strcat(str,str_buf);
   strcat(stro,".");
   strcat(stro,str_buf);
   strcat(stro,"ddc");
   //namef = ".mesh.epart."; //+str_buf;
   ifstream part_in;
   fstream part_out;
   part_out.open(stro, ios::out );
   // Output for gmsh
   strcpy(stro,fname);
   strcat(stro,"_gmsh.msh");
   fstream gmsh_out;
   gmsh_out.open(stro, ios::out );
   //gmsh_out<<"$NOD"<<endl;
   gmsh_out<<"$MeshFormat\n2 0 8\n$EndMeshFormat\n$Nodes"<<endl;
   gmsh_out<<NodesVector.size()<<endl;
   Node *node;
   for(i=0; i<(long)NodesVector.size(); i++)
   {
     gmsh_out<<i+1<<" ";
     node = NodesVector[i];  
     gmsh_out<<node->X()<<" ";
     gmsh_out<<node->Y()<<" ";
     gmsh_out<<node->Z()<<endl;
   }
   //gmsh_out<<"$ENDNOD"<<endl; 
   //gmsh_out<<"$ELM"<<endl; 
   gmsh_out<<"$EndNodes\n$Elements"<<endl;
   gmsh_out<<(long)ElementsVector.size()<<endl; 
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
   for(i=0; i<(long)ElementsVector.size(); i++)
   {
      part_in>>dom>>ws;
      elem = ElementsVector[i];
	  elem->setDomainIndex(dom);
//      ElementsVector[i]->AllocateLocalIndexVector(); 
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
   /*
   strcpy(str,fname);
   strcat(str,".npart");
   ifstream npart_in;

   npart_in.open(str);
   if(!npart_in.is_open())
   {
       cerr<<("Error: cannot open .npart file . It may not exist !");
       abort();
   } 
   for(i=0; i<(long)NodesVector.size(); i++)
   {
      npart_in>>dom>>ws;
	  node_dom.push_back(dom);
   }
   npart_in.close();
   */
   //Output ddc file
   // long *nod_dom = new long[max_dom];
   long *ele_dom = new long[max_dom];
   for(k=0; k<max_dom; k++)
   {
      ele_dom[k]=0;
      //nod_dom[k]=0;
      for(j=0; j<(long)ElementsVector.size(); j++)
      {
	     if(ElementsVector[j]->GetDomainIndex()==k)
             ele_dom[k] += 1;
      } 
      /*
      for(j=0; j<(long)NodesVector.size(); j++)
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
      for(j=0; j<(long)ElementsVector.size(); j++)
	  {
         ele = ElementsVector[j];
         for(kk=0; kk<ele->getNodesNumber(); kk++)
		 {
            ele->SetLocalNodeIndex(kk, -1);
            ele->AllocateLocalIndexVector();
            ele->setDomNodeIndex(kk, -1);
		 }
      }
      for(j=0; j<(long)ElementsVector.size(); j++)
      {
         ele = ElementsVector[j]; 
		 //ele->AllocateLocalIndexVector();
	     if(ele->GetDomainIndex()==k)   
		 {
            for(kk=0; kk<ele->getNodesNumber(); kk++)
			{  
               done = false;
               n_index = ele->GetLocalNodeIndex(kk);
               if(n_index>-1)
               {
                  ele->setDomNodeIndex(kk, n_index);                    
                  done = true;
               }
 			   if(!done)
			   { 
				   ele->setDomNodeIndex(kk, (long)nodes_dom.size()); //For test output
                   ele->SetLocalNodeIndex(kk, (long)nodes_dom.size());
				   nodes_dom.push_back(ele->getNodeIndex(kk));
			   } 
			}
            part_out<<ele->GetIndex()<<endl;
			eles_dom.push_back(ele->GetIndex()); //TEST OUT
		 }
      }        
	  part_out<<"$NODES_INNER "<<(long)nodes_dom.size()<<endl;
      for(j=0; j<(long)nodes_dom.size(); j++)
          part_out<<nodes_dom[j]<<endl;
      /*
	  part_out<<"$NODES_INNER "<<nod_dom[k]<<endl;
      for(j=0; j<(long)NodesVector.size(); j++)
      {
	     if(node_dom[j]==k)
            part_out<<NodesVector[j]->GetIndex()<<endl;
      } 
	  */
      //TEST OUT
      itoa(k,stro, 10);
      strcpy(str,fname);
      string aa = str;
      string bb= stro;
      string name_f = aa+bb+".msh";
      fstream test_out;
      test_out.open(name_f.c_str(), ios::out );

      //GMSH test_out<<"$NOD"<<endl;
 	  //GMSH test_out<<(long)nodes_dom.size()<<endl;
      test_out<<"#0#0#0#1#0.0#0#################################################################"<<endl;
      test_out<<"0 "<<(long)nodes_dom.size()<<" "<<(long)eles_dom.size()<<endl;
      for(j=0; j<(long)nodes_dom.size(); j++)
	  {
          nod = NodesVector[nodes_dom[j]];
        //GMSH  test_out<<j+1<<"  "
          test_out<<j<<"  "
                  << nod->X()<<"  "<< nod->Y()<<"  "<< nod->Z() <<endl;
	  }
      //GMSH test_out<<"$ENDNOD"<<endl;
      //GMSH test_out<<"$ELE"<<endl;
	  //GMSG test_out<<(long)eles_dom.size()<<endl;
      for(j=0; j<(long)eles_dom.size(); j++)
	  {
          ele = ElementsVector[eles_dom[j]]; 
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
             test_out<< ele->GetDomNodeIndex(kk)<<" ";
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

void Mesh::Write2METIS(ostream& os)
{
	 os<<(long)ElementsVector.size()<<" ";
     int e_type =0;
	 switch(ElementsVector[0]->getElementType())    
	 {
        case 1: cout<<"Not for 1D element"<<endl; abort(); 
        case 2: e_type =4; break; 
        case 3: e_type =3; break;  
        case 4: e_type =1; break;  
        case 5: e_type =2; break; 
        case 6: cout<<"Not for prismal element"<<endl; abort(); 
	 } 
     os<<e_type<<endl;
     for(long i=0; i<(long)ElementsVector.size(); i++)
        ElementsVector[i]->Write_index(os);
}

}//end namespace


