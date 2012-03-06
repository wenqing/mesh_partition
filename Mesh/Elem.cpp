#include "Elem.h"

#include <iomanip>

#include "Node.h"
#include "Edge.h"

//------------------------------------------------------ 
//   Topology definition of geometrical element.
//    WW. 10.01.2005
//------------------------------------------------------ 
namespace Mesh_Group
{
  
using namespace std;
using namespace Math_Group;
//-----------------------------------------------------
//2. Mesh
//    WW. 06.2005
Elem::Elem(const int Index):Grain(Index)
{

   nnodes = 0;
   nnodesHQ = 0;
   ele_dim = 1;         // Dimension of element
   PatchIndex = 0;
   //
   Volume = 0.0;
   
   Owner = NULL;
}
//    WW. 06.2005
Elem::~Elem()
{
   nodes_index.resize(0);
   locnodes_index.resize(0);
   nodes.resize(0);
   edges.resize(0);
   neighbors.resize(0);
   ghost_nodes.resize(0);
}

//    WW. 06.2005
Elem::  Elem( const int Index,  Elem* onwer, const int Face):
  Grain(Index), Owner(onwer)
{
   int i, j, k, n, ne; 
   static int faceIndex_loc[10];
   static int edgeIndex_loc[10];
 //  Owner = onwer;
   n = Owner->getElementFaceNodes(Face, faceIndex_loc);
   switch(Owner->ele_Type)
   {
       case 1:  // 1-D bar element
           break;          
       case 2: // 2-D quadrilateral element
           nnodes = 2;
           nnodesHQ = 3;    
           ele_dim = 1;
           ele_Type = 1;
           nfaces = 2;
           nedges = 0;
           break;           
       case 3: // 3-D hexahedral element 
           nnodes = 4;
           nnodesHQ = 8;      
           ele_dim = 2;
           ele_Type = 2;
           nfaces = 4;
           nedges = 4;
           break;           
       case 4:  // 2-D triagular element 
           nnodes = 2;
           nnodesHQ = 3;    
           ele_dim = 1;
           ele_Type = 1;
           nfaces = 2;
           nedges = 0;
           break;           
       case 5:  // 3-D tetrahedral element 
           nnodes = 3;
           nnodesHQ = 6;
           ele_dim = 2;
           ele_Type = 4;
           nfaces = 3;
           nedges = 3;
           break;           
       case 6: 
          if(Face<2) 
          {
              nnodes = 3;
              nnodesHQ = 6;
              ele_dim = 2;
              ele_Type = 4;
              nfaces = 3;
              nedges = 3;
           }
           else  
           {
              nnodes = 4;
              nnodesHQ = 8;      
              ele_dim = 2;
              ele_Type = 2;
              nfaces = 4;
              nedges = 4;              
           }
           break; // 3-D prismatic element 
    }

    PatchIndex =  Owner->PatchIndex;
    quadratic = Owner->quadratic;
    nodes_index.resize(n);
    nodes.resize(n);

    boundayC='B';
    for(i=0; i<n; i++)
    {
       nodes_index[i] =
                  Owner->nodes_index[faceIndex_loc[i]];
       nodes[i] = Owner->nodes[faceIndex_loc[i]];
							nodes[i]->boundayC = 'B';
    }
    // Face edges
    ne = Owner->getEdgesNumber();
    edges.resize(nnodes);
    edges_orientation.resize(nnodes);
    edges_orientation = 1;
    for(i=0; i<nnodes; i++)
    {
        k = (i+1)%nnodes; 
        for(j=0; j<ne; j++)
        {
           Owner->getLocalIndices_EdgeNodes(j, edgeIndex_loc);    
           if( (faceIndex_loc[i]==edgeIndex_loc[0]&&
                faceIndex_loc[k]==edgeIndex_loc[1])||															
               (faceIndex_loc[i]==edgeIndex_loc[1]&&
                faceIndex_loc[k]==edgeIndex_loc[0]) )
           {
               edges[i] = Owner->edges[j];
               if(faceIndex_loc[i]==edgeIndex_loc[1]&&
                      faceIndex_loc[k]==edgeIndex_loc[0] )
               edges_orientation[i] = -1; 
               edges[i]->boundayC = 'B';
               break;
            } 
         }
    } 
}

//    WW. 06.2005
string Elem::getName() const
{
  
   switch(ele_Type)
   {
     case 1:
       return "line";
       break;
     case 2:
       return "quad";
       break;
     case 3:
       return "hex";
       break;
     case 4:
       return "tri";
       break;
     case 5:
       return "tet";
       break;
     case 6:
       return "pris";
       break;
     default:
       return "none";
       break;
   }
}


void Elem::setLocalNodeIndex(const int li, const long n_lindex) 
{ 
	nodes[li]->local_index = n_lindex; 
}

long Elem::getLocalNodeIndex(const int li) const 
{ 
	return nodes[li]->local_index; 
}

//    WW. 06.2005
void Elem::Read(istream& is, int fileType)
{
   //fileType=0: msh
   //fileType=1: rfi
   //fileType=2: gmsh
   //fileType=3: GMS
   //fileType=4: SOL
   int idummy, et;
   string buffer, name;
   idummy=et=-1;
//   is.ignore(numeric_limits<int>::max(), '\n');  
  //----------------------------------------------------------------------
  // 1 Reading element type data
  switch(fileType){
    //....................................................................
    case 0: // msh
      is>>index>>PatchIndex;
      is>>buffer;
	  if(buffer.find("-1")!=string::npos)
         is>>name;
	  else
	    name = buffer;
      if(name.find("line")!=string::npos)
         ele_Type = 1;
      else if(name.find("quad")!=string::npos)
         ele_Type = 2;
      else if(name.find("hex")!=string::npos)
         ele_Type = 3;
      else if(name.find("tri")!=string::npos)
         ele_Type = 4;
      else if(name.find("tet")!=string::npos)
         ele_Type = 5;
      else if(name.find("pri")!=string::npos)
         ele_Type = 6;
      break;
    //....................................................................
    case 1: // rfi
	  is>>index>>PatchIndex>>name;
      if(name.find("line")!=string::npos)
         ele_Type = 1;
      else if(name.find("quad")!=string::npos)
         ele_Type = 2;
      else if(name.find("hex")!=string::npos)
         ele_Type = 3;
      else if(name.find("tri")!=string::npos)
         ele_Type = 4;
      else if(name.find("tet")!=string::npos)
         ele_Type = 5;
      else if(name.find("pri")!=string::npos)
         ele_Type = 6;
      break;
    //....................................................................
    case 2: // gmsh
      int gmsh_patch_index; 
      is>>index>>et>>gmsh_patch_index>>idummy>>nnodes;
      PatchIndex = gmsh_patch_index-1; //OK
      switch(et)
      {
         case 1: ele_Type = 1; break;
         case 2: ele_Type = 4; break;
         case 3: ele_Type = 2; break;
         case 4: ele_Type = 5; break;
         case 5: ele_Type = 3; break;
         case 6: ele_Type = 6; break;
      }
	  index--;
      break;
    //....................................................................
    case 3: // GMS
      ele_Type = 4;
      break;
    //....................................................................
    case 4: // gmsh
      ele_Type = 4;
      break;
  }
  //----------------------------------------------------------------------
  // 2 Element configuration
  switch(ele_Type)
   {
      case 1:
         nnodes = 2;
         nnodesHQ = 3;    
         ele_dim = 1;
         ele_Type = 1;
         nfaces = 2;
         nedges = 0;
         break;
      case 2:
         nnodes = 4;
         nnodesHQ = 9;      
         ele_dim = 2;
         ele_Type = 2;
         nfaces = 4;
         nedges = 4;
         break;
      case 3:
         nnodes = 8;
         nnodesHQ = 20;
         ele_dim = 3;
         nfaces = 6;
         nedges = 12;
         ele_Type = 3;
         break;
      case 4:
         nnodes = 3;
         nnodesHQ = 6;
         ele_dim = 2;
         ele_Type = 4;
         nfaces = 3;
         nedges = 3;
         break;
      case 5:
         nnodes = 4;
         nnodesHQ = 10;
         ele_dim = 3;
         ele_Type = 5;
         nfaces = 4;
         nedges = 6;
         break;
      case 6:
         nnodes = 6;
         nnodesHQ = 15;
         ele_dim = 3;
         ele_Type = 6;
         nfaces = 5;
         nedges = 9;
		 break;
   }
   nodes_index.resize(nnodes);
  //----------------------------------------------------------------------
  // 3 Reading element node data
  switch(fileType){
    //....................................................................
    case 0: // msh
      for(int i=0; i<nnodes; i++)
        is>>nodes_index[i];
      break;
    //....................................................................
    case 1: // rfi
      for(int i=0; i<nnodes; i++)
        is>>nodes_index[i];
      break;
    //....................................................................
    case 2: // gmsh
      for(int i=0; i<nnodes; i++){
        is>>nodes_index[i];
        nodes_index[i] -= 1;
      }
      break;
    //....................................................................
    case 3: // GMS
      for(int i=0; i<nnodes; i++){
        is>>nodes_index[i];
        nodes_index[i] -= 1;
      }
      break;
    //....................................................................
    case 4: // SOL
      for(int i=0; i<nnodes; i++){
        is>>nodes_index[i];
        nodes_index[i] -= 1;
      }
      is >> PatchIndex;
      break;
  }
  is>>ws;
  //----------------------------------------------------------------------
  // Initialize topological properties
  neighbors.resize(nfaces);
  for(int i=0; i<nfaces; i++)
    neighbors[i] = NULL;
  edges.resize(nedges);    
  edges_orientation.resize(nedges);
  for(int i=0; i<nedges; i++)
  {
    edges[i] = NULL;
	edges_orientation[i] = 1;
  }
}
//  WW. 03.2009
void Elem::WriteGmsh(ostream& os,  const int sdom_idx) const
{
   //int igeo=14;
 
   int et=1;
   int ntags=3;
   string deli = " ";

   int nn = nnodes;
   if(quadratic)
   {
      nn = nnodesHQ;
      switch(ele_Type)
      {
        case 1: et = 8; break;    //Line
        case 2: et = 10; break;    //Quad
        case 3: et = 12; break;    //Hex
        case 4: et = 9; break;    //Tri
        case 5: et = 11; break;    //Tet
        case 6: et = 18; break;    //Pris
      }
   }
   else
   {
      switch(ele_Type)
      {
        case 1: et = 1; break;    //Line
        case 2: et = 3; break;    //Quad
        case 3: et = 5; break;    //Hex
        case 4: et = 2; break;    //Tri
        case 5: et = 4; break;    //Tet
        case 6: et = 6; break;    //Pris
      }
   }
   os<<index+1<<deli<<et<<deli<<ntags<<deli<<PatchIndex+1<<deli<<PatchIndex+1<<deli<<sdom_idx<<deli;
   for(int i=0; i<nn; i++)
   os<<nodes_index[i]+1<<deli;
   os<<endl;
}

//  WW. 03.2009
void Elem::WriteGSmsh(ostream& os, bool quad) const
{
   string ename;
   string deli = " ";
 
   int nn = getNodesNumber(quad);

   switch(ele_Type)
   {
     case 1: ename = "line"; break;
     case 2: ename = "quad"; break;
     case 3: ename = "hex"; break;
     case 4: ename = "tri"; break;
     case 5: ename = "tet"; break;
     case 6: ename = "pris"; break;
   }
   os<<index<<deli<<PatchIndex<<deli<<ename<<deli;
   for(int i=0; i<nn; i++)
   {
//      nodes_index[i] = nodes[i]->getIndex();
      os<<nodes[i]->getIndex()<<deli;
   }
   os<<endl;
}

//  WW. 02.2012
void Elem::WriteVTK_Type(ostream& os,  bool isquad) const
{
   if(!isquad) 
   {
      switch(ele_Type)
      {
         case 1:  os<< "3  "<<endl;    break;
         case 2:  os<< "9  "<<endl;    break;
         case 3:  os<< "12 "<<endl;    break;
         case 4:  os<< "5  "<<endl;    break;
         case 5:  os<< "10 "<<endl;    break;
         case 6:  os<< "13 "<<endl;    break;
     }
   }
   else
   {
      switch(ele_Type)
      {
         case 1:  os<< "21  "<<endl;    break;
         case 2:  os<< "23  "<<endl;    break;
         case 3:  os<< "25 "<<endl;    break;
         case 4:  os<< "22  "<<endl;    break;
         case 5:  os<< "24 "<<endl;    break;
         default:   break;
     }
   }
}

//    WW. 06.2005
void Elem::WriteIndex(ostream& os) const
{
    os<<index<<deli<<PatchIndex<<deli<<getName()<<deli;
    for(int i=0; i<nnodes; i++)
      os<<nodes[i]->index<<deli;
    os<<endl;
}
void Elem::Write_index(ostream& os) const
{
	if(nodes.Size()>0)
	{
      for(int i=0; i<nnodes; i++)
         os<<nodes[i]->index+1<<deli;
	}
	else
	{
      for(int i=0; i<nnodes; i++)
         os<<nodes_index[i]+1<<deli;
	}
    os<<endl;
}
//    WW. 06.2005
void Elem::WriteAll(ostream& os) const
{
    os<<index<<deli<<PatchIndex<<deli<<getName()<<deli;
    //if(index==0) 
    os<<"Index X Y Z: "<<endl;
				for(int i=0; i<nodes.Size(); i++)
    {
       os<<nodes_index[i]
       <<deli<<nodes[i]->X()
       <<deli<<nodes[i]->Y()
       <<deli<<nodes[i]->Z()<<endl;
    }
}

void Elem::WriteNeighbors(ostream& os) const
{
    os<<"Neighbors of "<<index<<endl;
    for(int i=0; i<nfaces; i++)
      neighbors[i]->WriteAll(os);
    os<<"End neighbors of "<<index<<endl<<endl;;
}

void Elem::MarkingNodes(bool maker)
{
    int SizeV = nnodes;
    if(quadratic) SizeV = nnodesHQ;
  
    for (int i=0; i< SizeV;i++)
    {
		nodes[i]->Marking(false);
    }
}


//    WW. 06.2005
void Elem::setNodes(vec<Node*>&  ele_nodes, const bool ReSize)
{ 
    int SizeV = nnodes;
    if(quadratic) SizeV = nnodesHQ;
    if(ReSize)
    {
        nodes.resize(SizeV);
       	nodes_index.resize(SizeV);
    }
    for (int i=0; i< SizeV;i++)
    {
       nodes[i] = ele_nodes[i];
       nodes_index[i] = nodes[i]->getIndex();       
    }
}




//    WW. 06.2005
void  Elem::getLocalIndices_EdgeNodes(const int Edge, int *EdgeNodes)
{
	switch(ele_Type)
	{
       case 1: 
           break; // 1-D bar element 
       case 2: // 2-D quadrilateral element 
          EdgeNodes[0] = Edge;
          EdgeNodes[1] = (Edge+1)%4;
          break;             
       case 3: // 3-D hexahedral element
          if(Edge<8)
          {
             EdgeNodes[0] = Edge;
             EdgeNodes[1] = (Edge+1)%4+4*(int)(Edge/4);
          }
          else 
          {
             EdgeNodes[0] = Edge%4;
             EdgeNodes[1] = Edge%4+4;
          }
          break;  
       case 4:  // 2-D triagular element 
          EdgeNodes[0] = Edge;
          EdgeNodes[1] = (Edge+1)%3;
          break;
       case 5:  // 3-D tetrahedra
          if(Edge<3)
          {
             EdgeNodes[0] = Edge;
             EdgeNodes[1] = (Edge+1)%3;
          }
          else
          {
             EdgeNodes[0] = 3;
             EdgeNodes[1] = (Edge+1)%3;
          }

          break;  
       case 6: // 3-D prismatic element
          if(Edge<6)
          {
             EdgeNodes[0] = Edge;
             EdgeNodes[1] = (Edge+1)%3+3*(int)(Edge/3);
          }
          else 
          {
             EdgeNodes[0] = Edge%3;
             EdgeNodes[1] = Edge%3+3;
          }
          break;  
	}
}


/**************************************************************************
GetElementFaceNodes
Task: Get local indeces of an element face nodes
Return: number of nodes of a face
Programing:
06/2004 WW  
**************************************************************************/
int Elem::getElementFaces1D(int *FaceNode)
{
    FaceNode[0] = 0;
    FaceNode[1] = 1;
    return 2;
}
/**************************************************************************
GetElementFaceNodesTri
Task: Get local indeces of a traingle element face nodes
Augs.:
        const int Face :  Local index of element face 
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
Return: number of nodes of a face
Programing:
09/2004 WW  
**************************************************************************/
int Elem::getElementFacesTri(const int Face, int *FaceNode)
{
	if(!quadratic)
    {
       FaceNode[0] = Face;
       FaceNode[1] = (Face+1)%3;
       return 2;
    }
    else
    {
       FaceNode[0] = Face;
       FaceNode[1] = Face+3;
       FaceNode[2] = (Face+1)%3;
       return 3;
    }
}

/**************************************************************************
GetElementFaceNodesQuad
Task: Get local indeces of a quadralateral element face nodes
Augs.:
        const int Face :  Local index of element face 
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
Return: number of nodes of a face
Programing:
09/2004 WW  
**************************************************************************/
int Elem::getElementFacesQuad(const int Face, int *FaceNode)
{
    if(!quadratic)
    {
       FaceNode[0] = Face;
       FaceNode[1] = (Face+1)%4;
       return 2;
    }
    else
    {
       FaceNode[0] = Face;
       FaceNode[1] = Face+4;
       FaceNode[2] = (Face+1)%4;
       return 3;
    }
}

/**************************************************************************
GetElementFaceNodesHex
Task: Get local indeces of a hexahedra element face nodes
Augs.:
        const int Face :  Local index of element face 
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
Return: number of nodes of a face
Programing:
09/2004 WW  
**************************************************************************/
int Elem::getElementFacesHex(const int Face, int *FaceNode)
{
   int nn=4, k = 0;
   if(quadratic) nn = 8;
   switch(Face)
   {
      case 0:
         for(k=0; k<4; k++)
            FaceNode[k] = k;
         if(quadratic)
         {
            for(k=0; k<4; k++)
                FaceNode[k+4] = k+8;
         }
         break;
      case 1:
         for(k=0; k<4; k++)
            FaceNode[k] = k+4;
         if(quadratic)
         {
            for(k=0; k<4; k++)
                FaceNode[k+4] = k+12;
         }
         break;
      case 2:
         FaceNode[0] = 0;
         FaceNode[1] = 4;
         FaceNode[2] = 5;
         FaceNode[3] = 1;
         if(quadratic)
         {
            FaceNode[4] = 16;
            FaceNode[5] = 12;
            FaceNode[6] = 17;
            FaceNode[7] = 8;
         }
         break;
      case 3:
         FaceNode[0] = 1;
         FaceNode[1] = 5;
         FaceNode[2] = 6;
         FaceNode[3] = 2;
         if(quadratic)
         {
            FaceNode[4] = 17;
            FaceNode[5] = 13;
            FaceNode[6] = 18;
            FaceNode[7] = 9;
         }      

         break;
      case 4:
         FaceNode[0] = 2;
         FaceNode[1] = 6;
         FaceNode[2] = 7;
         FaceNode[3] = 3;
         if(quadratic)
         {
            FaceNode[4] = 18;
            FaceNode[5] = 14;
            FaceNode[6] = 19;
            FaceNode[7] = 10;
         }
         break;
      case 5:
         FaceNode[0] = 0;
         FaceNode[1] = 3;
         FaceNode[2] = 7;
         FaceNode[3] = 4;
         if(quadratic)
         {
            FaceNode[4] = 11;
            FaceNode[5] = 19;
            FaceNode[6] = 15;
            FaceNode[7] = 16;
         }
         break;
   }
   return nn;
}


/**************************************************************************
GetElementFaceNodesTet
Task: Get local indeces of a Tedrahedra element face nodes
Augs.:
        const int Face :  Local index of element face 
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
Return: number of nodes of a face
Programing:
09/2004 WW  
**************************************************************************/
int Elem::getElementFacesTet(const int Face, int *FaceNode)
{
   int nn=3;
   if(quadratic) nn =6;
   switch(Face)
   {
      case 0:
         FaceNode[0] = 1;
         FaceNode[1] = 2;
         FaceNode[2] = 3;
         if(quadratic)
         {
            FaceNode[3] = 5 ;
            FaceNode[4] = 8;
            FaceNode[5] = 7;
         }
         break;
      case 1:
         FaceNode[0] = 3;
         FaceNode[1] = 2;
         FaceNode[2] = 0;
         if(quadratic)
         {
            FaceNode[3] = 8 ;
            FaceNode[4] = 6;
            FaceNode[5] = 9;
         }
         break;
      case 2:
         FaceNode[0] = 1;
         FaceNode[1] = 3;
         FaceNode[2] = 0;
         if(quadratic)
         {
            FaceNode[3] = 7 ;
            FaceNode[4] = 9;
            FaceNode[5] = 4;
         }
         break;
      case 3:
         FaceNode[0] = 0;
         FaceNode[1] = 2;
         FaceNode[2] = 1;
         if(quadratic)
         {
            FaceNode[3] = 6 ;
            FaceNode[4] = 5;
            FaceNode[5] = 4;
         }
         break;

   }
   return nn;
}


/**************************************************************************
GetElementFaceNodesPri
Task: Get local indeces of a prismal element face nodes
Augs.:
        const int Face :  Local index of element face 
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes
Return: number of nodes of a face
Programing:
09/2004 WW  
**************************************************************************/
int Elem::getElementFacesPri(const int Face, int *FaceNode)
{
   int nn=3, k = 0;
   switch(Face)
   {
      case 0:
         nn = 3;
         for(k=0; k<3; k++)
            FaceNode[k] = k;
         if(quadratic)
         {
            for(k=0; k<3; k++)
                FaceNode[k+3] = k+6;
            nn = 6;
         }
         break;
      case 1:
         for(k=0; k<3; k++)
            FaceNode[k] = k+3;
         nn = 3;
         if(quadratic)
         {
            for(k=0; k<3; k++)
                FaceNode[k+3] = k+9;
            nn = 6;
         }
         break;
      case 2:
         FaceNode[0] = 1;
         FaceNode[1] = 2;
         FaceNode[2] = 5;
         FaceNode[3] = 4;
         nn = 4;
         if(quadratic)
         {
            FaceNode[4] = 7 ;
            FaceNode[5] = 14;
            FaceNode[6] = 10;
            FaceNode[7] = 13;
            nn = 8;
         }
         break;
      case 3:
         FaceNode[0] = 5;
         FaceNode[1] = 2;
         FaceNode[2] = 0;
         FaceNode[3] = 3;
         nn = 4;
         if(quadratic)
         {
            FaceNode[4] = 14 ;
            FaceNode[5] =  8;
            FaceNode[6] = 12;
            FaceNode[7] = 10;
            nn = 8;
         }
         break;
      case 4:
         FaceNode[0] = 0;
         FaceNode[1] = 1;
         FaceNode[2] = 4;
         FaceNode[3] = 3;
         nn = 4;
         if(quadratic)
         {
            FaceNode[4] = 6 ;
            FaceNode[5] = 13;
            FaceNode[6] = 9;
            FaceNode[7] = 12;
            nn = 8;
         }
         break;

   }
   return nn;
}



/**************************************************************************
GetElementFaces
Task: set element faces (Geometry)
Augs.:
        const int Face :  Local index of element face 
        const int order:  1 Linear. 2, quadratic
        int *FaceNode  :  Local index of face nodes

Programing:
09/2004 WW  
**************************************************************************/
int Elem::getElementFaceNodes(const int Face, int *FacesNode)
{
   switch(ele_Type)
   {
       case 1:  // 1-D bar element
           return getElementFaces1D(FacesNode);
           break;          
       case 2: // 2-D quadrilateral element
           return getElementFacesQuad(Face, FacesNode);
           break;           
       case 3: // 3-D hexahedral element 
           return getElementFacesHex(Face, FacesNode);
           break;           
       case 4:  // 2-D triagular element 
           return getElementFacesTri(Face, FacesNode);
           break;           
       case 5:  // 3-D tetrahedral element 
           return getElementFacesTet(Face, FacesNode);
           break;           
       case 6: 
           return getElementFacesPri(Face, FacesNode);
           break; // 3-D prismatic element 
    }
    return 0;
}



}//end namespace

