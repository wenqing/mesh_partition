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

//2.  Node declaration
class Node:public Grain
{
   public:
      Node(const int Index):Grain(Index)
	  {
         Coordinate = new double[3];
	  }
      Node(const int Index, const double x, 
           const double y, const double z=0.0);
      ~Node()
	   {
          delete [] Coordinate;
		  Coordinate = NULL;
		  ElementsBelonged.resize(0);
	   }

      // Operator
      void operator = (const Node& n);
      bool operator == (const Node & n);

      // Change members;
      // By component
      void setX(const double argX)  { Coordinate[0] = argX;}
      void setY(const double argY)  { Coordinate[1] = argY;}
      void setZ(const double argZ)  { Coordinate[2] = argZ;}
      void SetCoordinates(const double* argCoord); 

      // Access to members
      double X() const {return Coordinate[0];}
      double Y() const {return Coordinate[1];}
      double Z() const {return Coordinate[2];}
      void Coordinates(double *xyz) const
      { 
         for(int i=0; i<3; i++)  xyz[i] = Coordinate[i];
	  } 

	  void SetLocalIndex(const long l_index) {local_index = l_index; } 
	  long GetLocalIndex() const {return local_index;}
     
      // Output
      void Write(std::ostream& os = std::cout) const;

      std::vector<long>  ElementsBelonged;

   private:
      double *Coordinate;
      long local_index; // For domain decomposition
      friend class Edge;
      friend class Elem;
};

} 

#endif
