#ifndef Grain_INC
#define Grain_INC

#include<string>
#include<iostream>
 
//------------------------------------------------------ 
//   Topology declartion of geometrical element.
//   WW. 06.2005
//   WW  02.2012
//------------------------------------------------------ 

namespace Mesh_Group
{
 enum BCType  { INTERIOR, DIRCHLET, NEUMANN, CAUCHY, BOUNDARY };
 //enum Element_Type  {Line, Quad, Hex, Tri, Tet, Pris };

//1.  Mesh declaration
class Grain
{
   public:
      Grain(const int id);
      virtual  ~Grain() {}
      // Operator
      virtual void operator = (const Grain & g) {}
      virtual bool operator == (const Grain & g) {return false;}

      // Set members
      void secBC(const char BC_type) {boundayC = BC_type;}
      void SetOrder(const bool order) {quadratic = order;}
      void Marking(const bool state) {mark = state;}
      // Get members
      long GetIndex() const {return index;} 
      bool GetStatus() const {return mark;}
      bool Dirichlet() const  { return (BCType(boundayC) == DIRCHLET); }
      bool Neumann()   const  { return (BCType(boundayC) == NEUMANN);   }
      bool Cauchy ()   const  { return (BCType(boundayC) == CAUCHY);    }
      bool onBoundary() const { return (BCType(boundayC) == BOUNDARY);  }
      bool Interior() const { return (BCType(boundayC) == INTERIOR);  }

      // Output
      virtual void output(std::ostream& os = std::cout) const {};
   protected:
      long index;
      char boundayC;
      // Towards special purpose, 
      // e.g. marked to be refined or active 
      bool mark;
      // High order
      bool quadratic;
    
      // delimitor
      std::string deli;

};
    
} //end namespace

#endif
