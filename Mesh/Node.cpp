#include "Node.h"

#include <iomanip>
//------------------------------------------------------
//   Topology definition of geometrical element.
//    WW. 10.01.2005
//------------------------------------------------------
namespace Mesh_Group
{

using namespace std;

//1.  Node declaration
//    WW. 06.2005
Node:: Node(const MyInt Index, const double x,
            const double y, const double z):Grain(Index),global_index(Index)
{
   Coordinate = new double[3];
   Coordinate[0] = x;
   Coordinate[1] =y;
   Coordinate[2] =z;
   local_index = -1;
   index_org = Index;
}
//    WW. 06.2005
void Node::operator = (const Node& n)
{
   index = n.index;
   mark = n.mark;
   Coordinate[0] = n.Coordinate[0];
   Coordinate[1] = n.Coordinate[1];
   Coordinate[2] = n.Coordinate[2];
}
//    WW. 06.2005
bool Node::operator == (const Node& n)
{
   if(index == n.index)
      return true;
   else
      return false;
}
//    WW. 06.2005
// Output
void Node::Write(ostream& os) const
{
   os<<setw(14)<<index<<" ";
   os<<Coordinate[0]<<" "
     <<Coordinate[1]<<" "
     <<Coordinate[2]<<endl;
}

// WW 07.2013
void Node::WriteBIN(ostream& os) const
{

   Node_Str nd;
   nd.id = static_cast<MyInt> (index);
   nd.x = Coordinate[0];
   nd.y = Coordinate[1];
   nd.z = Coordinate[2];

   os.write( reinterpret_cast <const char*> (&nd), sizeof(Node_Str));
}

void Node::WriteCoordinates(std::ostream& os) const
{
   os <<Coordinate[0]<<" "
      <<Coordinate[1]<<" "
      <<Coordinate[2]<<endl;
}

//    WW. 06.2005
// Set
void Node::SetCoordinates(const double* argCoord)
{
   Coordinate[0] = argCoord[0];
   Coordinate[1] = argCoord[1];
   Coordinate[2] = argCoord[2];
}

}

