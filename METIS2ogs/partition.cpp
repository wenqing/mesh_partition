// expre_new_Operator.cpp
// compile with: /EHsc
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <time.h>
#include <vector>

#ifdef USE_METIS_SOURCE
extern "C" {
//#include "../metis-5.0.2/programs/metis_main.h"
#include "metis_main.h"
}
#endif

#include "Mesh.h"

using namespace std;
using namespace Mesh_Group;


/// returns used heap size in bytes or negative if heap is corrupted.
#ifdef WIN32
long HeapUsed()
{
   _HEAPINFO info = { 0, 0, 0 };
   long used = 0;
   int rc;

   while ((rc=_heapwalk(&info)) == _HEAPOK)
   {
      if (info._useflag == _USEDENTRY)
         used += (long)info._size;
   }
   if (rc != _HEAPEND && rc != _HEAPEMPTY)
      used = (used?-used:-1);

   return used;
}
#endif



#define ver "V3.1. 05.2014"

void Version()
{
   cout<<"\nOpenGeoSys interface to partitioning tools"<<endl;
   cout<<"\nCopyleft. If compiled with option USE_METIS_SOURCE, please refer to the license of METIS (or ohter tools)."<<endl;
   cout<<ver<<endl;
   cout<<"Written by wenqing.wang@ufz.de."<<endl<<endl;
}
void OptionList()
{
   string s_intro = "The task of this tool is twofold:"
                    " to convert ogs mesh file into the partitioning tool input file for domain decomposition,"
                    " and to use the paritioning tool's results to partition"
                    "the ogs finite element meshes for parallel computing.\n"
                    "Note: input mesh file must be given with its absolute path if option -mat is used.";
   cout << s_intro<<endl<<endl;
   cout << "Tasks:\n  --version\n  --help\n  --ogs2metis\n  --metis2ogs\n"<<endl;
   cout << "Option for --metis2ogs task:"<<endl;
   cout << "  -q                : generate quadratic elements. It can be ommitted if quadratic element is not used."<<endl;
   cout << "  -np [number]      : define the number of partitions."<<endl;
   cout << "  -e                : partition by element (non overlapped subdomain)"<<endl;
   cout << "  -n                : partition by node (overlapped subdomain)"<<endl;
   cout << "  -mat [file name without path]  : specify a file that contains file name of element-wsie material data."<<endl;
   cout << "  -asci             : ASCI output the partitioned mesh."<<endl;
   cout << "  -nvtk             : do not output the vtk file for the node index renumbered whole mesh."<<endl;
   cout << "  -ogsmsh           : output the renumbered ogs mesh."<<endl;
   cout << "  -odom             : output subdomain mesh."<<endl;
   cout << "  -cct              : output a CCT file for FCT (overlapped subdomain)"<<endl;
}

void FindFileNameInCommand(stringstream &ss, string &fname)
{
   while(!ss.eof())
   {
      ss>>fname;
      if(   !(fname.find("-")!=string::npos || fname.find("--")!=string::npos) )
      {
         return;
      }
   }
}

enum Task {metis2ogs, ogs2metis};
enum PartType {by_element, by_node};

int main(int argc, char* argv[])
{
   ifstream infile;
   fstream ofile;
   stringstream ss;

   Task this_task;
   PartType part_type = by_node;

   bool quad = false;


   Mesh_Group::MeshPartConfig mpc;
   mpc.is_vtk_out = true;
   mpc.osdom = false;
   mpc.out_renum_gsmsh = false;
   mpc.binary_output = true;
   mpc.out_cct = false;

   //ios::pos_type position;

   string s_buff;
   string fname;
   string fpath;
   string mat_file_name = "";

   int nparts = 1;
   string str_nparts;

   bool task_opt_given = false;
   bool part_opt_given = false;
   bool part_num_given = false;

   if(argc>1)
   {
      for(int i=1; i<argc; i++)
      {
         s_buff = argv[i];
         if(s_buff.compare("-e") == 0)
         {
            part_type = by_element;
            part_opt_given = true;
         }
         else if(s_buff.compare("-n") == 0)
         {
            part_type = by_node;
            part_opt_given = true;
         }

         if(s_buff.compare("-q") == 0)
            quad = true;

         if(s_buff.compare("-odom") == 0)
         {
            mpc.osdom = true;
         }

         if(s_buff.compare("-cct") == 0)
            mpc.out_cct = true;

         if(s_buff.find("-asci")!=string::npos)
         {
            mpc.binary_output = false;
         }
         if(s_buff.compare("-ogsmsh") == 0)
         {
            mpc.out_renum_gsmsh = true;
         }
         if(s_buff.compare("-nvtk") == 0)
         {
            mpc.is_vtk_out = false;
         }

         // Number of partitions
         if(s_buff.find("-np")!=string::npos)
         {
            //size_t pos;
            //pos = s_buff.find_first_of("p");
            //s_buff = s_buff.substr(pos+1);
            //str_nparts = s_buff;
            nparts = atoi( argv[i+1]);
            str_nparts =  argv[i+1];
            part_num_given = true;
         }

         // Number of partitions
         if(s_buff.find("-mat")!=string::npos)
         {
            mat_file_name = argv[i+1];
         }

         if(s_buff.find("ogs2metis")!=string::npos)
         {
            this_task = ogs2metis;
            task_opt_given = true;
         }
         else if(s_buff.find("metis2ogs")!=string::npos)
         {
            this_task = metis2ogs;
            task_opt_given = true;
         }
         else if(s_buff.find("--help")!=string::npos)
         {
            Version();
            OptionList();
            exit(0);
         }
         else if(s_buff.find("--version")!=string::npos)
         {
            cout<<ver;
            exit(0);
         }

         if(  s_buff[0] != '-')
         {
            fname = s_buff;
         }

      }
   }
   else //terminal
   {
      OptionList();
      Version();
      cout<<"\nInput task, options and file name (non extension):\n ";

      getline(cin, s_buff);
      ss.str(s_buff);
      if(s_buff.find("ogs2metis")!=string::npos)
      {
         this_task = ogs2metis;
         if(s_buff.find("-e")!=string::npos || s_buff.find("-n")!=string::npos||s_buff.find("-q")!=string::npos)
         {
            cout<<"Warning: option is not needed for this task"<<endl;
         }

         FindFileNameInCommand(ss, fname);
         task_opt_given = true;
      }
      else if(s_buff.find("metis2ogs")!=string::npos)
      {
         task_opt_given = true;
         this_task = metis2ogs;
         while(!ss.eof())
         {
            ss>>s_buff;

            if(s_buff.compare("-e") == 0)
            {
               part_type = by_element;
               part_opt_given = true;
            }
            else if(s_buff.compare("-n") == 0)
            {
               part_type = by_node;
               part_opt_given = true;
            }

            if(s_buff.compare("-q") == 0)
            {
               quad = true;
            }
            if(s_buff.compare("-odom") == 0)
            {
               mpc.osdom = true;
            }
            if(s_buff.compare("-cct") == 0)
            {
               mpc.out_cct = true;
            }
            if(s_buff.compare("-asci") == 0)
            {
               mpc.binary_output = false;
            }
            if(s_buff.compare("-ogsmsh") == 0)
            {
               mpc.out_renum_gsmsh = true;
            }

            if(s_buff.compare("-nvtk") == 0)
            {
               mpc.is_vtk_out = false;
            }
            if(s_buff.find("-np")!=string::npos)
            {
               //size_t pos;
               //pos = s_buff.find_first_of("p");
               //s_buff = s_buff.substr(pos+1);
               //nparts = atoi(s_buff.c_str());
               //str_nparts = s_buff;
               ss >> str_nparts;
               nparts = atoi(str_nparts.c_str());

               part_num_given = true;
            }

            if(s_buff.find("-mat")!=string::npos)
            {
               ss >> mat_file_name;
            }

            if(  s_buff[0] != '-' )
            {
               fname = s_buff;
            }
         }
      }
      ss.clear();

   }

   if( !task_opt_given )
   {
      cout<<"Task option (--metis2ogs or --ogs2metis) is not given. Stop now.\n";
      exit(EXIT_FAILURE);
   }

   if( this_task == metis2ogs)
   {
      if( !part_opt_given )
      {
         cout<<"Partitioning option (-n or -e) is not given. Stop now.\n";
         exit(EXIT_FAILURE);
      }
      if( !part_num_given )
      {
         cout<<"Partitioning number (e.g. -np 3) is not given. Stop now.\n";
         exit(EXIT_FAILURE);
      }
   }

   //Get the path to the folder where the input file is.
   size_t pos_end;
   pos_end = fname.find_last_of('\\');
   //
   if(pos_end != std::string::npos)
   {
      fpath = fname.substr(0,pos_end) + "\\";
   }
   else
   {
      pos_end = fname.find_last_of('/');
      if(pos_end != std::string::npos)
         fpath = fname.substr(0,pos_end) + "/";
   }


   s_buff = fname+".msh";
   infile.open(s_buff.c_str());
   if(!infile.is_open())
   {
      cerr<<("Error: cannot open msh file . It may not exist !");
      exit(1);
   }

   cout<<"File name is: "<<fname<<endl;
   if(fpath.size()>0)
      cout<<"File path is: "<<fpath<<endl;
   else
      cout<<"File path is: ./ "<<endl;


   clock_t elp_time;
   elp_time = -clock();

   Mesh_Group::Mesh *a_mesh = new Mesh(quad);

   bool rfiMesh = true;
   string line_string;
   getline(infile,line_string); // The first line
   if(line_string.find("#FEM_MSH")!=string::npos)
      rfiMesh = false;
   if(line_string.find("GeoSys-MSH")!=string::npos)
      rfiMesh = false;
   infile.seekg(0L,ios::beg);

   if(rfiMesh)
      a_mesh->ReadGrid(infile, quad);
   else
      a_mesh->ReadGridGeoSys(infile, quad);


   switch(this_task)
   {
      case ogs2metis:
         s_buff = fname+".mesh";
         ofile.open(s_buff.c_str(), ios::out | ios::trunc );
         a_mesh->Write2METIS(ofile);

         break;
      case metis2ogs:
         cout<<"\n***Compute mesh topology"<<endl;
         a_mesh->ConstructGrid();

         /// Partition mesh if metis source is include
         if(nparts>1)
         {
#ifdef USE_METIS_SOURCE
            int argc_m;
            argc_m = 3;
            char *argv_m[3];
            argv_m[0] = "-";
            s_buff = fname + ".mesh";

            argv_m[1] = &s_buff[0];
            argv_m[2] = &str_nparts[0];

            metis_main(argc_m, argv_m);
//#else
//         s_buff = fpath+"mpmetis "  + fname + ".mesh " + str_nparts;

//         if(!system(s_buff.c_str()))
//         {
//            cout<<"METIS executable file may not be found "<<endl;
//            exit(1);
//         }
#endif
         }

         cout<<"\n***Prepare subdomain mesh"<<endl;
         if(part_type == by_element)
            a_mesh->ConstructSubDomain_by_Elements(fname.c_str(), nparts, mpc.osdom);
         else if(part_type == by_node)
         {
            if(quad)
            {
               a_mesh->GenerateHighOrderNodes();
               a_mesh->setOrder(quad);
            }

            mpc.fname = fname;
            mpc.fpath = fpath;
            mpc.mat_fname = mat_file_name;
            mpc.num_parts = nparts;

            a_mesh->ConstructSubDomain_by_Nodes(mpc);
         }
         break;
      default:
         break;
   }


#ifdef WIN32
   cout<<"\n\tMemory usage: "<< HeapUsed()/1024./1024.<<"MB"<<endl;
#endif

   delete a_mesh;

   elp_time += clock();
   cout<<"\n***Total CPU time elapsed: "
       <<(double)elp_time / CLOCKS_PER_SEC<<"s"<<endl;

   return 0;

}
