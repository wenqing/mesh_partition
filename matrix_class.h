/*========================================================================
 GeoSys - class Matrix (Declaration)
 Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation. 
 Function:   See the declaration below
 programming:
  22/08/2004  WW  
==========================================================================*/
#ifndef matrix_class_INC

#define matrix_class_INC

#include<iostream>
#include<fstream>

//#define OverLoadNEW_DELETE

namespace Math_Group{
  using namespace std;

class Matrix
{
   public:
     Matrix(const int rows, const int cols=1);
     Matrix();
     explicit Matrix(const Matrix& m);
     
     void resize(const int rows, const int cols=1);


     virtual ~Matrix();

//----------------------------------------------
#ifdef OverLoadNEW_DELETE
     // Allocate memory
     void* operator new(size_t sz);
     // Release memory
     void operator delete(void* m);
#endif     
//----------------------------------------------
     // Operators
     virtual void operator = (const double a);
     virtual void operator *= (const double a);
     virtual void operator += (const double a);
     void operator = (const Matrix& m);
     void operator += (const Matrix& m);
     void operator -= (const Matrix& m);

     void GetTranspose(Matrix& m);

	 // vec_result = This*vec. vec_result must be initialized
     void multi(const double *vec, double *vec_result, const double fac=1.0);
     // m_result = this*m. m_result must be initialized
	 void multi(const Matrix& m, Matrix& m_result, const double fac=1.0);
	 // m_result = this*m1*m2. m_result must be initialized
     void multi(const Matrix& m1, const Matrix& m2, Matrix& m_result);

     // Access to members     
     virtual double& operator() (const int i, const int j=0) const;  
     void LimitSize(const int nRows, const int nCols=1);  

     int Rows() const {return nrows;}
     int Cols() const {return ncols;}
     int Size() const {return size;}

     // Print
     void Write(ostream& os=cout);
     void Write_BIN(fstream& os);
     void Read_BIN(fstream& is);
   protected:
     double *data;     
     int nrows, nrows0;
     int ncols, ncols0;  
     int size; 
     bool Sym;
};

// Symmetrical matrix. 12-01-2005. WW
class SymMatrix:public Matrix
{
   public:
     SymMatrix(const int dim);
     SymMatrix();
     explicit SymMatrix(const SymMatrix& m);
     
     void resize(const int dim);

     ~SymMatrix() {}

//----------------------------------------------
#ifdef OverLoadNEW_DELETE
     // Allocate memory
     void* operator new(size_t sz);
#endif     
//----------------------------------------------

     // Operators
     void operator = (const double a);
     void operator *= (const double a);
     void operator += (const double a);
     void operator = (const SymMatrix& m);
     void operator += (const SymMatrix& m);
     void operator -= (const SymMatrix& m);
     void LimitSize(const int dim);  

     // Access to members     
     double& operator() (const int i, const int j) const;  
};
typedef Matrix Vec;

/*========================================================================
 GeoSys - class my_vector (Declaration)
 Task:       Carry out vector operation 
 Function:   See the declaration below
 programming:
 05/2005  WW  
==========================================================================*/
template<class T> class vec
{									 
    public:
      vec(const int argSize);
	  vec() {size=0;}
      explicit vec(const vec<T>& v);

      virtual ~vec();
      // Operator
      virtual void operator = (T v) { for (int i=0; i<size; i++) entry[i] = v; }
      virtual void operator = (const vec<T>&);
	  virtual void resize(const int newh);
      virtual T& operator[] (int i) { return (T&) entry[i]; } 
      virtual const T& operator[] (int i) const {return (const T&) entry[i]; }
      virtual int Size() const { return size; } 
 
      T* Entry()        { return entry; }
      T* Entry()  const { return entry; }
 
      virtual void Write(ostream& os=cout) const;
    protected:	      
      T* entry;						 
      int size;

};
template<> class vec<void*> 
 {		
    public:
      vec(const int argSize);
	  vec() {size=0;}
      explicit vec(const vec<void*>& v);

      virtual ~vec();
      // Operator
      void operator = (void* v)  
	      { for (int i=0; i<size; i++) entry[i] = v; }
      void operator = (const vec<void*>& v);
      void*& operator[] (int i) { return entry[i]; } 
      const void*& operator[] (int i) const { return (const void*&) entry[i]; }
    
      // Access to memebers
      void** Entry()          { return entry; }
      const void** Entry()  const { return (const void**)entry; }

	  virtual void resize(const int newh);
      virtual int Size() const { return size; }           
      virtual void Write(ostream& os=cout) const;

    protected:	      
      void** entry;						 
      int size;
};
  
template<class T> class vec<T*> : public vec<void*>
{									 
   public:
      vec(const int Size) : vec<void*>(Size) { }
      vec()               : vec<void*>()     { }
      explicit vec(const vec<T*>& v) : vec<void*>(v) { }

      ~vec() { }
      
	  // Operator
      void operator = (T* v) { for (int i=0; i<size;i++) entry[i] = v; }
      void operator = (const vec<T*>& v);            
      T*& operator[] (int i) { return (T*&) entry[i]; } 
      const T*& operator[] (int i) const {return (const T*&) entry[i]; }

      T** Entry()        { return entry; }
      T** Entry()  const { return (const T**)entry; }


};
  
// Cross production x^y. WW 12.01.2005
//const Vec& operator ^ (Vec& x,  Vec& y);	 

#ifdef NewSparseMatrix
//WW
// Class definition
namespace Mesh_Group {class CFEMesh;}
using Mesh_Group::CFEMesh;
class SparseMatrix
{
   public:
     SparseMatrix(CFEMesh *a_mesh);
     ~SparseMatrix();
     virtual void operator = (const double a);
     virtual void operator *= (const double a);
     virtual void operator += (const double a);
     virtual void operator = (const SparseMatrix& m);
     virtual void operator += (const SparseMatrix& m);
     virtual void operator -= (const SparseMatrix& m);

     void GetTranspose(SparseMatrix& m);

	 // vec_result = This*vec. vec_result must be initialized
     void multi(const double *vec, double *vec_result, const double fac=1.0);
     // m_result = this*m. m_result must be initialized
	 void multi(const SparseMatrix& m, SparseMatrix& m_result, const double fac=1.0);
	 // m_result = this*m1*m2. m_result must be initialized
     void multi(const SparseMatrix& m1, const SparseMatrix& m2, SparseMatrix& m_result);

     // Access to members     
     virtual double& operator() (const int i, const int j=0) const;  

     int Rows() const {return nrows;}
     int Cols() const {return ncols;}
     int Size() const {return size;}

     // Print
     void Write(ostream& os=cout);
     
   private:
     // Data
     double* upper_triangle;       
     double* low_triangle;       
     double* diagnal;       
     // Index
     long* tree;
     long* row_in_tree;
     long* column_in_tree;
}
#endif
// End of class Matrix
}
//==========================================================================

#endif
