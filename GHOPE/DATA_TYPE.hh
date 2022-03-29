#ifndef DATA_TYPE_HH_
#define DATA_TYPE_HH_

#define _CHECK_BOUNDS_ 1

#include <cstddef>
#include <initializer_list>
#include <iterator>     // std::iterator_traits
#include <typeinfo>     // typeid
#include <iostream>
#include <ostream>
#include <vector>

//Contraction of vector by multiplication
inline long int inline_contract(const int* a, const int n) {
	long int A(1);
	for(int i = 0; i < n; ++i) {A *= a[i];}
	return A;
}

// ############## INDEXING EXPLENATION ###################
/*
int DIMS[];

Each dimension is given a size that is equal to the number of copies of the remaing dimension, e.g for a matrix
    DIMS[0] = Row size = Number of columns 
    DIMS[1] = Column size = Number of rows

Or in a hypercube
	DIMS[0] = Row size = Number of columns in a sheeth
    DIMS[1] = Column size = Number of rows in a sheeth
    DIMS[2] = Depth size = number of sheeth's

and so on.

The indexing to a vector is done by the cumulatuive chunking dependant on size

Consider a M1 x M2 x M3 x ... arbitrary tensor then the tensor index
T-index = (n1,n2,n3...)
maps to the storage vector as
S-index = n1*1 + M1*n2 + M1*M2*n3 + .... + product(i=1...x-1, Mi)*nx + ....

Thus the data will be stored in one consistent data allocation of configuration
| ŕow 1 of matrix 1 of hypercube 1 ...|, | row 2 of matrix 1 of .... |
and so on

#### FETCHING DATA

This directly translates to the fetching of data from this data structure.
To get data from row 5 column 3

one should think about it as accessing row 5 and picking the 3'd element
thus one needs to enter the element located at row INDEX 3 (not row, but index IN row) and ćolumn INDEX 5
as X(3,5) and NOT X(5,3)

This simplifies for higher dimensional structures since we have no real 
semantic convention for rows and columns in hyper geometry, e.g. to select the same element as above but in the
5th matrix sheeth of a cube one would simply
X(3,5,5)

As long as one keeps track of the difference between size and index no problems are expected.

FOR A MORE MATRIX OR MATLAB LIKE FETCHING SYSTEM USE .get(row, col) where the indecies have been reversed so that
To get data from row 5 column 3 in matrix X is
X.get(5,3)


*/

// #################### EXAMPLES BELOW ###########################
/*

tensor<double> a; // Unassigned tensor of doubles

int n[] = {100,100};
tensor<double> a1(n,2); // 100 X 100 matrix of unassigned doubles
tensor<double> a2(n,2,6.5); // 100 X 100 matrix of doubles filled with 6.5
tensor<double> b(a2); // exact copy of a2

int n2[] = {100,100,100};
tensor<double> bt(n2,3,0.0);

int getn[] = {6,4,4};
long int getZ;
double X = bt[bt.vec2ind(getn)]; //Get data direct from storage vector
X = bt[6 + 4*n2[0] + 4*n2[1]]; //equivalent expression
X = bt[6*bt.dim_C(0) + 4*bt.dim_C(1) + 4*bt.dim_C(2)]; //equivalent expression

X = bt[getn]; //Get data by tensor indecies, slightly slower

//Iterate trough storage vector with index
for(long int K=0; K < bt.size(); K++) {
    bt[K] = rand();
    X = bt[K];
}


//Iterate trough storage vector with iterator
double * T_ITER;
for(T_ITER = bt.begin(); T_ITER <= bt.end(); T_ITER++) {
    *T_ITER = (double)rand();
    X = *T_ITER;
}


//Iterate trough tensor with indecies
for(int k0 = 0; k0 < bt.size(0); k0++) {
    getn[0] = k0;
    for(int k1 = 0; k1 < bt.size(1); k1++) {
        getn[1] = k1;
        for(int k2 = 0; k2 < bt.size(2); k2++) {
            getn[2] = k2;
            bt[getn] = rand();
        }
    }
}

//Iterate trough tensor with indecies and conversion
for(int k0 = 0; k0 < bt.size(0); k0++) {
    for(int k1 = 0; k1 < bt.size(1); k1++) {
        for(int k2 = 0; k2 < bt.size(2); k2++) {
        	//bt.dim_C(1) = 1;
			getZ = k0 + k1*bt.dim_C(1) + k2*bt.dim_C(2);
            bt[getZ] = rand();
        }
    }
}

//Iterate trough part of tensor with iterator
        T_ITER = bt.begin() + 5 + bt.dim_C(1)*6;
        for(int k2 = 0; k2 < bt.size(2); k2++) {
            *T_ITER = (double)rand();
            T_ITER+=bt.dim_C(2);
        }


*/
// #################### EXAMPLES ABOVE ###########################



//We here define a multi-dimensional class tensor optmized for HPC usage

template <class T>
class tensor {
private:
	int m;		// Number if indices
	int *n;		// Size of each index
	int *cn;	// Cumulative size of each index
	long int Z;	// Total size if the vector storage space
	T *v;		// The actual storage space is stored as a pointer to a vector
				// as sequential memory is often preferable in performance applications
public:
	//Conversion of index vector to index
	inline long int vec2ind(const int* a);
	inline void ind2vec(int *a,const long int I);
	//Iterator
  	T* begin() { return Z>0 ? &v[0] : NULL; } //Pointer to first element
  	T* end() { return Z>0 ? &v[Z-1] : NULL; } //Pointer to last element

  	int dim_C(int D) const { return cn[D]; } //Cumulative multiplicative index transform from N-D tensor to 1-D vector
  	int dim_N() const { return m; }			 //The number of dimensions of tensor

	inline T back() { return Z>0 ? v[Z] : NULL; } //Last element in container (hypercube end corner of diagonal)

	tensor();										// Empty tensor without size
	explicit tensor(int N);							// Empty vector of size N
	explicit tensor(int N, int M);					// Empty matrix of size N x M
	explicit tensor(int N, int M, int L);			// Empty cube of size N x M x L
	explicit tensor(int N, const T &a);							// Filled vector of size N
	explicit tensor(int N, int M, const T &a);					// Filled matrix of size N x M
	explicit tensor(int N, int M, int L, const T &a);			// Filled cube of size N x M x L
	explicit tensor(int* N, int M);					// Empty tensor with size
	tensor(int* N, int M, const T &a);				// Initialize to constant value
	tensor(const tensor &rhs);						// Copy constructor
	typedef T value_type; 							// Make the type available externally
	tensor & operator=(const tensor &rhs);			// Assignment operator
	//tensor & slice(const int i, const int D);		// Take the i'th slice in dimension D
	//tensor & vector(const int i, const int P, const int D);		// Take the i'th vector in dimension P along dimension D
	
	void print(std::ostream& out);

	tensor row(const int i); //If m = 2, get row i
	tensor col(const int i); //If m = 2, get col i

	//THIS FUNCTION DOES NOT RESCALE CONTAINERS
	inline void copy(const tensor &rhs); 			// Copy data (same size)

	inline T & operator[](const long int i);			// Get the i'th element 
	inline const T & operator[](const long int i) const;// directly from the storage space

	//ASSUMES CORRECT INDEX DIMENSION OF 1,2 or 3 respectivly
	inline T & operator()(const int i);			// Get the (i,j)'th element
	inline const T & operator()(const int i) const;
	inline T & operator()(const int i,const int j);			// Get the (i,j)'th element
	inline const T & operator()(const int i,const int j) const;
	inline T & operator()(const int i,const int j,const int k);			// Get the (i,j,k)'th element
	inline const T & operator()(const int i,const int j,const int k) const;

	inline T & get(const int i,const int j);			// Get the (i,j)'th element
	inline const T & get(const int i,const int j) const;

	inline T & operator[](const int* i);			// Get the i'th element 
	inline const T & operator[](const int* i) const;// where i is a m-dimensional vector

	inline long int size() const;					// Get the size vector
	inline int size(int i) const;					// Get the size of i'th dimension
	
	// Resize (contents not preserved)
	void resize(int newN1); //to a vector
	void resize(int newN1, int newN2); //to a matrix
	void resize(int* newn, int newm); //to a general tensor
	void resize(const tensor &rhs); //to copy the size of another tensor
	// Resize and fill with a constant value
	void assign(int* newn, int newm, const T &a);
	~tensor();
};

// tensor definitions

template <class T>
void tensor<T>::print(std::ostream& out) {
	if(Z > 0) {
	    int a[m];

	    for(int i = 0; i < Z; i++) {
	    	if(m>2) {
		        if(i % cn[2] == 0) {
			        ind2vec(a,i);
			        out << "Matrix slice: (:,:,";
			        for(int j=2; j < m; j++) {
			        	out << a[j];
			        	if(j < m-1) {out << ",";}
			        	else {out << ")";}
			        }
			        out << " = ";
			        out << std::endl;
		        }
	    	}
	        out << v[i];
	        if((i+1) % n[0] != 0) {
	            out << " ";
	        }
	        else {
	        	out << std::endl;
	        }
	        
	    }	
	}
}


template <class T>
inline long int tensor<T>::vec2ind(const int* a) {
	long int IND=0;
	for(int i=0; i < m; i++) IND += a[i]*cn[i];
	return IND;
}

//Assumes a has been allocated to be int[m]
template <class T>
inline void tensor<T>::ind2vec(int *a,const long int I) {
	long int Itmp=I;
	for(int i=m-1; i >= 0; i--) {
		a[i] = Itmp/cn[i];
		Itmp -= a[i]*cn[i];
	}
}


template <class T>
tensor<T>::tensor() :  m(0), n(NULL), cn(NULL), Z(0), v(NULL) {}

template <class T>
tensor<T>::tensor(int* N, int M) : m(M), n(m>0 ? new int[m] : NULL), cn(m>0 ? new int[m] : NULL), Z(m>0 ? inline_contract(N,M) : 0) ,v(Z>0 ? new T[Z] : NULL) {
	n[0] = N[0];
	cn[0] = 1;
	for(int i=1; i<m; i++) {
		n[i] = N[i];
		cn[i] = N[i-1]*cn[i-1];
	}
}

template <class T>
tensor<T>::tensor(int N) : m(1), n(new int[m]), cn(new int[m]), Z(N) ,v(new T[Z]) {
	n[0] = N;
	cn[0] = 1;
}

template <class T>
tensor<T>::tensor(int N, int M) : m(2), n(new int[m]), cn(new int[m]), Z(N*M) ,v(new T[Z]) {
	n[0] = N;
	cn[0] = 1;
	
	n[1] = M;
	cn[1] = N;
}

template <class T>
tensor<T>::tensor(int N, int M, int L) : m(3), n(new int[m]), cn(new int[m]), Z(N*M*L) ,v(new T[Z]) {
	n[0] = N;
	cn[0] = 1;
	
	n[1] = M;
	cn[1] = N;
	
	n[2] = L;
	cn[2] = N*M;
}

template <class T>
tensor<T>::tensor(int N, const T &a) : m(1), n(new int[m]), cn(new int[m]), Z(N) ,v(new T[Z]) {
	n[0] = N;
	cn[0] = 1;
	for(long int i=0; i < Z; i++) v[i] = a;
}

template <class T>
tensor<T>::tensor(int N, int M, const T &a) : m(2), n(new int[m]), cn(new int[m]), Z(N*M) ,v(new T[Z]) {
	n[0] = N;
	cn[0] = 1;
	
	n[1] = M;
	cn[1] = N;
	for(long int i=0; i < Z; i++) v[i] = a;
}

template <class T>
tensor<T>::tensor(int N, int M, int L, const T &a) : m(3), n(new int[m]), cn(new int[m]), Z(N*M*L) ,v(new T[Z]) {
	n[0] = N;
	cn[0] = 1;
	
	n[1] = M;
	cn[1] = N;
	
	n[2] = L;
	cn[2] = N*M;
	for(long int i=0; i < Z; i++) v[i] = a;
}


template <class T>
tensor<T>::tensor(int* N, int M, const T &a) : m(M), n(m>0 ? new int[m] : NULL), cn(m>0 ? new int[m] : NULL), Z(m>0 ? inline_contract(N,M) : 0) ,v(Z>0 ? new T[Z] : NULL) {
	n[0] = N[0];
	cn[0] = 1;
	long int i;
	for(i=1; i<m; i++) {
		n[i] = N[i];
		cn[i] = N[i-1]*cn[i-1];
	}
	for(i=0; i<Z; i++) v[i] = a;
}

template <class T>
tensor<T>::tensor(const tensor &rhs) : m(rhs.m), n(m>0 ? new int[m] : NULL), cn(m>0 ? new int[m] : NULL), Z(rhs.Z>0 ? rhs.Z : 0) ,v(Z>0 ? new T[Z] : NULL) {
	long int i;
	for(i=0; i<m; i++) {
		n[i] = rhs.n[i];
		cn[i] = rhs.cn[i];
	}
	for(i=0; i<Z; i++) v[i] = rhs.v[i];
}

template <class T>
inline void tensor<T>::copy(const tensor &rhs) {
	//THIS FUNCTION DOES NOT RESCALE CONTAINERS
	#ifdef _CHECK_BOUNDS_
	if(m != rhs.m || Z != rhs.Z) {throw("Assignemt tensor: bounds does not match");}
	for(int i=0; i<m; i++) {
		if(n[i] != rhs.n[i]) {throw("Assignemt tensor: bounds does not match");}
	}
	#endif
	for(long int i=0; i<Z; i++) v[i] = rhs.v[i];
}

template <class T>
tensor<T> & tensor<T>::operator=(const tensor<T> &rhs) {
	
	if (this != &rhs) {
		long int i;
		if(m != rhs.m || Z != rhs.Z) {
			if (v != NULL) {delete [] (v);}
			if (n != NULL) {delete [] (n); delete [] (cn);}
			m = rhs.m;
			Z = rhs.Z;
			v= Z>0 ? new T[Z] : NULL;
			n= m>0 ? new int[m] : NULL;
			cn= m>0 ? new int[m] : NULL;
		}
		for(i=0; i<m; i++) {
			n[i] = rhs.n[i];
			cn[i] = rhs.cn[i];
		}
		for(i=0; i<Z; i++) v[i] = rhs.v[i];
	}
	return *this;
}


template <class T>
inline T & tensor<T>::operator[](const long int i) {
#ifdef _CHECK_BOUNDS_
	if(i >= Z) {throw("Tensor get element (direct): out of bounds");}
#endif
	return v[i];
}

template <class T>
inline const T & tensor<T>::operator[](const long int i) const {
#ifdef _CHECK_BOUNDS_
	if(i >= Z) {throw("Tensor get element (direct): out of bounds");}
#endif
	return v[i];
}

template <class T>
inline T & tensor<T>::operator()(const int i,const int j,const int k) {
#ifdef _CHECK_BOUNDS_
	if(i + cn[1]*j + cn[2]*k >= Z) {throw("Tensor get element (i,j,k): out of bounds");}
#endif
	return v[i + cn[1]*j + cn[2]*k];
}

template <class T>
inline const T & tensor<T>::operator()(const int i,const int j,const int k) const {
#ifdef _CHECK_BOUNDS_
	if(i + cn[1]*j + cn[2]*k >= Z) {throw("Tensor get element (i,j,k): out of bounds");}
#endif
	return v[i + cn[1]*j + cn[2]*k];
}

template <class T>
inline T & tensor<T>::operator()(const int i,const int j) {
#ifdef _CHECK_BOUNDS_
	if(i + cn[1]*j >= Z) {throw("Tensor get element (i,j): out of bounds");}
#endif
	return v[i + cn[1]*j];
}

template <class T>
inline const T & tensor<T>::operator()(const int i,const int j) const {
#ifdef _CHECK_BOUNDS_
	if(i + cn[1]*j >= Z) {throw("Tensor get element (i,j): out of bounds");}
#endif
	return v[i + cn[1]*j];
}

template <class T>
inline T & tensor<T>::get(const int i,const int j) {
#ifdef _CHECK_BOUNDS_
	if(j + cn[1]*i >= Z) {throw("Tensor get element (i,j): out of bounds");}
#endif
	return v[j + cn[1]*i];
}

template <class T>
inline const T & tensor<T>::get(const int i,const int j) const {
#ifdef _CHECK_BOUNDS_
	if(j + cn[1]*i >= Z) {throw("Tensor get element (i,j): out of bounds");}
#endif
	return v[j + cn[1]*i];
}

template <class T>
inline T & tensor<T>::operator()(const int i) {
#ifdef _CHECK_BOUNDS_
	if(i >= Z) {throw("Tensor get element (i): out of bounds");}
#endif
	return v[i];
}

template <class T>
inline const T & tensor<T>::operator()(const int i) const {
#ifdef _CHECK_BOUNDS_
	if(i >= Z) {throw("Tensor get element (i): out of bounds");}
#endif
	return v[i];
}

//ASSUMES CORRECT INDEX DIMENSION
template <class T>
inline T & tensor<T>::operator[](const int* i) {
	int j; 
	long int I = 0;
#ifdef _CHECK_BOUNDS_
	for(j=0; j < m; j++) {
		if(i[j] >= n[j]) {throw("Tensor get element (direct): out of bounds");}
	}
#endif
	for(j=0; j < m; j++) {
		I += i[j]*cn[j];
	}

	return v[I];
}

//ASSUMES CORRECT INDEX DIMENSION
template <class T>
inline const T & tensor<T>::operator[](const int* i) const {
	int j;
	long int I = 0;
#ifdef _CHECK_BOUNDS_
	for(int j=0; j < m; j++) {
		if(i[j] >= n[j]) {throw("Tensor get element (direct): out of bounds");}
	}
#endif
	for(j=0; j < m; j++) {
		I += i[j]*cn[j];
	}

	return v[I];
}


template <class T>
tensor<T> tensor<T>::row(const int i) {
	tensor<T> R(n[0]);

#ifdef _CHECK_BOUNDS_
	if(i >= n[1]) {throw("Trying to extract row out of bounds");}
#endif

	for(int k=0; k < n[0]; k++) {
		R[k] = v[k + cn[1]*i];
	}
	return R;
}
template <class T>
tensor<T> tensor<T>::col(const int i) {
	tensor<T> R(n[1]);

#ifdef _CHECK_BOUNDS_
	if(i >= n[0]) {throw("Trying to extract col out of bounds");}
#endif

	for(int k=0; k < n[1]; k++) {
		R(k) = v[i + cn[1]*k];
	}
	return R;
}


template <class T>
long int tensor<T>::size() const {
	return Z;
}

template <class T>
int tensor<T>::size(const int i) const {
#ifdef _CHECK_BOUNDS_
	if(i >= m) {throw("Tensor get size: out of bounds");}
#endif
	return n[i];
}

template <class T>
void tensor<T>::resize(int newN1) {
	if(newN1 != Z || m != 1) {
		delete [] (v);
		delete [] (n);
		delete [] (cn);
		m = newN1 > 0 ? 1 : 0;
		Z = newN1;

		v= Z>0 ? new T[Z] : NULL;
		n= m>0 ? new int[m] : NULL;
		cn= m>0 ? new int[m] : NULL;

		n[0] = newN1;
		cn[0] = 1;
	}
}

template <class T>
void tensor<T>::resize(int newN1, int newN2) {
	if(newN1 != n[0] || newN2 != n[1] || m != 2) {
		delete [] (v);
		delete [] (n);
		delete [] (cn);

		Z = newN1*newN2;
		m = Z > 0 ? 2 : 0;
		
		v= Z>0 ? new T[Z] : NULL;
		n= m>0 ? new int[m] : NULL;
		cn= m>0 ? new int[m] : NULL;

		n[0] = newN1;
		cn[0] = 1;
		n[0] = newN2;
		cn[0] = newN1;
	}
}



template <class T>
void tensor<T>::resize(const tensor &rhs) {
		delete [] (v);
		delete [] (n); 
		delete [] (cn);
		m = rhs.dim_N();
		Z = rhs.size();

		v= Z>0 ? new T[Z] : NULL;
		n= m>0 ? new int[m] : NULL;
		cn= m>0 ? new int[m] : NULL;

		if(m > 0) {
			n[0] = rhs.size(0);
			cn[0] = 1;
			for(int i=1; i<m; i++) {
				n[i] = rhs.size(i);
				cn[i] = rhs.size(i-1)*cn[i-1];
			}
		}
}

template <class T>
void tensor<T>::resize(int* newn, int newm) {
	if (newn != NULL && newm > 0) {
		delete [] (v);
		delete [] (n); 
		delete [] (cn);
		m = newm;
		Z = inline_contract(newn,newm);

		v= Z>0 ? new T[Z] : NULL;
		n= m>0 ? new int[m] : NULL;
		cn= m>0 ? new int[m] : NULL;

		n[0] = newn[0];
		cn[0] = 1;
		for(int i=1; i<m; i++) {
			n[i] = newn[i];
			cn[i] = newn[i-1]*cn[i-1];
		}
	}
}

template <class T>
void tensor<T>::assign(int* newn, int newm, const T &a) {
	if (newn != NULL && newm > 0) {
		delete [] (v);
		delete [] (n); 
		delete [] (cn);
		m = newm;
		Z = inline_contract(newn,newm);

		v= Z>0 ? new T[Z] : NULL;
		n= m>0 ? new int[m] : NULL;
		cn= m>0 ? new int[m] : NULL;

		n[0] = newn[0];
		cn[0] = 1;
		for(int i=1; i<m; i++) {
			n[i] = newn[i];
			cn[i] = newn[i-1]*cn[i-1];
		}

		for(long int i=0; i<Z; i++) v[i] = a;
	}
}


template <class T>
tensor<T>::~tensor()
{
	if (v != NULL) delete[] (v);
	if (n != NULL) delete[] (n);
	if (cn != NULL) delete[] (cn);
}



// ############################# BELOW ARE CONVERSION FROM vector<> TO tensor<>


template <class T>
void vector2tensor(tensor<T> &X,const std::vector<T> &Y) {
	int n0 = Y.size();
	X.resize(n0);
	for(int i=0; i < n0; i++) X[i] = Y[i];
}

template <class T>
void vvector2tensor(tensor<T> &X,const std::vector<std::vector<T> > &Y) {
	int n0 = Y.size();
	int n1 = Y[0].size();
	X.resize(n0,n1);
	for(int i=0; i < n0; i++) {
		for(int j=0; j < n1; j++) {
			X(j,i) = Y[i][j];
		}
	}
}

template <class T>
void tensor2vector(std::vector<T> &X, const tensor<T> &Y) {
	int n0 = Y.size(0);
	X.resize(n0);
	for(int i=0; i < n0; i++) X[i] = Y[i];
}

template <class T>
void tensor2vvector(std::vector<std::vector<T> > &X,const tensor<T> &Y) {
	int n0 = Y.size(0);
	int n1 = Y.size(1);
	X.resize(n0);
	for(int i=0; i < n0; i++) {
		X[i].resize(n1);
		for(int j=0; j < n1; j++) {
			X[i][j] = Y(j,i);
		}
	}
}


// ############################# BELOW ARE SORTING FUNCTIONS

template <class T>
tensor<T> sort_tensor_ret(const tensor<T> X) {
	tensor<T> Y(X);
	sort(Y.begin(),Y.end());
	return Y;
}

template <class T>
void sort_tensor(const tensor<T> &X) {
	sort(X.begin(),X.end());
}





#endif /* DATA_TYPE_HH_ */
