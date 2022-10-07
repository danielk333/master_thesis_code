#ifndef DATA_TYPE_HH_
#define DATA_TYPE_HH_

//#define _CHECK_BOUNDS_ 1

#include <cstddef>
#include <initializer_list>
#include <iterator>     // std::iterator_traits
#include <typeinfo>     // typeid

//Contraction of vector by multiplication
inline long int inline_contract(const int* a, const int n) {
	long int A(1);
	for(int i = 0; i < n; ++i) {A *= a[i];}
	return A;
}

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
	//Iterator
  	T* begin() { return Z>0 ? &v[0] : NULL; }
  	T* end() { return Z>0 ? &v[Z-1] : NULL; }
  	const int dim_C(int D) { return cn[D]; }
  	const int dim_N() { return m; }

	inline T back() { return Z>0 ? v[Z] : NULL; }
	tensor();										// Empty tensor without size
	explicit tensor(int* N, int M);					// Empty tensor with size
	tensor(int* N, int M, const T &a);				// Initialize to constant value
	tensor(const tensor &rhs);						// Copy constructor
	typedef T value_type; 							// Make the type available externally
	tensor & operator=(const tensor &rhs);			// Assignment operator
	//tensor & slice(const int i, const int D);		// Take the i'th slice in dimension D
	//tensor & vector(const int i, const int P, const int D);		// Take the i'th vector in dimension P along dimension D
	
	tensor & row(const int i); //If m = 2, get row i
	tensor & col(const int i); //If m = 2, get col i

	inline void copy(const tensor &rhs); 			// Copy data (same size)

	inline T & operator[](const long int i);			// Get the i'th element 
	inline const T & operator[](const long int i) const;// directly from the storage space

	//ASSUMES CORRECT INDEX DIMENSION

	inline T & operator()(const int i);			// Get the (i,j)'th element
	inline const T & operator()(const int i) const;

	inline T & operator()(const int i,const int j);			// Get the (i,j)'th element
	inline const T & operator()(const int i,const int j) const;

	inline T & operator()(const int i,const int j,const int k);			// Get the (i,j,k)'th element
	inline const T & operator()(const int i,const int j,const int k) const;

	inline T & operator[](const int* i);			// Get the i'th element 
	inline const T & operator[](const int* i) const;// where i is a m-dimensional vector
//	inline T & operator[](std::initializer_list<int> i);			// Get the i'th element
//	inline const T & operator[](std::initializer_list<int> i) const;// where i is a m-dimensional vector

	inline long int size() const;						// Get the size vector
	inline int size(int i) const;					// Get the size of i'th dimension
	void resize(int* newn, int newm); 				// Resize (contents not preserved)
	void assign(int* newn, int newm, const T &a); 	// Resize and fill with a constant value
	~tensor();
};

//Itterator
/*
template <class T>
struct tensor_iterator {
	typedef T* iterator;
  	typedef const T* const_iterator;
    iterator I;
    int *cn; // Cumulative size of each index
    tensor_iterator(): I(NULL), cn(NULL) {}
    tensor_iterator(tensor<T> X): I(X.begin()), cn(X.dim_N()>0 ? new int[X.dim_N()] : NULL) {
    	for(int i=0; i < X.dim_N(); i++) {
    		cn[i] = X.dim_C(i);
    	}
    }

    tensor_iterator & operator=(iterator Q) { I=Q; return *this; }
    tensor_iterator & operator=(tensor<T> X) { 
    	I=X.begin(); 
    	for(int i=0; i < X.m; i++) {
    		cn[i] = X.cn[i];
    	}
    	return *this; 
    }
    const tensor_iterator& operator++() { ++I; return *this; }
    const tensor_iterator& operator++(int) { I++; return *this; }
    void inc(int D) { I+=cn[D]; }

    const bool operator<(const_iterator Q) { return I<Q; }
    const bool operator>(const_iterator Q) { return I>Q; }
    const bool operator==(const_iterator Q) { return I==Q; }

    T& operator*() { return *I; }
};*/


// tensor definitions

template <class T>
inline long int tensor<T>::vec2ind(const int* a) {
	long int IND=0;
	for(int i=0; i < m; i++) IND += a[i]*cn[i];
	return IND;
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
	for(i=1; i<m; i++) {
		n[i] = rhs.n[i];
		cn[i] = rhs.cn[i];
	}
	for(i=0; i<Z; i++) v[i] = rhs.v[i];
}

template <class T>
inline void tensor<T>::copy(const tensor &rhs) {
	//THIS FUNCTION ASSUMES ALL DIMENSIONS ARE EQUAL
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
	long int I;
	I = i + cn[1]*j + cn[2]*k;
#ifdef _CHECK_BOUNDS_
	if(I >= Z) {throw("Tensor get element (i,j,k): out of bounds");}
#endif
	return v[I];
}

template <class T>
inline const T & tensor<T>::operator()(const int i,const int j,const int k) const {
	long int I;
	I = i + cn[1]*j + cn[2]*k;
#ifdef _CHECK_BOUNDS_
	if(I >= Z) {throw("Tensor get element (i,j,k): out of bounds");}
#endif
	return v[I];
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
tensor<T> & tensor<T>::row(const int i) {
	int D[] = {n[0]};
	tensor<T> R(D,1);

#ifdef _CHECK_BOUNDS_
	if(i >= n[1]) {throw("Trying to extract row out of bounds");}
#endif

	for(int k=0; k < n[0]; k++) {
		R(k) = v[k + cn[1]*i];
	}
	return R;
}
template <class T>
tensor<T> & tensor<T>::col(const int i) {
	int D[] = {n[1]};
	tensor<T> R(D,1);

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


/*
//ASSUMES CORRECT INDEX DIMENSION
template <class T>
inline T & tensor<T>::operator[](std::initializer_list<int> i) {
	int j, I = 0;
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
inline const T & tensor<T>::operator[](std::initializer_list<int> i) const {
	int j, I = 0;
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

*/
template <class T>
tensor<T>::~tensor()
{
	if (v != NULL) delete[] (v);
	if (n != NULL) delete[] (n);
	if (cn != NULL) delete[] (cn);
}


#endif /* DATA_TYPE_HH_ */
