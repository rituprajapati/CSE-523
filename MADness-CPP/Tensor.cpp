
class Vector {

	vector<ldouble> a;

public:
	Vector () {

		/* Do nothing */
	}

	/*********************************************/
	/* if arg is an int							 */
	/*********************************************/
	Vector (int arg ) {

		for ( int i = 0; i < arg ; i++ ){
			a.push_back( 0.0 );
		}
	}

	/*********************************************/
	/* if arg is an array						 */
	/*********************************************/
	Vector (vector<ldouble> arg ) {

		for ( int i = 0; i < arg.size(); i++ ){
			a.push_back( arg[i] );
		}
	}

	/*********************************************/
	/* copy constructor							 */
	/*********************************************/
	Vector ( const Vector& obj ){

		int n = obj.a.size();

		for ( int i = 0; i < n; i++ )
			a.push_back( obj.a[i] );

	}

	/*********************************************/
	/* new Vector from a slice of obj		 	 */
	/*********************************************/
	Vector ( const Vector& obj, int lo, int hi ){

		for ( int i = lo; i < hi; i++ ){
			this->a.push_back( obj.a[i] );
		}

	}
	/*********************************************/
	/* combining two Vectors 					 */
	/*********************************************/
	Vector ( Vector &obj1, Vector &obj2 ) {

		int k = 0;
		for ( int i = 0; i < obj1.len(); i++ )
			a.push_back( obj1.a[i] );

		for ( int i = 0; i < obj2.len(); i++ )
			a.push_back(obj2.a[i]);
	}

	/*********************************************/
	/* sum of square of all the values in a		 */
	/*********************************************/
	ldouble normf() {

		ldouble sum = 0.0;

		for (int i = 0; i < a.size(); i++ ){
			sum += ( a[i] * a[i]);
		}
		return sqrt(sum);
	}

	/*********************************************/
	/* sum of two vector object a				 */
	/*********************************************/
	ldouble inner( Vector& other ){

		ldouble sum = 0.0;

		for (int i = 0; i < a.size(); i++ ) {
			sum += a[i] * other.a[i];
		}

		return sum;
	}

	/*********************************************/
	/* update a[] acc to alpha and beta			 */
	/*********************************************/
	Vector gaxpy( ldouble alpha, Vector& other, ldouble beta) {
		
		for( int i = 0; i < this->a.size(); i++) {
			this->a[i] = alpha * this->a[i] + beta * other.a[i];
		}
		return *this;
	}

	/*********************************************/
	/* multiply each value of a with const s	 */
	/*********************************************/
	Vector scale( ldouble s ) {

		for (int i = 0; i < a.size(); i++ )
			a[i] *= s;

		return *this;
	}

	/*********************************************/
	/* if arg is an array						 */
	/*********************************************/
	Vector emul( Vector &other) {
		for (int i = 0 ; i < a.size(); i++ )
			a[i] *= other.a[i];
		return *this;
	} 

	/*********************************************/
	/* Vector v[ind]							 */
	/*********************************************/
	ldouble& operator[] ( int ind ) {
		return a[ind];
	}

	/*********************************************/
	/* return a Vector object defined over 		 */
	/* a[lo to hi]					 			 */
	/*********************************************/
	vector<ldouble> getSlice (int lo, int hi) {
		
		vector<ldouble>  tmp( a.begin()+lo, a.begin()+hi);
		return tmp;
	}

	/*********************************************/
	/* set a[lo to hi] to value					 */
	/*********************************************/
	void setSlice ( int lo, int hi, ldouble value ) {
		
		for(int i = lo; i < hi; i++ )
			a[i] = value;
	}

	/*********************************************/
	/* Print the object by a 					 */
	/*********************************************/
	void toStr() {
		
		for (int i = 0; i < a.size(); i++) {
			cout << a[i] << "  ";
		}
	}
   
	/*********************************************/
	/* obj.len() will return length of a 		 */
	/*********************************************/
    int len(){
        
        return this->a.size();
    }

    /*********************************************/
	/* obj3 = obj1 + obj2 				 		 */
	/*********************************************/
    Vector operator+ ( Vector other ){

    	Vector r( *this );
    	r.gaxpy (1.0, other, 1.0 );
    	return r;
    }

};



class Matrix {

	vector< Vector > b;

public:

	Matrix () {
		/* do nothing */
	}

	/***************************************************/
	/* Generate a 0 valued matrix given the dimensions */
	/***************************************************/
	Matrix ( int n, int m ){

		for (int i = 0; i < n; i++ ){ 
			Vector tmp(m);
			b.push_back ( tmp );
		}
	}

	/***************************************************/
	/* Copy constructor 							   */
	/***************************************************/	
	Matrix ( Matrix &obj ){

		for ( int i = 0; i < obj.len() ; i ++){
			Vector tmp = obj.b[i];
			b.push_back ( tmp );
		}
	}
	/***************************************************/
	/* Get number of rows of the matrix 			   */
	/***************************************************/
	int len() {

		return b.size();
	}

	/***************************************************/
	/* [] Overloading 								   */
	/***************************************************/
	Vector& operator[] ( int idx ){

		return b[idx];
	}

	/***************************************************/
	/* Multiply a vector and matrix 				   */
	/***************************************************/

	friend Vector operator* ( Matrix M, Vector V );

	friend Vector operator* ( Vector V, Matrix M );

};


/***************************************************/
/* Matrix * Vector 				 				   */
/***************************************************/
Vector operator* ( Matrix M, Vector V ) {

	int n, m = 0;
 	n = M.len();
 	if( n!= 0 )
 		m = M[0].len();

 	if( m != V.len() )
 		assert ( "Invalid dimensions\n" );

 	Vector r (n);

 	for ( int i = 0; i < n; i++ ) {
 		for ( int j = 0; j < m; j++ ){
 			r[i] += M[i][j] * V[j];
 		}
 	}

    return r;
}


/***************************************************/
/* Vector * Matrix 				 				   */
/***************************************************/
Vector operator* ( Vector V, Matrix M ) {
	
	int n, m = 0;
 	n = M.len();
 	if( n!= 0 )
 		m = M[0].len();

 	if( n != V.len() )
 		assert ( "Invalid dimensions\n" );

 	Vector r ( m );

 	for ( int i = 0; i < m; i++ ) {
 		for ( int j = 0; j < n; j++ ){
 			r[i] += V[j] * M[j][i];
 		}
 	}
    return r;
}




