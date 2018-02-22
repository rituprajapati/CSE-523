
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
// #include "Quad.cpp"
using namespace std;
#define ldouble long double

class CVector {

	public:
	/***************************************************/
	/* return the string representation of a vactor	   */
	/***************************************************/
	string toStr( vector< ldouble > tmp ) {

		string str = "";
		for (int i = 0; i < tmp.size(); i++) 
			str += to_string(tmp[i]);
		return str;
	}

	/***************************************************/
	/* change the value of vec from lo to hi to val	   */
	/***************************************************/
	void setSlice ( vector< ldouble > &vec, int lo, int hi, ldouble value ) {
		for(int i = lo; i < hi; i++ )
			vec[i] = value;
	}


   ldouble normf( vector< ldouble > a ) {

		ldouble sum = 0.0;

		for (int i = 0; i < a.size(); i++ )
			sum += a[i] * a[i];

		return sqrt(sum);
	}


	/*********************************************/
	/* sum of two vector object a				 */
	/*********************************************/
	ldouble inner( vector< ldouble > a, vector< ldouble > b ){
		ldouble sum = 0.0;

		for (int i = 0; i < a.size(); i++ ) {
			sum += a[i] * b[i];
		}

		return sum;
	}

	/*********************************************/
	/* if arg is an array						 */
	/*********************************************/
	vector< ldouble > emul( vector< ldouble > a1, vector< ldouble > a2 ) {

		for (int i = 0 ; i < a1.size(); i++ )
			a1[i] *= a2[i];
		return a1;
	} 

	/*********************************************/
	/* multiply each value of a with const s	 */
	/*********************************************/
	vector< ldouble > scale( vector< ldouble > a, ldouble s ) {

		for (int i = 0; i < a.size(); i++ )
			a[i] *= s;

		return a;
	}


	/*********************************************/
	/* update a[] acc to alpha and beta			 */
	/*********************************************/
	void gaxpy( vector< ldouble > &a1, vector< ldouble >a2, ldouble alpha, ldouble beta) {
		
		for( int i = 0; i < a1.size(); i++) {
			a1[i] = alpha * a1[i] + beta * a2[i];
		}
	}
};






class CMatrix {

	public:
	/***************************************************/
	/* Generate a 0 valued matrix given the dimensions */
	/***************************************************/
	vector<vector< ldouble > > Matrix ( int n, int m ){
		vector<vector< ldouble > > tmp(n);

		for (int i = 0; i < n; i++ ) 
			tmp.push_back(vector< ldouble > (m, 0));

		return tmp;
	}

	/***************************************************/
	/* Multiply a vector and matrix 				   */
	/***************************************************/
	vector< ldouble > mul ( vector< vector < ldouble > > M,  vector< ldouble > V ){

	 	int n = M.size();
	 	int m = M[0].size();

	 	// assert( m != V.size() ) 
	 	vector < ldouble > r(n, 0.0);
	 	

	 	for ( int i = 0; i < n; i++ ) {
	 		for ( int j = 0; j < m; j++ ){
	 			r[i] += M[i][j] * V[j];
	 		}
	 	}

        return r;
    }

};








class Function {

	/***************************************************/
	/* Private data								       */
	/***************************************************/
	int autorefine = 1;	
	int k;
	int thresh;
	int f;
	int max_level;
	unordered_map<int, unordered_map<int, vector< ldouble > > > d;
	unordered_map<int, unordered_map<int, vector< ldouble > > > s;
	int rm, r0, rp;
	int compressed;
	vector< vector< ldouble > > hg, hg0, hg1, hgT;
	vector< ldouble > quad_w, quad_x ;
	vector<vector< ldouble >> quad_phi, quad_phiT, quad_phiw;
	int quad_npt;

	CVector V_op;
	CMatrix M_op;
	Gauss G;

	/***************************************************/
	/* Constructor function for initializing data      */
	/***************************************************/
	Function( int k, int thr, int f = 0, int initial_level=2) {

		k = k;
		thresh = thr;
		f = f;
        max_level = 30;

        for (int i = 0; i <= max_level; i++ ) {
        	d[i] = {};
        	s[i] = {};
        }

        init_twoscale(k);
        init_quadrature(k);

        vector<int> m = make_dc_periodic(k);

        rm = m[0];
        r0 = m[1];
        rp = m[2];

        compressed = 0;

        if(f) {

        	int max_level = pow ( 2, initial_level );
        	for (int i = 0; i < max_level; i++) 
				refine(initial_level, i);
        }

	}

	/***************************************************/
	/* Convert the s[n][l] vector to string and print  */
	/***************************************************/
	void print_tree( int n = 0, int l = 0){

		if ( s[n].find(l) != s[n].end() ){
			string str = "";

			for (int i = 0; i < n; i++ )
				str = str + "   ";

			cout << str << "[" << n << "," << l << "] Leaf Node with Coefficients ";
			cout << str << V_op.toStr( s[n][l] );
			cout << " ";
		}
		else {
			print_tree( n+1, 2*l );
			print_tree( n+1, 2*l+1 );
		}
	}

	/***************************************************/
	/* Deep copy from one object to another			   */
	/***************************************************/
	Function copy(){
       
		Function result = Function(this->k, this->thresh);
		result.compressed = this->compressed;
		result.f = this->f;

		result.s = this->s;
		result.d = this->d;

		return result;

    }

    /***************************************************/
	/* Initialize the class Matrices				   */
	/***************************************************/
	void init_twoscale ( int k ) {

		vector<vector<double>> hg_ = twoscalecoeffs( k );

		hg 	= M_op.Matrix(2*k, 2*k);
		hg0 = M_op.Matrix(2*k, k);
		hg1 = M_op.Matrix(2*k, k);
		hgT = M_op.Matrix(2*k, 2*k);

		for (int i = 0; i < 2*k ; i++ ) {
			for (int j = 0; j < 2*k; j++ ) {
				hg[i][j]  = hg_[i][j];
				hgT[i][j] = hg_[j][i];
			}
		}

		for (int i = 0; i < 2*k; i++ ){
			for (int j = 0; j < k; j++ ) {
				hg0[i][j]  = hg_[i][j];
				hg1[i][j]  = hg_[i][j+k];
			}
		}
	}
	 
	/***************************************************/
	/* Initialize class matrices 					   */
	/***************************************************/
	void init_quadrature( int order) {

		int npt;
		G.gauss_legendre(order, this->quad_x, this->quad_w);
		this->quad_npt = npt = quad_w.size();

		quad_phi  = M_op.Matrix(npt, k);
		quad_phiT = M_op.Matrix( k, npt);
		quad_phiw = M_op.Matrix(npt, k);

		for ( int i = 0; i < npt; i++ ) {
			
			vector< ldouble > p = phi ( quad_x[i], k);

			for ( int m = 0; m < k; m++ ) {
				quad_phi[i][m] = p[m];
				quad_phiT[m][i] = p[m];
				quad_phiw[i][m] = quad_w[i] * p[m];
			}
		}
	}

	/***************************************************/
	/*  											   */
	/***************************************************/
	vector< ldouble > project ( int n, int l ){
        
		vector< ldouble > s( k, 0.0 );
		ldouble h = pow( 0.5, n );
		ldouble scale = sqrt( h );

		for ( int mu = 0; mu < quad_npt; mu++ ) { 

			ldouble x = ( l + quad_x[ mu ] ) * h ;
			ldouble f = f(x);

			for ( int i = 0; i < k; i++ )
				s[i] += scale * f * quad_phiw[mu][i];
		}

		return s;
	}


	/***************************************************/
	/*   refine numerical representation of f(x) 	   */
	/*   to desired tolerance				 		   */
	/***************************************************/
	void refine(int n, int l) {

		vector< ldouble > s0, s1, s, d;
		int k = this->k;

		s0 = project(n + 1, 2 * l);
		s1 = project(n + 1, 2 * l + 1);
		
		/* size of s would be 2*k */
		s.insert ( s0.end(), s1.begin(), s1.end() );
		
		d = M_op.mul( this->hgT, s );

		if ( V_op.normf( vector< ldouble >( d.begin()+k, d.end() ) ) < thresh 
			 || n >= max_level - 1 ) {
			this->s[n + 1][2 * l] = s0;
			this->s[n + 1][2 * l + 1] = s1;
		}
		else {

			refine(n + 1, 2 * l);
			refine(n + 1, 2 * l + 1);
		}
	}


    ldouble evaluate( int n, int l, ldouble x){
        
		vector< ldouble > p ;

		if ( s[n].find(l) != s[n].end() ) {
			p = phi ( x, this->k );
			return  ( V_op.inner( s[n][l], p ) * sqrt( pow ( 2.0, n)) );
		}

		else {
			n = n + 1;
			l = 2 * l;
			x = 2 * x;
			if ( x >= 1 ) {
				l = l+1;
				x = x-1;
			}
			return evaluate(n, l, x);
		}
	}


// def __call__(self, x):
// '''
// Evaluate function at x ... scaling function basis only
// call to self after creation
// looks like a Function evaluation
// say g = Function(5, 1e-3, f) so g(1.0) should be similar to f(1.0)
// '''
// if self.compressed: self.reconstruct()
// return self.__evaluate(0, 0, x)


	ldouble norm2 ( ){
		 
		if ( this->compressed )
			this->reconstruct();

		ldouble sum = 0.0;

		for ( auto itr = s.begin(); itr != s.end(); itr++ ){
			for ( auto itr1 = (itr->second).begin(); itr1 != (itr->second).end(); itr1++ )
				sum += pow ( V_op.normf ( itr1->second ), 2 );
		}

		return sqrt(sum);
	}

	void compress( int n = 0, int l = 0 ) {
		
		if ( this->compressed ) 
			return;

		if ( s[n+1].find( 2*l ) == (s[n+1].end()) )
			compress(n + 1, 2 * l);

		if ( s[n+1].find( 2*l + 1 ) == s[n+1].end() )
			compress(n + 1, 2 * l + 1);

		vector < ldouble > s, d;
		s.insert ( this->s[n+1][2*l].end(), this->s[n + 1][2 * l + 1].begin(), this->s[n + 1][2 * l + 1].end());

		d = M_op.mul ( hgT, s );

		this->s[n][l] = vector< ldouble >(d.begin(), d.begin()+k );
		this->d[n][l] = vector< ldouble >(d.begin() + k, d.end() );

		this->s[n+1].erase( 2*l );
		this->s[n+1].erase( 2*l + 1);

		if ( n==0 )
			this->compressed = 1;
	}

	void reconstruct( int n=0, int l=0) {
		
		if ( ! this->compressed )
			return;

		if ( this->d[n].find(l) != this->d[n].end() ) {

			vector< ldouble > d, s;
			d.insert ( this->s[n][l].end(), this->d[n][l].begin(), this->d[n][l].end() );
			this->d.erase(l);
			this->s.erase(l);

			/* apply the two scale relationship to get difference coeff */
            /* in 1d this is O(k^2) flops (in 3d this is O(k^4) flops). */
			s = M_op.mul ( this->hg, d );

			this->s[n + 1][2 * l] = vector< ldouble > (s.begin(), s.begin()+k );
			this->s[n + 1][2 * l + 1] = vector< ldouble > (s.begin()+k, s.end() );

			/* sub-trees can be done in parallel */
			reconstruct(n + 1, 2 * l);
			reconstruct(n + 1, 2 * l + 1);
		}

		if ( n == 0 )
			this->compressed = 0;
	}


	void mul_iter( Function f1, Function f2, int n=0, int l=0) {
		
		if ( f1.s[n].find(l) != f1.s[n].end() && f2.s[n].find(l) != f2.s[n].end() ) {
			if ( autorefine && n+1 <= this->max_level ) {   /* ?? Function.autorefine */
			
				/* refine both one more level */
				f1.recur_down( n, l, f1.s[n][l] );
				f2.recur_down( n, l, f2.s[n][l] );

				/* scale factor for this level = sqrt((2^d)^(n+1)) ?? */
				ldouble scale_factor = sqrt( pow ( 2.0, n+1) );

				/* multiply f1.s[n+1][2*l] and f2.s[n+1][2*l] */
				vector < ldouble > f, g;

				f = M_op.mul ( this->quad_phiT,  f1.s[n+1][2*l]);
				g = M_op.mul ( this->quad_phiT,  f2.s[n+1][2*l]);

				f = V_op.emul ( f, g);
				this->s[n+1][2*l] = V_op.scale( M_op.mul ( this->quad_phiw, f ), scale_factor);

				/* multiply f1.s[n+1][2*l+1] and f2.s[n+1][2*l+1] */
				f = M_op.mul ( this->quad_phiT, f1.s[n+1][2*l+1]);
				g = M_op.mul ( this->quad_phiT, f2.s[n+1][2*l+1]); 
				f = V_op.emul( f, g );
				this->s[n+1][2*l+1] = V_op.scale( M_op.mul ( this->quad_phiw, f ), scale_factor);
			}

			else {

				/* if autorefine is not set or we are at the max_level */
				/* live with what you get */
				vector < ldouble > f, g;
				f = M_op.mul ( this->quad_phiT, f1.s[n][l] ); 
				g = M_op.mul ( this->quad_phiT, f2.s[n][l] );
				f = V_op.emul(f, g );

				/* scale factor for this level = sqrt((2^d)^(n+1)) */
				this->s[n][l] = V_op.scale( M_op.mul ( this->quad_phiw, f ), sqrt( pow (2.0, n)));
			}
		}
		else {
			if ( f1.s[n].find(l) != f1.s[n].end() && f2.s[n].find(l) == f2.s[n].end()) {
				/* refine this box down to next level in f1 */
				f1.recur_down(n, l, f1.s[n][l]);
			}
			else if ( f1.s[n].find(l) == f1.s[n].end() && f2.s[n].find(l) != f2.s[n].end() ) {
				/* refine this box down to next level in f2 */
				f2.recur_down(n, l, f2.s[n][l]);

			}
				
			/* calls on sub-trees can go in parallel */
			this->mul_iter(f1, f2, n+1, 2*l);
			this->mul_iter(f1, f2, n+1, 2*l+1);
		}
	}

	/***************************************************/
	/* multiplication 								   */
    /* For multiply both operands need to be in the    */
    /* scaling function basis so possibly call 		   */
    /* reconstruct on one or both of them first  	   */
	/* to desired tolerance				 		       */
	/***************************************************/

	Function mul(  Function other ){

		if ( this->compressed )
			this->reconstruct();

		if ( this->compressed )
			other.reconstruct();

		Function result = Function(this->k, this->thresh);
		result.mul_iter( *this, other );

		this->sclean();
		other.sclean();

		return result;
	}

	/***************************************************/
	/*   * Overloading  					       	   */
	/***************************************************/
	Function operator* ( Function other ) {
		return this->mul( other );
	}


	/***************************************************/
	/*   recursive "iteration" for gaxpy.       	   */
	/***************************************************/
	void gaxpy_iter ( ldouble alpha, Function other, ldouble beta, int n=0, int l=0 ) {
	
		if ( this->d[n].find(l) != this->d[n].end() || other.d[n].find(l) != other.d[n].end() ){

			if ( this->d[n].find(l) != this->d[n].end() && other.d[n].find(l) != other.d[n].end() )
				V_op.gaxpy( this->d[n][l], other.d[n][l], alpha, beta);
			
			else if ( this->d[n].find(l) == this->d[n].end() && other.d[n].find(l) != other.d[n].end() )
				this->d[n][l] = V_op.scale ( other.d[n][l], beta);

			else if ( this->d[n].find(l) != this->d[n].end() && other.d[n].find(l) == other.d[n].end() )
				this->d[n][l] = V_op.scale( this->d[n][l], alpha);

			/* calls on sub-trees can go in parallel */
			gaxpy_iter(alpha, other, beta, n+1, 2*l );
			gaxpy_iter(alpha, other, beta, n+1, 2*l+1 );
		}
	}

	/******************************************************/
	/*  only in multi-wavelet basis					      */
	/*	i.e. run compress first on one or both operands.  */ 
	/*	self = alpha*self + beta*other                    */
	/*	Other is not changed.      	   					  */
	/******************************************************/
	Function gaxpy( ldouble alpha, Function other, ldouble beta) {

		if ( !this->compressed )
			this->compress();
		if ( !other.compressed )
			other.compress();

		V_op.gaxpy ( this->s[0][0], other.s[0][0], alpha, beta );
		this->gaxpy_iter(alpha, other, beta);

		/* return self so operations can be chained */
		return *this;
	}


	/***************************************************/
	/*   basic addition 					       	   */
	/***************************************************/
	Function add ( Function other ){

		return this->copy().gaxpy( 1.0, other, 1.0 );
	}	


	/***************************************************/
	/*   + overloading 	 					       	   */
	/***************************************************/
	Function operator+ ( Function other ) {
		return this->add ( other );
	}


	/***************************************************/
	/*   basic subtraction 					       	   */
	/***************************************************/

	Function sub ( Function other ){

		return this->copy().gaxpy(1.0, other, -1.0);
	}

	/***************************************************/
	/*   - overloading for Function  	 			   */
	/***************************************************/	
	Function operator- ( Function other ){

		return this->sub(other);
	}



	/***************************************************/
	/* Mostly for debugging, print summary of  		   */
	/* coefficients, optionally printing the norm of   */
	/* each block 	 			  					   */
	/***************************************************/
	void summarize( int printcoeff = 0 ) {

		cout << "sum coefficients";

		for (auto itr = this->s.begin(); itr != this->s.end(); itr++) {

			ldouble sum = 0.0;
			int n = itr->first;

			for ( auto itr2 = this->s[itr->first].begin(); itr2 != this->s[itr->first].end(); itr2++) {
				int l = itr2->first;
				if ( printcoeff )
					// cout << "%3d %6d %.2e" << (n, l, self.s[n][l].normf())
					cout << n << l << V_op.normf( this->s[n][l]);
				else
					sum +=  pow( V_op.normf( this->s[n][l]), 2);
			}
			if ( !printcoeff ) {
				if ( this->s[n].size() != 0 ) 
					cout << "level=" << n << "boxes=" << this->s[n].size() << "norm=" << sqrt(sum);
			}
		}

		cout << "difference coefficients";
		for (auto itr = this->d.begin(); itr != this->d.end(); itr++) {
			int n = itr->first;
			ldouble sum = 0.0;

			for ( auto itr2 = this->d[itr->first].begin(); itr2 != this->d[itr->first].end(); itr2++) {
				int l = itr2->first;

				if ( printcoeff )
					// cout << "%3d %6d %.2e" << (n, l, self.s[n][l].normf())
					cout << n << l << V_op.normf( this->d[n][l]);
				else
					sum +=  pow( V_op.normf( this->d[n][l]), 2);
			}

			if ( !printcoeff ) {
				if ( this->d[n].size() != 0 ) 
					cout << "level=" << n << "boxes=" << this->d[n].size() << "norm=" << sqrt(sum);
			}
		}
	}




	/*************************************************************/
	/* In s are scaling coefficients for box n,l ... apply		 */ 
	/* twoscale to generate the corresponding coefficients on    */
	/* level n+1 and insert the results into the tree of scaling */ 
	/* function coefficients.  					       	 	     */
	/*************************************************************/

	void recur_down ( int n, int l, vector< ldouble > s ) {
		
		int k = this->k;
		vector< ldouble > d (2*k, 0.0);
		for (int i = 0; i < k; i++ )
			d[i] = s[i];

		s = M_op.mul( this->hg, d );
		this->s[n + 1][2 * l] = vector<ldouble> (s.begin(), s.begin()+k);
		this->s[n + 1][2 * l + 1] = vector<ldouble> (s.begin()+k, s.end());
	}



	/***************************************************/
	/* get_coeffs 				    	 			   */
	/***************************************************/	

	vector< ldouble > get_coeffs ( int n, int l ) {

		vector< ldouble > s;
		if ( l < 0 || l >= pow( 2, n) ) 
			return vector< ldouble > (this->k, 0.0 );

		if ( this->s[n].find(l) != this->s[n].end() ) 
			return this->s[n][l];

		if ( n > 0 ){
			s = this->get_coeffs(n-1, l/2);
			if ( s.size() == 0)
				return {};
		}
		else 
			return {};

		this->recur_down(n-1, l/2, s);
		return this->s[n][l];
	}



	/***************************************************/
	/*  sclean 					    	 			   */
	/***************************************************/
	void sclean ( int n = 0, int l = 0, int cleaning = 0 ) {

		if ( cleaning )
			this->s[n].erase(l);
		else
			cleaning = ( this->s[n].find(l) != this->s[n].end() ) ? 1: 0;

		/* Sub trees can run in parallel */

		if ( n < this->max_level ) {

			if ( !cleaning || this->s[n + 1].find(2 * l) != this->s[n+1].end() )
				this->sclean(n + 1, 2 * l, cleaning);

			if ( !cleaning || this->s[n + 1].find(2 * l + 1) != this->s[n+1].end() )
				this->sclean(n + 1, 2 * l + 1, cleaning);
		}
	}


	

};























