
#include "header.h"

/**************************************************************************/
/**************************************************************************/
/* CLASS - FUNCTION 													  */
/* Multiresolution representation of 1-d functions using a multiwavelet   */
/* basis (similar to discontinuous spectral element with a hierarchal 	  */
/* decomposition). 														  */
/**************************************************************************/
/**************************************************************************/

class Function {

public:
	/***************************************************/
	/* class data								       */
	/***************************************************/
	int k;
	ldouble thresh;
	ldouble ( *f ) ( ldouble) ;
	int max_level;
	unordered_map<int, unordered_map< int, Vector > > d;
	unordered_map<int, unordered_map< int, Vector > > s;
	Matrix rm, r0, rp;
	int compressed;

	Matrix quad_phi, quad_phiT, quad_phiw;
	Matrix hg, hg0, hg1, hgT;

	Vector quad_w, quad_x ;
	int quad_npt;

	Gauss G;

	Function () {
		/* do_nothing */
	}

	/***************************************************/
	/* Constructor function for initializing data      */
	/***************************************************/
	Function( int k_, double thr, ldouble( *func_ptr )( ldouble) = NULL, int initial_level = 2) {

		this->k = k_;
		this->thresh = thr;
		this->f = func_ptr;
        this->max_level = 30;

        for (int i = 0; i <= max_level; i++ ) {
        	this->d[i] = {};
        	this->s[i] = {};
        }


        this->init_twoscale(k);
        this->init_quadrature(k);


        vector< Matrix > m = make_dc_periodic(this->k);

        if ( m.size() != 3 )
        	assert ( "Vector size issue\n" );
        this->rm = m[0];
        this->r0 = m[1];
        this->rp = m[2];

       
        this->compressed = 0;

        if(this->f != NULL) {
        	int max_level_ = pow ( 2, initial_level );

        	for (int i = 0; i < max_level_; i++) 
				refine(initial_level, i);
        }
	}

	/***************************************************/
	/* Copy constructor 						       */
	/***************************************************/
	Function ( const Function &other ){

		k = other.k;
		thresh = other.thresh;
		f = other.f;
		max_level = other.max_level;
		d = other.d;
		s = other.s;
		rm = other.rm;
		r0 = other.rm;
		rp = other.rp;;
		compressed = other.compressed;

		quad_phi = other.quad_phi;
		quad_phiT = other.quad_phiT;
		quad_phiw = other.quad_phiw;

		hg  = other.hg;
		hg0 = other.hg0;
		hg1 = other.hg1;
		hgT = other.hgT;

		quad_w = other.quad_w;
		quad_x = other.quad_x;
		quad_npt = other. quad_npt;
	}

	/***************************************************/
	/* Convert the s[n][l] Vector to string and print  */
	/***************************************************/
	void print_tree( int n = 0, int l = 0){

		if ( s.find(n) == s.end() ){
			cout << " print_tree -- 364 \n";
			return;
		}

		if ( s[n].find(l) != s[n].end() ){
			string str = "";

			for (int i = 0; i < n; i++ )
				str = str + "   ";

			cout << str << " [ " << n << " , " << l << " ] Leaf Node with Coefficients \n";
			cout << str << "  ";
			s[n][l].toStr();
			cout << "\n";
		}
		else {
			print_tree( n+1, 2*l );
			print_tree( n+1, 2*l+1 );
		}
	}

	/***************************************************/
	/* Deep copy from one object to another			   */
	/***************************************************/
	void copy( Function *other){
       
		this->compressed = other->compressed;
		this->f = other->f;
		this->s = other->s;
		this->d = other->d;

		for ( auto itr = other->s.begin(); itr != other->s.end(); itr ++) {
			for ( auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); itr2++)
				addpair ( this->s, itr->first, itr2->first, itr2->second);
		}

		for ( auto itr = other->d.begin(); itr != other->d.end(); itr ++) {
			for ( auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); itr2++)
				addpair ( this->d, itr->first, itr2->first, itr2->second);
		}
    }

    /***************************************************/
	/* Initialize the class Matrices				   */
	/***************************************************/
	void init_twoscale ( int k ) {

		vector< vector<ldouble>> hg_ = twoscalecoeffs( k );

		Matrix tmp1( 2*k, 2*k );
		Matrix tmp2( 2*k, k);
		Matrix tmp3( 2*k, k);
		Matrix tmp4( 2*k, 2*k);

		hg 	= tmp1;
		hg0 = tmp2;
		hg1 = tmp3;
		hgT = tmp4;

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
		G.gauss_legendre(order, quad_x, quad_w);
		this->quad_npt = npt = quad_w.len();

		Matrix tmp1(npt, k);
		Matrix tmp2( k, npt);
		Matrix tmp3(npt, k);

		this->quad_phi  = tmp1;
		this->quad_phiT = tmp2;
		this->quad_phiw = tmp3;

		for ( int i = 0; i < npt; i++ ) {

			Vector p( phi ( quad_x[i], k) );
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
	Vector project ( int n, int l ){
		Vector s( this->k );
		ldouble h = pow( 0.5, n );
		ldouble scale = sqrt( h );

		for ( int mu = 0; mu < this->quad_npt; mu++ ) {

			ldouble x = ( l + this->quad_x[ mu ] ) * h ;
			ldouble f_ = this->f(x);

			for ( int i = 0; i < k; i++ ){

				s[i] += scale * f_ * this->quad_phiw[mu][i];
			}
		}

		return s;
	}


	/***************************************************/
	/*   refine numerical representation of f(x) 	   */
	/*   to desired tolerance				 		   */
	/***************************************************/
	void refine(int n, int l) {

		int k = this->k;

		Vector s0 = project(n + 1, 2 * l);
		Vector s1 = project(n + 1, 2 * l + 1);

		/* size of s would be 2*k */
		Vector s(s0, s1);
		Vector d = (s * hgT);
		Vector tmp ( d, k, d.len() );

		if ( tmp.normf() < this->thresh  || n >= this->max_level - 1 ) {
			addpair ( this->s, n+1, 2*l, s0);
			addpair ( this->s, n+1, 2*l+1, s1);
		}
		else {
			this->refine(n + 1, 2 * l);
			this->refine(n + 1, 2 * l + 1);
		}
	}

	/***************************************************/
	/* eval f(x) using adaptively refined numerical    */
	/* representation of f(x) 						   */
	/***************************************************/

	ldouble evaluate( int n, int l, ldouble x){
        
		if ( this->s[n].find(l) != this->s[n].end() ) {
			Vector p = ( phi ( x, this->k ));
			return  ( p.inner( this->s[n][l] ) * sqrt( pow ( 2.0, n)) );
		}
		else {
			n = n + 1;
			l = 2 * l;
			x = 2.0 * x;
			if ( x >= 1 ) {
				l = l+1;
				x = x-1;
			}
			return evaluate(n, l, x);
		}
	}

	/*********************************************************/
	/* Evaluate function at x ... scaling function basis only*/
	/* call to self after creation 							 */
	/* looks like a Function evaluation 					 */
	/* say g = Function(5, 1e-3, f) so g(1.0) should be  	 */
	/* similar to f(1.0) 									 */
	/*********************************************************/

	ldouble __f( ldouble x){

		if ( this->compressed ){
			this->reconstruct();
		}
		return this->evaluate(0, 0, x);

	}

	/***************************************************/
	/*  Return sqrt(integral(f(x)**2))			 	   */
	/***************************************************/
	ldouble norm2 ( ){
		 
		if ( this->compressed )
			this->reconstruct();

		ldouble sum = 0.0;

		for ( auto itr = this->s.begin(); itr != this->s.end(); itr++ ){
			for ( auto itr1 = (itr->second).begin(); itr1 != (itr->second).end(); itr1++ )
				sum += pow ( (itr1->second).normf(), 2 );
		}

		return sqrt(sum);
	}


	/*********************************************************/
	/*change from scaling function basis to multi-wavelet.   */
	/* basis (s -> d) tree is filled out with s[0][0] and d  */
    /* n is level in tree l is box index			         */
	/*********************************************************/

	void compress( int n = 0, int l = 0 ) {

		if ( this->compressed ) 
			return;

		if ( this->s[n+1].find(2*l) == this->s[n+1].end() )
			this->compress(n + 1, 2 * l);

		if ( this->s[n+1].find( 2*l + 1 ) == this->s[n+1].end() )
			this->compress(n + 1, 2 * l + 1);

		int k = this->k;

		Vector s( this->s[n+1][2*l], this->s[n + 1][2 * l + 1] );

		Vector d = ( s * this->hgT );

		Vector tmp1( d, 0, k );
		Vector tmp2( d, k, d.len());

		addpair( this->s, n, l, tmp1);
		addpair( this->d, n, l, tmp2);

		this->s[n+1].erase( 2*l );
		this->s[n+1].erase( 2*l + 1);

		if ( n==0 )
			this->compressed = 1;
	
	}


	/**********************************************************************/
	/* change from multi-wavelet basis to scaling function basis (d -> s) */
	/* tree just has s at leaves 										  */
	/* n is level in tree 												  */
	/* l is box index 													  */
	/**********************************************************************/
	
	void reconstruct( int n=0, int l=0 ) {
		
		if ( ! this->compressed )
			return;

		if ( this->d[n].find(l) != this->d[n].end() ) {

			Vector d( this->s[n][l], this->d[n][l] );

			this->d[n].erase(l);
			this->s[n].erase(l);

			Vector  s = ( d * this->hg) ;

			Vector tmp1( s, 0, k );
			Vector tmp2( s, k, s.len());

			addpair( this->s, n+1, 2*l, tmp1);
			addpair( this->s, n+1, 2*l+1, tmp2);

			reconstruct(n + 1, 2 * l);
			reconstruct(n + 1, 2 * l + 1);
		}

		if ( n == 0 )
			this->compressed = 0;
	}


	/**********************************************************************/
	/* recursive "iteration" for mul 							          */
	/* multiply f1 and f2 put result into self 							  */
	/**********************************************************************/

	void mul_iter( Function &f1, Function &f2, int n=0, int l=0) {
		
		if ( f1.s[n].find(l) != f1.s[n].end() && f2.s[n].find(l) != f2.s[n].end() ) {

			if ( autorefine &&  ( n+1 <= this->max_level) ) {   
			
				/* refine both one more level */
				f1.recur_down( n, l, f1.s[n][l] );
				f2.recur_down( n, l, f2.s[n][l] );

				/* scale factor for this level = sqrt((2^d)^(n+1)) ?? */
				ldouble scale_factor = sqrt( pow ( 2.0, n+1) );

				/* multiply f1.s[n+1][2*l] and f2.s[n+1][2*l] */
				Vector f = f1.s[n+1][2*l] * this->quad_phiT;
				Vector g = f2.s[n+1][2*l] * this->quad_phiT;
				f.emul(g);

				Vector tmp =  f * this->quad_phiw;
				addpair ( this->s, n+1, 2*l, tmp.scale (scale_factor));

				/* multiply f1.s[n+1][2*l+1] and f2.s[n+1][2*l+1] */
				Vector f_ = f1.s[n+1][2*l+1] * this->quad_phiT;
				Vector g_ = f2.s[n+1][2*l+1] * this->quad_phiT; 

				f_.emul( g_ );
				Vector tmp_ = f_ * this->quad_phiw;
				addpair ( this->s, n+1, 2*l+1, tmp_.scale(scale_factor));
			}
			else {

				/* if autorefine is not set or we are at the max_level */
				/* live with what you get */
				Vector f = f1.s[n][l] * this->quad_phiT ; 
				Vector g = f2.s[n][l] * this->quad_phiT;
				f.emul( g );
				ldouble scale_factor = sqrt (pow (2.0, n));

				/* scale factor for this level = sqrt((2^d)^(n+1)) */
				Vector tmp = f * this->quad_phiw;
				addpair ( this->s, n, l, tmp.scale(  scale_factor ));
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

	Function operator* (  Function other ){

		if ( this->compressed )
			this->reconstruct();

		if ( other.compressed )
			other.reconstruct();

		Function result(this->k, this->thresh);
		result.mul_iter( *this, other );

		this->sclean();
		other.sclean();

		return result;
	}



	/***************************************************/
	/*   recursive "iteration" for gaxpy.       	   */
	/***************************************************/
	void gaxpy_iter ( ldouble alpha, Function other, ldouble beta, int n=0, int l=0 ) {
		
		if ( this->d[n].find(l) != this->d[n].end() || other.d[n].find(l) != other.d[n].end() ){

			if ( this->d[n].find(l) != this->d[n].end() && other.d[n].find(l) != other.d[n].end() )
				this->d[n][l].gaxpy( alpha, other.d[n][l], beta);
			
			else if ( this->d[n].find(l) == this->d[n].end() && other.d[n].find(l) != other.d[n].end() ){
				Vector tmp =  other.d[n][l];
				addpair ( d,  n , l, tmp.scale(beta));
			}

			else if ( this->d[n].find(l) != this->d[n].end() && other.d[n].find(l) == other.d[n].end() )
				this->d[n][l].scale( alpha);

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

		this->s[0][0].gaxpy ( alpha, other.s[0][0], beta );
		this->gaxpy_iter(alpha, other, beta);

		/* return self so operations can be chained */
		return *this;
	}


	/***************************************************/
	/*   basic addition 					       	   */
	/***************************************************/
	Function operator+ ( Function other ){

		Function result( this->k, this->thresh );
		result.copy ( this );
		result.gaxpy( 1.0, other, 1.0 );

		return result;
	}	


	/***************************************************/
	/*   basic subtraction 					       	   */
	/***************************************************/
	Function operator- ( Function other ){

		Function result( this->k, this->thresh );
		result.copy ( this );
		result.gaxpy( 1.0, other, -1.0 );

		return result;
	}

	/***************************************************/
	/* Mostly for debugging, print summary of  		   */
	/* coefficients, optionally printing the norm of   */
	/* each block 	 			  					   */
	/***************************************************/
	void summarize( int printcoeff = 0 ) {

		cout << "\nsum coefficients\n";

		for (auto itr = this->s.begin(); itr != this->s.end(); itr++) {

			ldouble sum = 0.0;
			int n = itr->first;

			for ( auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); itr2++) {
				int l = itr2->first;
				if ( printcoeff ){
					cout << left << setw(15) << n;
					cout << left << setw(15) << l;
					cout << left << setw(15) << this->s[n][l].normf( );
				}
				else
					sum +=  pow( this->s[n][l].normf( ), 2);
			}
			if ( !printcoeff ) {
				if ( this->s[n].size() != 0 ) {
					cout << left << setw (10 ) << "level=";
					cout << left << setw (15 ) << n;
					cout << left << setw (10 ) << "boxes=";
					cout << left << setw (15 ) << this->s[n].size();
					cout << left << setw (10 ) << "norm=";
					cout << left << setw (15 ) << sqrt(sum);
					cout << "\n";
				}
			}
		}

		cout << "\ndifference coefficients\n";
		for (auto itr = this->s.begin(); itr != this->s.end(); itr++) {
			int n = itr->first;
			ldouble sum = 0.0;

			for ( auto itr2 = this->d[n].begin(); itr2 != this->d[n].end(); itr2++) {
				int l = itr2->first;

				if ( printcoeff )
					cout << setw(10) << n << setw(10) << l << this->d[n][l].normf( );
				else
					sum +=  pow( this->d[n][l].normf( ), 2);
			}

			if ( !printcoeff ) {
				if ( this->d[n].size() != 0 ) {
					cout << " level=   " << setw(10) << n << "   boxes=  "<< setw(10) << this->d[n].size() << "  norm=. " << sqrt(sum);
					cout  << "\n";
				}	
			}
		}
	}

	/*************************************************************/
	/* In s are scaling coefficients for box n,l ... apply		 */ 
	/* twoscale to generate the corresponding coefficients on    */
	/* level n+1 and insert the results into the tree of scaling */ 
	/* function coefficients.  					       	 	     */
	/*************************************************************/

	void recur_down ( int n, int l, Vector & s ) {

		int k = this->k;
		Vector d (2*k );

		for (int i = 0; i < k; i++ )
			d[i] = s[i]; 

		Vector s__ = ( d * this->hg );

		Vector tmp1 ( s__, 0, k );
		Vector tmp2 ( s__, k, s__.len());
		addpair( this->s, n + 1, 2 * l, tmp1 );
		addpair( this->s, n + 1, 2 * l + 1, tmp2 );

	}


	/***************************************************/
	/* get_coeffs 				    	 			   */
	/***************************************************/	
	Vector get_coeffs ( int n, int l ) {

		Vector s;
		if ( l < 0 || l >= pow( 2, n) ) {
			Vector tmp( this->k );
			return tmp;
		}
		
		if ( this->s.find(n) != this->s.end() && this->s[n].find(l) != this->s[n].end() ) {
			return this->s[n][l];
		}

		if ( n > 0 ){
			Vector s__ = ( this->get_coeffs(n - 1, l / 2 ) );
			if ( s__.len() == 0)
				return s__;
			else{
				s = s__;
			}
		}
		else{ 
			return s;
		}
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


	/****************************************************************/
	/* Differentiate the function, which corresponds to application */
	/* of a block triadiagonal matrix.  For an adaptively refined   */
	/* target function we may need to refine boxes down until three */
	/* boxes exist in the same scale.				    	 	    */
	/****************************************************************/	
	Function diff( Function &result, int n = 0, int l = 0 ){

		if ( n == 0){
			if ( this->compressed ) 
				this->reconstruct();
		}
		
		if ( ( this->s.find(n) != this->s.end()) && ((this->s[n]).find(l) == (this->s[n]).end()) ){
			
			this->diff( result, n+1, 2*l );
			this->diff( result, n+1, 2*l+1 );
		}
		else {
			Vector sm = ( this->get_coeffs(n,l-1) );
			Vector sp = ( this->get_coeffs(n,l+1) );
			Vector s0 = ( this->s[n][l] );

			if ( sm.len() && s0.len() && sp.len() ) {
				Vector r =  ( this->rp * sm ) +  ( this->r0 * s0 ) + ( this->rm * sp );

				Vector tmp( r.scale( pow( 2.0, n) ) );
				addpair( result.s, n, l, tmp);
			}
			else {
				this->recur_down( n, l, s0 );
				this->diff( result, n+1, 2*l );
				this->diff( result, n+1, 2*l+1 );
			}
		}

		if ( n == 0 ) {
			this->sclean();
			return result;
		}

		Function c;
		return c;
	}


	/*****************************************************************/
	/*  Return coefficients for box [n][l] given function values at  */
	/* quadrature points within same box.				   	         */
	/*****************************************************************/

	Vector quad_values_to_coeff( Vector values, int n, int l, Matrix transform_matrix) {

		Vector coeff = ( values * transform_matrix ).scale( sqrt( pow(2, -n) ) ); 
		return coeff;
	}

	/*****************************************************************/
	/* Return the largest occupied level (i.e. finest refinement)    */ 
	/* in the scaling function basis.					   	         */
	/*****************************************************************/
	int finest_level () {
	
		int n;
		if ( this->compressed )
			 this->reconstruct();

		for ( n = this->max_level; n >= 0; n-- ){
			if( this->s[n].size() > 0 )
				break;
		}
		return n;
	}


	/*****************************************************************/
	/*  Return a list of occupied levels (i.e. all n such that 		 */
	/*	self.s[n] != {}).				   	         				 */
	/*****************************************************************/
	vector<int> occupied_levels( ){

		if ( this->compressed )
			this->reconstruct();

		vector<int> result;
		int n;

		for ( n = this->max_level; n >= 0 ; n-- ){
			if( this->s[n].size() != 0 )
				result.push_back(n);
		}
		return result;
	}


	/*****************************************************************/
	/* Return a list of tuples of all leaves    					 */
	/* (i.e. finest level n & l)				   	         		 */
	/*****************************************************************/
	vector<pair<int, int>> get_leaves( ) {

		vector<int> nrange = this->occupied_levels();
		vector<pair<int, int>> result;

		for ( int i = 0; i < nrange.size(); i++ ){
			int n = nrange[i];
			for ( auto itr = this->s[n].begin(); itr != this->s[n].end(); itr++ )
				result.push_back( { n, itr->first } );
		}

		return result;
	}

	/*****************************************************************/
	/* refine numerical representation of f(x) to desired tolerance  */
	/* but don't refine any more than the finest level of other 	 */
	/* n is level in tree 											 */
	/* l is box index 							   					 */
	/*****************************************************************/
	void refine_limited ( Function other, int n, int l) {
	
		if ( other.compressed )
			other.reconstruct();

		// project f(x) at next level

		Vector s0 = this-> project(n + 1, 2 * l);
		Vector s1 = this->project(n + 1, 2 * l + 1);

		int k = this->k;

		if( s0.len() != k || s1.len() != k )
			cout << "Wrong initialization\n";

		Vector s( s0, s1);

		// apply the two scale relationship to get difference coeff
		// in 1d this is O(k^2) flops (in 3d this is O(k^4) flops)

		Vector d = ( s * this->hgT );

		// check to see if within tolerance
		// normf() is Frobenius norm == 2-norm for vectors
		Vector V( d, k, d.len() );

		if ( ( V.normf() < this->thresh) ||
			 ( other.s[n + 1].find(2 * l) !=  other.s[n + 1].end() && 
			   other.s[n + 1].find(2 * l + 1) != other.s[n + 1].end() ) )
		{
			// put into tree at level n+1
			addpair( this->s, n+1, 2*l, s0);
			addpair( this->s, n+1, 2*l+1, s1);

		}	
		else if ( other.s[n + 1].find(2 * l) != other.s[n + 1].end() )
		{
			addpair ( this->s, n+1, 2*l, s0);
			this->refine_limited( other, n + 1, 2 * l + 1);
		}
		else if ( other.s[n + 1].find(2 * l + 1) != other.s[n + 1].end() )
		{
			this->refine_limited(other, n + 1, 2 * l);
			addpair ( this->s, n+1, 2*l+1, s1);
		}
		else
		{
			// these recursive calls on sub-trees can go in parallel
			this->refine_limited(other, n + 1, 2 * l);
			this->refine_limited(other, n + 1, 2 * l + 1);
		}
	}



	/*****************************************************************/
	/* Take the inner product of the function with another 			 */
	/* madness function. 											 */
	/*****************************************************************/
	ldouble inner( Function other){

		if ( ! this->compressed )
			this->compress();

		if ( !other.compressed )
			other.compress();

		ldouble result = this->s[0][0].inner( other.s[0][0] ); 

		int n = 0;

		while  ( this->d.find(n) != this->d.end() &&
			     other.d.find(n) != other.d.end() ) {

			for ( auto itr = this->d[n].begin(); itr != this->d[n].end(); itr++ ){
				int l = itr->first;
				if( other.d[n].find(l) != other.d[n].end() ){
					result += this->d[n][l].inner( other.d[n][l] );
				}
			}
			n += 1;
		}
		return result;
	}


	/*****************************************************************/
	/* Take the inner product in the box [n][l] with an external  	 */
	/* analytic expression (i.e. not a madness function)			 */
	/*****************************************************************/

	ldouble box_quad( Function other, int n, int l){
	
		if ( this->compressed )
			this->reconstruct();

		Vector x(this->quad_npt);
		Vector g(this->quad_npt);

		ldouble h = pow( 0.5, n);
		ldouble scale = sqrt(h);

		for (int mu = 0; mu < this->quad_npt; mu++){

			x[mu] = ( this->quad_x[mu] + l ) * h;
			g[mu] = other.__f(x[mu]);
		}

		Vector tmp = this->quad_phiw * this->s[n][l];
		return  ( tmp.inner(g) * scale );
	}


	/*****************************************************************/
	/* Call box_quad iteratively until convergence 					 */
	/* Take the inner product in the box [n][l] with an external  	 */
	/* analytic expression (i.e. not a madness function)			 */
	/*****************************************************************/

	ldouble box_quad_iter( Function other, int n, int l, ldouble old=0.0, ldouble thresh = 0, bool debug = false){ 
	
		ldouble result;
		if (old == 0.0)
			old = this->box_quad(other, n, l);

		if ( !thresh )
			thresh = this->thresh;
		    
		this->get_coeffs(n + 1, 2 * l);
		this->get_coeffs(n + 1, 2 * l + 1);
		ldouble i1 = this->box_quad(other, n + 1, 2 * l);
		ldouble i2 = this->box_quad(other, n + 1, 2 * l + 1);
		ldouble new_ = i1 + i2;

		if ( debug ){
			cout << "--- in box [" << n << "][" << l << "]:" << "\n";
			cout << "\t      new = " << new_ << "n";
			cout << "\t      old = " << old  << "\n";
			cout << "\t rel_diff = " << abs((new_ - old)/new_) << "\n";
		}

		if ( abs((new_ - old)/new_) <= thresh ){
			result = new_;
		}
		else{
			result = this->box_quad_iter(other, n + 1, 2 * l, i1, thresh, debug);
			result += this->box_quad_iter(other, n + 1, 2 * l + 1, i2, thresh, debug);
		}
		return result;
	}


	/*****************************************************************/
	/* Take the inner product of the function with an external 		 */
	/* analytic function (i.e. not a madness function). 			 */
	/*****************************************************************/

	ldouble inner_ext( Function other) {

		if ( this->compressed )
			this->reconstruct();

		this->sclean();
		vector< pair< int, int > > leaves = this->get_leaves();
		ldouble result = 0.0;

		for ( int i = 0; i < leaves.size(); i++ ){
			result += this->box_quad_iter(other, leaves[i].first, leaves[i].second);
		}

		this->sclean();
		return result;
	}

};



int main () {
	
	int npt = 20 ;
    int k = 5;
    ldouble thresh = 1e-5 ;

    vector< ldouble > range;
	for (ldouble i = 0.0; i < (ldouble)npt+1.0; i++)
		range.push_back(i);

    /*****************************************************************/
	/* Test 1	 		 											 */
	/*****************************************************************/
	
	for( int j = 0; j < 3; j++ ){

		Function F( k,  thresh, test[j] );
		cout << "\nNorm of function is " << F.norm2() << "\n";

		for ( int i = 0; i < range.size(); i++ ){
			ldouble x = (ldouble) range[i];
			x = x/(ldouble)npt;
			string s1 = "f(" + to_string(x) + ")=";
			string s2 = "Exact(" + to_string(x) + ")=";

			cout << left << setw(20) << s1;
			cout << left << setw(20) << F.__f( x );
			cout << left << setw(20) << s2;
			cout << left << setw(20) << test[j](x);
			cout << left << setw(10) << "Err=";
			cout << left << setw(20) << F.__f(x) - test[j](x);
			cout << "\n";
		}

		cout << "\ncoefficients before compressing";
	    F.summarize();

	    F.compress();
	    cout <<  "\ncoefficients after compressing";
	    F.summarize();

	    F.reconstruct();
	    cout << "\ncoefficients after reconstructing";
	    F.summarize();

	  	cout << "\n";
	    Function fd( F.k, F.thresh);
	    Function df ( F.diff( fd ) );

	    for ( int i = 0; i < range.size(); i++ ){
			ldouble x = (ldouble) range[i];
			x = x/(ldouble)npt;
			string s1 = "f(" + to_string(x) + ")=";
			string s2 = "Exact(" + to_string(x) + ")=";
			cout << left << setw(20) << s1;
			cout << left << setw(20) << df.__f( x );
			cout << left << setw(20) << s2;
			cout << left << setw(20) << dtest[j](x);
			cout << left << setw(10) << "Err=";
			cout << left << setw(20) << df.__f(x) - dtest[j](x);

			cout << "\n";
		}
	}

	/*****************************************************************/
    /* Addition test which in turn tests gaxpy 						 */
    /*****************************************************************/
	for( int j = 0; j < 3; j++ ){
    
	    cout << "\n";

	    Function f1(k,thresh,test[0]);
	    cout << "norm of f1 is  " << f1.norm2() << "\n";

	    Function f2(k,thresh,test[j]);
	    cout << "norm of f2 is  " << f2.norm2() << "\n";

	    Function f3 = f1 + f2;
	    cout <<  "norm of f3 = f1 + f2 is  " << f3.norm2() << "\n";

	    f3.summarize();

		for ( int i = 0 ; i < range.size(); i++){
			ldouble x = range[i];
			x = x/ (ldouble) npt;

	        ldouble f3_x = f3.__f(x);
	        ldouble exact_x = test[0](x) + test[j](x);
	        ldouble err_x = f3_x - exact_x;
	        string s1 = "f3(" + to_string(x) + ")=";
	        string s2 = "Exact(" + to_string(x) + ")=";

			cout << left << setw(20) << s1;
			cout << left << setw(20) << f3_x;
			cout << left << setw(20) << s2;
			cout << left << setw(20) << exact_x;
			cout << left << setw(10) << "Err=";
			cout << left << setw(20) << err_x;
			cout << "\n";
	        if ( err_x > thresh )
	            cout << left << setw(20) << "outside thresh" << thresh - err_x << "\n";
		}
	}

    /*****************************************************************/
    /* multiplication test which in turn tests gaxpy 				 */
    /*****************************************************************/

	autorefine = 1;
	for ( int j = 0; j < 3; j++ ){
	   
	    Function f1(k,thresh,test[0]);
	    cout <<  "\nnorm of f1 is  " << f1.norm2() << "\n";

	    Function f2(k,thresh,test[j]);
	    cout << "norm of f2 is  " << f2.norm2() << "\n";

	    Function f3 = f1 * f2;
	    cout << "norm of f3 = f1 * f2 is  " << f3.norm2() << "\n";
	    f3.summarize();

		for ( int i = 0; i < range.size(); i++ ){
			ldouble x = range[i];
			x /= ( ldouble )npt;
			ldouble f3_x = f3.__f(x);
			ldouble exact_x = test[0](x) * test[j](x);
			ldouble err_x = f3_x - exact_x;
			string s1 = "f3(" + to_string(x) + ")=";
	        string s2 = "Exact(" + to_string(x) + ")=";

			cout << left << setw(20) << s1;
			cout << left << setw(20) << f3_x;
			cout << left << setw(20) << s2;
			cout << left << setw(20) << exact_x;
			cout << left << setw(10) << "Err=";
			cout << left << setw(20) << err_x << "\n";

			if( err_x > thresh )
				cout << "outside thresh  " <<  thresh - err_x << "\n";
		}
	}

}




