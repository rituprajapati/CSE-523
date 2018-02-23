#include<math.h>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <assert.h>

using namespace std;
#define ldouble  long double


class Vector {

	vector<ldouble> a;

public:
	Vector () {
		cout << "Do nothing\n";
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
	/* combining two Vectors 					 */
	/*********************************************/

	Vector ( Vector obj1, Vector obj2 ) {


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

		for (int i = 0; i < a.size(); i++ )
			sum += a[i] * a[i];

		return sqrt(sum);
	}

	/*********************************************/
	/* sum of two vector object a				 */
	/*********************************************/
	ldouble inner( Vector other ){

		ldouble sum = 0.0;

		for (int i = 0; i < a.size(); i++ ) {
			sum += a[i] * other.a[i];
		}

		return sum;
	}

	/*********************************************/
	/* update a[] acc to alpha and beta			 */
	/*********************************************/
	Vector gaxpy( ldouble alpha, Vector other, ldouble beta) {
		
		for( int i = 0; i < a.size(); i++) {
			a[i] = alpha * a[i] + beta * other.a[i];
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
	Vector emul( Vector other) {
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

	Vector getSlice (int lo, int hi) {
		vector<ldouble>  tmp( a.begin()+lo, a.begin()+hi);
		Vector obj(tmp);
		return obj;
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
	string toStr() {
		string str = "";
		for (int i = 0; i < a.size(); i++) {
			str += to_string( a[i]);
			str += "  ";
		}
		return str;
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
		cout << "do nothing\n\n";
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
	
	cout << "Matrix * Vector\n";

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

	cout << "Vector * Matrix\n";
 	n = M.len();
 	if( n!= 0 )
 		m = M[0].len();

 	if( n != V.len() )
 		assert ( "Invalid dimensions\n" );

 	Vector r (m);

 	for ( int i = 0; i < m; i++ ) {
 		for ( int j = 0; j < n; j++ ){
 			r[i] += V[j] * M[j][i];
 		}
 	}
    return r;
}




class Gauss {

	vector<vector< ldouble > > points;
	vector<vector< ldouble > > weights;

public:
	Gauss () {

		points.push_back ( {} );
		weights.push_back ( {} );

		points.push_back ( {0.50000000000000000} );
		weights.push_back ( { 1.00000000000000000} );

		points.push_back (  {0.78867513459481287, 0.21132486540518713} );
		weights.push_back ( {0.50000000000000011, 0.50000000000000011} );

		points.push_back ( {0.88729833462074170, 0.50000000000000000, 0.11270166537925830} );
		weights.push_back ( { 0.27777777777777751, 0.44444444444444442, 0.27777777777777751} );

		points.push_back ( {0.93056815579702623, 0.66999052179242813, 0.33000947820757187, 0.06943184420297371} );
		weights.push_back ( { 0.17392742256872701, 0.32607257743127305, 0.32607257743127305, 0.17392742256872701} );
		
		points.push_back ( {0.95308992296933193, 0.76923465505284150, 0.50000000000000000, 0.23076534494715845, 0.04691007703066802} );
		weights.push_back ( { 0.11846344252809465, 0.23931433524968321, 0.28444444444444444, 0.23931433524968321, 0.11846344252809465} );
		
		points.push_back ( {0.96623475710157603, 0.83060469323313224, 0.61930959304159849, 0.38069040695840156, 0.16939530676686776, 0.03376524289842397} );
		weights.push_back ( { 0.08566224618958508, 0.18038078652406936, 0.23395696728634546, 0.23395696728634546, 0.18038078652406936, 0.08566224618958508} );

		points.push_back ( {0.97455395617137930, 0.87076559279969723, 0.70292257568869854, 0.50000000000000000, 0.29707742431130141, 0.12923440720030277, 0.02544604382862070} );
		weights.push_back ( { 0.06474248308443417, 0.13985269574463832, 0.19091502525255946, 0.20897959183673470, 0.19091502525255946, 0.13985269574463832, 0.06474248308443417} );
		
		points.push_back ( {0.98014492824876809, 0.89833323870681336, 0.76276620495816450, 0.59171732124782495, 0.40828267875217511, 0.23723379504183550, 0.10166676129318664, 0.01985507175123191} );
		weights.push_back ( { 0.05061426814518921, 0.11119051722668717, 0.15685332293894369, 0.18134189168918102, 0.18134189168918102, 0.15685332293894369, 0.11119051722668717, 0.05061426814518921} );
		
		points.push_back ( {0.98408011975381304, 0.91801555366331788, 0.80668571635029518, 0.66212671170190451, 0.50000000000000000, 0.33787328829809554, 0.19331428364970482, 0.08198444633668206, 0.01591988024618696} );
		weights.push_back (  { 0.04063719418078738, 0.09032408034742866, 0.13030534820146775, 0.15617353852000135, 0.16511967750062989, 0.15617353852000135, 0.13030534820146775, 0.09032408034742841, 0.04063719418078738});

		points.push_back ( {0.98695326425858587, 0.93253168334449232, 0.83970478414951222, 0.71669769706462361, 0.57443716949081558, 0.42556283050918442, 0.28330230293537639, 0.16029521585048778, 0.06746831665550773, 0.01304673574141413} );
		weights.push_back ( { 0.03333567215434403, 0.07472567457529021, 0.10954318125799088, 0.13463335965499817, 0.14776211235737641, 0.14776211235737641, 0.13463335965499817, 0.10954318125799088, 0.07472567457529021, 0.03333567215434403} );
		
		points.push_back ( {0.98911432907302843, 0.94353129988404771, 0.86507600278702468, 0.75954806460340585, 0.63477157797617245, 0.50000000000000000, 0.36522842202382755, 0.24045193539659410, 0.13492399721297532, 0.05646870011595234, 0.01088567092697151});
		weights.push_back ( { 0.02783428355808731, 0.06279018473245226, 0.09314510546386703, 0.11659688229599521, 0.13140227225512346, 0.13646254338895031, 0.13140227225512346, 0.11659688229599521, 0.09314510546386703, 0.06279018473245226, 0.02783428355808731} );
	}

public:
	void gauss_legendre (int order, Vector &p, Vector &w ) {
		Vector newP( points[order] );
		Vector newW( weights[order] );
		
		p = newP;
		w = newW;
	}
};




/*****************************************************************/
/*  Return the level-0 blocks rm, r0, rp of the central			 */
/*  difference derivative operator with periodic boundary 		 */
/*  conditions on either side.									 */
/*****************************************************************/
vector< Matrix >  make_dc_periodic ( int k){
    
    Matrix r0(k,k);
    Matrix rp(k,k);
    Matrix rm(k,k);
    
    vector< Matrix > res(3);

    ldouble iphase = 1.0, jphase, gammaij, Kij;

    for ( int i = 0; i < k; i++ ){
        
        jphase = 1.0;

        for ( int j = 0; j < k; j ++ ){

            gammaij = sqrt( (2*i+1)*(2*j+1) );

            if ( (i-j) > 0  &&  ((i-j) %2 ) == 1 )
                Kij = 2.0;
            else
                Kij = 0.0;

            r0[i][j] = 0.5 * (1.0 - iphase*jphase - 2.0*Kij) * gammaij;
            rm[i][j] = 0.5 * jphase * gammaij;
            rp[i][j] = -0.5 * iphase * gammaij;
            jphase = -jphase;
        }
        iphase = -iphase;
    }
    res[0] = rm;;
    res[1] = r0;
    res[2] = rp;
    return res;
}


class Function {

	/***************************************************/
	/* Private data								       */
	/***************************************************/
	int autorefine = 1;	
	int k;
	int thresh;
	int f;
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

public:

	Function () {
		cout << "do_nothing\n";
	}
	/***************************************************/
	/* Constructor function for initializing data      */
	/***************************************************/
	Function( int k, int thr, int f = 0, int initial_level = 2) {

		k = k;
		thresh = thr;
		f = f;
        max_level = 30;

        for (int i = 0; i <= max_level; i++ ) {
        	this->d[i] = {};
        	this->s[i] = {};
        }

        init_twoscale(k);
        init_quadrature(k);

        vector< Matrix > m = make_dc_periodic(k);

        if ( m.size() != 3 )
        	assert ( "Vector size issue\n" );
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

	Function ( Function &other ){

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

			cout << str << "[" << n << "," << l << "] Leaf Node with Coefficients ";
			cout << str << s[n][l].toStr();
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
       
		Function result( k, thresh );
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

		Matrix hg_ = twoscalecoeffs( k );

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
		quad_npt = npt = quad_w.len();

		Matrix tmp1(npt, k);
		Matrix tmp2( k, npt);
		Matrix tmp3(npt, k);

		quad_phi  = tmp1;
		quad_phiT = tmp2;
		quad_phiw = tmp3;

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
        
		Vector s( k );
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

		Vector s0, s1;
		int k = k;

		s0 = project(n + 1, 2 * l);
		s1 = project(n + 1, 2 * l + 1);
		
		/* size of s would be 2*k */
		Vector s(s0, s1);
		
		Vector d = s * hgT;

		Vector tmp ( d.getSlice(k, d.len()) );

		if ( tmp.normf() < thresh  || n >= max_level - 1 ) {
			this->s[n + 1][2 * l] = s0;
			this->s[n + 1][2 * l + 1] = s1;
		}
		else {

			refine(n + 1, 2 * l);
			refine(n + 1, 2 * l + 1);
		}
	}




	ldouble evaluate( int n, int l, ldouble x){
        
		if ( s[n].find(l) != s[n].end() ) {
			Vector p = phi ( x, this->k );
			return  ( p.inner( s[n][l] ) * sqrt( pow ( 2.0, n)) );
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

		if ( s[n+1].find(2*l) == s[n+1].end() )
			compress(n + 1, 2 * l);

		if ( s[n+1].find( 2*l + 1 ) == s[n+1].end() )
			compress(n + 1, 2 * l + 1);

		int k = this->k;

		Vector s ( this->s[n+1][2*l], this->s[n + 1][2 * l + 1] );

		Vector d = s * hgT;

		this->s[n][l] = d.getSlice(0, k);
		this->d[n][l] = d.getSlice(k, d.len());

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

			/* apply the two scale relationship to get difference coeff */
            /* in 1d this is O(k^2) flops (in 3d this is O(k^4) flops). */
			Vector  s = d * this->hg;

			this->s[n + 1][2 * l] = s.getSlice( 0, k); 
			this->s[n + 1][2 * l + 1] = s.getSlice( k, s.len() );

			/* sub-trees can be done in parallel */
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
			if ( autorefine && n+1 <= this->max_level ) {   /* ?? Function.autorefine */
			
				/* refine both one more level */
				f1.recur_down( n, l, f1.s[n][l] );
				f2.recur_down( n, l, f2.s[n][l] );

				/* scale factor for this level = sqrt((2^d)^(n+1)) ?? */
				ldouble scale_factor = sqrt( pow ( 2.0, n+1) );

				/* multiply f1.s[n+1][2*l] and f2.s[n+1][2*l] */
				Vector f = f1.s[n+1][2*l] * this->quad_phiT;
				Vector g = f2.s[n+1][2*l] * this->quad_phiT;
				f.emul(g);

				this->s[n+1][2*l] = ( f * this->quad_phiw ).scale (scale_factor); 

				/* multiply f1.s[n+1][2*l+1] and f2.s[n+1][2*l+1] */
				Vector f_ = f1.s[n+1][2*l+1] * this->quad_phiT;
				Vector g_ = f2.s[n+1][2*l+1] * this->quad_phiT; 

				f_.emul( g_ );

				this->s[n+1][2*l+1] = (f * this->quad_phiw ).scale(scale_factor);
			}

			else {

				/* if autorefine is not set or we are at the max_level */
				/* live with what you get */
				
				Vector f = f1.s[n][l] * this->quad_phiT ; 
				Vector g = f2.s[n][l] * this->quad_phiT;
				f.emul( g );

				/* scale factor for this level = sqrt((2^d)^(n+1)) */
				this->s[n][l] =  ( f * this->quad_phiw).scale( sqrt( pow (2.0, n)) );
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
	/*   * Overloading  					       	   */
	/***************************************************/
	// Function operator* ( Function other ) {
	// 	Function f = mul(other);
	// 	return f;
	// }

	/***************************************************/
	/*   recursive "iteration" for gaxpy.       	   */
	/***************************************************/
	void gaxpy_iter ( ldouble alpha, Function other, ldouble beta, int n=0, int l=0 ) {
	
		if ( this->d[n].find(l) != this->d[n].end() || other.d[n].find(l) != other.d[n].end() ){

			if ( this->d[n].find(l) != this->d[n].end() && other.d[n].find(l) != other.d[n].end() )
				this->d[n][l].gaxpy( alpha, other.d[n][l], beta);
			
			else if ( this->d[n].find(l) == this->d[n].end() && other.d[n].find(l) != other.d[n].end() ){
				Vector tmp =  other.d[n][l];
				this->d[n][l] = tmp.scale (beta);
			}

			else if ( this->d[n].find(l) != this->d[n].end() && other.d[n].find(l) == other.d[n].end() )
				this->d[n][l].scale( alpha);

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

		this->s[0][0].gaxpy ( alpha, other.s[0][0], beta );
		this->gaxpy_iter(alpha, other, beta);

		/* return self so operations can be chained */
		return *this;
	}

	/***************************************************/
	/*   basic addition 					       	   */
	/***************************************************/
	Function operator+ ( Function other ){

		Function result( k, thresh );
		result.compressed = this->compressed;
		result.f = this->f;

		result.s = this->s;
		result.d = this->d;
		result.gaxpy( 1.0, other, 1.0 );

		return result;
	}	


	/***************************************************/
	/*   basic subtraction 					       	   */
	/***************************************************/
	Function operator- ( Function other ){

		Function result( k, thresh );
		result.compressed = this->compressed;
		result.f = this->f;

		result.s = this->s;
		result.d = this->d;
		result.gaxpy( 1.0, other, -1.0 );

		return result;
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

			for ( auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); itr2++) {
				int l = itr2->first;
				if ( printcoeff )
					// cout << "%3d %6d %.2e" << (n, l, self.s[n][l].normf())
					cout << n << l << this->s[n][l].normf( );
				else
					sum +=  pow( this->s[n][l].normf( ), 2);
			}
			if ( !printcoeff ) {
				if ( this->s[n].size() != 0 ) 
					cout << "level=" << n << "boxes=" << this->s[n].size() << "norm=" << sqrt(sum);
			}
		}

		cout << "difference coefficients";
		for (auto itr = this->s.begin(); itr != this->s.end(); itr++) {
			int n = itr->first;
			ldouble sum = 0.0;

			for ( auto itr2 = this->d.begin(); itr2 != this->d.end(); itr2++) {
				int l = itr2->first;

				if ( printcoeff )
					// cout << "%3d %6d %.2e" << (n, l, self.s[n][l].normf())
					cout << n << l << this->d[n][l].normf( );
				else
					sum +=  pow( this->d[n][l].normf( ), 2);
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

	void recur_down ( int n, int l, Vector & s ) {
		
		int k = this->k;
		Vector d (2*k );

		for (int i = 0; i < k; i++ )
			d[i] = s[i]; 

		Vector s__ = d * this->hg;

		for ( int i = 0; i < s__.len(); i++ )
			s[i] = s__[i]; 

		this->s[n + 1][2 * l] = s.getSlice(0,k);
		this->s[n + 1][2 * l + 1] = s.getSlice(k, s.len());
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

		if ( this->s[n].find(l) != this->s[n].end() ) 
			return this->s[n][l];

		if ( n > 0 ){
			s = this->get_coeffs(n-1, l/2);
			if ( s.len() == 0)
				return s;
		}
		else 
			return s;

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
	Function diff( Function result, int n = 0, int l = 0 ){
	
		if ( n == 0){
			if ( this->compressed ) 
				this->reconstruct();
			result = Function( this->k, this->thresh );
		}
		
		if ( this->s[n].find(l) == this->s[n].end() ){
			/* Sub trees can run in parallel */
			/* Run down tree until we hit scaling function coefficients */
			this->diff( result, n+1, 2*l );
			this->diff( result, n+1, 2*l+1 );
		}
		else {
			// These can also go in parallel since may involve
			// recurring up & down the tree.
			Vector sm = this->get_coeffs(n,l-1);
			Vector sp = this->get_coeffs(n,l+1);
			Vector s0 = this->s[n][l];

			if ( sm.len() && s0.len() && sp.len() ) {
				Vector r =  ( this->rp * sm ) +  ( this->r0 * s0 ) + ( this->rm * sp );
				result.s[n][l] = r.scale( pow( 2.0, n) ); 
			}
			else {
				this->recur_down( n, l, s0 );
				// Sub trees can run in parallel
				this->diff( result, n+1, 2*l );
				this->diff( result, n+1, 2*l+1 );
			}
		}

		if ( n == 0 ) {
			this->sclean();
			return result;
		}
		

	}

	/***************************************************/
	/* evaluate_function  	????				   	   */
	/***************************************************/
	// ldouble evaluate_function ( x,n,l):
	
	// coordinate = (x[i]+l)*(2.0**(n))
	// return self.f(coordinate)


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

		for ( n = this->max_level; n >=0 ; n-- ){
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

		if ( ( d.getSlice(k, d.len() ).normf() <this->thresh) ||
			 ( other.s[n + 1].find(2 * l) !=  other.s[n + 1].end() && 
			   other.s[n + 1].find(2 * l + 1) != other.s[n + 1].end() ) )
		{
			// put into tree at level n+1
			this->s[n + 1][2 * l] = s0;
			this->s[n + 1][2 * l + 1] = s1;

		}	
		else if ( other.s[n + 1].find(2 * l) != other.s[n + 1].end() )
		{
			this->s[n + 1][2 * l] = s0;
			this->refine_limited( other, n + 1, 2 * l + 1);
		}
		else if ( other.s[n + 1].find(2 * l + 1) != other.s[n + 1].end() )
		{
			this->refine_limited(other, n + 1, 2 * l);
			this->s[n + 1][2 * l + 1] = s1;
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

};


int main () {
	Vector tmp( 10 );
	Matrix M ;
	// M = Matrix ( 10, 20);
	// for ( int i = 0; i < )
	cout << "Hello";
}















