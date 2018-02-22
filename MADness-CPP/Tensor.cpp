#include<math.h>
#include <vector>
using namespace std;


class Vector {

	vector<double> a;

	/*********************************************/
	/* if arg is an int							 */
	/*********************************************/
	Vector (int arg, int wrap = 0) {

		a = vector<double> (n, 0.0);
	}

	/*********************************************/
	/* if arg is an array						 */
	/*********************************************/
	Vector (vector<double> arg, int wrap = 0) {

		a = vector<double> (arg.begin(), arg.end());
	}

	/*********************************************/
	/* copy constructor							 */
	/*********************************************/
	Vector ( Vector& obj ){

		a = obj.a;
	}
	/*********************************************/
	/* sum of square of all the values in a		 */
	/*********************************************/
	long double normf() {

		long double sum = 0.0;

		for (int i = 0; i < a.length(); i++ )
			sum += a[i] * a[i];

		return sqrt(sum);
	}

	/*********************************************/
	/* sum of two vector object a				 */
	/*********************************************/
	long double inner( Vector other ){
		long double sum = 0.0;

		for (int i = 0; i < a.length(); i++ ) {
			sum += a[i] * other.a[i];
		}

		return sum;
	}

	/*********************************************/
	/* update a[] acc to alpha and beta			 */
	/*********************************************/
	Vector gaxpy( double alpha, Vector other, double beta) {
		
		for( int i = 0; i < a.length(); i++) {
			this->a[i] = alpha * this->a[i] + beta * other.a[i];
		}

		return *this;
	}

	/*********************************************/
	/* multiply each value of a with const s	 */
	/*********************************************/
	Vector scale( double s ) {

		for (int i = 0; i < a.length(); i++ )
			a[i] *= s;

		return *this;
	}

	/*********************************************/
	/* if arg is an array						 */
	/*********************************************/
	Vector emul( Vector other) {
		for (int i = 0 ; i < a.length(); i++ )
			a[i] *= other.a[i];
		return *this;
	} 

	/*********************************************/
	/* Vector v[ind]							 */
	/*********************************************/
	double operator[] ( int ind ) {
		return a[ind];
	}

	/*********************************************/
	/* Vector v[ind] = value					 */
	/*********************************************/
	// double operator[] ( int ind, double value ) {
	// 	a[ind] = value;
	// }

	/*********************************************/
	/* return a Vector object defined over 		 */
	/* a[lo to hi]					 			 */
	/*********************************************/
	Vector getSlice (int lo, int hi) {
		vector<double>  tmp(a.begin()+lo, a.begin()+hi);
		Vector obj = Vector(tmp, wrap=1);
	}

	/*********************************************/
	/* set a[lo to hi] to value					 */
	/*********************************************/
	void setSlice ( int lo, int hi, double value ) {
		for(int i = lo; i < hi; i++ )
			a[i] = value;
	}

	/*********************************************/
	/* Print the object by a 					 */
	/*********************************************/
	string toStr() {
		string str = "";
		for (int i = 0; i < a.length(); i++) 
			str += toString(a[i]);
		return str;
	}
   
	/*********************************************/
	/* obj.len() will return length of a 		 */
	/*********************************************/
    int len(){
        return a.length();
    }

    /*********************************************/
	/* obj3 = obj1 + obj2 				 		 */
	/*********************************************/
    Vector operator+ ( Vector other ){

    	Vector r = Vector( &this );
    	return r.gaxpy (1.0, other, 1.0 );

    }
};




















