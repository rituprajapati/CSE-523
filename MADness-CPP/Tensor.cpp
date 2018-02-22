#include<math.h>


class Vector() {

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
		for(int i = 0; i < arg.length(); i++) {
			a.push_back(arg[i]);
		}
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
	/* sume of two vector object a				 */
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
			a[i] = alpha * a[i] + beta * other.a[i];
		}

		return *this;
	}

	/*********************************************/
	/* multiply each value of a with a const	 */
	/*********************************************/
	Vector scale( double s ) {

		for (int i = 0; i < a.length(); i++ )
			a[i] *= s;

		return *this;
	}

	/*********************************************/
	/* if arg is an array						 */
	/*********************************************/
	def emul(self, other):
        for i in xrange(len(self.a)):
            self.a[i] *= other.a[i]
        return self
};




















