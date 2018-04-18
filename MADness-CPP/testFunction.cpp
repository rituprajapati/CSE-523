
/*****************************************************************/
/* gaussian with square normalized to 1							 */
/*****************************************************************/

ldouble test1 ( ldouble x ){
	
	// ldouble a = 500.0;
	// ldouble tmp = 2 * a / pi;
	// ldouble tmp2 = -a * pow(( x-0.5 ), 2);
	// return pow( tmp , 0.25 ) * exp( tmp2 );
	return x;

}

/*****************************************************************/
/* superposition of multiple gaussians							 */
/*****************************************************************/
ldouble test2 (  ldouble x ) { 

	return ( test1( x-0.3) + test1 ( x ) + test1 ( x+0.3 ) );
}

/*****************************************************************/
/* a more interesting (singular and oscillating) function 		 */
/* ... note it is never computed exactly at the singularity 	 */
/* except by the test code below.							     */
/*****************************************************************/

ldouble test3 ( ldouble x ){
	ldouble a = 100.0 * pi;

	if ( x == 0.5 ){
		return 0.0;
	}
	else{
		return cos ( a * ( x-0.5) ) / sqrt ( abs ( x-0.5 ));
	}
}

/*****************************************************************/
/* derivative of test1 											 */
/*****************************************************************/

ldouble dtest1 ( ldouble x ) {

	// ldouble a = 500.0;
	// ldouble tmp1 =  -2.0 * a * (x-0.5);
	// ldouble tmp2 =  pow( 2 * a /pi, 0.25 );
	// ldouble tmp3 =  exp( -a * pow( (x-0.5), 2 ));
	// return tmp1 * tmp2 * tmp3;
	return 1;
}


/*****************************************************************/
/* derivative of test2 											 */
/*****************************************************************/

ldouble dtest2 ( ldouble x ){

	return  ( dtest1( x-0.3 ) + dtest1( x )  + dtest1( x+0.3) );
}


/*****************************************************************/
/* derivative of test3 											 */
/*****************************************************************/

ldouble dtest3 ( ldouble x ){

	ldouble a = 100.0 * pi;

	if ( x == 0.5 ){
		return 0.0;
	}
	else{

		ldouble s = 1.0;
		if ( x < 0.5 ){
			s = -1.0;
		}
		ldouble tmp1 = -a * sin( a * (x-0.5) );
		ldouble tmp2 = sqrt( abs(x-0.5) );
		ldouble tmp3 = s * 0.5 * cos( a * (x-0.5) );
		ldouble tmp4 = pow( abs( x-0.5 ), 1.5 );

		return ( tmp1 / tmp2 - tmp3 / tmp4 );
	}
} 

ldouble ( *test[3] ) ( ldouble x ) = { test1, test2, test3 };
ldouble ( *dtest[3] ) ( ldouble x ) = { dtest1, dtest2, dtest3 };


