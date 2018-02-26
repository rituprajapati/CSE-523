
/*****************************************************************/
/* Takes input a map, n, l and a vector                          */
/* insert the combination { n, { l, Vector}} in map              */
/*****************************************************************/

void addpair ( unordered_map<int, unordered_map< int, Vector > >& mp, int n, int l, Vector v){

	if( mp.find( n ) != mp.end() ){
		mp[n].insert( make_pair( l, v) );
	}
	else{
		unordered_map< int, Vector > tmp;
		tmp.insert ( make_pair ( l, v));
		mp.insert ( make_pair( n, tmp));
	}

}
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
