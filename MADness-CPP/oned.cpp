
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

class Function {

	/***************************************************/
	/* Private data								       */
	/***************************************************/
	int autorefine = 1;	
	int k;
	int threh;
	int f;
	int max_level;
	unordered_map<int, unordered_map<int, int> > d;
	unordered_map<int, unordered_map<int, int> > s;
	int rm, r0, rp;
	int compressed;

	/***************************************************/
	/* Constructor function for initializing data      */
	/***************************************************/
	Function( int k, int thresh, int f = 0, int initial_level=2) {

		
		k = k;
		thresh = thresh;
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

	void print_tree( int n = 0, int l = 0){

		if ( s[n].find(l) != s[n].end() ){
			string str = "";

			for (int i = 0; i < n; i++ )
				str = str + "   ";

			cout << str << "[" << n << "," << l << "] Leaf Node with Coefficients ";
			cout << str << s[n][l];
			cout << " ";
		}
		else {
			print_tree(n+1,2*l);
			print_tree(n+1,2*l+1);
		}
	}


	void init_twoscale ( int k) {

		vector<vector<int>> hg_ = twoscalecoeffs( k );

		hg 	= Matrix(2*k, 2*k);
		hg0 = Matrix(2*k, k);
		hg1 = Matrix(2*k, k);
		hgT = Matrix(2*k, 2*k);

		for (int i = 0; i < 2*k ; i++ ) {
			for (int j = 0; j < 2*k; j++ ) {
				hg[i, j]  = hg[i][j]
				hgT[i,j] = hg[j][i]
			}
		}

		for (int i = 0; i < 2*k; i++ ){
			for (int j = 0; j < k; j++ ) {
				hg0[i,j]  = hg[i][j]
				hg1[i,j]  = hg[i][j+k]
			}
		}
	}
	 // Function copy(){
       
		// Function result = Function(this.k, this.thresh);
		// result.compressed = this.compressed;
		// result.f = this.f;

		// for ( auto itr = this.s.begin(); itr != this.s.end(); itr++ )
		// 	for l in self.s[n].keys():
		// 		result.s[n][l] = Vector(self.s[n][l])


		// for n in self.d.keys():
		// for l in self.d[n].keys():
		// result.d[n][l] = Vector(self.d[n][l])
		// return result

  //   }
};























