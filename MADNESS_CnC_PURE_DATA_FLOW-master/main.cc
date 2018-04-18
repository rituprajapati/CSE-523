#include "oned.h"


void Multiplication( int k, int max_level, double thresh ){

  int npt = 20;

  vector< double > range;
  for (double i = 0.0; i < (double)npt+1.0; i++)
    range.push_back(i);

  /*****************************************************************/
  /* Multiplication test 
  /*****************************************************************/
  for( int j = 0; j < 3; j++ ){
    
      cout << "\n";

      multiplication_test mul_test_obj( k, thresh, max_level, test[0], test[j]);
      mul_test_obj.projectA_tag.put(std::make_pair(0, 0));
      mul_test_obj.projectB_tag.put(std::make_pair(0, 0)); 
      mul_test_obj.wait();

      mul_test_obj.norm2_f1_tag.put(std::make_pair(0, 0));
      mul_test_obj.wait();
      cout << "norm of f1 is  " << sqrt( mul_test_obj.norm2_result[0]) << "\n";
      mul_test_obj.norm2_result[0] = 0.0;


      mul_test_obj.norm2_f2_tag.put(std::make_pair(0, 0));
      mul_test_obj.wait();
      cout << "norm of f2 is  " << sqrt(mul_test_obj.norm2_result[0]) << "\n";
      mul_test_obj.norm2_result[0] = 0.0;


      mul_test_obj.norm2_mul_tag.put(std::make_pair(0, 0));
      mul_test_obj.wait();
      cout << "norm of Multiplication is  " << sqrt(mul_test_obj.norm2_result[0]) << "\n";
      mul_test_obj.norm2_result[0] = 0.0;

      for ( int i = 0 ; i < range.size(); i++){

        double x = range[i];
        x = x/ (double) npt;

        double f3_x = mul_test_obj.__evaluate(0,0, x, mul_test_obj.mul_item);
        double exact_x = test[0](x) * test[j](x);
        double err_x = f3_x - exact_x;
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

}





void Addition( int k, int max_level, double thresh ){

   int npt = 20;

    vector< double > range;
    for (double i = 0.0; i < (double)npt+1.0; i++)
      range.push_back(i);

  /*****************************************************************/
  /* Addition test which in turn tests gaxpy             */
  /*****************************************************************/
  for( int j = 0; j < 3; j++ ){
    
      cout << "\n";

      addition_test add_obj(k, thresh, max_level, test[0], test[j]);
      add_obj.projectA_tag.put(std::make_pair(0, 0));
      add_obj.projectB_tag.put(std::make_pair(0, 0));   
      add_obj.wait();

      add_obj.norm2_f1_tag.put(std::make_pair(0, 0));
      add_obj.wait();
      cout << "norm of f1 is  " << sqrt(add_obj.norm2_result[0]) << "\n";
      add_obj.norm2_result[0] = 0.0;


      add_obj.norm2_f2_tag.put(std::make_pair(0, 0));
      add_obj.wait();
      cout << "norm of f2 is  " << sqrt(add_obj.norm2_result[0]) << "\n";
      add_obj.norm2_result[0] = 0.0;


      add_obj.norm2_add_tag.put(std::make_pair(0, 0));
      add_obj.wait();
      cout << "norm of addition is  " << sqrt(add_obj.norm2_result[0]) << "\n";
      add_obj.norm2_result[0] = 0.0;

    for ( int i = 0 ; i < range.size(); i++){

      double x = range[i];
      x = x/ (double) npt;

      double f3_x = add_obj.__evaluate(0,0,x, add_obj.reconstruct_result_item);
      double exact_x = test[0](x) + test[j](x);
      double err_x = f3_x - exact_x;
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

}



void Differentiation( int k, int max_level, double thresh ){

  int npt = 20;

  vector< double > range;
  for (double i = 0.0; i < (double)npt+1.0; i++)
    range.push_back(i);


  for( int j = 0; j < 3; j++ ){

    diff_test diff_test_obj( k, thresh, max_level, test[j]);

    //Initialize the 0th and 1st level with 
  
    diff_test_obj.project_tag.put( std::make_pair(0, 0) );
    diff_test_obj.wait();

    diff_test_obj.norm2_tag.put(std::make_pair(0, 0));
    diff_test_obj.wait();

    cout << "\nNorm of function is " << sqrt(diff_test_obj.norm2_result[0]) << "\n";

    for ( int i = 0; i < range.size(); i++ ){

      double x = (double) range[i];
      x = x/(double)npt;
      string s1 = "f(" + to_string(x) + ")=";
      string s2 = "Exact(" + to_string(x) + ")=";

      double F_eval =  diff_test_obj.__evaluate(0,0,x, diff_test_obj.project_item);
      double F_test = test[j](x);


      cout << left << setw(20) << s1;
      cout << left << setw(20) << F_eval;
      cout << left << setw(20) << s2;
      cout << left << setw(20) << F_test;
      cout << left << setw(10) << "Err=";
      cout << left << setw(20) << F_eval - F_test;
      cout << "\n";
    }

    cout << "\nDifferentiation Evaluation\n";

    for ( int i = 0; i < range.size(); i++ ){

      double x = (double) range[i];
      x = x/(double)npt;

      string s1 = "f(" + to_string(x) + ")=";
      string s2 = "Exact(" + to_string(x) + ")=";

      double dF_eval =  diff_test_obj.__evaluate(0,0,x, diff_test_obj.diff_result_item);
      double dF_test = dtest[j](x);

      cout << left << setw(20) << s1;
      cout << left << setw(20) << dF_eval;
      cout << left << setw(20) << s2;
      cout << left << setw(20) << dF_test;
      cout << left << setw(10) << "Err=";
      cout << left << setw(20) << dF_eval - dF_test;

      cout << "\n";
    }
  }

}



void Gaxpy_Sub( int k, int max_level, double thresh ){

  int npt = 20;
  vector< double > range;
  for (double i = 0.0; i < (double)npt+1.0; i++)
    range.push_back(i);

  /*****************************************************************/
  /* Multiplication test 
  /*****************************************************************/
  for( int j = 0; j < 3; j++ ){
    
      // cout << "Gaxpy Sub test on test[0] and test[" <<j << "]\n";  
      gaxpy_sub_test test_obj( k, thresh, max_level, test[0], test[j]);
      test_obj.projectA_tag.put(std::make_pair(0, 0));
      test_obj.projectB_tag.put(std::make_pair(0, 0)); 
      test_obj.wait();

      test_obj.norm2_tag.put(std::make_pair(0, 0));
      test_obj.wait();
      cout << "norm of (test[0]-test[" << j <<  "])- (test[0]-test[" <<j << "]) is  " << sqrt( test_obj.norm2_result[0]) << "\n";
  }

}



int main(int argc, char *argv[]) {
   int k = 5;
   int max_level;
   double thresh;
   int num_threads;

   std::istringstream iss1( argv[1] );
   if( ! (iss1 >> max_level)){
    cout << "Unable to assign max_levels\n";
    return -1;
   }

   std::istringstream iss2( argv[2] );
   if( ! (iss2 >> thresh)){
    cout << "Unable to assign thresh\n";
    return -1;
   }

   std::istringstream iss3( argv[3] );
   if( ! (iss3 >> num_threads)){
    cout << "Unable to assign thresh\n";
    return -1;
   }

   // cout << "MAX_levels : " << max_level << " Thresh: " << thresh << "\n";
   high_resolution_clock::time_point t1 = high_resolution_clock::now();

   // cout <<"\n\n---------------------";
   // cout << "\nDIFFERENTIATION TEST";
   // cout <<"\n---------------------";

   // Differentiation( k, max_level, thresh );

   // cout <<"\n\n--------------";
   // cout << "\nADDITION TEST";
   // cout << "\n--------------";

   // Addition(k, max_level, thresh);

   // cout <<"\n\n-------------------";
   // cout << "\nMULTIPLICATION TEST";
   // cout <<   "\n-------------------";

   // Multiplication( k, max_level, thresh );

   cout <<"\n\n-----------------------------------------------------------------------";
   cout << "\nGAXPY SUBTRACT with max_level:" << max_level << ", Thresh:" << thresh << ", Num of CNC_THREADS: " << num_threads;
   cout <<   "\n-----------------------------------------------------------------------\n";

   Gaxpy_Sub( k, max_level, thresh );

   high_resolution_clock::time_point t2 = high_resolution_clock::now();
   auto duration = duration_cast<microseconds>( t2 - t1 ).count();
   std::cout << "Time taken during computation: " << duration/1000000.0 << std::endl;
   return 0;
}



