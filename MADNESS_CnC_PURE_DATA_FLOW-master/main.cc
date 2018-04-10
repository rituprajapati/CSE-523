#include <cnc/cnc.h>
#include <utility>
#include <iostream>

#include <cmath>
#include <stdlib.h>

#include "Vector.h"
#include "Matrix.h"
#include "Node.h"
#include "Twoscalecoeffs.h"
#include "Quadrature.h"

const double pi = 3.1415926535897;
// http://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
#include <chrono>


using namespace std;
// http://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
using namespace std::chrono;


/*******************************************************************************************/
/* This program will have pure-data flow model of CnC which gives the potential that
/* the implementations (i.e., codes) inside the structs such as Project, Mult, Sub, etc
/* can be used in other computations as well.
/*******************************************************************************************/


/**********************************************************************/
/**********************************************************************/
/* This struct is used by different steps to put the 
/* output in item collection
/**********************************************************************/
/**********************************************************************/

template<typename K, typename V>
struct OutputTerminal {
   CnC::item_collection<K, V> *item_collection;
   std::vector<CnC::tag_collection<K> *> next_op_tags;

   /*----------------------------------------------------------------*/
   /* Output terminal constructor
   /*----------------------------------------------------------------*/
   OutputTerminal(CnC::item_collection<K, V> *item_collection,
                  std::vector<CnC::tag_collection<K> *> next_op_tags)
   : item_collection(item_collection), next_op_tags(next_op_tags) {}

   /*----------------------------------------------------------------*/
   /* This method is invoked to put item in the output_item_collection 
   /* as well as putting the tag to the next tag_collections
   /*----------------------------------------------------------------*/
   void put(K key, V value) const {
      item_collection->put(key, value);
      for (CnC::tag_collection<K> *next_op_tag: next_op_tags) {
         next_op_tag->put(key);
      }
   }

   /*----------------------------------------------------------------*/
   /* This method is invoked to put tag to the next tag_collections
   /* --> mainly used for recursion purposes
   /*----------------------------------------------------------------*/
   void put(K key) const {
      for (CnC::tag_collection<K> *next_op_tag: next_op_tags) {
         next_op_tag->put(key);
      }
   }
};


/**********************************************************************/
/**********************************************************************/
/*  Forward declaration
/**********************************************************************/
/**********************************************************************/

struct CnCContext;
double funcA(double x);
double funcB(double x);
Vector sub(const Vector &v1, const Vector &v2);
double sub_scale_factor(double n);
void init_twoscale(int );


struct step_Base {

public:
   std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals;
   std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals;

   step_Base( 
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          :
          input_terminals(input_terminals),
          output_terminals(output_terminals) {}
};

/**********************************************************************/
/**********************************************************************/
/* Struct  - Project
/* Used By - projectA_step, projectB_step   
/**********************************************************************/
/**********************************************************************/
struct Project : step_Base {

   double (*func)(double); 

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   Vector sValue(int n, int l, CnCContext &context) const;

  /*----------------------------------------------------------------*/
  /* Project Constructor
  /*----------------------------------------------------------------*/
   Project( double (*func)(double), 
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          :
          func(func),
          step_Base(input_terminals, output_terminals)
          {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - BinaryOp
/* Used By - subtract_1_step
/**********************************************************************/
/**********************************************************************/
struct BinaryOp : step_Base{


   using funcT = Vector (*)(const Vector &, const Vector&);
   funcT func;   

   double (*scale_factor)(double);

   Vector unfilter(Vector inputVector, int k, Matrix * hg) const;
   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   /*----------------------------------------------------------------*/
   /* BinaryOp Constructor
   /*----------------------------------------------------------------*/
   BinaryOp(  const funcT &func, 
              double (*scale_factor)(double), 
              std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
              std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
            : 
            func(func), 
            scale_factor(scale_factor), 
            step_Base( input_terminals, output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - Printer
/* Used By - printer_step
/**********************************************************************/
/**********************************************************************/
struct Printer : step_Base{

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   /*----------------------------------------------------------------*/
   /* Printer Constructor
   /*----------------------------------------------------------------*/
   Printer( std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
         step_Base(input_terminals, output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - Compress_Prolog
/* Used By - compress_prolog_FuncA_step, compress_prolog_FuncB_step
/**********************************************************************/
/**********************************************************************/
struct Compress_Prolog :step_Base{

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   /*----------------------------------------------------------------*/
   /* Compress_Prolog Constructor
   /*----------------------------------------------------------------*/
   Compress_Prolog( 
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
          step_Base(input_terminals, output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - Compress_doIt
/* Used By - compress_doIt_funcA_step, compress_doIt_funcB_step
/**********************************************************************/
/**********************************************************************/
struct Compress_doIt :step_Base{

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   /*----------------------------------------------------------------*/
   /* Compress_doIt Constructor
   /*----------------------------------------------------------------*/
   Compress_doIt( 
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
          step_Base(input_terminals, output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - GaxpyOp
/* Used By - gaxpyOp_step
/**********************************************************************/
/**********************************************************************/
struct GaxpyOp :step_Base{

   double alpha;
   double beta;

   int execute(const std::pair<int, int> &node, CnCContext &context) const;
   /*----------------------------------------------------------------*/
   /* GaxpyOp Constructor
   /*----------------------------------------------------------------*/
   GaxpyOp( 
            double alpha,
            double beta,
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
          alpha( alpha ),
          beta( beta),
          step_Base(input_terminals, output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - Reconstruct_Prolog
/* Used By - reconstruct_prolog_step
/**********************************************************************/
/**********************************************************************/
struct Reconstruct_Prolog : step_Base{

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   /*----------------------------------------------------------------*/
   /* Reconstruct_Prolog Constructor
   /*----------------------------------------------------------------*/
   Reconstruct_Prolog( 
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
         step_Base(input_terminals, output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - Reconstruct_doIt
/* Used By - reconstruct_doIt_step
/**********************************************************************/
/**********************************************************************/
struct Reconstruct_doIt : step_Base{

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   /*----------------------------------------------------------------*/
   /* Reconstruct_doIt Constructor
   /*----------------------------------------------------------------*/
   Reconstruct_doIt( 
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
          step_Base(input_terminals, output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - Norm2
/* Used By - norm2_step
/**********************************************************************/
/**********************************************************************/
struct Norm2 : step_Base{

   int execute(const std::pair<int, int> &node, CnCContext &context) const;
   /*----------------------------------------------------------------*/
   /* Norm2 Constructor
   /*----------------------------------------------------------------*/
   Norm2( 
          std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
          std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
          step_Base(input_terminals, output_terminals) {}
};



/**********************************************************************/
/**********************************************************************/
/* Struct  - Diff_Prolog
/* Used By - diff_prolog_step
/**********************************************************************/
/**********************************************************************/
struct Diff_Prolog :step_Base{

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   /*----------------------------------------------------------------*/
   /* Diff_Prolog Constructor
   /*----------------------------------------------------------------*/
   Diff_Prolog( 
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
          step_Base(input_terminals, output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - Diff_doIt
/* Used By - Diff_doIt_step
/**********************************************************************/
/**********************************************************************/
struct Diff_doIt : step_Base {

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   Vector unfilter(const Vector &inputVector, int k, const Matrix * hg) const;

  /*----------------------------------------------------------------*/
  /* Diff_doIt Constructor
  /*----------------------------------------------------------------*/
   Diff_doIt(  
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          :
          step_Base(input_terminals, output_terminals)
          {}
};

/**********************************************************************/
/**********************************************************************/
/* Struct  - InnerProduct
/* Used By - inner_step
/**********************************************************************/
/**********************************************************************/
// input_terminals - 2
// output terminals - 2, 1 with tag
// 
struct InnerProduct : step_Base {

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

  /*----------------------------------------------------------------*/
  /* InnerProduct Constructor
  /*----------------------------------------------------------------*/
   InnerProduct(  
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          :
          step_Base(input_terminals, output_terminals)
          {}
};



/**********************************************************************/
/**********************************************************************/
/* Struct  -  evaluate
/* Used By -  evaluate_step
/**********************************************************************/
/**********************************************************************/
// struct  Evaluate : step_Base{

//    int execute(const std::pair<int, int> &node, CnCContext &context) const;
//    double x;
//    /*----------------------------------------------------------------*/
//    /* Reconstruct_doIt Constructor
//    /*----------------------------------------------------------------*/
//    Evaluate( 
//             double x,
//             std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
//             std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
//           : 
//           x(x),
//           step_Base(input_terminals, output_terminals) {}
// };



/**********************************************************************/
/**********************************************************************/
/* Struct  - CnCContext
/**********************************************************************/
/**********************************************************************/
struct CnCContext : public CnC::context<CnCContext> {

public:
   int k, quad_npt, max_level;

   double thresh;

   double (*a_function)(double);
   double (*b_function)(double);

   Matrix *hg, *hg0, *hg1, *hgT;
   Matrix *rm, *r0, *rp;

   Vector *quad_w, *quad_x;
   Matrix *quad_phi, *quad_phiT, *quad_phiw;

   tbb::concurrent_vector<double> inner_results;
   tbb::concurrent_vector<double> norm2_results;
   CnC::item_collection<std::pair<int, int>, Node> evaluate_item;

   CnCContext(int k, double thresh, int max_level)
   : 
    CnC::context<CnCContext>(), 
    evaluate_item(*this),
    k(k), 
    thresh(thresh), 
    max_level(max_level)
    {
      // cout << "Inside CnCContext " << k << " " << thresh << " " << max_level << "\n";
      init_twoscale(k);
      init_quadrature(k);
      make_dc_periodic();
 
   }

   void init_twoscale(int k) {
      double  (*hgInput)[22] = twoscalecoeffs(k);

      hg = new Matrix(2*k, 2*k);
      hg0 = new Matrix(2*k, k);
      hg1 = new Matrix(2*k, k);
      hgT = new Matrix(2*k, 2*k);

      for (int i = 0; i < 2 * k; ++i) {
         for (int j = 0; j < 2 * k; ++j) {
            hg->set_item(i, j, hgInput[i][j]);
            hgT->set_item(i, j, hgInput[j][i]);
         }
      }

      for (int i = 0; i < 2 * k; ++i) {
         for (int j = 0; j < k; ++j) {
            hg0->set_item(i, j, hgInput[i][j]);
            hg1->set_item(i, j, hgInput[i][j+k]);
         }
      }
   }

   void init_quadrature(int order) {
      double *x = gauss_legendre_point(order);
      double *w = gauss_legendre_weight(order);

      quad_w = new Vector(w, 0, order);
      quad_x = new Vector(x, 0, order);

      int npt = order;
      quad_npt = npt;

      quad_phi = new Matrix(npt, k);
      quad_phiT = new Matrix(k, npt);
      quad_phiw = new Matrix(npt, k);

      for (int i = 0; i < npt; ++i) {
         double * p = phi((*quad_x)[i], k);
         for (int m = 0; m < k; ++m) {
            quad_phi->set_item(i, m, p[m]);
            quad_phiT->set_item(m, i, p[m]);
            quad_phiw->set_item(i, m, w[i] * p[m]);
         }
      }
   }


   void make_dc_periodic() {
      rm = new Matrix(k, k);
      r0 = new Matrix(k, k);
      rp = new Matrix(k, k);

      double iphase = 1.0;
      for (int i = 0; i < k; ++i) {
         double jphase = 1.0;

         for (int j = 0; j < k; ++j) {
            double gammaij = sqrt(( 2 * i + 1) * ( 2 * j + 1));
            double Kij;
            if ((( i -  j ) > 0) && (((i - j ) % 2) == 1 )) {
               Kij = 2.0;
            } else {
               Kij = 0.0;
            }

            r0->set_item(i, j, (0.5 * (1.0 - iphase * jphase - 2.0 * Kij) * gammaij));
            rm->set_item(i, j, (0.5 * jphase * gammaij));
            rp->set_item(i, j, (-0.5 * iphase * gammaij));

            jphase = -1 * jphase;
         }
         iphase = -1 * iphase;
      }
   }


  double __evaluate (  int n, int l, double x ) {

    Node nodeInfo;

    this->evaluate_item.get(make_pair(n,l), nodeInfo);
    int k = this->k;

    if( !nodeInfo.has_children ){
      double *p = phi(x, k);
      return ( nodeInfo.s.inner( *p) * sqrt( pow ( 2.0, n)) ) ; // doubt!! Inner needs a Vector bu p here is a pointer???
    } 
    else{
      n += 1; l *= 2; x = 2.0 * x;

      if ( x >= 1 ) {
        l = l+1; x = x-1;
      }
      return __evaluate(n, l, x );
    }
  } 

};


// This method is used to instantiate the general struct BinaryOp to be specifically subtraction
// of two mathematical functions
Vector sub(const Vector &v1, const Vector &v2) {
   Vector result(v1);
   for (unsigned int i = 0; i < v2.length(); ++i) {
      result.data[i] -= v2.data[i];
   }
   return result;
}

// used in the method execute of class BinaryOp for the subtraction instance 
double sub_scale_factor(double n) {
   return 1.0;
}

// Mathematical test functions test1 and test2 would use the following guassian method
double gaussian(double x, double a, double coeff) {
    return coeff*exp(-a*x*x);
}

// Mathematical test function test1
double funcA(double x) {
    static const int N = 100;
    static double a[N], X[N], c[N];
    static bool initialized = false;

    if (!initialized) {
        for (int i=0; i<N; i++) {
            a[i] = 1000*drand48();
            X[i] = drand48();
            c[i] = pow(2*a[i]/M_PI,0.25);
        }
        initialized = true;
    }

    double sum = 0.0;
    for (int i=0; i<N; i++) sum += gaussian(x-X[i], a[i], c[i]);
    return sum;
}

// Mathematical test function test2
double funcB(double x) {
    static const int N = 100;
    static double a[N], X[N], c[N];
    static bool initialized = false;

    if (!initialized) {
        for (int i=0; i<N; i++) {
            a[i] = 1000*drand48();
            X[i] = drand48();
            c[i] = pow(2*a[i]/M_PI,0.25);
        }
        initialized = true;
    }

    double sum = 0.0;
    for (int i=0; i<N; i++) sum += gaussian(x-X[i], a[i], c[i]);
    return sum;
}

Vector Project::sValue(int n, int l, CnCContext &context) const {
   Vector s(context.k);
   // cout << "inside Project::sValue\n";
   Vector &quad_x_ref = *(context.quad_x);
   Matrix &quad_phiw_ref = *(context.quad_phiw);

   double h = pow(0.5, n);
   double scale = sqrt(h);
   for (int mu = 0; mu < context.quad_npt; ++mu) {
      double x = (l + quad_x_ref[mu]) * h;
      double fValue = func(x);
      for (int i = 0; i < (context.k); ++i) {
         s[i] = s[i] + (scale * fValue * (quad_phiw_ref.get_item(mu, i)));
      }
   }
   return s;
}

int Project::execute(const std::pair<int, int> &node, CnCContext &context) const {

   // cout << "inside Project::execute\n";

   Vector s0 = sValue(node.first + 1, 2 * node.second, context);
   Vector s1 = sValue(node.first + 1, 2 * node.second + 1, context);

   int k = context.k;
   Vector s(s0 | s1); // concatenation of s0 and s1
   Vector d(s * (*(context.hgT)));


   // Node(int n, int l, int k, const Vector &s_input, const Vector &d_input, bool has_children)

   if (d.normf(k, 2 * k) < context.thresh || node.first >= context.max_level - 1) {
      output_terminals[1].put(node, Node(node.first, node.second, k, Vector(), Vector(), true));
      output_terminals[1].put(std::make_pair(node.first + 1, 2 * node.second), Node(node.first + 1, 2 * node.second, k, s0, Vector(), false));
      output_terminals[1].put(std::make_pair(node.first + 1, 2 * node.second + 1), Node(node.first + 1, 2 * node.second + 1, k, s1, Vector(), false));
   }
   else {
      output_terminals[1].put(node, Node(node.first, node.second, k, Vector(), Vector(), true));
      output_terminals[0].put(std::make_pair(node.first + 1, 2 * node.second));
      output_terminals[0].put(std::make_pair(node.first + 1, 2 * node.second + 1));   
   }
   return CnC::CNC_Success;
}

Vector BinaryOp::unfilter(Vector inputVector, int k, Matrix * hg) const {
   Vector inputVector_copy(inputVector);
   Vector vector_d(2 * k);
   vector_d.set_slice_from_another_vector(0, k, inputVector);
   Vector vector_s2 = (vector_d * (*hg));
   return vector_s2;
}

int BinaryOp::execute(const std::pair<int, int> &node, CnCContext &context) const {
   Node left;
   Node right;

   int k = context.k;

   input_terminals[0]->get(node, left);
   input_terminals[1]->get(node, right);

   if (left.s.length() != 0 && right.s.length() != 0) { // If both of them are at the leaf level
      double scale_fact = scale_factor(node.first);
      Vector f_vector(left.s * (*(context.quad_phiT)));
      Vector g_vector(right.s * (*(context.quad_phiT)));

      Vector temp = func(f_vector, g_vector);
      Vector resultVector((temp * (*context.quad_phiw)).scale(scale_fact));
      
      output_terminals[1].put(node, Node(node.first, node.second, k, resultVector, Vector(), false));
   }
   else {
      if (left.s.length() != 0) {
         Vector left_unfiltered = unfilter(left.s, k, context.hg);
         output_terminals[0].put(std::make_pair(node.first + 1, 2 * node.second), Node(node.first + 1, 2 * node.second, k, left_unfiltered.get_slice(0, k), Vector(), false));
         output_terminals[0].put(std::make_pair(node.first + 1, 2 * node.second + 1), Node(node.first + 1, 2 * node.second + 1, k, left_unfiltered.get_slice(k, 2 * k), Vector(), false));
      }
      else if (right.s.length() != 0) {
         Vector right_unfiltered = unfilter(right.s, k, context.hg);
         output_terminals[2].put(std::make_pair(node.first + 1, 2 * node.second), Node(node.first + 1, 2 * node.second, k, right_unfiltered.get_slice(0, k), Vector(), false));
         output_terminals[2].put(std::make_pair(node.first + 1, 2 * node.second + 1), Node(node.first + 1, 2 * node.second + 1, k, right_unfiltered.get_slice(k, 2 * k), Vector(), false));
      }
      output_terminals[1].put(node, Node(node.first, node.second, k, Vector(), Vector(), true));
   }
   
   return CnC::CNC_Success;
}


int Compress_Prolog::execute( const std::pair<int, int> &node, CnCContext &context ) const {

    Node input;
    input_terminals[0]->get(node, input);
    int k = context.k;

    if ( !input.has_children ) { //if the node is a leaf

      if( node.first == 0){
        output_terminals[1].put( node, input);
      }
      else{

        output_terminals[1].put( node, Node( node.first, node.second, k, Vector(), Vector(k), false ));

        if( node.second & 0x1uL)
           output_terminals[2].put( std::make_pair( node.first-1, node.second/2), input);
        else
           output_terminals[0].put( std::make_pair( node.first-1, node.second/2), input);
      }

    }
    return CnC::CNC_Success;
}


int Compress_doIt::execute( const std::pair<int, int> &node, CnCContext &context ) const {


    Node left;
    Node right;
    int k = context.k;

    input_terminals[0]->get(node, left);
    input_terminals[1]->get(node, right);

    Vector s( left.s | right.s );
    Vector d(s * (*(context.hgT)));

    Vector sValue(d.data, 0, k);
    Vector dValue(d.data, k, 2 * k);

    if( node.first == 0){
      output_terminals[1].put( node, Node( node.first, node.second, k, sValue, dValue, true) );
    }
    else{
      output_terminals[1].put( node, Node( node.first, node.second, k, Vector(), dValue, true) );

      if( node.second & 0x1uL){
        output_terminals[2].put( std::make_pair( node.first-1, node.second/2), Node( node.first-1, node.second/2, k, sValue, Vector(), false ));
      }
      else{
        output_terminals[0].put( std::make_pair( node.first-1, node.second/2), Node( node.first-1, node.second/2, k, sValue, Vector(), false ));
      }
    }

    return CnC::CNC_Success;
}



int GaxpyOp::execute( const std::pair<int, int> &node, CnCContext &context ) const {

    Node left;
    Node right;

    int k = context.k;

    input_terminals[0]->get(node, left);
    input_terminals[1]->get(node, right);

    Vector tempD(left.d);
    tempD.gaxpy( alpha, right.d, beta);

    Vector tempS;

    if( node.first == 0 && node.second == 0){
      tempS = left.s;
      tempS.gaxpy( alpha, right.s, beta);
    }

    output_terminals[1].put( node, Node( node.first, node.second, k, tempS, tempD, left.has_children || right.has_children ) );

    if( left.has_children && !right.has_children ){

      output_terminals[2].put( make_pair(node.first + 1, node.second * 2), Node(node.first+1, node.second * 2, k, Vector(), Vector(k), false));
      output_terminals[2].put( make_pair(node.first + 1, node.second * 2 + 1), Node(node.first+1, node.second*2+1, k, Vector(), Vector(k), false));

    }

    if( !left.has_children && right.has_children ){

      output_terminals[0].put( make_pair(node.first + 1, node.second * 2), Node(node.first+1, node.second*2, k, Vector(), Vector(k), false));
      output_terminals[0].put( make_pair(node.first + 1, node.second * 2 + 1), Node(node.first+1, node.second*2+1, k, Vector(), Vector(k), false));

    }

    return CnC::CNC_Success;
}


int Reconstruct_Prolog::execute( const std::pair<int, int> &node, CnCContext &context ) const {

    Node input;

    input_terminals[0]->get( node, input);

    if( node.first == 0) {
        output_terminals[0].put( node, input);
    }
    return CnC::CNC_Success;
}


int Reconstruct_doIt::execute( const std::pair<int, int> &node, CnCContext &context ) const {

    Node s_coeff;
    Node node_information;

    int k = context.k;
    input_terminals[0]->get( node, s_coeff);
    input_terminals[1]->get( node, node_information);

    Vector s = s_coeff.s;

    if( node_information.has_children ){

        Vector v1( s| node_information.d );
        Vector v2( v1 * (*context.hg) );

        Vector leftChildS(v2.data, 0, k);
        Vector rightChildS(v2.data, k, 2 * k);

        output_terminals[0].put(make_pair(node.first + 1, node.second * 2), Node( node.first + 1, node.second * 2, k, leftChildS, Vector(), false));
        output_terminals[0].put(make_pair(node.first + 1, node.second * 2 + 1), Node( node.first + 1, node.second * 2 + 1, k, rightChildS, Vector(), false));

        output_terminals[1].put( node, Node( node.first, node.second, k, Vector(), Vector(), true));
    }
    else{
        output_terminals[1].put( node, Node( node.first, node.second, k, s, Vector(), true));
    }

    return CnC::CNC_Success;
}


int Printer::execute(const std::pair<int, int> &node, CnCContext &context) const {
   Node nodeInfo;
   input_terminals[0]->get(node, nodeInfo);
   std::cout << "Printer:: Node with info: (Key: (" << node.first << ", " << node.second << "), " << nodeInfo.toString() << ")" << std::endl;
   return CnC::CNC_Success;
}



int Norm2::execute( const std::pair<int, int> &node, CnCContext &context ) const {

   Node nodeInfo;
   input_terminals[0]->get(node, nodeInfo);

   context.norm2_results.push_back( pow( nodeInfo.s.normf(), 2) );

   if( nodeInfo.has_children ){

      output_terminals[0].put( make_pair(node.first + 1, node.second * 2), Node());
      output_terminals[0].put( make_pair(node.first + 1, node.second * 2 + 1), Node());

   }
   return CnC::CNC_Success;
}



int Diff_Prolog::execute( const std::pair<int, int> &node, CnCContext &context ) const {

    Node nodeInfo;
    // cout << "Diff_Prolog::execute\n";

    input_terminals[0]->get(node, nodeInfo);

    output_terminals[0].put( make_pair( node.first, node.second==0ul ? (1ul<<node.first)-1 : node.second-1), nodeInfo);
    output_terminals[1].put( node, nodeInfo);
    output_terminals[2].put( make_pair( node.first, node.second==((1ul<<node.first)-1) ? 0 : node.second+1), nodeInfo);

   return CnC::CNC_Success;
}


Vector Diff_doIt::unfilter(const Vector &inputVector, int k, const Matrix * hg) const {
   
  // cout << "inside Diff_doIt::unfilter\n";
  Vector vector_d(2 * k);
  vector_d.set_slice_from_another_vector(0, k, inputVector);
  Vector vector_s2 = (vector_d * (*hg));
  return vector_s2;
}


int Diff_doIt::execute( const std::pair<int, int> &node, CnCContext &context ) const {

  Node left, center, right;

  // cout << "inside Diff_doIt::execute\n";

  input_terminals[0]->get(node, left);
  input_terminals[1]->get(node, center);
  input_terminals[2]->get(node, right);

  int k = context.k;

  if ((left.s.length() != 0 && left.has_children) || (left.s.length() == 0 && !left.has_children)) {
      std::cout << "ERROR at left" << std::endl;
  }

  if ((center.s.length() != 0 && center.has_children)|| (center.s.length() == 0 && !center.has_children)) {
      std::cout << "ERROR at center" << std::endl;
  }

  if ((right.s.length() != 0 && right.has_children) ||(right.s.length() == 0 && !right.has_children)) {
      std::cout << "ERROR at right" << std::endl;
  }

  if (left.s.length() != 0 && center.s.length() != 0 && right.s.length() != 0) {
     Vector r = ((*context.rp) * left.s) + ((*context.r0) * center.s) + ((*context.rm) * right.s);
     output_terminals[3].put(node, Node( node.first, node.second, k, r.scale( pow(2.0, node.first)), Vector(), false));
  }
  else {

     output_terminals[3].put(node, Node( node.first, node.second, k, Vector(), Vector(), true));

     if (left.s.length() != 0) {   
        Vector unfiltered = unfilter(left.s, k, context.hg);
        output_terminals[0].put(make_pair( node.first+1, node.second*2), Node( left.n +1, left.l*2 +1, k, unfiltered.get_slice(k, 2 * k), Vector(), false));
     }
       
     if (center.s.length() != 0) {

       Vector unfiltered = unfilter(center.s, k, context.hg);
 
       output_terminals[2].put(make_pair( node.first+1, node.second*2), Node( node.first +1, node.second*2 +1, k, unfiltered.get_slice(k, 2 * k), Vector(), false));
       output_terminals[0].put(make_pair( node.first+1, node.second*2+1), Node( node.first +1, node.second*2, k, unfiltered.get_slice(0,k), Vector(), false));
       output_terminals[1].put(make_pair( node.first+1, node.second*2), Node( node.first +1, node.second*2, k, unfiltered.get_slice(0,k), Vector(), false));
       output_terminals[1].put(make_pair( node.first+1, node.second*2+1), Node( node.first +1, node.second*2+1, k, unfiltered.get_slice(k,2*k), Vector(), false));
     }

     if (right.s.length() != 0) {
        Vector unfiltered = unfilter(right.s, k, context.hg);
        output_terminals[2].put(make_pair( node.first+1, node.second*2+1), Node( right.n +1, right.l*2 , k, unfiltered.get_slice(0,k), Vector(), false));
     }
  }
   
  return CnC::CNC_Success;
}



int InnerProduct::execute( const std::pair<int, int> &node, CnCContext &context ) const{

   Node left;
   Node right;

   input_terminals[0]->get(node, left);
   input_terminals[1]->get(node, right);

   int k = context.k;

   if (left.has_children && right.has_children) {


    context.inner_results.push_back(left.d.inner(right.d));

    if (node.first == 0) { // It is the root
      context.inner_results.push_back(left.s.inner(right.s));
    }

    //Trigger the tags for children
    output_terminals[0].put( make_pair(node.first + 1, node.second * 2), Node());
    output_terminals[0].put( make_pair(node.first + 1, node.second * 2 + 1), Node());

   }
}




/*****************************************************************/
/* gaussian with square normalized to 1              */
/*****************************************************************/

double test1 ( double x ){
  
  double a = 500.0;
  double tmp = 2 * a / pi;
  double tmp2 = -a * pow(( x-0.5 ), 2);
  return pow( tmp , 0.25 ) * exp( tmp2 );

}

/*****************************************************************/
/* superposition of multiple gaussians               */
/*****************************************************************/
double test2 (  double x ) { 

  return ( test1( x-0.3) + test1 ( x ) + test1 ( x+0.3 ) );
}

/*****************************************************************/
/* a more interesting (singular and oscillating) function      */
/* ... note it is never computed exactly at the singularity    */
/* except by the test code below.                  */
/*****************************************************************/

double test3 ( double x ){
  double a = 100.0 * pi;

  if ( x == 0.5 ){
    return 0.0;
  }
  else{
    return cos ( a * ( x-0.5) ) / sqrt ( abs ( x-0.5 ));
  }
}

/*****************************************************************/
/* derivative of test1                       */
/*****************************************************************/

double dtest1 ( double x ) {

  double a = 500.0;
  double tmp1 =  -2.0 * a * (x-0.5);
  double tmp2 =  pow( 2 * a /pi, 0.25 );
  double tmp3 =  exp( -a * pow( (x-0.5), 2 ));
  return tmp1 * tmp2 * tmp3;
}


/*****************************************************************/
/* derivative of test2                       */
/*****************************************************************/

double dtest2 ( double x ){

  return  ( dtest1( x-0.3 ) + dtest1( x )  + dtest1( x+0.3) );
}


/*****************************************************************/
/* derivative of test3                       */
/*****************************************************************/

double dtest3 ( double x ){

  double a = 100.0 * pi;

  if ( x == 0.5 ){
    return 0.0;
  }
  else{

    double s = 1.0;
    if ( x < 0.5 ){
      s = -1.0;
    }
    double tmp1 = -a * sin( a * (x-0.5) );
    double tmp2 = sqrt( abs(x-0.5) );
    double tmp3 = s * 0.5 * cos( a * (x-0.5) );
    double tmp4 = pow( abs( x-0.5 ), 1.5 );

    return ( tmp1 / tmp2 - tmp3 / tmp4 );
  }
} 

double ( *test[3] ) ( double x ) = { test1, test2, test3 };
double ( *dtest[3] ) ( double x ) = { dtest1, dtest2, dtest3 };





struct diff_test: CnCContext{

  /*----------------------------------------------------------------*/
  /* Item Collections
  /*----------------------------------------------------------------*/
   CnC::item_collection<std::pair<int, int>, Node> project_item;

   CnC::item_collection<std::pair<int, int>, Node> diff_left_item;
   CnC::item_collection<std::pair<int, int>, Node> diff_right_item;
   CnC::item_collection<std::pair<int, int>, Node> diff_centre_item;

   CnC::item_collection<std::pair<int, int>, Node> diff_result_item;

   CnC::tag_collection<std::pair<int, int>> project_tag;
   CnC::tag_collection<std::pair<int, int>> diff_prolog_tag;
   CnC::tag_collection<std::pair<int, int>> diff_doIt_tag;
   CnC::tag_collection<std::pair<int, int>> printer_tag;

   using OutputTerminalType = OutputTerminal<std::pair<int, int>, Node>;

   CnC::step_collection<Project> project_step;
   CnC::step_collection<Diff_Prolog> diff_prolog_step;
   CnC::step_collection<Diff_doIt> diff_doIt_step;
   CnC::step_collection<Printer> printer_step;

   diff_test( int k, double thresh, int max_level )
    :
    CnCContext( k, thresh, max_level),
    project_item(*this),
    diff_left_item(*this),
    diff_right_item(*this),
    diff_centre_item(*this),
    diff_result_item(*this),
    project_tag(*this),
    diff_prolog_tag(*this),
    diff_doIt_tag(*this),
    printer_tag(*this),

    /*----------------------------------------------------------------*/
    /* Declare project_step
    /*----------------------------------------------------------------*/

    project_step(
                  *this, 
                  "project_step", 
                  Project(
                          &funcA, 
                          std::vector<CnC::item_collection<std::pair<int, int>, Node> *>{},
                          std::vector<OutputTerminalType> {
                              OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&project_tag}),
                              OutputTerminalType(&project_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&diff_prolog_tag})}
                          )
                  ),

    /*----------------------------------------------------------------*/
    /* Declare diff_prolog_step
    /*----------------------------------------------------------------*/

    diff_prolog_step(
                  *this, 
                  "diff_prolog_step", 
                  Diff_Prolog(
                          std::vector<CnC::item_collection<std::pair<int, int>, Node> *>{ &project_item },
                          std::vector<OutputTerminalType> {
                              OutputTerminalType(&diff_left_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {}),
                              OutputTerminalType(&diff_centre_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&diff_doIt_tag}),
                              OutputTerminalType(&diff_right_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {})}
                          )
                  ),

    /*----------------------------------------------------------------*/
    /* Declare diff_doIt_step
    /*----------------------------------------------------------------*/

    diff_doIt_step(
                  *this, 
                  "diff_doIt_step", 
                  Diff_doIt(
                          std::vector<CnC::item_collection<std::pair<int, int>, Node> *>{ &diff_left_item, &diff_centre_item, &diff_right_item },
                          std::vector<OutputTerminalType> {
                              OutputTerminalType(&diff_left_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {}),
                              OutputTerminalType(&diff_centre_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&diff_doIt_tag}),
                              OutputTerminalType(&diff_right_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {}),
                              OutputTerminalType(&diff_result_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&printer_tag})}
                          )
                  ),
 
    /*----------------------------------------------------------------*/
    /* Declare printer_step
    /*----------------------------------------------------------------*/
    printer_step(
                *this, 
                "printer_step", 
                Printer( 
                  std::vector<CnC::item_collection<std::pair<int, int>, Node> *>{&diff_result_item}, 
                  std::vector<OutputTerminalType>{})
                )
    {

      // cout << "Inside diff_test constructor\n";


      /*----------------------------------------------------------------*/
      /* Tag Prescription
      /*----------------------------------------------------------------*/
      project_tag.prescribes(project_step, *this);
      diff_prolog_tag.prescribes(diff_prolog_step, *this);
      diff_doIt_tag.prescribes(diff_doIt_step, *this);
      printer_tag.prescribes(printer_step, *this);

      /*----------------------------------------------------------------*/
      /* Steps Produce Consume
      /*----------------------------------------------------------------*/

      project_step.produces(project_item);

      diff_prolog_step.consumes(project_item);
      diff_prolog_step.produces(diff_left_item);
      diff_prolog_step.produces(diff_centre_item);
      diff_prolog_step.produces(diff_right_item);

      diff_doIt_step.consumes(diff_left_item);
      diff_doIt_step.consumes(diff_centre_item);
      diff_doIt_step.consumes(diff_right_item);

      diff_doIt_step.produces(diff_left_item);
      diff_doIt_step.produces(diff_centre_item);
      diff_doIt_step.produces(diff_right_item);
      diff_doIt_step.produces(diff_result_item);

      printer_step.consumes(diff_result_item);

    }

};




struct addition_test : CnCContext {
   
  /*----------------------------------------------------------------*/
  /* Item Collections
  /*----------------------------------------------------------------*/
   CnC::item_collection<std::pair<int, int>, Node> projectA_item;
   CnC::item_collection<std::pair<int, int>, Node> projectB_item; 
  

   CnC::item_collection<std::pair<int, int>, Node> compress_prologA_left_item;
   CnC::item_collection<std::pair<int, int>, Node> compress_prologA_right_item;
   CnC::item_collection<std::pair<int, int>, Node> funcA_coeff_compressed_item;
   CnC::item_collection<std::pair<int, int>, Node> compress_prologB_left_item;
   CnC::item_collection<std::pair<int, int>, Node> compress_prologB_right_item;
   CnC::item_collection<std::pair<int, int>, Node> funcB_coeff_compressed_item;
   CnC::item_collection<std::pair<int, int>, Node> gaxpy_result_item;
   CnC::item_collection<std::pair<int, int>, Node> reconstruct_result_item;
   CnC::item_collection<std::pair<int, int>, Node> s_coeff_item;


  /*----------------------------------------------------------------*/
  /* Tag Collections
  /*----------------------------------------------------------------*/
   CnC::tag_collection<std::pair<int, int>> projectA_tag;
   CnC::tag_collection<std::pair<int, int>> projectB_tag;
   CnC::tag_collection<std::pair<int, int>> printer_tag;


   CnC::tag_collection<std::pair<int, int>> compress_prolog_FuncA_tag;
   CnC::tag_collection<std::pair<int, int>> compress_prolog_FuncB_tag;
   CnC::tag_collection<std::pair<int, int>> compress_doIt_funcA_tag;
   CnC::tag_collection<std::pair<int, int>> compress_doIt_funcB_tag;
   CnC::tag_collection<std::pair<int, int>> gaxpyOP_tag;
   CnC::tag_collection<std::pair<int, int>> reconstruct_prolog_tag;
   CnC::tag_collection<std::pair<int, int>> reconstruct_doIt_tag;

  /*----------------------------------------------------------------*/
  /* Step Collections
  /*----------------------------------------------------------------*/
   using OutputTerminalType = OutputTerminal<std::pair<int, int>, Node>;

   CnC::step_collection<Project> projectA_step;
   CnC::step_collection<Project> projectB_step;
   CnC::step_collection<Printer> printer_step;

   CnC::step_collection<Compress_Prolog> compress_prolog_FuncA_step;
   CnC::step_collection<Compress_Prolog> compress_prolog_FuncB_step;
   CnC::step_collection<Compress_doIt> compress_doIt_funcA_step;
   CnC::step_collection<Compress_doIt> compress_doIt_funcB_step;
   CnC::step_collection<GaxpyOp> gaxpyOp_step;
   CnC::step_collection<Reconstruct_Prolog> reconstruct_prolog_step;
   CnC::step_collection<Reconstruct_doIt>  reconstruct_doIt_step;

  /*----------------------------------------------------------------*/
  /* addition_test Constructor
  /*----------------------------------------------------------------*/
   addition_test(int k, double thresh, int max_level)
   : 
     CnCContext( k, thresh, max_level),
     projectA_item(*this),
     projectB_item(*this), 
     projectA_tag(*this), 
     projectB_tag(*this), 
     printer_tag(*this), 
     
     compress_prologA_left_item(*this),
     compress_prologA_right_item(*this),
     funcA_coeff_compressed_item(*this),
     compress_prologB_left_item(*this),
     compress_prologB_right_item(*this),
     funcB_coeff_compressed_item(*this),
     gaxpy_result_item(*this),
     reconstruct_result_item(*this),
     s_coeff_item(*this),

     compress_prolog_FuncA_tag(*this),
     compress_prolog_FuncB_tag(*this),
     compress_doIt_funcA_tag(*this),
     compress_doIt_funcB_tag(*this),
     gaxpyOP_tag(*this),
     reconstruct_prolog_tag(*this),
     reconstruct_doIt_tag(*this),

     /*----------------------------------------------------------------*/
     /* Declare projectA_step
     /*----------------------------------------------------------------*/
     projectA_step(
                  *this, 
                  "projectA_step", 
                  Project(
                          &funcA, std::vector<CnC::item_collection<std::pair<int, int>, Node> *>{},
                          std::vector<OutputTerminalType> {
                              OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&projectA_tag,}),
                              OutputTerminalType(&projectA_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&compress_prolog_FuncA_tag})}
                          )
                  ),
 
     /*----------------------------------------------------------------*/
     /* Declare projectB_step
     /*----------------------------------------------------------------*/
     projectB_step(
                  *this, 
                  "projectB_step", 
                  Project(
                          &funcB, 
                          std::vector<CnC::item_collection<std::pair<int, int>, Node> *>{},
                          std::vector<OutputTerminalType> {
                              OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&projectB_tag}),
                              OutputTerminalType(&projectB_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&compress_prolog_FuncB_tag})}
                          )
                  ),
    
  
    /*----------------------------------------------------------------*/
    /* Declare printer_step
    /*----------------------------------------------------------------*/
    printer_step(
                *this, 
                "printer_step", 
                Printer( 
                  std::vector<CnC::item_collection<std::pair<int, int>, Node> *>{&reconstruct_result_item}, 
                  std::vector<OutputTerminalType>{})
                ),

    /*----------------------------------------------------------------*/
    /* Declare compress_prolog_FuncA_step
    /*----------------------------------------------------------------*/
    compress_prolog_FuncA_step ( 
                                *this, 
                                "compress_prolog_FuncA_step", 
                                Compress_Prolog( 
                                    std::vector<CnC::item_collection<std::pair<int, int>, Node> *> {&projectA_item}, 
                                    std::vector<OutputTerminalType>{
                                    OutputTerminalType(&compress_prologA_left_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&compress_doIt_funcA_tag}),
                                    OutputTerminalType(&funcA_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> { &gaxpyOP_tag}),
                                    OutputTerminalType(&compress_prologA_right_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {})})
                              ),


    /*----------------------------------------------------------------*/
    /* Declare compress_prolog_FuncB_step
    /*----------------------------------------------------------------*/
    compress_prolog_FuncB_step ( 
                                *this, 
                                "compress_prolog_FuncB_step", 
                                Compress_Prolog( 
                                    std::vector<CnC::item_collection<std::pair<int, int>, Node> *> {&projectB_item}, 
                                    std::vector<OutputTerminalType>{
                                    OutputTerminalType(&compress_prologB_left_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&compress_doIt_funcB_tag}),
                                    OutputTerminalType(&funcB_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {}),
                                    OutputTerminalType(&compress_prologB_right_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {})})
                              ),

    /*----------------------------------------------------------------*/
    /* Declare compress_doIt_funcA_step
    /*----------------------------------------------------------------*/
    compress_doIt_funcA_step ( 
                                *this, 
                                "compress_doIt_funcA_step", 
                                Compress_doIt( 
                                    std::vector<CnC::item_collection<std::pair<int, int>, Node> *> {&compress_prologA_left_item, &compress_prologA_right_item}, 
                                    std::vector<OutputTerminalType>{
                                    OutputTerminalType(&compress_prologA_left_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&compress_doIt_funcA_tag}),
                                    OutputTerminalType(&funcA_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&gaxpyOP_tag}),
                                    OutputTerminalType(&compress_prologA_right_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {})})
                              ),


    /*----------------------------------------------------------------*/
    /* Declare compress_doIt_funcB_step
    /*----------------------------------------------------------------*/
    compress_doIt_funcB_step ( 
                                *this, 
                                "compress_doIt_funcB_step", 
                                Compress_doIt( 
                                    std::vector<CnC::item_collection<std::pair<int, int>, Node> *> {&compress_prologB_left_item, &compress_prologB_right_item}, 
                                    std::vector<OutputTerminalType>{
                                    OutputTerminalType(&compress_prologB_left_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&compress_doIt_funcB_tag}),
                                    OutputTerminalType(&funcB_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {}),
                                    OutputTerminalType(&compress_prologB_right_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {})})
                              ),


    /*----------------------------------------------------------------*/
    /* Declare gaxpyOp_step
    /*----------------------------------------------------------------*/
    gaxpyOp_step ( 
                  *this, 
                  "gaxpyOp_step", 
                  GaxpyOp( 
                      1.0,
                      1.0,
                      std::vector<CnC::item_collection<std::pair<int, int>, Node> *> {&funcA_coeff_compressed_item, &funcB_coeff_compressed_item}, 
                      std::vector<OutputTerminalType>{
                      OutputTerminalType(&funcA_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&gaxpyOP_tag}),
                      OutputTerminalType(&gaxpy_result_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&reconstruct_prolog_tag}),
                      OutputTerminalType(&funcB_coeff_compressed_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {})})
                ),

    /*----------------------------------------------------------------*/
    /* Declare reconstruct_prolog_step
    /*----------------------------------------------------------------*/
    reconstruct_prolog_step ( 
                            *this, 
                            "reconstruct_prolog_step", 
                            Reconstruct_Prolog( 
                                std::vector<CnC::item_collection<std::pair<int, int>, Node> *> {&gaxpy_result_item}, 
                                std::vector<OutputTerminalType>{
                                OutputTerminalType(&s_coeff_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&reconstruct_doIt_tag})})
                          ),

  
    /*----------------------------------------------------------------*/
    /* Declare reconstruct_doIt_step
    /*----------------------------------------------------------------*/
    reconstruct_doIt_step( 
                          *this, 
                          "reconstruct_doIt_step", 
                          Reconstruct_doIt( 
                              std::vector<CnC::item_collection<std::pair<int, int>, Node> *> {&s_coeff_item, &gaxpy_result_item}, 
                              std::vector<OutputTerminalType>{
                                OutputTerminalType(&s_coeff_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&reconstruct_doIt_tag}),
                                OutputTerminalType(&reconstruct_result_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&printer_tag})})
                         )

    {
      
      /*----------------------------------------------------------------*/
      /* Tag Prescription
      /*----------------------------------------------------------------*/
      projectA_tag.prescribes(projectA_step, *this);
      projectB_tag.prescribes(projectB_step, *this);
      printer_tag.prescribes(printer_step, *this);
      compress_prolog_FuncA_tag.prescribes( compress_prolog_FuncA_step, *this);
      compress_prolog_FuncB_tag.prescribes( compress_prolog_FuncB_step, *this);
      compress_doIt_funcA_tag.prescribes( compress_doIt_funcA_step, *this);
      compress_doIt_funcB_tag.prescribes( compress_doIt_funcB_step, *this);
      gaxpyOP_tag.prescribes( gaxpyOp_step, *this);
      reconstruct_prolog_tag.prescribes( reconstruct_prolog_step, *this);
      reconstruct_doIt_tag.prescribes( reconstruct_doIt_step, *this);

      /*----------------------------------------------------------------*/
      /* Steps Produce Consume
      /*----------------------------------------------------------------*/
      projectA_step.produces(projectA_item);

      projectB_step.produces(projectB_item);

      compress_prolog_FuncA_step.consumes( projectA_item);
      compress_prolog_FuncA_step.produces( compress_prologA_left_item );
      compress_prolog_FuncA_step.produces( compress_prologA_right_item );
      compress_prolog_FuncA_step.produces( funcA_coeff_compressed_item );
      
      compress_prolog_FuncB_step.consumes( projectB_item);
      compress_prolog_FuncB_step.produces( compress_prologB_left_item );
      compress_prolog_FuncB_step.produces( compress_prologB_right_item );
      compress_prolog_FuncB_step.produces( funcB_coeff_compressed_item );

      compress_doIt_funcA_step.consumes( compress_prologA_left_item );
      compress_doIt_funcA_step.consumes( compress_prologA_right_item );
      compress_doIt_funcA_step.produces( funcA_coeff_compressed_item);
      compress_doIt_funcA_step.produces( compress_prologA_left_item);
      compress_doIt_funcA_step.produces( compress_prologA_right_item);

      compress_doIt_funcB_step.consumes( compress_prologB_left_item );
      compress_doIt_funcB_step.consumes( compress_prologB_right_item );
      compress_doIt_funcB_step.produces( funcB_coeff_compressed_item);
      compress_doIt_funcB_step.produces( compress_prologB_left_item );
      compress_doIt_funcB_step.produces( compress_prologB_right_item );

      gaxpyOp_step.consumes( funcA_coeff_compressed_item);
      gaxpyOp_step.consumes( funcB_coeff_compressed_item);
      gaxpyOp_step.produces( funcA_coeff_compressed_item);
      gaxpyOp_step.produces( funcB_coeff_compressed_item);
      gaxpyOp_step.produces( gaxpy_result_item);

      reconstruct_prolog_step.consumes( gaxpy_result_item );
      reconstruct_prolog_step.produces( s_coeff_item);

      reconstruct_doIt_step.consumes( s_coeff_item);
      reconstruct_doIt_step.consumes( gaxpy_result_item);
      reconstruct_doIt_step.produces( s_coeff_item);
      reconstruct_doIt_step.produces( reconstruct_result_item);

      printer_step.consumes(reconstruct_result_item);

   }

};


int main(int argc, char *argv[]) {
   int k = 5;
   int npt = 20;
   int max_level = 30;
   double thresh = atof(argv[2]);

   high_resolution_clock::time_point t1 = high_resolution_clock::now();


   //Addition test
   addition_test add_obj(k, thresh, max_level);
   add_obj.projectA_tag.put(std::make_pair(0, 0));
   add_obj.projectB_tag.put(std::make_pair(0, 0));   
   add_obj.wait();

   //Diff test
   // diff_test diff_test_obj( k, thresh, max_level);
   // diff_test_obj.project_tag.put( std::make_pair(0, 0) );
   // diff_test_obj.wait();


   high_resolution_clock::time_point t2 = high_resolution_clock::now();
   auto duration = duration_cast<microseconds>( t2 - t1 ).count();
   std::cout << duration/1000000.0 << std::endl;


   return 0;
}
