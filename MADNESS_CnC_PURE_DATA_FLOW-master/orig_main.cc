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

// http://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
#include <chrono>


using namespace std;
// http://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
using namespace std::chrono;


// This program will have pure-data flow model of CnC which gives the potential that
// the implementations (i.e., codes) inside the structs such as Project, Mult, Sub, etc
// can be used in other computations as well.


template<typename K, typename V>
struct OutputTerminal {
   CnC::item_collection<K, V> *item_collection;
   std::vector<CnC::tag_collection<K> *> next_op_tags;

   // Constructor
   OutputTerminal(CnC::item_collection<K, V> *item_collection,
                  std::vector<CnC::tag_collection<K> *> next_op_tags)
   : item_collection(item_collection), next_op_tags(next_op_tags) {}

   // This method is invoked to put item in the output_item_collection as well as
   // putting the tag to the next tag_collections
   void put(K key, V value) const {
      item_collection->put(key, value);
      for (CnC::tag_collection<K> *next_op_tag: next_op_tags) {
         next_op_tag->put(key);
      }
   }

   // This method is invoked to put tag to the next tag_collections --> mainly used for
   // recursion purposes
   void put(K key) const {
      for (CnC::tag_collection<K> *next_op_tag: next_op_tags) {
         next_op_tag->put(key);
      }
   }
};


// forward declaration of the CnCContext which has been used in different parts of the structurs Project, BinaryOp and Printer
struct CnCContext;

struct Project {
   std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals;
   std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals;
   double (*func)(double); 

   Vector sValue(int n, int l, CnCContext &context) const;
   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   // Project Constructor
   Project(double (*func)(double), std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
           std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
   : func(func),
     input_terminals(input_terminals),
     output_terminals(output_terminals) {}
};

struct BinaryOp {
   std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals;
   std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals;

   using funcT = Vector (*)(const Vector &, const Vector&);
   funcT func;   

   double (*scale_factor)(double);

   Vector unfilter(Vector inputVector, int k, Matrix * hg) const;
   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   // BinaryOp Constructor
   BinaryOp(const funcT &func, double (*scale_factor)(double), std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
   : func(func), scale_factor(scale_factor), input_terminals(input_terminals), output_terminals(output_terminals) {}
};

struct Printer {
   std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals;
   std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals;

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   // Printer Constructor
   Printer(std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
           std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
   : input_terminals(input_terminals), output_terminals(output_terminals) {}
};



// forward declaration of two mathematical functions which we use in our MADNESS computation
double funcA(double x);
double funcB(double x);
Vector sub(const Vector &v1, const Vector &v2);
double sub_scale_factor(double n);

struct CnCContext : public CnC::context<CnCContext> {
   int k,
       quad_npt,
       max_level;

   double thresh;

   double (*a_function)(double);
   double (*b_function)(double);

   Matrix *hg, *hg0, *hg1, *hgT;
   Matrix *rm, *r0, *rp;

   Vector *quad_w, *quad_x;
   Matrix *quad_phi, *quad_phiT, *quad_phiw;

   
   // item_collections
   CnC::item_collection<std::pair<int, int>, Node> projectA_item;
   CnC::item_collection<std::pair<int, int>, Node> projectB_item; 
   CnC::item_collection<std::pair<int, int>, Node> subtract_item;
  
   // tag_collections
   CnC::tag_collection<std::pair<int, int>> projectA_tag;
   CnC::tag_collection<std::pair<int, int>> projectB_tag;
   CnC::tag_collection<std::pair<int, int>> subtract_tag;
   CnC::tag_collection<std::pair<int, int>> printer_tag;

   // stpe_collections   
   using OutputTerminalType = OutputTerminal<std::pair<int, int>, Node>;

   CnC::step_collection<Project> projectA_step;
   CnC::step_collection<Project> projectB_step;
   CnC::step_collection<BinaryOp> subtract_step;
   CnC::step_collection<Printer> printer_step;


   // Function Definitions
   
   // CnCContext Constructor
   CnCContext(int k, double thresh, int max_level)
   : CnC::context<CnCContext>(), projectA_item(*this),
     projectB_item(*this), subtract_item(*this), projectA_tag(*this), projectB_tag(*this), subtract_tag(*this),
     printer_tag(*this), 
     
     projectA_step(*this, "projectA_step", Project(&funcA, std::vector<CnC::item_collection<std::pair<int, int>, Node> *>{},
        std::vector<OutputTerminalType> {
           OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&projectA_tag}),
           OutputTerminalType(&projectA_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&subtract_tag})})),
 
     projectB_step(*this, "projectB_step", Project(&funcB, std::vector<CnC::item_collection<std::pair<int, int>, Node> *>{},
        std::vector<OutputTerminalType> {
           OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&projectB_tag}),
           OutputTerminalType(&projectB_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {})})),
    
    subtract_step(*this, "subtract_step", BinaryOp(&sub, &sub_scale_factor, std::vector<CnC::item_collection<std::pair<int, int>, Node> *> {&projectA_item, &projectB_item},
       std::vector<OutputTerminalType>{
          OutputTerminalType(&projectA_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&subtract_tag}),
          OutputTerminalType(&subtract_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&printer_tag}),
          OutputTerminalType(&projectB_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&subtract_tag})})),

    printer_step(*this, "printer_step", Printer(std::vector<CnC::item_collection<std::pair<int, int>, Node> *>{&subtract_item}, std::vector<OutputTerminalType>{})),


    k(k), thresh(thresh), max_level(max_level) {
      
      projectA_tag.prescribes(projectA_step, *this);
      projectB_tag.prescribes(projectB_step, *this);
      subtract_tag.prescribes(subtract_step, *this);
      printer_tag.prescribes(printer_step, *this);

      projectA_step.produces(projectA_item);
      projectB_step.produces(projectB_item);

      subtract_step.consumes(projectA_item);
      subtract_step.consumes(projectB_item);
      subtract_step.produces(projectA_item);
      subtract_step.produces(projectB_item);
      subtract_step.produces(subtract_item);

      printer_step.consumes(subtract_item);

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
};

//Vector sub(const Vector &v1, const Vector &v2);
//double sub_scale_factor(double n);


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

int Printer::execute(const std::pair<int, int> &node, CnCContext &context) const {
   Node nodeInfo;
   input_terminals[0]->get(node, nodeInfo);
   std::cout << "Printer:: Node with info: (Key: (" << node.first << ", " << node.second << "), " << nodeInfo.toString() << ")" << std::endl;
   return CnC::CNC_Success;
}


int main(int argc, char *argv[]) {
   int k = 6;
   int npt = 20;
   int max_level = atoi(argv[1]);
   double thresh = atof(argv[2]);

   high_resolution_clock::time_point t1 = high_resolution_clock::now();

   CnCContext computation(k, thresh, max_level);
   computation.projectA_tag.put(std::make_pair(0, 0));
   computation.projectB_tag.put(std::make_pair(0, 0));   
   computation.wait();

   high_resolution_clock::time_point t2 = high_resolution_clock::now();
   auto duration = duration_cast<microseconds>( t2 - t1 ).count();
   std::cout << duration/1000000.0 << std::endl;


   return 0;
}
