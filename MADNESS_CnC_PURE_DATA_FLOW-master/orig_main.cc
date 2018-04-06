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
double alpha = 1.0;
double beta = -1.0;

/**********************************************************************/
/**********************************************************************/
/* Struct  - Project
/* Used By - projectA_step, projectB_step   
/**********************************************************************/
/**********************************************************************/
struct Project {

   std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals;
   std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals;
   double (*func)(double); 

   Vector sValue(int n, int l, CnCContext &context) const;
   int execute(const std::pair<int, int> &node, CnCContext &context) const;

  /*----------------------------------------------------------------*/
  /* Project Constructor
  /*----------------------------------------------------------------*/
   Project( double (*func)(double), 
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          :
          func(func),
          input_terminals(input_terminals),
          output_terminals(output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - BinaryOp
/* Used By - subtract_1_step
/**********************************************************************/
/**********************************************************************/
struct BinaryOp {

   std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals;
   std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals;

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
            input_terminals(input_terminals), 
            output_terminals(output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - Printer
/* Used By - printer_step
/**********************************************************************/
/**********************************************************************/
struct Printer {

   std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals;
   std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals;

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   /*----------------------------------------------------------------*/
   /* Printer Constructor
   /*----------------------------------------------------------------*/
   Printer( std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
          input_terminals(input_terminals), 
          output_terminals(output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - Compress_Prolog
/* Used By - compress_prolog_FuncA_step, compress_prolog_FuncB_step
/**********************************************************************/
/**********************************************************************/
struct Compress_Prolog {

   std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals;
   std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals;

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   /*----------------------------------------------------------------*/
   /* Compress_Prolog Constructor
   /*----------------------------------------------------------------*/
   Compress_Prolog( 
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
          input_terminals(input_terminals), 
          output_terminals(output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - Compress_doIt
/* Used By - compress_doIt_funcA_step, compress_doIt_funcB_step
/**********************************************************************/
/**********************************************************************/
struct Compress_doIt {

   std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals;
   std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals;

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   /*----------------------------------------------------------------*/
   /* Compress_doIt Constructor
   /*----------------------------------------------------------------*/
   Compress_doIt( 
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
          input_terminals(input_terminals), 
          output_terminals(output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - GaxpyOp
/* Used By - gaxpyOp_step
/**********************************************************************/
/**********************************************************************/
struct GaxpyOp {

   std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals;
   std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals;

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
          input_terminals(input_terminals), 
          output_terminals(output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - Reconstruct_Prolog
/* Used By - reconstruct_prolog_step
/**********************************************************************/
/**********************************************************************/
struct Reconstruct_Prolog {

   std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals;
   std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals;

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   /*----------------------------------------------------------------*/
   /* Reconstruct_Prolog Constructor
   /*----------------------------------------------------------------*/
   Reconstruct_Prolog( 
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
          input_terminals(input_terminals), 
          output_terminals(output_terminals) {}
};


/**********************************************************************/
/**********************************************************************/
/* Struct  - Reconstruct_doIt
/* Used By - reconstruct_doIt_step
/**********************************************************************/
/**********************************************************************/
struct Reconstruct_doIt {

   std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals;
   std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals;

   int execute(const std::pair<int, int> &node, CnCContext &context) const;

   /*----------------------------------------------------------------*/
   /* Reconstruct_doIt Constructor
   /*----------------------------------------------------------------*/
   Reconstruct_doIt( 
            std::vector<CnC::item_collection<std::pair<int, int>, Node> *> input_terminals,
            std::vector<OutputTerminal<std::pair<int, int>, Node>> output_terminals)
          : 
          input_terminals(input_terminals), 
          output_terminals(output_terminals) {}
};




/**********************************************************************/
/**********************************************************************/
/* Struct  - CnCContext
/**********************************************************************/
/**********************************************************************/
struct CnCContext : public CnC::context<CnCContext> {

   int k, quad_npt, max_level;

   double thresh;

   double (*a_function)(double);
   double (*b_function)(double);

   Matrix *hg, *hg0, *hg1, *hgT;
   Matrix *rm, *r0, *rp;

   Vector *quad_w, *quad_x;
   Matrix *quad_phi, *quad_phiT, *quad_phiw;

   
  /*----------------------------------------------------------------*/
  /* Item Collections
  /*----------------------------------------------------------------*/
   CnC::item_collection<std::pair<int, int>, Node> projectA_item;
   CnC::item_collection<std::pair<int, int>, Node> projectB_item; 
   CnC::item_collection<std::pair<int, int>, Node> subtract_1_item;
  

   CnC::item_collection<std::pair<int, int>, Node> compress_prologA_left_item;
   CnC::item_collection<std::pair<int, int>, Node> compress_prologA_right_item;
   CnC::item_collection<std::pair<int, int>, Node> funcA_coeff_compressed_item;
   CnC::item_collection<std::pair<int, int>, Node> compress_prologB_left_item;
   CnC::item_collection<std::pair<int, int>, Node> compress_prologB_right_item;
   CnC::item_collection<std::pair<int, int>, Node> funcB_coeff_compressed_item;
   CnC::item_collection<std::pair<int, int>, Node> gaxpy_result_item;
   CnC::item_collection<std::pair<int, int>, Node> reconstruct_result_item;
   CnC::item_collection<std::pair<int, int>, Node> s_coeff_item;
   CnC::item_collection<std::pair<int, int>, Node> subtract_2_item;


  /*----------------------------------------------------------------*/
  /* Tag Collections
  /*----------------------------------------------------------------*/
   CnC::tag_collection<std::pair<int, int>> projectA_tag;
   CnC::tag_collection<std::pair<int, int>> projectB_tag;
   CnC::tag_collection<std::pair<int, int>> subtract_1_tag;
   CnC::tag_collection<std::pair<int, int>> printer_tag;


   CnC::tag_collection<std::pair<int, int>> compress_prolog_FuncA_tag;
   CnC::tag_collection<std::pair<int, int>> compress_prolog_FuncB_tag;
   CnC::tag_collection<std::pair<int, int>> compress_doIt_funcA_tag;
   CnC::tag_collection<std::pair<int, int>> compress_doIt_funcB_tag;
   CnC::tag_collection<std::pair<int, int>> gaxpyOP_tag;
   CnC::tag_collection<std::pair<int, int>> subtract_2_tag;
   CnC::tag_collection<std::pair<int, int>> reconstruct_prolog_tag;
   CnC::tag_collection<std::pair<int, int>> reconstruct_doIt_tag;


  /*----------------------------------------------------------------*/
  /* Step Collections
  /*----------------------------------------------------------------*/
   using OutputTerminalType = OutputTerminal<std::pair<int, int>, Node>;

   CnC::step_collection<Project> projectA_step;
   CnC::step_collection<Project> projectB_step;
   CnC::step_collection<BinaryOp> subtract_1_step;
   CnC::step_collection<Printer> printer_step;


   CnC::step_collection<Compress_Prolog> compress_prolog_FuncA_step;
   CnC::step_collection<Compress_Prolog> compress_prolog_FuncB_step;
   CnC::step_collection<Compress_doIt> compress_doIt_funcA_step;
   CnC::step_collection<Compress_doIt> compress_doIt_funcB_step;
   CnC::step_collection<GaxpyOp> gaxpyOp_step;
   CnC::step_collection<Reconstruct_Prolog> reconstruct_prolog_step;
   CnC::step_collection<Reconstruct_doIt>  reconstruct_doIt_step;
   CnC::step_collection<BinaryOp> subtract_2_step;

  /*----------------------------------------------------------------*/
  /* CnCContext Constructor
  /*----------------------------------------------------------------*/
   CnCContext(int k, double thresh, int max_level)
   : 
     CnC::context<CnCContext>(), 
     projectA_item(*this),
     projectB_item(*this), 
     subtract_1_item(*this), 
     projectA_tag(*this), 
     projectB_tag(*this), 
     subtract_1_tag(*this),
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
     subtract_2_item(*this),

     compress_prolog_FuncA_tag(*this),
     compress_prolog_FuncB_tag(*this),
     compress_doIt_funcA_tag(*this),
     compress_doIt_funcB_tag(*this),
     gaxpyOP_tag(*this),
     subtract_2_tag(*this),
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
                              OutputTerminalType(&projectA_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&subtract_1_tag, &compress_prolog_FuncA_tag})}
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
    /* Declare subtract_1_step
    /*----------------------------------------------------------------*/
    subtract_1_step(
                  *this, 
                  "subtract_1_step", 
                  BinaryOp(
                          &sub, 
                          &sub_scale_factor, 
                          std::vector<CnC::item_collection<std::pair<int, int>, Node> *> {&projectA_item, &projectB_item},
                          std::vector<OutputTerminalType>{
                              OutputTerminalType(&projectA_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&subtract_1_tag}),
                              OutputTerminalType(&subtract_1_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {}),
                              OutputTerminalType(&projectB_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {})}
                          )
                  ),

    /*----------------------------------------------------------------*/
    /* Declare printer_step
    /*----------------------------------------------------------------*/
    printer_step(
                *this, 
                "printer_step", 
                Printer( 
                  std::vector<CnC::item_collection<std::pair<int, int>, Node> *>{&subtract_2_item}, 
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
                      alpha,
                      beta,
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
                                OutputTerminalType(&reconstruct_result_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&subtract_2_tag})})
                         ),



    /*----------------------------------------------------------------*/
    /* Declare subtract_2_step
    /*----------------------------------------------------------------*/
    subtract_2_step( 
                    *this, 
                    "subtract_2_step", 
                    BinaryOp( 
                        &sub, 
                        &sub_scale_factor, 
                        std::vector<CnC::item_collection<std::pair<int, int>, Node> *> {&reconstruct_result_item, &subtract_1_item}, 
                        std::vector<OutputTerminalType>{
                          OutputTerminalType(&reconstruct_result_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&subtract_2_tag}),
                          OutputTerminalType(&subtract_2_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&printer_tag}),
                          OutputTerminalType(&subtract_1_item, std::vector<CnC::tag_collection<std::pair<int, int>> *> {})})
                   ),

    k(k), 
    thresh(thresh), 
    max_level(max_level) {
      
      /*----------------------------------------------------------------*/
      /* Tag Prescription
      /*----------------------------------------------------------------*/
      projectA_tag.prescribes(projectA_step, *this);
      projectB_tag.prescribes(projectB_step, *this);
      subtract_1_tag.prescribes(subtract_1_step, *this);
      printer_tag.prescribes(printer_step, *this);

      compress_prolog_FuncA_tag.prescribes( compress_prolog_FuncA_step, *this);
      compress_prolog_FuncB_tag.prescribes( compress_prolog_FuncB_step, *this);
      compress_doIt_funcA_tag.prescribes( compress_doIt_funcA_step, *this);
      compress_doIt_funcB_tag.prescribes( compress_doIt_funcB_step, *this);
      gaxpyOP_tag.prescribes( gaxpyOp_step, *this);
      subtract_2_tag.prescribes( subtract_2_step, *this);
      reconstruct_prolog_tag.prescribes( reconstruct_prolog_step, *this);
      reconstruct_doIt_tag.prescribes( reconstruct_doIt_step, *this);

      /*----------------------------------------------------------------*/
      /* Steps Produce Consume
      /*----------------------------------------------------------------*/
      projectA_step.produces(projectA_item);

      projectB_step.produces(projectB_item);

      subtract_1_step.consumes(projectA_item);
      subtract_1_step.consumes(projectB_item);
      subtract_1_step.produces(projectA_item);
      subtract_1_step.produces(projectB_item);
      subtract_1_step.produces(subtract_1_item);

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

      subtract_2_step.consumes( reconstruct_result_item);
      subtract_2_step.consumes( subtract_1_item);
      subtract_2_step.produces( reconstruct_result_item);
      subtract_2_step.produces( subtract_1_item);
      subtract_2_step.produces( subtract_2_item);

      printer_step.consumes(subtract_2_item);

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
    tempD.gaxpy( 1.0, right.d, -1.0);

    Vector tempS;

    if( node.first == 0 && node.second == 0){
      tempS = left.s;
      tempS.gaxpy(1.0, right.s, -1.0);
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
