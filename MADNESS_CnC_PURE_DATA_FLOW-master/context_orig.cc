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

   tbb::concurrent_vector<double> inner_results;
   tbb::concurrent_vector<double> norm2_results;
   
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


   CnC::item_collection<std::pair<int, int>, Node> evaluate_item; // Have to decide who is gonna put item here
   CnC::item_collection<std::pair<int, int>, Node> norm2_item; // Have to decide who is gonna put item here

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

   // CnC::tag_collection<std::pair<int, int>> evaluate_tag;
   CnC::tag_collection<std::pair<int, int>> norm2_tag;

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

   CnC::step_collection<Norm2> norm2_step;
   // CnC::step_collection<Evaluate> evaluate_step;

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

     norm2_item(*this),
     norm2_tag(*this),
     evaluate_item( *this),


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


    /*----------------------------------------------------------------*/
    /* Declare subtract_2_step
    /*----------------------------------------------------------------*/
    norm2_step( 
              *this, 
              "norm2_step", 
              Norm2( 
                  std::vector<CnC::item_collection<std::pair<int, int>, Node> *> {&norm2_item}, 
                  std::vector<OutputTerminalType>{
                    OutputTerminalType(nullptr, std::vector<CnC::tag_collection<std::pair<int, int>> *> {&norm2_tag})})
             ),

    /*----------------------------------------------------------------*/
    /* Declare printer_step
    /*----------------------------------------------------------------*/
    // evaluate_step(
    //             *this, 
    //             "evaluate_step", 
    //             Evaluate( 
    //               x,
    //               std::vector<CnC::item_collection<std::pair<int, int>, Node> *>{&evaluate_item}, 
    //               std::vector<OutputTerminalType>{})
    //             ),


    k(k), 
    thresh(thresh), 
    max_level(max_level)
    // norm2_results(max_level , 0.0),
    // inner_results( max_level, 0.0) 
    {
      
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

      norm2_step.consumes(norm2_item);

      init_twoscale(k);
      init_quadrature(k);
      make_dc_periodic();
 
   }
