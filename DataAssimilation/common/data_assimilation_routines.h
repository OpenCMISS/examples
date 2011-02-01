
//command: mpiexec -n 1 ./01  -measurement-noise-scale 0.1 -global-noise-scale  0.1 -initialization-error 100 -measurement-noise-amount 100
static char help[] = "rUKF approach.\n\
  ... \n\n";

#undef __FUNCT__
#define __FUNCT__ "main"
#define PETSC_OPERATIONS_DEBUG 1 

#include "slepcsvd.h" 
#include "petscksp.h"
#include "petscmat.h"
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "data_assimilation_utility_routines.h"

#define ERROR_TYPE_COMMON 1 
#define INFO_PETSC_INT  1
#define INFO_PETSC_REAL 2
#define INFO_PETSC_VEC  3
#define INFO_PETSC_MAT  4
#define INFO_PETSC_INTARRAY  5
#define INFO_PETSC_REALARRAY  6

#define FILTER_TYPE_RUKF  7

#define MEAUSREMENT_TYPE_FROM_HEADER_AND_DAT_FILES 1

#define DIAGNOSIS_DISP_EVERYTHING   0x8000
#define DIAGNOSIS_DISP_MEDIUM       0x4000
#define DIAGNOSIS_DISP_MINIMUM      0x2000
#define DIAGNOSIS_DISP_INITIAL_CONF 0x0001

#define app_context_assert(x,msg,app_context) {if(!(x)) app_context_report_error(ERROR_TYPE_COMMON, msg, app_context, __LINE__,__FUNCT__,__FILE__);}

#ifndef MEASUREMENT_DIR
  #define MEASUREMENT_DIR  "./" //where to find the measurement files 
#endif

#define MAX_SIZE_PATH_NAME 256
#define MAX_SIZE_FILE_NAME 128

//-------------------------------------------------------------------------------------------------
//---structures definitions---
struct _APP_CONTEXT_TYPE{
  MPI_Comm              mpi_comm_global;//process in the work group               
  MPI_Comm              mpi_comm_local; //this process               
  PetscViewer           pv_global;             
  PetscViewer           pv_local;              
  PetscErrorCode        ierr;
  PetscInt              diagnosisMode;
  
  PetscReal measurement_noise_scale; // for R
  PetscReal global_noise_scale;      // for both R and P_0
//  PetscReal initialization_error;    // for theta_0
  PetscReal measurement_noise_amount; // for y_k

  PetscInt n_x;//ndof 
  PetscInt n_theta; //npdof
  PetscReal* x_0; 
  PetscReal  x_scale_constant; 
  PetscReal* x_scale; //currently initialized from x_scale_constant;
  PetscReal* theta_truth; 
  PetscReal* theta_0; 
  PetscReal* theta_scale; 
  
  //temporary storage for computation
  //-SVD related
  Mat S_ensemble; //ensemble matrix, only for storage purpose (for doing SVD), app_context->n_x+app_context->n_theta by model->n_ensem
  Mat S_ensemble2; //ensemble matrix, only for storage purpose (for S_ensemble*S_ensemble2^T)
  Mat P_z; //BIG! error covariance for the augmented state vecotr, used by rUKF for doing SVD,app_context->n_x+app_context->n_theta by app_context->n_x+app_context->n_theta
  Mat S_z; //SVD result of P_z, app_context->n_x+app_context->n_theta by filter->subfilter->r
  Mat S_bar; //used in the measurement update by rUKF,  measurement->n_y by filter->subfilter->r
  Mat S_bar2; //temperory storage for a matrix of the size measurement->n_y by filter->subfilter->r
  //-general temperory usage
  Mat tmp_r; // r by r matrix 
  Mat tmp_nz_by_r;
  Mat tmp_r_by_ny;
  Vec tmp_z; //n_x+n_theta by 1 vector  
  Vec tmp_x; //n_x by 1 vector  
  Vec tmp_theta; //n_theta by 1 vector
  Vec tmp_y; //measurement->n_y by 1 vector
  PetscInt* tmp_x_indices;
  PetscInt* tmp_theta_indices;
  
  //measurement info
  char* measurement_dir_name;
  char* measurement_header_file_name;
  char* measurement_data_file_name;
  char* measurement_matrix_file_name;
  char* measurement_error_covariance_matrix_file_name;
  int measurementType;
};
typedef struct _APP_CONTEXT_TYPE APP_CONTEXT_TYPE;


struct _MODEL_TYPE{
  struct _MODEL_TYPE* (*f)(struct _MODEL_TYPE*, APP_CONTEXT_TYPE*);    //f: model forward operator : n (number of ensembles) states at k --> n states at k+1 (n>=1)

  PetscInt n_x;//ndof 
  PetscInt n_theta;
  PetscInt n_z;
  
  Vec x_k_hat;
  Vec theta_k_hat;
  Vec x_k_minus_1_hat;
  Vec theta_k_minus_1_hat;
  
  PetscInt k; //time step
  PetscReal t; //time 
 
  PetscInt n_ensem; //number of ensembles
  Vec* x_k_minus_1;
  Vec* theta_k_minus_1;
  Vec* x_k;
  Vec* theta_k;
//  measurementDisplacementScale=[];
//ipfibreNums=[];
//fixed_contractility_factor=[];
//phase=[];
//Cai_i_max=[];
//parameter_permutation_matrx =[];
//repeatScaleFacotrs =[];
//repeatNumber=[]; 
};
typedef struct _MODEL_TYPE MODEL_TYPE;

struct _FILTER_SUBTYPE{
  PetscInt filter_type;
  PetscInt r; //rank ( r>=(n_theat>0)?n_theat:n_x, r<= n_theta+n_x;
  Mat P_x; //BIG! error covariance for the state vecotr, not used by rUKF
  Mat P_theta; //error covariance for the parameter vector, not used by rUKF
  
  Mat S_x; //square-root of P_x
  Mat S_theta;   //square-root of P_theta
};

typedef struct _FILTER_SUBTYPE FILTER_SUBTYPE;

struct _FILTER_TYPE{
  PetscInt n_x;            //ndof, n_x>=0 
  PetscInt n_theta;        //ndof, n_theta>=0 
  PetscInt n_z;
  Mat K; //Kalman gain matrix, n_x + n_theta by measurement->n_y
//  Vec x;
//  Vec theta;
//  PetscInt filter_type;
  FILTER_SUBTYPE* subfilter; //FILTER_SUBTYPE
};
typedef struct _FILTER_TYPE FILTER_TYPE;

struct _OPTIMIZER_TYPE{
  PetscInt n_theta; //ndof 
};
typedef struct _OPTIMIZER_TYPE OPTIMIZER_TYPE;


struct _MEAUSREMENT_TYPE{
  //Vec k_y;
  //Vec t_y;
  Vec y; //measurement_petsc_vector, current measurement in use(maybe interpolated), n_y by 1
  Mat H; //measurement matrix, n_y by n_z
  Mat R_d; //measurement error covariance, n_y by n_y
  Mat inv_R_d; //measurement error covariance inverse, n_y by n_y
  
  int n_y;//length of single measurement vector
  int n_z;//length of state vector
  int K; //number of measurements
  int nMeta; //default 2

  double* map_index2time;     //nMeta(2, index, time) by K 
  double* all_measurements;   //n_y by K
  double* single_measurement; //n_y by 1
  double* H_values; //n_y by n_z;
  double* R_d_values;//n_y by n_y
//  Vec theta_y; //ground-truth parameter 
//xiCoordsFileName=[];
//sampleType=[];
//nSampleOrder=[];
//elementsID=[];
//facesID=[];
//tempDirID =[];
//experimentID =[];


};
typedef struct _MEAUSREMENT_TYPE MEAUSREMENT_TYPE;


//-------------------------------------------------------------------------------------------------
//---all routines---

PetscErrorCode 
app_context_report_error(PetscErrorCode errorcode, const char* msg, APP_CONTEXT_TYPE* app_context, int lineNum, const char* funct,  const char* file);

PetscErrorCode 
app_context_report_info(PetscInt info_type, PetscInt info_size, void* info, const char* msg, APP_CONTEXT_TYPE* app_context);

int 
app_context_destroy_and_finalize(APP_CONTEXT_TYPE* app_context);

APP_CONTEXT_TYPE* 
app_context_create_and_initialize(int argc,char **args);

MODEL_TYPE* 
model_create_and_initialize(MODEL_TYPE* (*model_opeartor)(MODEL_TYPE*, APP_CONTEXT_TYPE*), APP_CONTEXT_TYPE* app_context);

MODEL_TYPE* 
model_set_opeartor_and_intial_condition(MODEL_TYPE* model, MODEL_TYPE* (*model_opeartor)(MODEL_TYPE*, APP_CONTEXT_TYPE*), APP_CONTEXT_TYPE* app_context);

FILTER_TYPE* 
filter_create_and_initialize(APP_CONTEXT_TYPE* app_context,MODEL_TYPE* model);

FILTER_TYPE* 
filter_set_intial_condition(FILTER_TYPE* filter, MODEL_TYPE* model, APP_CONTEXT_TYPE* app_context);

OPTIMIZER_TYPE* 
optimizer_create_and_initialize(APP_CONTEXT_TYPE* app_context);

OPTIMIZER_TYPE* 
optimizer_set_intial_condition(OPTIMIZER_TYPE*  optimizer, APP_CONTEXT_TYPE* app_context);

//create and initialize the filter
FILTER_SUBTYPE* 
filter_subtype_create_and_initialize(FILTER_TYPE* filter, MODEL_TYPE* model, APP_CONTEXT_TYPE* app_context);

// Initialize the filter subtype
//  PetscInt filter_type;
//  PetscInt r; //rank ( r>=(n_theat>0)?n_theat:1, r<= n_theta+n_x;
//  Mat S_x;
//  Mat S_theta;  
FILTER_SUBTYPE* 
filter_subtype_set_intial_condition(FILTER_TYPE* filter, MODEL_TYPE* model, APP_CONTEXT_TYPE* app_context);

//-------------------------------------------------------------------------------------------------
//---app_context routines---
PetscErrorCode 
app_context_report_error(PetscErrorCode errorcode, const char* msg, APP_CONTEXT_TYPE* app_context, int lineNum, const char* funct,  const char* file)
{
  static int previous_lineNum=0; 
  if(errorcode){
    PetscPrintf(app_context->mpi_comm_global, "== fatal error: %s: at line %d, function %s(), file %s. prevous line number: %d \n", msg, lineNum, funct, file,previous_lineNum);
    //CHKERRQ(app_context->ierr); //if errorcode is non-zero,     
    SETERRQ(errorcode,msg);
  }else{
    previous_lineNum = lineNum;
    return 0;
  }
}

PetscErrorCode 
app_context_report_info(PetscInt info_type, PetscInt info_size, void* info, const char* msg, APP_CONTEXT_TYPE* app_context)
{
  app_context->ierr = PetscPrintf(app_context->mpi_comm_global, "== info: %s: ", msg);
  //app_context_report_error(app_context->ierr,"app_context_report_info: PetscPrintf",app_context,__LINE__,__FUNCT__,__FILE__);

  switch(info_type)
  {
    case INFO_PETSC_INT:
    case INFO_PETSC_INTARRAY:
      app_context->ierr = PetscIntView(info_size, (PetscInt*)info, app_context->pv_global);
      app_context_report_error(app_context->ierr,"app_context_report_info: INFO_PETSC_INT",app_context,__LINE__,__FUNCT__,__FILE__);
      break;
    case INFO_PETSC_REAL:
    case INFO_PETSC_REALARRAY:
      app_context->ierr = PetscRealView(info_size, (PetscReal*)info, app_context->pv_global);
      app_context_report_error(app_context->ierr,"app_context_report_info: INFO_PETSC_REAL",app_context,__LINE__,__FUNCT__,__FILE__);
      break;
    case INFO_PETSC_VEC:
      app_context->ierr = PetscPrintf(app_context->mpi_comm_global, "\n", msg);
      app_context_report_error(app_context->ierr,"app_context_report_info: INFO_PETSC_VEC",app_context,__LINE__,__FUNCT__,__FILE__);
      app_context->ierr = VecView( *((Vec*)info), app_context->pv_global);
      app_context_report_error(app_context->ierr,"app_context_report_info: INFO_PETSC_VEC",app_context,__LINE__,__FUNCT__,__FILE__);
      break;
    case INFO_PETSC_MAT:
      app_context->ierr = PetscPrintf(app_context->mpi_comm_global, "\n", msg);
      app_context_report_error(app_context->ierr,"app_context_report_info: INFO_PETSC_MAT",app_context,__LINE__,__FUNCT__,__FILE__);
      app_context->ierr = MatView( *((Mat*)info), app_context->pv_global);    
      app_context_report_error(app_context->ierr,"app_context_report_info: INFO_PETSC_MAT",app_context,__LINE__,__FUNCT__,__FILE__);
      break;
    default:
      app_context->ierr = PetscPrintf(app_context->mpi_comm_global, "== wrong info_type: %d\n", info_type);
      app_context_report_error(app_context->ierr,"app_context_report_info: default",app_context,__LINE__,__FUNCT__,__FILE__);
      break;
  }
  app_context->ierr = PetscPrintf(app_context->mpi_comm_global, "\n");
  app_context_report_error(app_context->ierr,"app_context_report_info: PetscPrintf",app_context,__LINE__,__FUNCT__,__FILE__);
  return 0;
}

APP_CONTEXT_TYPE* 
app_context_create_and_initialize(int argc,char **args)
{
  APP_CONTEXT_TYPE* app_context = NULL;
  PetscTruth flg = PETSC_FALSE, flg_1=PETSC_FALSE;
  int i=0;
//  char options[2][PETSC_MAX_PATH_LEN];     /* input file name */
  
 //app_context->ierr = PetscMalloc(sizeof(APP_CONTEXT_TYPE), (void **) &app_context);
 app_context = (APP_CONTEXT_TYPE*) malloc(sizeof(APP_CONTEXT_TYPE));
  if(app_context==NULL){
    PetscPrintf(MPI_COMM_WORLD, "== fatal error: app_context_create_and_initialize: app_context==NULL");
    PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,ERROR_TYPE_COMMON,1,"app_context_create_and_initialize: app_context==NULL");
    return NULL;
  }
  
 //app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscMalloc",app_context);

  //program initialization 

//  PetscInitialize(&argc,&args,(char *)0,help);
  SlepcInitialize(&argc,&args,(char*)0,(char*) 0);

  app_context->mpi_comm_local=MPI_COMM_SELF;
  app_context->mpi_comm_global=MPI_COMM_WORLD;
  app_context->pv_global=PETSC_VIEWER_STDOUT_WORLD;
  app_context->pv_local=PETSC_VIEWER_STDOUT_SELF;
  PetscViewerSetFormat(app_context->pv_local, PETSC_VIEWER_ASCII_MATLAB) ; //PETSC_VIEWER_ASCII_INDEX);
  PetscViewerSetFormat(app_context->pv_global, PETSC_VIEWER_ASCII_MATLAB) ;//PETSC_VIEWER_ASCII_INDEX);
    
  // --- program parameter initialization  ----
  
  //default values
  app_context->diagnosisMode = DIAGNOSIS_DISP_INITIAL_CONF | DIAGNOSIS_DISP_EVERYTHING;
  app_context->measurement_noise_scale = 0.25;
  app_context->global_noise_scale = 0.0001;
//  app_context->initialization_error = 0.5;
  app_context->measurement_noise_amount = 0.0;

  app_context->n_x = 2;
  
  app_context->x_scale_constant=1.0;
  
  PetscMalloc(sizeof(PetscReal)*app_context->n_x, (void**) &app_context->x_0);
  app_context->x_0[0] = 1.0; 
  app_context->x_0[1] = 1.0; 

  app_context->n_theta = 2 ;
  PetscMalloc(sizeof(PetscReal)*app_context->n_theta, (void**) &app_context->theta_truth);
  app_context->theta_truth[0] = 1.0; 
  app_context->theta_truth[1] = 1.0; 
  
  PetscMalloc(sizeof(PetscReal)*app_context->n_theta, (void**) &app_context->theta_0);
  app_context->theta_0[0] = 1.5; 
  app_context->theta_0[1] = 0.5; 
  
  PetscMalloc(sizeof(PetscReal)*app_context->n_theta, (void**) &app_context->theta_scale);
  app_context->theta_scale[0] = 1.0; 
  app_context->theta_scale[1] = 1.0; 
  
  //user-defined values
  app_context->ierr = PetscOptionsGetReal(PETSC_NULL,"-measurement-noise-scale",&app_context->measurement_noise_scale,&flg);
  app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetReal",app_context,__LINE__,__FUNCT__,__FILE__);
  
  app_context->ierr = PetscOptionsGetReal(PETSC_NULL,"-global-noise-scale",&app_context->global_noise_scale,&flg);
  app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetReal",app_context,__LINE__,__FUNCT__,__FILE__);
  
//  app_context->ierr = PetscOptionsGetReal(PETSC_NULL,"-initialization-error",&app_context->initialization_error,&flg);
//  app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetReal",app_context,__LINE__,__FUNCT__,__FILE__);
  
  app_context->ierr = PetscOptionsGetReal(PETSC_NULL,"-measurement-noise-amount",&app_context->measurement_noise_amount,&flg);
  app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetReal",app_context,__LINE__,__FUNCT__,__FILE__);
  
  //app_context->x_0
  app_context->ierr = PetscOptionsGetInt(PETSC_NULL,"-n-x",&app_context->n_x,&flg);
  app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetInt",app_context,__LINE__,__FUNCT__,__FILE__);
  if(flg){ //resize the array
    PetscFree(app_context->x_0);
    app_context->ierr=PetscMalloc(sizeof(PetscReal)*app_context->n_x, (void **) &app_context->x_0);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscMalloc",app_context,__LINE__,__FUNCT__,__FILE__);
  }
  flg_1 = flg;
  PetscOptionsHasName(PETSC_NULL,"-x-0",&flg);
  if( flg_1 | flg ){
    app_context_report_error(ERROR_TYPE_COMMON,"app_context_create_and_initialize: -n-x and -x_0 should be given together",app_context,__LINE__,__FUNCT__,__FILE__);
  }
  PetscOptionsHasName(PETSC_NULL,"-x-0",&flg);
  if(flg){
    app_context->ierr = PetscOptionsGetRealArray(PETSC_NULL,"-x-0",app_context->x_0, &app_context->n_x,&flg);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetInt",app_context,__LINE__,__FUNCT__,__FILE__);
  }
  
  //x_scale_constant
  app_context->ierr = PetscOptionsGetReal(PETSC_NULL,"-x-scale",&app_context->x_scale_constant, &flg);
  app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetReal x_scale",app_context,__LINE__,__FUNCT__,__FILE__);

  //x-scale
  app_context->ierr = PetscMalloc(sizeof(PetscReal)*app_context->n_x, (void **) &app_context->x_scale);
  app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscMalloc",app_context,__LINE__,__FUNCT__,__FILE__);
  for(i=0; i<app_context->n_x; i++)
    app_context->x_scale[i]=app_context->x_scale_constant;
  
  //app_context->theta_truth
  app_context->ierr = PetscOptionsGetInt(PETSC_NULL,"-n-theta-truth",&app_context->n_theta,&flg);
  app_context_assert(app_context->n_theta>=0, "", app_context);
  app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetInt",app_context,__LINE__,__FUNCT__,__FILE__);
  if(flg){ //resize the array
    PetscFree(app_context->x_0);
    PetscMalloc(sizeof(PetscReal)*app_context->n_theta, (void **) &app_context->theta_truth);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscMalloc",app_context,__LINE__,__FUNCT__,__FILE__);
  }
  flg_1 = flg;
  PetscOptionsHasName(PETSC_NULL,"-theta-truth",&flg);
  if( app_context->n_theta>0 && (flg_1 | flg) ){
    app_context_report_error(ERROR_TYPE_COMMON,"app_context_create_and_initialize: -n-theta-truth and -theta_truth should be given together",app_context,__LINE__,__FUNCT__,__FILE__);
  }
  if(app_context->n_theta>0 && flg){
    app_context->ierr = PetscOptionsGetRealArray(PETSC_NULL,"-theta-truth",app_context->theta_truth,&app_context->n_theta, &flg);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetRealArray",app_context,__LINE__,__FUNCT__,__FILE__);
  }  
  
  //app_context->theta_0
  app_context->ierr = PetscOptionsGetInt(PETSC_NULL,"-n-theta-0",&app_context->n_theta,&flg);
  app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetInt",app_context,__LINE__,__FUNCT__,__FILE__);
  if(flg){ //resize the array
    if(app_context->n_theta != sizeof(app_context->theta_truth)/sizeof(PetscReal))
      app_context_report_error(app_context->ierr,"app_context_create_and_initialize: -n-theta-0 is different from -n-theta-truth",app_context,__LINE__,__FUNCT__,__FILE__);
    PetscFree(app_context->x_0);
    PetscMalloc(sizeof(PetscReal)*app_context->n_theta, (void **) &app_context->theta_0);
  }
  flg_1 = flg;
  PetscOptionsHasName(PETSC_NULL,"-theta-0",&flg);
  if( flg_1 | flg ){
    app_context_report_error(ERROR_TYPE_COMMON,"app_context_create_and_initialize: -n-theta-0 and -theta_0 should be given together",app_context,__LINE__,__FUNCT__,__FILE__);
  }
  if(flg){
    app_context->ierr = PetscOptionsGetRealArray(PETSC_NULL,"-theta-0",app_context->theta_0,&app_context->n_theta, &flg);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetIntArray",app_context,__LINE__,__FUNCT__,__FILE__);
  }
  
  //app_context->theta_scale
  app_context->ierr = PetscOptionsGetInt(PETSC_NULL,"-n-theta-scale",&app_context->n_theta,&flg);
  app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetInt",app_context,__LINE__,__FUNCT__,__FILE__);
  if(flg){ //resize the array
    if(app_context->n_theta != sizeof(app_context->theta_truth)/sizeof(PetscReal))
      app_context_report_error(app_context->ierr,"app_context_create_and_initialize: -n-theta-scale is different from -n-theta",app_context,__LINE__,__FUNCT__,__FILE__);
    PetscFree(app_context->theta_scale);
    PetscMalloc(sizeof(PetscReal)*app_context->n_theta, (void **) &app_context->theta_scale);
  }
  flg_1 = flg;
  PetscOptionsHasName(PETSC_NULL,"-theta-scale",&flg);
  if( flg_1 | flg ){
    app_context_report_error(ERROR_TYPE_COMMON,"app_context_create_and_initialize: n_theta_scale and -theta-scale should be given together",app_context,__LINE__,__FUNCT__,__FILE__);
  }
  if(flg){
    app_context->ierr = PetscOptionsGetRealArray(PETSC_NULL,"-theta-scale",app_context->theta_scale,&app_context->n_theta, &flg);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscOptionsGetRealArray",app_context,__LINE__,__FUNCT__,__FILE__);
  }  
  
  //echo configuration information 
  if(app_context->diagnosisMode & DIAGNOSIS_DISP_INITIAL_CONF){
    app_context->ierr = app_context_report_info(INFO_PETSC_REAL,1, &app_context->measurement_noise_scale,"measurement_noise_scale", app_context);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: app_context_report_info",app_context,__LINE__,__FUNCT__,__FILE__);
    
    app_context->ierr = app_context_report_info(INFO_PETSC_REAL,1, &app_context->global_noise_scale,"global_noise_scale", app_context);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: app_context_report_info",app_context,__LINE__,__FUNCT__,__FILE__);

    //app_context->ierr = app_context_report_info(INFO_PETSC_REAL,1, &app_context->initialization_error,"initialization_error", app_context);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: app_context_report_info",app_context,__LINE__,__FUNCT__,__FILE__);
    
    app_context->ierr = app_context_report_info(INFO_PETSC_INT,1, &app_context->n_x,"n_x", app_context);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: app_context_report_info",app_context,__LINE__,__FUNCT__,__FILE__);
    
    app_context->ierr = app_context_report_info(INFO_PETSC_REALARRAY, app_context->n_x, app_context->x_0,"x_0", app_context);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: app_context_report_info",app_context,__LINE__,__FUNCT__,__FILE__);
    
    app_context->ierr = app_context_report_info(INFO_PETSC_REALARRAY, app_context->n_x, app_context->x_scale,"x_scale", app_context);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: app_context_report_info",app_context,__LINE__,__FUNCT__,__FILE__);
    
    app_context->ierr = app_context_report_info(INFO_PETSC_INT,1, &app_context->n_theta,"n_theta", app_context);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: app_context_report_info",app_context,__LINE__,__FUNCT__,__FILE__);
    
    app_context->ierr = app_context_report_info(INFO_PETSC_REALARRAY, app_context->n_theta, app_context->theta_truth,"theta_truth", app_context);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: app_context_report_info",app_context,__LINE__,__FUNCT__,__FILE__);
    
    app_context->ierr = app_context_report_info(INFO_PETSC_REALARRAY, app_context->n_theta, app_context->theta_0,"theta_0", app_context);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: app_context_report_info",app_context,__LINE__,__FUNCT__,__FILE__);
    
    app_context->ierr = app_context_report_info(INFO_PETSC_REALARRAY, app_context->n_theta, app_context->theta_scale,"theta_scale", app_context);
    app_context_report_error(app_context->ierr,"app_context_create_and_initialize: app_context_report_info",app_context,__LINE__,__FUNCT__,__FILE__);    
  }  
  
  //--allocate temp space
  VecCreate(app_context->mpi_comm_global, &app_context->tmp_x);
  VecCreate(app_context->mpi_comm_global, &app_context->tmp_theta);
  VecCreate(app_context->mpi_comm_global, &app_context->tmp_z);

  VecSetSizes(app_context->tmp_x, PETSC_DECIDE,app_context->n_x);
  VecSetSizes(app_context->tmp_theta, PETSC_DECIDE,app_context->n_theta);
  VecSetSizes(app_context->tmp_z, PETSC_DECIDE,app_context->n_x+app_context->n_theta);
  app_context->ierr=PetscMalloc(sizeof(PetscInt)*app_context->n_x, &app_context->tmp_x_indices);
  app_context_report_error(app_context->ierr,"PetscMalloc: app_context->tmp_x_indices",app_context,__LINE__,__FUNCT__,__FILE__);    
  app_context->ierr=PetscMalloc(sizeof(PetscInt)*app_context->n_theta, &app_context->tmp_theta_indices);
  app_context_report_error(app_context->ierr,"PetscMalloc: app_context->tmp_theta_indices",app_context,__LINE__,__FUNCT__,__FILE__);    

  VecSetFromOptions(app_context->tmp_x);
  VecSetFromOptions(app_context->tmp_theta);
  VecSetFromOptions(app_context->tmp_z);
  
  //measurement
  app_context->measurement_dir_name  = new char[MAX_SIZE_PATH_NAME+1];
  app_context->measurement_data_file_name = new char[MAX_SIZE_FILE_NAME+MAX_SIZE_PATH_NAME+1];
  app_context->measurement_header_file_name = new char[MAX_SIZE_FILE_NAME+MAX_SIZE_PATH_NAME+1];
  app_context->measurement_matrix_file_name = new char[MAX_SIZE_FILE_NAME+MAX_SIZE_PATH_NAME+1];
  app_context->measurement_error_covariance_matrix_file_name = new char[MAX_SIZE_FILE_NAME+MAX_SIZE_PATH_NAME+1];
  
  PetscOptionsGetString(PETSC_NULL,"-measurement_dir_name", app_context->measurement_dir_name,MAX_SIZE_PATH_NAME,&flg);
  if(flg==PETSC_FALSE) sprintf(app_context->measurement_dir_name, "%smeasurement/", MEASUREMENT_DIR);
  
  PetscOptionsGetString(PETSC_NULL,"-measurement_header_file_name", app_context->measurement_header_file_name,MAX_SIZE_FILE_NAME+MAX_SIZE_PATH_NAME,&flg);
  if(flg==PETSC_FALSE) 
    sprintf(app_context->measurement_header_file_name, "%sheader", app_context->measurement_dir_name);
  else if(!strchr(app_context->measurement_header_file_name,'/'))
    sprintf(app_context->measurement_header_file_name, "%s%s", app_context->measurement_dir_name,app_context->measurement_header_file_name);
  else{}
  
  PetscOptionsGetString(PETSC_NULL,"-measurement_matrix_file_name", app_context->measurement_matrix_file_name,MAX_SIZE_FILE_NAME+MAX_SIZE_PATH_NAME,&flg);
  if(flg==PETSC_FALSE) 
  sprintf(app_context->measurement_matrix_file_name, "%sH.dat", app_context->measurement_dir_name);
  else if(!strchr(app_context->measurement_matrix_file_name,'/'))
    sprintf(app_context->measurement_matrix_file_name, "%s%s", app_context->measurement_dir_name,app_context->measurement_matrix_file_name);
  else{}

  PetscOptionsGetString(PETSC_NULL,"-measurement_error_covariance_matrix_file_name", app_context->measurement_error_covariance_matrix_file_name,MAX_SIZE_FILE_NAME+MAX_SIZE_PATH_NAME,&flg);
  if(flg==PETSC_FALSE) 
  sprintf(app_context->measurement_error_covariance_matrix_file_name, "%sR_d.dat", app_context->measurement_dir_name);
  else if(!strchr(app_context->measurement_error_covariance_matrix_file_name,'/'))
    sprintf(app_context->measurement_error_covariance_matrix_file_name, "%s%s", app_context->measurement_dir_name,app_context->measurement_error_covariance_matrix_file_name);
  else{}

  PetscOptionsGetInt(PETSC_NULL,"-measurementType",&app_context->measurementType,&flg);
  if(flg==PETSC_FALSE) app_context->measurementType=MEAUSREMENT_TYPE_FROM_HEADER_AND_DAT_FILES;
    
  return app_context;
}

int 
app_context_destroy_and_finalize(APP_CONTEXT_TYPE* app_context)
{
  if(app_context==NULL){
    PetscPrintf(MPI_COMM_WORLD, "== fatal error: app_context_create_and_initialize: app_context==NULL");
    PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,ERROR_TYPE_COMMON,1,"app_context_finalize: app_context==NULL");
    return 1;
  }
  
//  app_context->ierr = PetscFinalize(); CHKERRQ(app_context->ierr);
  app_context->ierr = SlepcFinalize(); CHKERRQ(app_context->ierr);

  return 0;
}

//-------------------------------------------------------------------------------------------------
//---measurement routines---

//load single mesurement from .dat file with simple format: % dimension length(V) 1\n V \n
// into measurement->single_measurement 
int
measurement_io_read_single_measurement_dat(char* filename, MEAUSREMENT_TYPE* measurement, APP_CONTEXT_TYPE* app_context)
{
  int measurement_size_1 = 0, measurement_size_2 = 0;
  
  if(measurement->single_measurement)    free(measurement->single_measurement);
  measurement->single_measurement = load_matrix_2D_ascii(filename, measurement_size_1, measurement_size_2);  
  assert(measurement_size_2==1); assert(measurement_size_1 > 1);
  measurement->n_y = measurement_size_1;
  
  return 0;
}

//write .dat file with simple format: length(V) \n V \n
int
measurement_io_write_single_measurement_dat(char* filename, MEAUSREMENT_TYPE* measurement, APP_CONTEXT_TYPE* app_context)
{
  int ret = 0;
  
  if(measurement->single_measurement) ret = save_matrix_2D_ascii(filename, measurement->single_measurement, measurement->n_y, 1, NULL);  
  assert(ret==0); 
  
  return 0;
}

//load all meausrment into measurement data structure
int
measurement_io_read_all_measurements(MEAUSREMENT_TYPE* measurement, APP_CONTEXT_TYPE* app_context)
{
    int meta_info_size_1 = 0, meta_info_size_2 = 0;

    //header + dat files
  if(app_context->measurementType==MEAUSREMENT_TYPE_FROM_HEADER_AND_DAT_FILES){ 
    
    //load header file
    measurement->map_index2time = load_matrix_2D_ascii(app_context->measurement_header_file_name, meta_info_size_1, meta_info_size_2);
    assert(meta_info_size_1 == 2); assert(meta_info_size_2 > 2);
    measurement->K = meta_info_size_2;
    
    //load each measurement
    measurement->all_measurements = new double[measurement->n_y*measurement->K];
    for (int i=0; i<measurement->K; i++){
      sprintf(app_context->measurement_data_file_name, "%s%04d.dat", app_context->measurement_dir_name, (int)measurement->map_index2time[i*2+0]);
      measurement_io_read_single_measurement_dat(app_context->measurement_data_file_name, measurement, app_context);
      for (int j=0; j<measurement->n_y; j++) measurement->all_measurements[j*measurement->K+i]=measurement->single_measurement[j];
    }
  }
  
  return 0;
}

//create (Vec) measurement->y from (double*) measurement->single_measurement
int
measurement_petsc_vector_assemble(MEAUSREMENT_TYPE* measurement, APP_CONTEXT_TYPE* app_context)
{  
  static int* index = NULL;
  if(index==NULL) {
    index = new int[measurement->n_y];
    VecCreate(app_context->mpi_comm_global, &measurement->y);
    VecSetSizes(measurement->y, PETSC_DECIDE, measurement->n_y);
    VecSetFromOptions(measurement->y);
    for (int i = 0; i<measurement->n_y; i++) index[i]=i;
  }
  app_context->ierr=VecSetValues(measurement->y,measurement->n_y, index, measurement->single_measurement,ADD_VALUES);
  app_context_report_error(app_context->ierr,"measurement_petsc_vector_assemble: VecSetValues",app_context,__LINE__,__FUNCT__,__FILE__);
  app_context->ierr=VecAssemblyBegin(measurement->y);
  app_context_report_error(app_context->ierr,"measurement_petsc_vector_assemble: VecAssemblyBegin",app_context,__LINE__,__FUNCT__,__FILE__);
  app_context->ierr=VecAssemblyEnd(measurement->y);
  app_context_report_error(app_context->ierr,"measurement_petsc_vector_assemble: VecAssemblyEnd",app_context,__LINE__,__FUNCT__,__FILE__);
  
  return 0;
}

//create (double*) measurement->single_measurement from (Vec) measurement->y
int
measurement_petsc_vector_disassemble(MEAUSREMENT_TYPE* measurement, APP_CONTEXT_TYPE* app_context)
{
  static int* index = NULL;
  if(index==NULL) {
    index = new int[measurement->n_y];
    for (int i = 0; i<measurement->n_y; i++) index[i]=i;
  }
  app_context->ierr=VecGetValues(measurement->y,measurement->n_y, index, measurement->single_measurement);
  app_context_report_error(app_context->ierr,"measurement_petsc_vector_disassemble: VecGetValues",app_context,__LINE__,__FUNCT__,__FILE__);
  app_context->ierr=VecAssemblyBegin(measurement->y);
  app_context_report_error(app_context->ierr,"measurement_petsc_vector_disassemble: VecAssemblyBegin",app_context,__LINE__,__FUNCT__,__FILE__);
  app_context->ierr=VecAssemblyEnd(measurement->y);  
  app_context_report_error(app_context->ierr,"measurement_petsc_vector_disassemble: VecAssemblyEnd",app_context,__LINE__,__FUNCT__,__FILE__);
  
  return 0;
}

//initialize the interpolation for all meausrements, from all_measurements and map_index2time
int
measurement_interpolation_initialize(MEAUSREMENT_TYPE* measurement, APP_CONTEXT_TYPE* app_context)
{
  interp_v(0, measurement->n_y, measurement->K, &measurement->map_index2time[1*measurement->K+0], measurement->all_measurements, -1, NULL);  
  return 0;
}

//finalize the interpolation for all meausrements
int
measurement_interpolation_finalize(MEAUSREMENT_TYPE* measurement, APP_CONTEXT_TYPE* app_context)
{
  interp_v(2, measurement->n_y, measurement->K, NULL, NULL, -1, NULL);
  return 0;
}

//get the meausrement at any given time point, stored in single_measurement
int
measurement_get(double time, MEAUSREMENT_TYPE* measurement, APP_CONTEXT_TYPE* app_context)
{
  if(1){ //interpolation
    interp_v(1, measurement->n_y, measurement->K, NULL, NULL,time, measurement->single_measurement);    
    return 0;
  }
}

//load measurement matrix from a data file into measurement->H_values;
int 
measurement_io_load_measurement_matrix_dat(char* filename, MEAUSREMENT_TYPE* measurement, APP_CONTEXT_TYPE* app_context)
{
  int measurement_matrix_size_1 = 0, measurement_matrix_size_2 = 0;
  
  if(measurement->H_values){PetscPrintf(app_context->mpi_comm_global, "measurement matrix was already loaded (measurement->H_values=%x), so no need to load it from %s again\n", measurement->H_values,filename); assert(0);}

  measurement->H_values = load_matrix_2D_ascii(filename, measurement_matrix_size_1, measurement_matrix_size_2);  
  assert(measurement_matrix_size_2==app_context->n_x+app_context->n_theta);
  assert(measurement_matrix_size_1==measurement->n_y);
  measurement->n_z = measurement_matrix_size_2;
  
  return 0;
}
//load measurement error covariance matrix from a data file into measurement->R_d_values;
int 
measurement_io_load_measurement_error_covariance_matrix_dat(char* filename, MEAUSREMENT_TYPE* measurement, APP_CONTEXT_TYPE* app_context)
{
  int measurement_matrix_size_1 = 0, measurement_matrix_size_2 = 0;
  
  if(measurement->R_d_values){PetscPrintf(app_context->mpi_comm_global, "measurement error covariance matrix was already loaded (measurement->R_d_values=%x), so no need to load it from %s again\n", measurement->R_d_values,filename); assert(0);}

  measurement->R_d_values = load_matrix_2D_ascii(filename, measurement_matrix_size_1, measurement_matrix_size_2);  
  assert(measurement_matrix_size_2==measurement->n_y);
  assert(measurement_matrix_size_1==measurement->n_y);
  return 0;
}
//prepare measurement (IO, create petsc data structure etc.)
MEAUSREMENT_TYPE* 
measurement_create_and_initialize(APP_CONTEXT_TYPE* app_context)
{
  MEAUSREMENT_TYPE* measurement = new MEAUSREMENT_TYPE;
  memset(measurement, 0, sizeof(MEAUSREMENT_TYPE));
  if(0){
    measurement->n_y = 2;
    measurement->n_z = 4;
    measurement->single_measurement=new double[measurement->n_y];
    return measurement;
  }
  
  //load all measurements
  measurement_io_read_all_measurements(measurement, app_context);
  
  //load measurement matrix
  measurement_io_load_measurement_matrix_dat(app_context->measurement_matrix_file_name, measurement, app_context);
  { //create measurement->H (Mat type)
    int* index_row = new int[measurement->n_y];
    int* index_col = new int[measurement->n_z];
    MatCreate(app_context->mpi_comm_global, &measurement->H);
    MatSetSizes(measurement->H, PETSC_DECIDE, PETSC_DECIDE, measurement->n_y, measurement->n_z);
    MatSetFromOptions(measurement->H);
    for (int i = 0; i<measurement->n_y; i++) index_row[i]=i;
    for (int i = 0; i<measurement->n_z; i++) index_col[i]=i;
  
    app_context->ierr=MatSetValues(measurement->H,measurement->n_y,index_row,measurement->n_z, index_col,measurement->H_values,ADD_VALUES);
    MatAssemblyBegin(measurement->H,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(measurement->H,MAT_FINAL_ASSEMBLY);
    app_context_report_error(app_context->ierr,"measurement_create_and_initialize: MatSetValues",app_context,__LINE__,__FUNCT__,__FILE__);
    free(index_row); free(index_col);
  }


  //load measurement error covariance matrix
  measurement_io_load_measurement_error_covariance_matrix_dat(app_context->measurement_error_covariance_matrix_file_name, measurement, app_context);
  { //create measurement->R_d (Mat type)
    int* index_row = new int[measurement->n_y];
    MatCreateMPIDense(app_context->mpi_comm_global, PETSC_DECIDE, PETSC_DECIDE, measurement->n_y, measurement->n_y, PETSC_NULL, &measurement->R_d);
    //MatCreate(app_context->mpi_comm_global, &measurement->R_d);
    //MatSetSizes(measurement->R_d, PETSC_DECIDE, PETSC_DECIDE, measurement->n_y, measurement->n_y);
    MatSetFromOptions(measurement->R_d);
    for (int i = 0; i<measurement->n_y; i++) index_row[i]=i;
  
    app_context->ierr=MatSetValues(measurement->R_d,measurement->n_y,index_row,measurement->n_y, index_row,measurement->R_d_values,ADD_VALUES);
    MatAssemblyBegin(measurement->R_d,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(measurement->R_d,MAT_FINAL_ASSEMBLY);
    app_context_report_error(app_context->ierr,"measurement_create_and_initialize: MatSetValues",app_context,__LINE__,__FUNCT__,__FILE__);
    free(index_row);
  }

   //Create measurement error covariance matrix inverse inv_R_d
   MatDuplicate(measurement->R_d, MAT_COPY_VALUES, &measurement->inv_R_d);

  //prepare the interpolation
  measurement_interpolation_initialize(measurement, app_context);
 
  //usage:
  //ret = measurement_get(0, measurement, app_context);
  //if(ret==0) measurement_petsc_vector_assemble(measurement,app_context);
  
  //load the measurement error covariance R_d
  
  
  return measurement;
}


//-------------------------------------------------------------------------------------------------
//---models routines---
MODEL_TYPE* 
model_create_and_initialize(MODEL_TYPE* (*model_opeartor)(MODEL_TYPE*, APP_CONTEXT_TYPE* ), APP_CONTEXT_TYPE* app_context)
{
  MODEL_TYPE* model = NULL;
  app_context->ierr = PetscMalloc(sizeof(MODEL_TYPE), (void**) &model);
  app_context_report_error(app_context->ierr,"app_context_create_and_initialize: PetscMalloc",app_context,__LINE__,__FUNCT__,__FILE__);
    
  model_set_opeartor_and_intial_condition(model, model_opeartor,app_context);
    
  return model;  
}


// Initialize 
//1) model operator (model->f), 
//2) dimension (model->n_x, model->n_theta), 
//3) time stamp (),
//4) state (model->x_k_hat,model->theta_k_hat) ,
//5) previous state (model->x_k_minus_1_hat,model->theta_k_minus_1_hat) 
MODEL_TYPE* 
model_set_opeartor_and_intial_condition(MODEL_TYPE* model,MODEL_TYPE* (*model_opeartor)(MODEL_TYPE*, APP_CONTEXT_TYPE*), APP_CONTEXT_TYPE* app_context)
{
  PetscInt i=0;
  
  if(1){  //model 01: 2 parameter, 2 state,  simplest model
    
    model->f = *model_opeartor;  
    model->k = 0;
    model->t = 0.0;
    
    model->n_x=app_context->n_x;
//    app_context->ierr=VecCreateMPIWithArray(app_context->mpi_comm_global, model->n_x, PETSC_DECIDE, app_context->x_0, &model->x_k_hat);
//    app_context->ierr=VecCreateMPIWithArray(app_context->mpi_comm_global, model->n_x, model->n_x, app_context->x_0, &model->x_k_hat);

    app_context->ierr=VecCreate(app_context->mpi_comm_global, &model->x_k_hat);
    app_context->ierr=VecSetSizes(model->x_k_hat, PETSC_DECIDE,app_context->n_x);
    app_context->ierr=VecSetFromOptions(model->x_k_hat);
    for(i=0;i<model->n_x; i++) app_context->tmp_x_indices[i]=i;
    VecSetValues(model->x_k_hat,model->n_x,app_context->tmp_x_indices,app_context->x_0,ADD_VALUES); 
    VecAssemblyBegin(model->x_k_hat);
    VecAssemblyEnd(model->x_k_hat);

    app_context_report_error(app_context->ierr,"model_set_opeartor_and_intial_condition: VecCreateMPIWithArray",app_context,__LINE__,__FUNCT__,__FILE__);

    model->n_theta=app_context->n_theta;
//app_context->ierr=VecCreateMPIWithArray(app_context->mpi_comm_global, model->n_theta, PETSC_DECIDE, app_context->theta_0, &model->theta_k_hat);
//    app_context->ierr=VecCreateMPIWithArray(app_context->mpi_comm_global, model->n_theta, model->n_theta, app_context->theta_0, &model->theta_k_hat);
    app_context->ierr=VecCreate(app_context->mpi_comm_global, &model->theta_k_hat);
    app_context->ierr=VecSetSizes(model->theta_k_hat, PETSC_DECIDE,app_context->n_theta);
    app_context->ierr=VecSetFromOptions(model->theta_k_hat);
    for(i=0;i<model->n_theta; i++) app_context->tmp_theta_indices[i]=i;
    VecSetValues(model->theta_k_hat,model->n_theta,app_context->tmp_theta_indices,app_context->theta_0,ADD_VALUES); 
    VecAssemblyBegin(model->theta_k_hat);
    VecAssemblyEnd(model->theta_k_hat);
    app_context_report_error(app_context->ierr,"model_set_opeartor_and_intial_condition: VecCreateMPIWithArray",app_context,__LINE__,__FUNCT__,__FILE__);
    
    if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
      app_context_report_info(INFO_PETSC_REALARRAY,app_context->n_x,  app_context->x_0,"model_set_opeartor_and_intial_condition, app_context->x_0", app_context);        
      app_context_report_info(INFO_PETSC_REALARRAY, app_context->n_theta, app_context->theta_0,"model_set_opeartor_and_intial_condition, app_context->theta_0", app_context);    
      app_context_report_info(INFO_PETSC_VEC,1,  &model->x_k_hat,"model_set_opeartor_and_intial_condition, model->x_0_hat", app_context);        
      app_context_report_info(INFO_PETSC_VEC, 1, &model->theta_k_hat,"model_set_opeartor_and_intial_condition, model->theta_0_hat", app_context);    
    }

    
  }else if(0){ //model 02: 
    
    app_context_report_error(app_context->ierr,"model_set_opeartor_and_intial_condition: not implemented yet",app_context,__LINE__,__FUNCT__,__FILE__);
  }else{
    
    app_context_report_error(app_context->ierr,"model_set_opeartor_and_intial_condition: not implemented yet",app_context,__LINE__,__FUNCT__,__FILE__);
  }
  
  //other standard setups
  VecDuplicate(model->x_k_hat, &model->x_k_minus_1_hat);
  VecDuplicate(model->theta_k_hat, &model->theta_k_minus_1_hat);
  
  //diagnosis
  if(app_context->diagnosisMode & DIAGNOSIS_DISP_INITIAL_CONF){
    app_context_report_info(INFO_PETSC_REAL, 1, &model->t,"model_set_opeartor_and_intial_condition, model->t", app_context);
    app_context_report_info(INFO_PETSC_INT, 1, &model->k,"model_set_opeartor_and_intial_condition, model->k", app_context);

    app_context_report_info(INFO_PETSC_INT, 1, &model->n_x,"model_set_opeartor_and_intial_condition, model->n_x", app_context);
    app_context_report_info(INFO_PETSC_INT, 1, &model->n_theta,"model_set_opeartor_and_intial_condition, model->n_theta", app_context);
    
    app_context_report_info(INFO_PETSC_VEC, 1, &model->theta_k_hat,"model_set_opeartor_and_intial_condition, model->theta_0_hat", app_context);    
    app_context_report_info(INFO_PETSC_VEC, 1, &model->theta_k_minus_1_hat,"model_set_opeartor_and_intial_condition, model->theta_0_minus_1_hat", app_context);    
  }
  if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
    app_context_report_info(INFO_PETSC_VEC,1,  &model->x_k_hat,"model_set_opeartor_and_intial_condition, model->x_0_hat", app_context);        
    app_context_report_info(INFO_PETSC_VEC, 1, &model->x_k_minus_1_hat,"model_set_opeartor_and_intial_condition, model->x_0_minus_1_hat", app_context);    
  }
  
  return model;  
}


//-------------------------------------------------------------------------------------------------
//---filter routines---

//create and initialize the filter
FILTER_TYPE* 
filter_create_and_initialize(APP_CONTEXT_TYPE* app_context,MODEL_TYPE* model)
{
  FILTER_TYPE* filter = NULL;
  app_context->ierr = PetscMalloc(sizeof(FILTER_TYPE), (void**) &filter);
  app_context_report_error(app_context->ierr,"filter_create_and_initialize: PetscMalloc",app_context,__LINE__,__FUNCT__,__FILE__);
  
  filter->n_x = app_context->n_x;
  filter->n_theta = app_context->n_theta;
  filter = filter_set_intial_condition(filter,model, app_context);

  return filter;
}


// Initialize 
//1) covariance for state, measurement (portion_afilter->), 
//2) dimension (filter->n_x, filter->n_theta), 
//3) time stamp (filter->),
//4) state (K) 
FILTER_TYPE* 
filter_set_intial_condition(FILTER_TYPE* filter, MODEL_TYPE* model, APP_CONTEXT_TYPE* app_context)
{
  int i=0;

  //allocate space and initialize covariance matrix (filter subtype)
  filter->subfilter = filter_subtype_create_and_initialize(filter,model, app_context);

  //model->n_ensem = 2;
  //set in filter_subtype_set_intial_condition()

  //allocate space for sigma points (in model struture)  
  PetscMalloc(sizeof(Vec*)*model->n_ensem, &model->x_k_minus_1);
  PetscMalloc(sizeof(Vec*)*model->n_ensem, &model->theta_k_minus_1);
  PetscMalloc(sizeof(Vec*)*model->n_ensem, &model->x_k);
  PetscMalloc(sizeof(Vec*)*model->n_ensem, &model->theta_k);

  for(i=0; i<model->n_ensem; i++){
    VecCreate(app_context->mpi_comm_global, &model->x_k_minus_1[i]);
    VecCreate(app_context->mpi_comm_global, &model->theta_k_minus_1[i]);
    VecCreate(app_context->mpi_comm_global, &model->x_k[i]);
    VecCreate(app_context->mpi_comm_global, &model->theta_k[i]);
    
    VecSetSizes(model->x_k_minus_1[i], PETSC_DECIDE,app_context->n_x);
    VecSetSizes(model->theta_k_minus_1[i], PETSC_DECIDE,app_context->n_theta);
    VecSetSizes(model->x_k[i], PETSC_DECIDE,app_context->n_x);
    VecSetSizes(model->theta_k[i], PETSC_DECIDE,app_context->n_theta);
    
    VecSetFromOptions(model->x_k_minus_1[i]);
    VecSetFromOptions(model->theta_k_minus_1[i]);
    VecSetFromOptions(model->x_k[i]);
    VecSetFromOptions(model->theta_k[i]);
  }
    
  return filter;
}

//create and initialize the filter
FILTER_SUBTYPE* 
filter_subtype_create_and_initialize(FILTER_TYPE* filter, MODEL_TYPE* model, APP_CONTEXT_TYPE* app_context)
{
  filter->subfilter = NULL;
  app_context->ierr = PetscMalloc(sizeof(FILTER_SUBTYPE), (void**) &filter->subfilter);
  app_context_report_error(app_context->ierr,"filter_set_intial_condition: PetscMalloc",app_context,__LINE__,__FUNCT__,__FILE__);
  
  filter->subfilter = filter_subtype_set_intial_condition(filter,model, app_context);
  return filter->subfilter;
}

// Initialize the filter subtype
//  PetscInt filter_type;
//  PetscInt r; //rank ( r>=(n_theat>0)?n_theat:1, r<= n_theta+n_x;
//  Mat S_x;
//  Mat S_theta;  
FILTER_SUBTYPE* 
filter_subtype_set_intial_condition(FILTER_TYPE* filter, MODEL_TYPE* model, APP_CONTEXT_TYPE* app_context)
{
  //local values
  PetscScalar S_x_value; //all 1  default
  PetscScalar* S_theta_values; //1,1,1... default
  PetscInt    i; //temperory usage
  
  //check pointer
  app_context_report_error(filter==NULL,"FILTER_TYPE* filter is null",app_context,__LINE__,__FUNCT__,__FILE__);
  app_context_report_error(model==NULL,"MODEL_TYPE* model is null",app_context,__LINE__,__FUNCT__,__FILE__);
    
  FILTER_SUBTYPE* subfilter=filter->subfilter;
  app_context_report_error(subfilter==NULL,"filter->subfilter is null",app_context,__LINE__,__FUNCT__,__FILE__);
  
  //initialize the subfilter
  subfilter->filter_type = FILTER_TYPE_RUKF; //hard-code rUKF
  
  // Initialize
  //  Mat S_x;
  //  Mat S_theta; 
  switch(subfilter->filter_type)
  {
    case FILTER_TYPE_RUKF: //reduced-order UKF
      //choose rank
      if(model->n_theta==0){ //state estimation only
        subfilter->r = model->n_x - 1;
      }//chanllenge it by setting it to be not equal to n_x
      else{      //joint estimation
        subfilter->r = model->n_theta + 0 ; //chanllenge it by setting it to be not equal to n_theta
      }
      model->n_ensem = subfilter->r * 2;
      
      //allocate space for S_x, S_theta; app_context->S_ensemble, app_context->P_z, app_context->S_z
      app_context->ierr=MatCreateMPIDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, model->n_x,subfilter->r, PETSC_NULL, &subfilter->S_x); 
      app_context_report_error(app_context->ierr,"subfilter->S_x     MatCreateMPIDense",app_context,__LINE__,__FUNCT__,__FILE__);
      app_context->ierr=MatCreateMPIDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, model->n_theta,subfilter->r, PETSC_NULL, &subfilter->S_theta); 
      app_context_report_error(app_context->ierr,"subfilter->S_theta MatCreateMPIDense",app_context,__LINE__,__FUNCT__,__FILE__);
      app_context->ierr=MatCreateMPIDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, model->n_x+model->n_theta,model->n_ensem, PETSC_NULL, &app_context->S_ensemble); 
      app_context_report_error(app_context->ierr,"app_context->S_ensemble     MatCreateMPIDense",app_context,__LINE__,__FUNCT__,__FILE__);
      app_context->ierr=MatCreateMPIDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, model->n_x+model->n_theta,model->n_ensem, PETSC_NULL, &app_context->S_ensemble2); 
      app_context_report_error(app_context->ierr,"app_context->S_ensemble2     MatCreateMPIDense",app_context,__LINE__,__FUNCT__,__FILE__);
      app_context->ierr=MatCreateMPIDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, model->n_x+model->n_theta,model->n_x+model->n_theta, PETSC_NULL, &app_context->P_z); 
      app_context_report_error(app_context->ierr,"app_context->P_z            MatCreateMPIDense",app_context,__LINE__,__FUNCT__,__FILE__);
      
      app_context->ierr=MatCreateMPIDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, model->n_x+model->n_theta,subfilter->r, PETSC_NULL, &app_context->S_z); 
      app_context_report_error(app_context->ierr,"app_context->S_z            MatCreateMPIDense",app_context,__LINE__,__FUNCT__,__FILE__);
      MatAssemblyBegin(app_context->S_z  ,MAT_FINAL_ASSEMBLY ); //for later svd usage
      MatAssemblyEnd  (app_context->S_z  ,MAT_FINAL_ASSEMBLY ); //for later svd usage

      MatSetFromOptions( subfilter->S_x );
      MatSetFromOptions( subfilter->S_theta );
      MatSetFromOptions( app_context->S_ensemble );
      MatSetFromOptions( app_context->P_z );
      MatSetFromOptions( app_context->S_z );
      

      //initialization S_x and S_theta
     //values
     app_context->ierr=PetscMalloc(sizeof(PetscScalar)*model->n_theta, &S_theta_values);
     for(i=0;i<model->n_theta; i++){
      S_theta_values[i] = 1.0;
     }
     S_x_value = 1.0;
     
     //set (set diagonal element of S_theta first)
      for(i=0; i<model->n_theta; i++){ 
          MatSetValue(subfilter->S_theta, i,i, S_theta_values[i], INSERT_VALUES);   
      }
      MatAssemblyBegin(subfilter->S_theta,MAT_FINAL_ASSEMBLY );
      MatAssemblyEnd  (subfilter->S_theta,MAT_FINAL_ASSEMBLY );

      //(then set diagonal element of S_x )
      for(i=0; i < subfilter->r-model->n_theta; i++){ 
        MatSetValue(subfilter->S_x, i, i+model->n_theta, S_x_value, INSERT_VALUES);   
      }
      MatAssemblyBegin(subfilter->S_x,MAT_FINAL_ASSEMBLY );
      MatAssemblyEnd  (subfilter->S_x,MAT_FINAL_ASSEMBLY );

      //diagnosis
      if(app_context->diagnosisMode & DIAGNOSIS_DISP_INITIAL_CONF){
        app_context_report_info(INFO_PETSC_MAT, 1, &subfilter->S_x,"subfilter->S_x_0", app_context);
        app_context_report_info(INFO_PETSC_MAT, 1, &subfilter->S_theta,"subfilter->S_theta_0", app_context);
      }
            
      break;
    default: //others
        app_context_report_error(1,"filter->filter_type not implemented",app_context,__LINE__,__FUNCT__,__FILE__);      
  }
  
  return subfilter;
}

void
app_context_create_and_initialize_dependently(APP_CONTEXT_TYPE* app_context, FILTER_TYPE* filter, MODEL_TYPE* model, MEAUSREMENT_TYPE* measurement)
{
    //S_bar
    app_context->S_bar = petsc_create_matrix(1, measurement->n_y, filter->subfilter->r,app_context->mpi_comm_global);
    //S_bar2
    MatDuplicate(app_context->S_bar, MAT_COPY_VALUES, &app_context->S_bar2);
    //tmp_r
    app_context->tmp_r = petsc_create_matrix(1, filter->subfilter->r, filter->subfilter->r,app_context->mpi_comm_global);
    //tmp_nz_by_r
    app_context->tmp_nz_by_r = petsc_create_matrix(1, app_context->n_x+app_context->n_theta, filter->subfilter->r,app_context->mpi_comm_global);
    //tmp_r_by_ny
    app_context->tmp_r_by_ny = petsc_create_matrix(1, filter->subfilter->r, measurement->n_y, app_context->mpi_comm_global);
    //tmp_y
    app_context->tmp_y = petsc_create_vector(1, measurement->n_y, app_context->mpi_comm_global);
    
     //Kalman gain matrix: K
    filter->K = petsc_create_matrix(1, filter->n_x+filter->n_theta,measurement->n_y,app_context->mpi_comm_global);

}

// generate sigma points from filter->covariance and model->x_k_hat
MODEL_TYPE*
filter_generate_sigma_points(FILTER_TYPE* filter, MODEL_TYPE* model, APP_CONTEXT_TYPE* app_context)
{
  int i;
  PetscScalar localScale;  //localization parameter
  
  app_context_assert(model->n_ensem>=1, "", app_context);
  if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
    app_context_report_info(INFO_PETSC_INT, 1, &model->k,"filter_generate_sigma_points, model->k", app_context);
    app_context_report_info(INFO_PETSC_REAL, 1, &model->t,"filter_generate_sigma_points, model->t", app_context);
    app_context_report_info(INFO_PETSC_INT, 1, &i,"i",  app_context);
    app_context_report_info(INFO_PETSC_VEC, 1, &model->x_k_hat,"filter_generate_sigma_points, model->x_k_hat", app_context);    
    app_context_report_info(INFO_PETSC_VEC, 1, &model->theta_k_hat,"filter_generate_sigma_points, model->theta_k_hat", app_context);    
    app_context_report_info(INFO_PETSC_MAT, 1, &filter->subfilter->S_x,"/*filter_generate_sigma_points, filter->subfilter->S_x*/", app_context);    
    app_context_report_info(INFO_PETSC_MAT, 1, &filter->subfilter->S_theta,"filter_generate_sigma_points, filter->subfilter->S_theta", app_context);    
  }

  if(0){ //temp: assume filter->covariance is 0
    for(i=0; i<model->n_ensem; i++){  //model->n_ensem>=1
      VecCopy(model->x_k_hat, model->x_k[i]);
      VecCopy(model->theta_k_hat, model->theta_k[i]);
    }
  }
      
  if(1){ //generate sigma points from filter->covariance and model->x_k_hat
    if(filter->subfilter->filter_type!=FILTER_TYPE_RUKF){
    //calculate the square-root matrix of filter->covariance
    //1. allocate memory
    //2. do svd, storge them in filter->subfilter->S_x, filter->subfilter->S_theta
    }
    
    //3. piling up first r eig vectors weighted by eig values
    localScale = 1; //hard coded 1, localization parameter
    for(i=0; i<filter->subfilter->r; i++){
      MatGetColumnVector(filter->subfilter->S_x,app_context->tmp_x,i);
      VecWAXPY(model->x_k[i], 1*localScale, app_context->tmp_x, model->x_k_hat);
      VecWAXPY(model->x_k[i+filter->subfilter->r], -1*localScale, app_context->tmp_x, model->x_k_hat);

      //TODO!!Arguments are incompatible!
      //[0]PETSC ERROR: Incompatible vector global lengths!
      //solution: generalize MatGetColumnVector-->get the global vector instead of only the local portion!
      
      MatGetColumnVector(filter->subfilter->S_theta,app_context->tmp_theta,i);
      VecWAXPY(model->theta_k[i], 1*localScale, app_context->tmp_theta, model->theta_k_hat);
      VecWAXPY(model->theta_k[i+filter->subfilter->r], -1*localScale, app_context->tmp_theta, model->theta_k_hat);
      
      //diagnosis
      if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
        app_context_report_info(INFO_PETSC_VEC, 1, &model->theta_k[i],"filter_generate_sigma_points, model->theta_k[i]", app_context);    
        app_context_report_info(INFO_PETSC_VEC,1,  &model->x_k[i],"filter_generate_sigma_points, model->x_k[i]", app_context);        
      }
    }
  }
      
  return model;
}

FILTER_TYPE* 
filter_time_update(FILTER_TYPE* filter, MODEL_TYPE* model, APP_CONTEXT_TYPE* app_context)
{
  int i=0;
  PetscErrorCode    ierr;
  PetscInt nEigs; //number of eig vectors to be retrieved

  //sigma points
  filter_generate_sigma_points(filter, model, app_context);
  
  //model evaluation
  model->f(model, app_context);
  
  //compute mean
  ierr=VecAXPY(model->x_k_hat,-1, model->x_k_hat); //nullize
  for(i=0; i<model->n_ensem; i++){
    ierr=VecAXPY(model->x_k_hat,1, model->x_k[i]);
  }
  ierr=VecScale(model->x_k_hat, 1.0/model->n_ensem );

  ierr=VecAXPY(model->theta_k_hat,-1, model->theta_k_hat);//nullize
  for(i=0; i<model->n_ensem; i++){
    ierr=VecAXPY(model->theta_k_hat,1, model->theta_k[i]);
  }
  ierr=VecScale(model->theta_k_hat, 1.0/model->n_ensem );


  //1. compute covariance (filter->subfilter->S_x and filter->subfilter->S_theta)
  //todo: add inflation parameter later
  //1.1. join them together  (use MatSetColumnVec) 
  for(i=0; i<model->n_ensem; i++){
    ierr=VecWAXPY(app_context->tmp_x, -1.0, model->x_k_hat, model->x_k[i]); app_context_report_error(ierr,"VecWAXPY",app_context,__LINE__,__FUNCT__,__FILE__);
    ierr=MatSetColumnVec(&app_context->S_ensemble, &app_context->tmp_x, app_context->tmp_x_indices, i, 0, 0, model->n_x); 
    ierr=VecWAXPY(app_context->tmp_theta, -1.0, model->theta_k_hat, model->theta_k[i]);
    app_context_report_error(ierr,"VecWAXPY",app_context,__LINE__,__FUNCT__,__FILE__);
    ierr=MatSetColumnVec(&app_context->S_ensemble, &app_context->tmp_theta, app_context->tmp_theta_indices, i, model->n_x, 0, model->n_theta); 
    if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
      app_context_report_info(INFO_PETSC_VEC,1,  &model->x_k[i],"compute covariance, filter_time_update, model->x_k[i]", app_context);        
      app_context_report_info(INFO_PETSC_VEC, 1, &model->theta_k[i],"compute covariance, filter_time_update, model->theta_k[i]", app_context);    
    }
  }
  MatAssemblyBegin(app_context->S_ensemble,MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd  (app_context->S_ensemble,MAT_FINAL_ASSEMBLY );
  if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
    app_context_report_info(INFO_PETSC_MAT, 1, &app_context->S_ensemble,"compute covariance, filter_time_update, app_context->S_ensemble", app_context);    
  }
      
  //1.2. P_z=S_ensemble*S_ensemble' 
  MatCopy(app_context->S_ensemble,app_context->S_ensemble2,SAME_NONZERO_PATTERN);
  if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
    app_context_report_info(INFO_PETSC_MAT, 1, &app_context->S_ensemble,"compute covariance, filter_time_update, app_context->S_ensemble", app_context);    
    app_context_report_info(INFO_PETSC_MAT, 1, &app_context->S_ensemble2,"compute covariance, filter_time_update, app_context->S_ensemble2", app_context);    
  }
  MatMultilXYT(app_context->S_ensemble, app_context->S_ensemble2, &app_context->P_z, app_context->mpi_comm_global);
  if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
    app_context_report_info(INFO_PETSC_MAT, 1, &app_context->P_z,"compute covariance, filter_time_update, app_context->P_z", app_context);    
  }
//   
  //2. do svd
  nEigs = filter->subfilter->r;
  MatSVD(&app_context->P_z, app_context->mpi_comm_global, &nEigs, NULL, &app_context->S_z, NULL,1);
  if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
    app_context_report_info(INFO_PETSC_INT, 1, &nEigs,"number of computed nEigs vectors, filter_time_update, nEigs", app_context);    
    app_context_report_info(INFO_PETSC_MAT, 1, &app_context->S_z,"result of SVD, filter_time_update, app_context->S_z", app_context);    
  }
  if(nEigs!=filter->subfilter->r){
    PetscPrintf(app_context->mpi_comm_global, "app_context->P_z only has %d eig vectors (filter->subfilter->r=%d)!!\n", nEigs, filter->subfilter->r);
    assert(0); //static rank
  }

  //3. separate them (use MatSetColumnVec or MatSetMatBlock)
  MatSetMatBlock(app_context->mpi_comm_global, &filter->subfilter->S_x,     app_context->S_z, 0, 0, 0,          0, model->n_x,    filter->subfilter->r);
  MatSetMatBlock(app_context->mpi_comm_global, &filter->subfilter->S_theta, app_context->S_z, 0, 0, model->n_x, 0, model->n_theta,filter->subfilter->r);
  if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
    app_context_report_info(INFO_PETSC_MAT, 1, &filter->subfilter->S_x,"S_x, after SVD, filter_time_update, filter->subfilter->S_x", app_context);    
    app_context_report_info(INFO_PETSC_MAT, 1, &filter->subfilter->S_theta,"S_theta, after SVD, filter_time_update, filter->subfilter->S_theta", app_context);    
    app_context_report_info(INFO_PETSC_INT, 1, &nEigs,"==exit of filter_time_update==", app_context);    
  }
  
  app_context_report_error(ierr,"exit of filter_time_update",app_context,__LINE__,__FUNCT__,__FILE__);
  return filter;
}
  
FILTER_TYPE* 
filter_measurement_update(FILTER_TYPE* filter, MODEL_TYPE* model, MEAUSREMENT_TYPE* measurement, APP_CONTEXT_TYPE* app_context)
{
  PetscErrorCode        ierr=0;
  app_context_report_error(ierr,"entry of filter_measurement_update",app_context,__LINE__,__FUNCT__,__FILE__);
  
  if(0){ //generate synthetic data
    VecDuplicate(model->x_k_hat,&measurement->y);
    measurement_petsc_vector_disassemble(measurement,app_context);
    sprintf(app_context->measurement_data_file_name, "%s%04d.dat", app_context->measurement_dir_name, model->k); //assume model->k==model->t
    measurement_io_write_single_measurement_dat(app_context->measurement_data_file_name, measurement, app_context); 
    return filter;
  }
  
  //2. covariance update and gain matrix caculation----------
  //        S_bar = H*S_minus;
  //        inv_R_d = inv(R_d);
  //        tmp = eye(npdof,npdof) + S_bar'*inv_R_d*S_bar;
  //        inv_tmp = inv(tmp);
  //        try
  //         S_plus = S_minus * chol(inv_tmp,'lower');
  //        catch ME
  //            disp('exception:Matrix must be positive definite with real diagonal.');           
  //            keyboard;
  //        end
  //        P_plus = S_plus * S_plus';
  //        K_x = S_minus * inv_tmp * S_bar'*inv_R_d;
  
//          S_minus (app_context->S_z) = [filter->subfilter->S_x; filter->subfilter->S_theta];
  for(int i=0; i<filter->subfilter->r; i++){
    ierr=MatGetColumnVector(filter->subfilter->S_x, app_context->tmp_x, i); 
    app_context_report_error(ierr,"",app_context,__LINE__,__FUNCT__,__FILE__);
    ierr=MatSetColumnVec(&app_context->S_z, &app_context->tmp_x, app_context->tmp_x_indices, i, 0, 0, model->n_x); 
    app_context_report_error(ierr,"",app_context,__LINE__,__FUNCT__,__FILE__);
    ierr=MatGetColumnVector(filter->subfilter->S_theta, app_context->tmp_theta, i); 
    app_context_report_error(ierr,"",app_context,__LINE__,__FUNCT__,__FILE__);
    ierr=MatSetColumnVec(&app_context->S_ensemble, &app_context->tmp_theta, app_context->tmp_theta_indices, i, model->n_x, 0, model->n_theta);  
    app_context_report_error(ierr,"",app_context,__LINE__,__FUNCT__,__FILE__);
    if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
      app_context_report_info(INFO_PETSC_VEC,1,  &app_context->tmp_x,"filter_measurement_update, filter->subfilter->S_x(i,:)", app_context);        
      app_context_report_info(INFO_PETSC_VEC, 1, &app_context->tmp_theta,"filter_measurement_update, filter->subfilter->S_theta(i,:)", app_context);    
    }
  }
  MatAssemblyBegin(app_context->S_z,MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd  (app_context->S_z,MAT_FINAL_ASSEMBLY );
  if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
    app_context_report_info(INFO_PETSC_MAT, 1, &app_context->S_z,"filter_measurement_update, app_context->S_z", app_context);    
  }

//          S_bar(app_context->S_bar) = H*S_minus;
  MatMultilXY(measurement->H, app_context->S_z, &app_context->S_bar, app_context->mpi_comm_global);

//        inv_R_d = inv(R_d);
  MatInverse(measurement->R_d, &measurement->inv_R_d, app_context->mpi_comm_global);
//        tmp = eye(npdof,npdof) + S_bar'*inv_R_d*S_bar;
  MatMultilXY(measurement->inv_R_d, app_context->S_bar, &app_context->S_bar2, app_context->mpi_comm_global);
  MatMultilXTY(app_context->S_bar, app_context->S_bar2, &app_context->tmp_r,  app_context->mpi_comm_global);
  MatShift(app_context->tmp_r,1);
//        inv_tmp = inv(tmp);
  Mat tmp = app_context->tmp_r;
  MatInverse(app_context->tmp_r, &app_context->tmp_r,  app_context->mpi_comm_global);
  MatDestroy(tmp);
//        K_x = S_minus * inv_tmp * S_bar'*inv_R_d;
  MatMultilXY(app_context->S_z, app_context->tmp_r, &app_context->tmp_nz_by_r, app_context->mpi_comm_global);
  MatMultilXTY(app_context->S_bar, measurement->inv_R_d, &app_context->tmp_r_by_ny, app_context->mpi_comm_global);
  MatMultilXY(app_context->tmp_nz_by_r,app_context->tmp_r_by_ny, &filter->K, app_context->mpi_comm_global);
//        S_plus(app_context->S_z) = S_minus * chol(inv_tmp,'lower');
  MatFactorize(app_context->tmp_r, app_context->mpi_comm_global, 1);
  MatMultilXY (app_context->S_z, app_context->tmp_r, &app_context->tmp_nz_by_r, app_context->mpi_comm_global);
  MatCopy(app_context->tmp_nz_by_r, app_context->S_z, SAME_NONZERO_PATTERN);
//        [filter->subfilter->S_x; filter->subfilter->S_theta]=S_minus (app_context->S_z)
  MatSetMatBlock(app_context->mpi_comm_global, &filter->subfilter->S_x    , app_context->S_z,0,0,0,         0, model->n_x,     filter->subfilter->r);
  MatSetMatBlock(app_context->mpi_comm_global, &filter->subfilter->S_theta, app_context->S_z,0,0,model->n_x,0, model->n_theta, filter->subfilter->r);
  
//(ignore)        P_plus = S_plus * S_plus';
  
//        clean up  
  MatDestroy(measurement->inv_R_d);

  
  //3. update state---------------------------------------
  //        y = model.u_t(:);                  //synthetic measurement
  //        correction =measurementAmplifyFactor * K_x*(H*z-y);
  //        tmp = model.parameters(:) + correction(ndof+1:nadof);

//        y = model.u_t(:);                  //synthetic measurement
  //get the measurement
  int ret = measurement_get(model->k, measurement, app_context);
  if(ret==0) measurement_petsc_vector_assemble(measurement,app_context);

//        z(app_context->tmp_z) = [model->x_k_hat;model->theta_k_hat];
  VecSetVecBlock(app_context->mpi_comm_global, &app_context->tmp_z, model->x_k_hat, 0, 0, model->n_x);
  VecSetVecBlock(app_context->mpi_comm_global, &app_context->tmp_z, model->theta_k_hat, model->n_x, 0, model->n_theta);
// H*z
  MatMultilXY(measurement->H, app_context->S_z, &app_context->S_bar, app_context->mpi_comm_global);

//        correction(app_context->tmp_z) =measurementAmplifyFactor * K_x*(H*z-y);
  VecScale(measurement->y, -1);
  MatMultAdd(measurement->H, app_context->tmp_z,measurement->y,app_context->tmp_y);
  MatMult   (filter->K, app_context->tmp_y, app_context->tmp_z);
//        [model->x_k_hat;model->theta_k_hat] = z(app_context->tmp_z);
  VecSetVecBlock(app_context->mpi_comm_global, &model->x_k_hat, app_context->tmp_z, 0, 0, model->n_x);
  VecSetVecBlock(app_context->mpi_comm_global, &model->theta_k_hat, app_context->tmp_z, 0, model->n_x, model->n_theta);
            
  //4. filter monitor---------------------------------------
  if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
    app_context_report_info(INFO_PETSC_VEC, 1, &model->x_k_hat,"x_k_hat, after filter_measurement_update, model->x_k_hat", app_context);    
    app_context_report_info(INFO_PETSC_VEC, 1, &model->theta_k_hat,"theta_k_hat, after filter_measurement_update, model->theta_k_hat", app_context); 
    app_context_report_info(INFO_PETSC_MAT, 1, &filter->subfilter->S_x,"S_x, after filter_measurement_update, filter->subfilter->S_x", app_context);    
    app_context_report_info(INFO_PETSC_MAT, 1, &filter->subfilter->S_theta,"S_theta, after filter_measurement_update, filter->subfilter->S_theta", app_context);    
    app_context_report_info(INFO_PETSC_MAT, 1, &filter->K,"K, after filter_measurement_update, filter->K", app_context);    
  }

  app_context_report_error(ierr,"exit of filter_measurement_update",app_context,__LINE__,__FUNCT__,__FILE__);
  return filter;
}

//-------------------------------------------------------------------------------------------------
//---optimizer routines (not implemented yet. TODO!!)---

OPTIMIZER_TYPE* 
optimizer_create_and_initialize(APP_CONTEXT_TYPE* app_context)
{
  OPTIMIZER_TYPE* optimizer = NULL;
  app_context->ierr = PetscMalloc(sizeof(OPTIMIZER_TYPE), (void**) &optimizer);
  app_context_report_error(app_context->ierr,"filter_create_and_initialize: PetscMalloc",app_context,__LINE__,__FUNCT__,__FILE__);
  
  optimizer_set_intial_condition(optimizer,app_context);

  return optimizer;
}

OPTIMIZER_TYPE* 
optimizer_set_intial_condition(OPTIMIZER_TYPE*  optimizer, APP_CONTEXT_TYPE* app_context)
{
  return optimizer;
}


//---- Pesudo-code (matlab) of reduced-order UKF ----------------------------
//--- global variable declaration and initialization

//configurationfile='';
//estimationFromRealData = 0; %defualt value
//estimator = 'seik';


//---configuration
//        load(configurationfile);
// save configuration state


// set up simulation folders and default variables


// -------------- filter variables---------
//nadof = ndof+npdof;
//r=npdof;
//u_0 = zeros(ndof,1);%only the size matter, no the values
//H = [eye(ndof,ndof) zeros(ndof,npdof)];
//R_d = eye(ndof,ndof)*R_scale;
//model.u_t = u_0;

//parameters_truth = normalize('normalize', parameters_truth,model.parametersRange);

//if(exist('pertubulation','var')==1)
//  model.parameters = parameters_truth + (rand(size(parameters_truth))-0.5)*2 .* pertubulation;
//else
//  model.parameters = normalize('normalize', model.parameters,model.parametersRange);
//end

//S_plus_theta = diag(min(model.parameters, 1-model.parameters)) * S_plus_theta_scale;
//S_plus = [zeros(ndof,npdof); S_plus_theta];


//io_save_variable([objective '/-----configuration-----'],[]); %mark this simulation

//io_dynamical_simulation_snapshot('initialize',t:Delta_t/nIter:N+Delta_t*(nIter-1)/nIter);


// looping
//while 1
//---time update
        
        // model simulation
        // filter equations
//---meausrment update
            //load measurement   
            //covariance update and gain matrix caculation
            //update state
            
            //post-processing (disp,save)
    
//save parameter history