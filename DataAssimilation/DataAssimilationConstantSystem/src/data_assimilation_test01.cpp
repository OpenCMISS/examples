#define MEASUREMENT_DIR  "./" //where to find the measurement files, relative to the path of running the exe file

#include "data_assimilation_routines.h"

//model operator 01: simplest
MODEL_TYPE* 
model_opeartor_01(MODEL_TYPE* model, APP_CONTEXT_TYPE* app_context);

//model operator 01: constant system
//evalute model with (model->n_ensem) sigma points, 
//forwad the system state (model->theta_k[i],model->x_k[i]) and previous states (model->x_k_minus_1[i],model->theta_k_minus_1[i])  
//by 1 (model->k) time step or 0.01s (model->t)
MODEL_TYPE* 
model_opeartor_01(MODEL_TYPE* model, APP_CONTEXT_TYPE* app_context)
{
  PetscInt i=0;

  for(i=0; i<model->n_ensem; i++){
    //1. update the system state vector
    VecCopy(model->x_k[i],model->x_k_minus_1[i]);
    
    //2. update the system parameter vector
    VecCopy(model->theta_k[i],model->theta_k_minus_1[i]);
    //VecCopy(model->theta_k[i],model->x_k[i]);
    
    //debug diagnosis
    if(app_context->diagnosisMode & DIAGNOSIS_DISP_INITIAL_CONF){
      app_context_report_info(INFO_PETSC_INT, 1, &model->k,"===simulation==== \nmodel_opeartor_01, model->k", app_context);
      app_context_report_info(INFO_PETSC_REAL, 1, &model->t,"model_opeartor_01, model->t", app_context);

      app_context_report_info(INFO_PETSC_INT, 1, &i,"i",  app_context);
      
      app_context_report_info(INFO_PETSC_INT, 1, &model->n_x,"model_opeartor_01, model->n_x", app_context);
      app_context_report_info(INFO_PETSC_INT, 1, &model->n_theta,"model_opeartor_01, model->n_theta", app_context);
      
      app_context_report_info(INFO_PETSC_VEC, 1, &model->theta_k[i],"model_opeartor_01, model->theta_k_hat", app_context);    
      app_context_report_info(INFO_PETSC_VEC, 1, &model->theta_k_minus_1[i],"model_opeartor_01, model->theta_k_minus_1_hat", app_context);    
    }
    if(app_context->diagnosisMode & DIAGNOSIS_DISP_EVERYTHING){
      app_context_report_info(INFO_PETSC_VEC,1,  &model->x_k[i],"model_opeartor_01, model->x_k_hat", app_context);        
      app_context_report_info(INFO_PETSC_VEC, 1, &model->x_k_minus_1[i],"model_opeartor_01, model->x_k_minus_1_hat", app_context);    
    }
  }
  
  //update the system time stamp
  model->k=model->k+1;
  model->t=model->t+0.01;
  
  return model;
}


//-------------------------------------------------------------------------------------------------
//---test routines---
int
test01(int argc, char **args)
{
  APP_CONTEXT_TYPE* app_context = app_context_create_and_initialize(argc, args);
  MEAUSREMENT_TYPE* measurement = measurement_create_and_initialize(app_context);
  MODEL_TYPE*       model       = model_create_and_initialize      (&model_opeartor_01, app_context);
  FILTER_TYPE*      filter      = filter_create_and_initialize     (app_context, model);
  
  app_context_create_and_initialize_dependently(app_context, filter, model, measurement);
  filter_set_intial_condition(filter, model, app_context);
  
  for(int k=0; k<3; k++){
    filter_time_update       (filter, model,              app_context);
    filter_measurement_update(filter, model, measurement, app_context);
  }
  
  app_context_destroy_and_finalize(app_context);

  return 0;
}

//-------------------------------------------------------------------------------------------------
//---driver routines---
int 
main(int argc,char **args)
{
  return test01(argc,args);
}
