extern "C" {
  void uniaxialextensionexample_initializeC();
  void uniaxialextensionexample_finalizeC();
  void uniaxialextensionexample_solveC(double* theta);
}

int 
main(int argc,char **args)
{
  double theta[2]={1.0, 1.0};
  uniaxialextensionexample_initializeC();
  uniaxialextensionexample_solveC(theta);
  uniaxialextensionexample_finalizeC();
  
  return 0;
}
