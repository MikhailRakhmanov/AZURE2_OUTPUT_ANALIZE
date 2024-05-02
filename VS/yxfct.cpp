//#include "stdafx.h"
#include <stdio.h>
#include <math.h>
//------------------------------------------------------------------------------
//              Subroutine computes N!/M! by short method
//
double YXFCT(int M,int N)
{
  long double yx,FCTOR;
  int NUMAX,ICTRL,NU;

  NUMAX=M-N;
  if(NUMAX<0)
  {
    ICTRL=1;
    NUMAX=-NUMAX;
    FCTOR=M;
  }
  else if(NUMAX>0)
  {
    ICTRL=0;
    FCTOR=N;
  }
  else return 1.;

  yx=1.0;
  for(NU=1;NU<=NUMAX;NU++)
  {
    FCTOR=FCTOR+1.0;
    yx=yx*FCTOR;
  }

  if(ICTRL==0) yx=1.0/yx;
  return (double) yx;
}

