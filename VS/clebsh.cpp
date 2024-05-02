//#include "stdafx.h"
#include <math.h>
#include <conio.h>

double YXFCT(int M, int N);
double PHASEF(int N);
double factorial(int n);
long double long_factorial(int n);
//------------------------------------------------------------------------------
//      Clebsch-Gordon coefficient, all momenta are doubled
//        ( jx1 , mx1, jx2, mx2 | jx3, mx1+mx2 )
double Clebsh_Gordon(int JX1,int JX2,int JX3,int MX1,int MX2)
{
  int J1,J2,J3,M1,M2,ICNTR,IT,JZ1,JZ2,JZ3,JT1,JT2,JT3,JT4,JT5,J4,J5,NUMIN,NUMAX,NU;
  double VCC=0.0,PHAS,FCTOR;
		J1=JX1;
		J2=JX2; 
		J3=JX3;
		M1=MX1;
		M2=MX2;
		if(J1<J2) goto l20;
		if(J3<J2) goto l30;
		ICNTR=0;
		goto l40;
 l20: if(J3<J1) goto l30;
		ICNTR=-1;
		IT=J1;
		J1=J2;
		J2=IT;
		IT=M1;
		M1=M2;
		M2=IT;
		goto l40;
 l30: ICNTR=1;
		IT=J2;
		J2=J3;
		J3=IT;
		M2=-M1-M2;
 l40:
		JZ1=(J1+J2-J3)/2;
		if(JZ1<0) goto l150;
		JZ2=(J1+J3-J2)/2;
		if(JZ2<0) goto l150;
		JZ3=(J2+J3-J1)/2;
		if(JZ3<0) goto l150;
		if(J1-abs(M1)<0) goto l150;
		if(J2-abs(M2)<0) goto l150;
		if(J3-abs(M1+M2)<0) goto l150;
		if((J1+J2+J3)%2 != 0) goto l150;
		JT1=(J1-J3+M2)/2;
		JT2=(J2-J3-M1)/2;
		if(JT1>JT2) NUMIN=JT1; else NUMIN=JT2;
		if(NUMIN<0) NUMIN=0;
		JT3=(J1-M1)/2;
		JT4=(J2+M2)/2;
		if(JT3<JT4) NUMAX=JT3; else NUMAX=JT4;
		if(NUMAX>JZ1) NUMAX=JZ1;
		JT5=(J2-M2)/2;
		if(NUMAX<NUMIN) goto l150;
		J4=J1/2;
		J5=J3/2;
		PHAS=PHASEF(NUMIN);
		for(NU=NUMIN;NU<=NUMAX;NU++){
		VCC=VCC
      	+ PHAS*(YXFCT(JT3-NU,J4)*YXFCT(NU-JT2,J5))/
           (factorial(JT4-NU)*factorial(NU-JT1)*factorial(JZ1-NU)*factorial(NU));
		PHAS=-PHAS;
											 }
		FCTOR=YXFCT(J4,(J1+M1)/2)*YXFCT(J4,JT3)*YXFCT((J1+J2+J3)/2+1,JZ2)*
				YXFCT(J5,(J3+M1+M2)/2)*YXFCT(J5,(J3-M1-M2)/2)*factorial(JZ1)*factorial(JZ3)*
				factorial(JT4)*factorial(JT5)*fabs(J3+1.0);
		VCC=sqrt(FCTOR)*VCC;
		if(ICNTR<0) goto l120;
		if(ICNTR==0)goto l150;
		VCC=VCC*sqrt(fabs(J2+1.0)/fabs(J3+1.0))*PHASEF(JT3);
		goto l150;
l120: VCC=VCC*PHASEF(JZ1);
l150: return(VCC);
}
//------------------------------------------------------------------------------
double CG_coeff(int l1,int m1,int l2,int m2,int l,int m)// only for int l,l1,l2
{
	if(m1+m2!=m) return 0;
   return Clebsh_Gordon(2*l1,2*l2,2*l,2*m1,2*m2);
}
//------------------------------------------------------------------------------
// Clebsh-Gordon coef in case of zerro magnetic numbers:
//                     (I0J0|K0)
//------------------------------------------------------------------------------
double CG0(int A, int B, int C)
{
  int G2 = A+B+C;
  if(G2%2) return 0;
  int G = G2/2,
      A2 = 2*A, B2 = 2*B, C2 = 2*C,
      parity = 1,
      GA = G-A, GB = G-B, GC = G-C,
      GA2 = G2-A2, GB2 = G2-B2, GC2 = G2-C2;
  long double
  		cg;
  //----------------------------------------------------------------------------
  if((G-C)%2) parity = -1;
  //----------------------------------------------------------------------------
  GA = G-A; GB = G-B; GC = G-C;
  GA2 = G2-A2; GB2 = G2-B2; GC2 = G2-C2;

  cg = parity*sqrtl(C2+1.0L);
  if(A<B && A<C)
  {
     cg = cg*YXFCT(GA,G)/(long_factorial(GB)*long_factorial(GC));
     cg = cg*sqrtl(YXFCT(G2+1,GA2)*long_factorial(GB2)*long_factorial(GC2));
  }
  else if(B<A && B<C)
  {
     cg = cg*YXFCT(GB,G)/(long_factorial(GA)*long_factorial(GC));
     cg = cg*sqrtl(YXFCT(G2+1,GB2)*long_factorial(GA2)*long_factorial(GC2));
  }
  else
  {
     cg = cg*YXFCT(GC,G)/(long_factorial(GA)*long_factorial(GB));
     cg = cg*sqrtl(YXFCT(G2+1,GC2)*long_factorial(GA2)*long_factorial(GB2));
  }
  return (double)cg;
}
