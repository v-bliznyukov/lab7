

#include <stdio.h>
#include <math.h>
void EulersMethod(int x0, int y0, int n, float h, double LTE[], double GTE[]){
  float x[n+1];
  for(int i=0; i<n+1; i++){
    x[i]=x0+h*i;
  }
  double yEuler[n+1];
  yEuler[0]=y0;
  for(int i=1; i<n+1; i++){
    yEuler[i]=yEuler[i-1]+h*(pow(yEuler[i-1]/x[i-1], 2)+yEuler[i-1]/x[i-1] -1);
  }
  double yExact[n+1];
  for(int i=0; i<n+1; i++){
    yExact[i]=(x[i]*(1+pow(x[i],2)/3))/(1-pow(x[i],2)/3);
  }
  GTE[0]=0;
  for(int i=1; i<n+1; i++){
    GTE[i]=yExact[i]-yEuler[i];
  }
  LTE[0]=0;
  for (int i=1; i<n+1; i++){
     LTE[i]=yExact[i]-yExact[i-1]-h*(pow(yExact[i-1]/x[i-1], 2)+yExact[i-1]/x[i-1] -1);
  }
  printf("          EULER METHOD:\n");
  printf("x    y(exact) y(approx) LTE   GTE\n");
  for (int i=0; i<n+1; i++){
    
    printf("%.3f  ", round(1000*x[i])/1000);
    printf("%.3f  ", round(1000*yExact[i])/1000);
    printf("%.3f   ",round(1000*yEuler[i])/1000);
    printf("%.3f  ", round(1000*LTE[i])/1000);
    printf("%.3f  ", round(1000*GTE[i])/1000);
    printf("\n");
  }

}

void ImprovedEuler(int x0, int y0, int n, float h, double LTE[], double GTE[]){
  float x[n+1];
  for(int i=0; i<n+1; i++){
    x[i]=x0+h*i;
  }
  float yImproved[n+1];
  
  yImproved[0]=y0;
  for(int i=1; i<n+1; i++){
    float y1=yImproved[i-1];
    float x1=x[i-1];
    float A1=(pow(y1/x1, 2)+y1/x1 -1);
    
    float x2=x1+h;
    float y2=y1+h*A1;
    float A2=(pow(y2/x2, 2)+y2/x2 -1);

    yImproved[i]=yImproved[i-1]+h*0.5*(A1+A2);
  }
  double yExact[n+1];
  for(int i=0; i<n+1; i++){
    yExact[i]=(x[i]*(1+pow(x[i],2)/3))/(1-pow(x[i],2)/3);
  }
  for(int i=0; i<n+1; i++){
    GTE[i]=yExact[i]-yImproved[i];
  }
  LTE[0]=0;
  for (int i=0; i<n; i++){
    float y1=yExact[i];
    float x1=x[i];
    float A1=(pow(y1/x1, 2)+y1/x1 -1);
    
    float x2=x1+h;
    float y2=y1+h*A1;
    float A2=(pow(y2/x2, 2)+y2/x2 -1);
    float ExactI=yExact[i]+h*0.5*(A1+A2);
    LTE[i+1]=yExact[i+1]-ExactI;
  }
  printf("      IMPROVED EULER METHOD:\n");
  printf("x    y(exact) y(approx) LTE   GTE\n");
  for (int i=0; i<n+1; i++){
    
    printf("%.3f  ", round(1000*x[i])/1000);
    printf("%.3f  ", round(1000*yExact[i])/1000);
    printf("%.3f   ",round(1000*yImproved[i])/1000);
    printf("%.3f  ", round(1000*LTE[i])/1000);
    printf("%.3f  ", round(1000*GTE[i])/1000);
    printf("\n");
  }

}

void Runge_Kutta(int x0, int y0, int n, float h, double LTE[], double GTE[]){
  float x[n+1];
  for(int i=0; i<n+1; i++){
    x[i]=x0+h*i;
  }
  double yRK[n+1];
  yRK[0]=y0;
  for(int i=1; i<n+1; i++){
    double y1=yRK[i-1];
    double x1=x[i-1];
    double A1=(pow(y1/x1, 2)+y1/x1 -1);
    double y2=y1+h*A1/2;
    double x2=x1+h/2;
    double A2=(pow(y2/x2, 2)+y2/x2 -1);
    double x3=x1+h/2;
    double y3=y1+A2*h/2;
    double A3=(pow(y3/x3, 2)+y3/x3 -1);
    double x4=x1+h;
    double y4=y1+h*A3;
    double A4=(pow(y4/x4, 2)+y4/x4 -1);;

    yRK[i]=yRK[i-1]+h*(A1+2*A2+2*A3+A4)/6;
  }
  double yExact[n+1];
  for(int i=0; i<n+1; i++){
    yExact[i]=(x[i]*(1+pow(x[i],2)/3))/(1-pow(x[i],2)/3);
  }
  for(int i=0; i<n+1; i++){
    GTE[i]=yExact[i]-yRK[i];
  }
  LTE[0]=0;
  for (int i=0; i<n; i++){
    double y1=yExact[i];
    double x1=x[i];
    double A1=(pow(y1/x1, 2)+y1/x1 -1);
    double y2=y1+h*A1/2;
    double x2=x1+h/2;
    double A2=(pow(y2/x2, 2)+y2/x2 -1);
    double x3=x1+h/2;
    double y3=y1+A2*h/2;
    double A3=(pow(y3/x3, 2)+y3/x3 -1);
    double x4=x1+h;
    double y4=y1+h*A3;
    double A4=(pow(y4/x4, 2)+y4/x4 -1);;
    double Exact;
    Exact=yExact[i]+h*(A1+2*A2+2*A3+A4)/6;
     LTE[i+1]=yExact[i+1]-Exact;
  }
  printf("        RUNGE-KUTTA METHOD:\n");
  printf("x    y(exact) y(approx) LTE   GTE\n");
  for (int i=0; i<n+1; i++){
    
    printf("%.3f  ", round(1000*x[i])/1000);
    printf("%.3f  ", round(1000*yExact[i])/1000);
    printf("%.3f   ",round(1000*yRK[i])/1000);
    printf("%.3f  ", round(1000*LTE[i])/1000);
    printf("%.3f  ", round(1000*GTE[i])/1000);
    printf("\n");
  }

}
int comparefloat(float f1, float f2){
  float precision = 0.000001;
  if (((f1-precision)<f2)&&((f1+precision)>f2)){
   return 1;
  }
  else{
    return 0;
  }
}
void compare (double LTE1[], double GTE1[], double LTE2[], double GTE2[], int steps1, int steps2,float h1, float h2, float x0, int x, int y){
  float steps;
  float x1;
  float x2;
  steps=steps2;
  x1=x0+h1;
  x2=x0+h2;
  printf("        RATE OF CHANGE TABLE:\n");
  printf("x  LTE%d/LTE%d  GTE%d/GTE%d\n",x,y,x,y);
  printf("%.3f  ", x0);
  printf("%.3f  ", 0.0);
  printf("%.3f   ",0.0);
  printf("\n");
  int j=1;
   for (int i=1; i<steps+1;i++){
    if(comparefloat(x1,x2)){

    printf("%.3f  ", round(1000*x2)/1000);
    
    printf("%.3f  ", LTE1[j]/LTE2[i]);
    printf("%.3f   ",GTE1[j]/GTE2[i]);
    printf("\n");
    x1+=h1;
    j++;
    }
    x2+=h2;
   }
}
int main(void) {

float beginning=1;
float ending=1.5;
int n;
float x0=1;
float y0=2;
//------------------------//
float h1=0.1;
n=(ending-beginning)/h1;
double LTE1[n+1];
double GTE1[n+1];
double RKLTE1[n+1];
double RKGTE1[n+1];
ImprovedEuler(x0,y0,n,h1,LTE1,GTE1);
Runge_Kutta(x0,y0,n,h1,RKLTE1,RKGTE1);
//------------------------//
float h2=0.05;
n=(ending-beginning)/h2;
double LTE2[n+1];
double GTE2[n+1];
double RKLTE2[n+1];
double RKGTE2[n+1];
ImprovedEuler(x0,y0,n,h2,LTE2,GTE2);
Runge_Kutta(x0,y0,n,h2,RKLTE2,RKGTE2);

//------------------------//
float h3=0.01;
n=(ending-beginning)/h3;
double LTE3[n+1];
double GTE3[n+1];
double RKLTE3[n+1];
double RKGTE3[n+1];
ImprovedEuler(x0,y0,n,h3,LTE3,GTE3);
Runge_Kutta(x0,y0,n,h3,RKLTE3,RKGTE3);
//------------------------//
printf("Rate of change table for Improved Euler(1&2)\n the rate of change for LTE is O(h^3), GTE is O(h^2)\n");
compare(LTE1,GTE1,LTE2,GTE2,(ending-beginning)/h1,(ending-beginning)/h2,h1,h2,beginning,1,2);
printf("Rate of change table for Improved Euler(1&3)\n the rate of change for LTE is O(h^3), GTE is O(h^2)\n");
compare(LTE1,GTE1,LTE3,GTE3,(ending-beginning)/h1,(ending-beginning)/h3,h1,h3,beginning,1,3);
//------------------------//
printf("Rate of change table for Rugne-Kutta(1&2)\n the rate of change for LTE is O(h^5), GTE is O(h^4)\n");
compare(RKLTE1,RKGTE1,RKLTE2,RKGTE2,(ending-beginning)/h1,(ending-beginning)/h2,h1,h2,beginning,1,2);
printf("Rate of change table for Rugne-Kutta(1&3)\n the rate of change for LTE is O(h^5), GTE is O(h^4)\n");

compare(RKLTE1,RKGTE1,RKLTE3,RKGTE3,(ending-beginning)/h1,(ending-beginning)/h3,h1,h3,beginning,1,3);
//-------------------//
double ELTE1[n+1];
double ELTE2[n+1];
double ELTE3[n+1];
double EGTE1[n+1];
double EGTE2[n+1];
double EGTE3[n+1];
n=(ending-beginning)/h1;
EulersMethod(x0,y0,n,h1, ELTE1, EGTE1);
n=(ending-beginning)/h2;
EulersMethod(x0,y0,n,h2, ELTE2, EGTE2);
n=(ending-beginning)/h3;
EulersMethod(x0,y0,n,h3, ELTE3, EGTE3);
printf("Rate of change table for Euler metod(1&2)\nthe rate of change for LTE is O(h^2), GTE is O(h)\n");
compare(ELTE1,EGTE1,ELTE2,EGTE2,(ending-beginning)/h1,(ending-beginning)/h2,h1,h2,beginning,1,2);
printf("Rate of change table for Euler metod(1&3)\nthe rate of change for LTE is O(h^2), GTE is O(h)\n");
compare(ELTE1,EGTE1,ELTE3,EGTE3,(ending-beginning)/h1,(ending-beginning)/h3,h1,h3,beginning,1,3);
  return 0;
}

