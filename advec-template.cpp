#include<iostream>
#include<algorithm>
#include<fstream>
#include<cmath>
#include<vector>

using namespace std;

float rkstep(float u, float dt, float alpha, float rhs)
{
  u = u + alpha*dt*rhs;
  return u;
}

float flux(float stn[5]) //calculate the right-hand side
{
  float fr;

  float gamma1, gamma2, gamma3;
  gamma1 = 0.1f;
  gamma2 = 0.6f;
  gamma3 = 0.3f;
    
  float p1, p2, p3;

  p1 = (1.0/3.0)*stn[0]+(-7.0/6.0)*stn[1]+(11.0/6.0)*stn[2];
  p2 = (-1.0/6.0)*stn[1]+(5.0/6.0)*stn[2]+(1.0/3.0)*stn[3];
  p3 = (1.0/3.0)*stn[2]+(5.0/6.0)*stn[3]+(-1.0/6.0)*stn[4];

  /*Calculating smoothness indicator */
  float b1 = (13.0/12.0)*(stn[0]-2.0*stn[1]+stn[2])*(stn[0]-2.0*stn[1]+stn[2])
    +0.25*(3.0*stn[0]-4.0*stn[1]+stn[2])*(3.0*stn[0]-4.0*stn[1]+stn[2]);
  
  float b2 = (13.0/12.0)*(stn[1]-2.0*stn[2]+stn[3])*(stn[1]-2.0*stn[2]+stn[3])
    +0.25*(3.0*stn[1]-4.0*stn[2]+stn[3])*(3.0*stn[1]-4.0*stn[2]+stn[3]);
  
  float b3 = (13.0/12.0)*(stn[2]-2.0*stn[3]+stn[4])*(stn[2]-2.0*stn[3]+stn[4])
    +0.25*(3.0*stn[2]-4.0*stn[3]+stn[4])*(3.0*stn[2]-4.0*stn[3]+stn[4]);
  
  
  float wtilde1 = gamma1/((1.0e-6 + b1)+(1.0e-6 + b2)+(1.0e-6 + b3) );
 float wtilde2 = gamma2/((1.0e-6 + b1)+(1.0e-6 + b2)+(1.0e-6 + b3) );
 float wtilde3 = gamma3/((1.0e-6 + b1)+(1.0e-6 + b2)+(1.0e-6 + b3) );
  

 float w1 = wtilde1/(wtilde1+wtilde2+wtilde3);
 float w2 = wtilde2/(wtilde1+wtilde2+wtilde3);
 float w3 = wtilde3/(wtilde1+wtilde2+wtilde3);
    
  fr = w1*p1 + w2*p2 + w3*p3; // flux at xi+1/2
  
  return fr;
}

int main()
{

  ofstream fp; //file output object
  float stn[5]; // stencil storage
  int imax;
  int i;
  float pi = 4.0*atan(1.0);

  cout<<"Enter imax: "<<endl;
  cin>>imax;
  
  float dx = 2.0*pi/(imax-1);
  float dt;
  cout<<"Enter dt (dx="<<dx<<"): \n";
  cin>>dt;

  float c=1.0; //advection speed
  
  vector<float> x, u, fx; //x->grid, u->qty, fx->flux
  
  x.resize(imax);
  fx.resize(imax);
  u.resize(imax);
  
  for(i=0; i<imax; ++i){
    x[i] = i*dx;
    u[i] = 0.0;
    if(x[i]>=0.25*pi and x[i]<=0.75*pi) u[i]=1.0; //function with jump
    fx[i] = c*u[i]; //flux
  }


  for(i=1; i<imax-2; ++i){

    //Boundary cell u[0],u[1],u[2]
    if(i==1){
      //special treatment
      
    }
    else{

    }
    
  }
  
}
