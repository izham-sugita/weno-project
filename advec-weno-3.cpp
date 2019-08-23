#include<iostream>
#include<algorithm> //for swap()
#include<fstream> //for file I/O
#include<cmath> // math function
#include<vector>

using namespace std;

float rkstep(float u, float dt,
	     float alpha, float source)
{
  u = u + alpha*dt*source;
  return u;
}

float rk1(float u, float dt, float rhs)
{
  u = u + dt*rhs;
  return u;
}

//------------------
float flux(float stn[5]) //calculate the right-hand side
{
  float fr;

  float gamma1, gamma2, gamma3;
  gamma1 = 1.0/16.0;
  gamma2 = 5.0/8.0;
  gamma3 = 5.0/16.0;
    
  float p1, p2, p3;

  p1 = (3.0/8.0)*stn[0]+(-5.0/4.0)*stn[1]+(15.0/8.0)*stn[2];
  p2 = (-1.0/8.0)*stn[1]+(3.0/4.0)*stn[2]+(3.0/8.0)*stn[3];
  p3 = (3.0/8.0)*stn[2]+(3.0/4.0)*stn[3]+(-1.0/8.0)*stn[4];

  /*Calculating smoothness indicator */
  float b1 = 0.3333*(4.0*stn[0]*stn[0]
	       -19.0*stn[0]*stn[1]
	       +25.0*stn[1]*stn[1]
	       +11.0*stn[0]*stn[2]
	       -31.0*stn[1]*stn[2]
	       +10.0*stn[2]*stn[2]);

  float b2 = 0.3333*(4.0*stn[1]*stn[1]
	       -13.0*stn[1]*stn[2]
	       +13.0*stn[2]*stn[2]
	       +5.0*stn[1]*stn[3]
	       -13.0*stn[2]*stn[3]
	       +4.0*stn[3]*stn[3]);

  float b3 = 0.3333*(10.0*stn[2]*stn[2]
	       -31.0*stn[2]*stn[3]
	       +25.0*stn[3]*stn[3]
	       +11.0*stn[2]*stn[4]
	       -19.0*stn[3]*stn[4]
	       +4.0*stn[4]*stn[4]);

  
  float wtilde1 = gamma1/((1.0e-6 + b1)*(1.0e-6 + b1));
  float wtilde2 = gamma2/((1.0e-6 + b2)*(1.0e-6 + b2));
  float wtilde3 = gamma3/((1.0e-6 + b3)*(1.0e-6 + b3));
  

 float w1 = wtilde1/(wtilde1+wtilde2+wtilde3);
 float w2 = wtilde2/(wtilde1+wtilde2+wtilde3);
 float w3 = wtilde3/(wtilde1+wtilde2+wtilde3);

  fr = w1*p1 + w2*p2 + w3*p3; // flux at xi+1/2
 
  return fr;
}

//------------------


int main()
{
  vector<float> uex, u, uinit;
  int imax;
  float dx, dt;
  float c=1.0; //advection velocity

  /*Input*/
  cout<<"Enter imax "<<endl;
  cin>>imax;

  dx = 12.0/(imax-1);
  cout<<"Enter dt (dx= "<<dx<<" ) "<<endl;;
  cin>>dt;

  int itermax;
  cout<<"Enter itermax "<<endl;
  cin>>itermax;
  
  /*Initial condition*/
  uex.resize(imax);
  u.resize(imax);
  uinit.resize(imax);

  /*
  for(int i=0; i<imax; ++i){
    u[i] = 0.0;
    uex[i] = 0.0;
    if(i*dx>=1.0f && i*dx<=2.0f){
      u[i] = 1.0;
    }
  }  
  */

  for(int i=0; i<imax; ++i){
    u[i] = 1.0;
    if(i*dx>1.0f){
      u[i] = 0.0;
    }

    uex[i] = u[i];
  }
  
  for(int i=0; i<imax; ++i) uinit[i] = u[i];
  
  /*Write to file for initial condition*/
  ofstream fp;
  fp.open("init.csv");
  fp<<"x, u\n";
  for(int i=0; i<imax; ++i){
    fp<<i*dx<<",\t"<<u[i]<<"\n";
  }
  fp.close();
  
  /*Creating analytical solution*/
  int iter=0;
  do{

    for(int i=2; i<imax-1; ++i){
      uex[i] = u[i-1];  
    }
    //update
    swap(u,uex);
    iter +=1;

  }while(iter<itermax+1);

  /*Analytical solution output*/
  fp.open("analytical.csv");
  fp<<"x, u\n";
  for(int i=0; i<imax; ++i){
    fp<<i*dx<<",\t"<<u[i]<<"\n";
  }
  fp.close();

  /*starting the WENO scheme*/
  cout<<"Solution marching time for analytical solution"<<endl;
  cout<<itermax*dx<<endl;
  cout<<"Enter new itermax to match"<<endl;
  cin>>itermax;
  
  
  vector<float> u1,u2;
  u1.resize(imax);
  u2.resize(imax);
  for(int i=0; i<imax; ++i){
    u1[i] = uinit[i];
    u2[i] = uinit[i];
    u[i] = uinit[i];
  }

  float source;
  float utemp;
  float uin;
  float u_1;
  float alpha;
  iter =0;
 do{

   float stn[5];
   float fl, fr;
   float fleft, fright;
   float uleft, uright;
   
   for(int i=3; i<imax-3; ++i){
   uin = u[i];
/*------------------------------------------------------*/      
   stn[0] = c*u[i-2];
    stn[1] = c*u[i-1];
    stn[2] = c*u[i]; //interface [i,i+1]
    stn[3] = c*u[i+1];
    stn[4] = c*u[i+2];
    fleft = flux(stn);
    stn[0] = u[i-2];
    stn[1] = u[i-1];
    stn[2] = u[i]; //interface [i,i+1]
    stn[3] = u[i+1];
    stn[4] = u[i+2];
    uleft = flux(stn);
    /*--------------------------------------*/    
    stn[0] = c*u[i+3];
    stn[1] = c*u[i+2];
    stn[2] = c*u[i+1]; // interface [i,i+1]
    stn[3] = c*u[i];
    stn[4] = c*u[i-1];
    fright = flux(stn);
    stn[0] = u[i+3];
    stn[1] = u[i+2];
    stn[2] = u[i+1]; // interface [i,i+1]
    stn[3] = u[i];
    stn[4] = u[i-1];
    uright = flux(stn);
    /*Roe flux*/
    fr = 0.5*(fleft + fright) - 0.5*abs(c)*(uright-uleft);
    /*---------------------------------------------------*/
    /*for interface [i,i-1] */
      stn[0] = c*u[i-3];
    stn[1] = c*u[i-2];
    stn[2] = c*u[i-1]; //interface [i,i-1]
    stn[3] = c*u[i];
    stn[4] = c*u[i+1];
    fleft = flux(stn);
    stn[0] = u[i-3];
    stn[1] = u[i-2];
    stn[2] = u[i-1]; //interface [i,i-1]
    stn[3] = u[i];
    stn[4] = u[i+1];
    uleft = flux(stn);
    /*--------------------------------------*/    
    stn[0] = c*u[i+2];
    stn[1] = c*u[i+1];
    stn[2] = c*u[i]; // interface [i-1,i]
    stn[3] = c*u[i-1];
    stn[4] = c*u[i-2];
    fright = flux(stn);
    stn[0] = u[i+2];
    stn[1] = u[i+1];
    stn[2] = u[i]; // interface [i-1,i]
    stn[3] = u[i-1];
    stn[4] = u[i-2];
    uright = flux(stn);
    /*Roe flux*/
    fl = 0.5*(fleft + fright) - 0.5*abs(c)*(uright-uleft);
    /*---------------------------------------------------*/
/*------------------------------------------------------*/      
    
   source = -(fr-fl)/dx;

   utemp = rk1(uin, dt, source);

   u1[i] = utemp;
   
   } //first step

   for(int i=3; i<imax-3; ++i){
      /*------------------------------------------------------*/      
   stn[0] = c*u1[i-2];
    stn[1] = c*u1[i-1];
    stn[2] = c*u1[i]; //interface [i,i+1]
    stn[3] = c*u1[i+1];
    stn[4] = c*u1[i+2];
    fleft = flux(stn);
    stn[0] = u1[i-2];
    stn[1] = u1[i-1];
    stn[2] = u1[i]; //interface [i,i+1]
    stn[3] = u1[i+1];
    stn[4] = u1[i+2];
    uleft = flux(stn);
    /*--------------------------------------*/    
    stn[0] = c*u1[i+3];
    stn[1] = c*u1[i+2];
    stn[2] = c*u1[i+1]; // interface [i,i+1]
    stn[3] = c*u1[i];
    stn[4] = c*u1[i-1];
    fright = flux(stn);
    stn[0] = u1[i+3];
    stn[1] = u1[i+2];
    stn[2] = u1[i+1]; // interface [i,i+1]
    stn[3] = u1[i];
    stn[4] = u1[i-1];
    uright = flux(stn);
    /*Roe flux*/
    fr = 0.5*(fleft + fright) - 0.5*abs(c)*(uright-uleft);
    /*---------------------------------------------------*/
    /*for interface [i,i-1] */
      stn[0] = c*u1[i-3];
    stn[1] = c*u1[i-2];
    stn[2] = c*u1[i-1]; //interface [i,i-1]
    stn[3] = c*u1[i];
    stn[4] = c*u1[i+1];
    fleft = flux(stn);
    stn[0] = u1[i-3];
    stn[1] = u1[i-2];
    stn[2] = u1[i-1]; //interface [i,i-1]
    stn[3] = u1[i];
    stn[4] = u1[i+1];
    uleft = flux(stn);
    /*--------------------------------------*/    
    stn[0] = c*u1[i+2];
    stn[1] = c*u1[i+1];
    stn[2] = c*u1[i]; // interface [i-1,i]
    stn[3] = c*u1[i-1];
    stn[4] = c*u1[i-2];
    fright = flux(stn);
    stn[0] = u1[i+2];
    stn[1] = u1[i+1];
    stn[2] = u1[i]; // interface [i-1,i]
    stn[3] = u1[i-1];
    stn[4] = u1[i-2];
    uright = flux(stn);
    /*Roe flux*/
    fl = 0.5*(fleft + fright) - 0.5*abs(c)*(uright-uleft);
    /*---------------------------------------------------*/
/*------------------------------------------------------*/ 
    //utemp = rk1(uin, dt, source);
    
   uin = u[i];
   u_1 = u1[i];
   source = -(fr-fl)/dx;
   
   utemp = rk1(u_1,dt,source);

   
   u2[i] = 0.75*uin + 0.25*utemp;

   } //second step

   for(int i=3; i<imax-3; ++i){
   uin = u[i];
   u_1 = u1[i];

/*------------------------------------------------------*/      
   stn[0] = c*u2[i-2];
    stn[1] = c*u2[i-1];
    stn[2] = c*u2[i]; //interface [i,i+1]
    stn[3] = c*u2[i+1];
    stn[4] = c*u2[i+2];
    fleft = flux(stn);
    stn[0] = u2[i-2];
    stn[1] = u2[i-1];
    stn[2] = u2[i]; //interface [i,i+1]
    stn[3] = u2[i+1];
    stn[4] = u2[i+2];
    uleft = flux(stn);
    /*--------------------------------------*/    
    stn[0] = c*u2[i+3];
    stn[1] = c*u2[i+2];
    stn[2] = c*u2[i+1]; // interface [i,i+1]
    stn[3] = c*u2[i];
    stn[4] = c*u2[i-1];
    fright = flux(stn);
    stn[0] = u2[i+3];
    stn[1] = u2[i+2];
    stn[2] = u2[i+1]; // interface [i,i+1]
    stn[3] = u2[i];
    stn[4] = u2[i-1];
    uright = flux(stn);
    /*Roe flux*/
    fr = 0.5*(fleft + fright) - 0.5*abs(c)*(uright-uleft);
    /*---------------------------------------------------*/
    /*for interface [i,i-1] */
      stn[0] = c*u2[i-3];
    stn[1] = c*u2[i-2];
    stn[2] = c*u2[i-1]; //interface [i,i-1]
    stn[3] = c*u2[i];
    stn[4] = c*u2[i+1];
    fleft = flux(stn);
    stn[0] = u2[i-3];
    stn[1] = u2[i-2];
    stn[2] = u2[i-1]; //interface [i,i-1]
    stn[3] = u2[i];
    stn[4] = u2[i+1];
    uleft = flux(stn);
    /*--------------------------------------*/    
    stn[0] = c*u2[i+2];
    stn[1] = c*u2[i+1];
    stn[2] = c*u2[i]; // interface [i-1,i]
    stn[3] = c*u2[i-1];
    stn[4] = c*u2[i-2];
    fright = flux(stn);
    stn[0] = u2[i+2];
    stn[1] = u2[i+1];
    stn[2] = u2[i]; // interface [i-1,i]
    stn[3] = u2[i-1];
    stn[4] = u2[i-2];
    uright = flux(stn);
    /*Roe flux*/
    fl = 0.5*(fleft + fright) - 0.5*abs(c)*(uright-uleft);
    /*---------------------------------------------------*/
/*------------------------------------------------------*/ 
   
   uin = u[i];
   u_1 = u2[i];
   source = -(fr-fl)/dx;
   utemp = rk1(u_1,dt,source);
   
   u1[i] = 0.33333f*uin + 0.66666f*utemp;
   
   } //third step
   
   swap(u,u1);
   
   
     iter +=1;
     
 }while(iter<itermax+1);

 //WENO output
 fp.open("weno-ssprk.csv");
 fp<<"x, fx\n";
 for(int i=0; i<imax; ++i){
   fp<<i*dx<<",\t"<<u[i]<<"\n";
 }
 fp.close();
 
}
