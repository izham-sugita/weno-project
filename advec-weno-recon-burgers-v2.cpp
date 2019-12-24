#include<iostream>
#include<algorithm> //for swap()
#include<fstream> //for file I/O
#include<cmath> // math function
#include<vector>
#include<string>

#include<sstream>
#include<iomanip>

#define pi 4.0*atan(1.0)

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

template<typename T>int sign(T val){
  return (T(0) < val) - (val < T(0));
}

/*File output*/
void output(int steps, vector<float> x, vector<float> fx)
{
  ofstream fp;
  stringstream buf;
  string filenumber;

  buf<<setfill('0');
  filenumber = to_string(steps);
  buf<<setw(5)<<filenumber;

  string filename="./data/f"+buf.str()+".csv";
  fp.open(filename, ios::out);
  fp<<"x, fx\n";
  for(int i=0; i<x.size(); ++i){
    fp<<x[i]<<", "
      <<fx[i]<<"\n";
  }
  fp.close();
  buf.str(string()); //clear buffer
}

//flux2 is a WENO reconstruction procedure.
float weno_recon(float stn[5]) //calculate the right-hand side
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
  
    
 float wtilde1 = gamma1/( (1.0e-6 + b1)*(1.0e-6 + b1) );
 float wtilde2 = gamma2/( (1.0e-6 + b2)*(1.0e-6 + b2) );
 float wtilde3 = gamma3/( (1.0e-6 + b3)*(1.0e-6 + b3) );
  

 float w1 = wtilde1/(wtilde1+wtilde2+wtilde3);
 float w2 = wtilde2/(wtilde1+wtilde2+wtilde3);
 float w3 = wtilde3/(wtilde1+wtilde2+wtilde3);
    
  fr = w1*p1 + w2*p2 + w3*p3; // flux at xi+1/2
  
  return fr;
}

float cdweno2(float stn[3]){
  float dfdx;
  float gamma1 = 0.5;
  float gamma2 =  0.5; 
  

  float p1 = stn[0];
  float p2 = stn[1];
  

  float eps = 1.0e-6;

  float b1 = stn[0]*stn[0];
  float b2 = stn[1]*stn[1];
  

  float wtilde1 = gamma1/( (eps + b1)*(eps + b1) );
  float wtilde2 = gamma2/( (eps + b2)*(eps + b2) );
  

  float w1 = wtilde1/(wtilde1 + wtilde2);
  float w2 = wtilde2/(wtilde1 + wtilde2);
    
  //reconstructed dfdx
  dfdx = w1*p1 + w2*p2;
    
  return dfdx;
}



int main()
{
  vector<float> uex, u, uinit;
  int imax;
  float dx, dt;
  float c = 1.0; //advection velocity

  cout<<"Speed sign "<<sign(c)<<endl;
  
  /*Input*/
  cout<<"Enter imax "<<endl;
  cin>>imax;

  dx = (2.0*pi)/(imax-1);
  cout<<"Enter dt (dx= "<<dx<<" ) "<<endl;;
  cin>>dt;

  int itermax;
  cout<<"Enter itermax "<<endl;
  cin>>itermax;
  
  /*Initial condition*/
  uex.resize(imax);
  u.resize(imax);
  uinit.resize(imax);
  
  vector<float> ccyclic;
  vector<float> xg;
  ccyclic.resize(imax);
  xg.resize(imax);
  
  for(int i=0; i<imax; ++i){
       xg[i] = i*dx;

       u[i] = sin(xg[i]);
    
  }  

  uinit = u; // i lov vector!
    
  
  /*Write to file for initial condition*/
  ofstream fp;
  fp.open("init.csv");
  fp<<"x, u\n";
  for(int i=0; i<imax; ++i){
    fp<<i*dx<<",\t"<<u[i]<<"\n";
  }
  fp.close();
  
  /*Creating analytical solution*/
  /*
  int iter=0;
  do{

    for(int i=2; i<imax-1; ++i){
      uex[i] = u[i-1];  
    }
    //update
    swap(u,uex);
    iter +=1;

  }while(iter<itermax+1);
  */
  
  /*Analytical solution output*/
  /*
  fp.open("analytical.csv");
  fp<<"x, u\n";
  for(int i=0; i<imax; ++i){
    fp<<i*dx<<",\t"<<u[i]<<"\n";
  }
  fp.close();
  */
  
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
    //u[i] = uinit[i];
  }

  float source;
  float utemp;
  float uin;
  float u_1;
  float alpha;
  int iter =0;

  int shift;
  if(sign(c) > 0) shift = 0;
  if(sign(c) < 0) shift = 1;
    
   float stn[5];
   float dfdx;
   
   do{
   
     //for(int i=3; i<imax-3; ++i){
   for(int i=2; i<imax-2; ++i){
     /*Handling boundary condition*/
     if(i==1){
       shift = 1;
        // WENO4
       stn[abs(4*shift-0)] = c*( u[i-2+shift] - u[i-3+shift] )/dx;
   stn[abs(4*shift-1)] = c*( u[i-1+shift] - u[i-2+shift] )/dx;
   stn[abs(4*shift-2)] = c*( u[i+shift] - u[i-1+shift] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = c*( u[i+1+shift] - u[i+shift] )/dx;
   stn[abs(4*shift-4)] = c*( u[i+2+shift] - u[i+1+shift] )/dx;
   dfdx = -u[i]*weno_recon(stn); //reconstructed derivative
       
     }else if(i == imax-2){
       shift = 0;
              // WENO4
       stn[abs(4*shift-0)] = c*( u[i-2+shift] - u[i-3+shift] )/dx;
   stn[abs(4*shift-1)] = c*( u[i-1+shift] - u[i-2+shift] )/dx;
   stn[abs(4*shift-2)] = c*( u[i+shift] - u[i-1+shift] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = c*( u[i+1+shift] - u[i+shift] )/dx;
   stn[abs(4*shift-4)] = c*( u[i+2+shift] - u[i+1+shift] )/dx;
   dfdx = -u[i]*weno_recon(stn); //reconstructed derivative
       
     }else{
     
     if(sign(u[i]) >= 0 ) shift =0;
     if(sign(u[i]) < 0 ) shift = 1;

   // WENO4
   stn[abs(4*shift-0)] = c*( u[i-2+shift] - u[i-3+shift] )/dx;
   stn[abs(4*shift-1)] = c*( u[i-1+shift] - u[i-2+shift] )/dx;
   stn[abs(4*shift-2)] = c*( u[i+shift] - u[i-1+shift] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = c*( u[i+1+shift] - u[i+shift] )/dx;
   stn[abs(4*shift-4)] = c*( u[i+2+shift] - u[i+1+shift] )/dx;
   dfdx = -u[i]*weno_recon(stn); //reconstructed derivative

     }

     uin = u[i];
     utemp = rk1(uin, dt, dfdx);
     u1[i] = utemp;
   
   } //first step

   for(int i=2; i<imax-1; ++i){
 /*Handling boundary condition*/
     if(i==1){
       shift = 1;
        // WENO4
       stn[abs(4*shift-0)] = c*( u1[i-2+shift] - u1[i-3+shift] )/dx;
   stn[abs(4*shift-1)] = c*( u1[i-1+shift] - u1[i-2+shift] )/dx;
   stn[abs(4*shift-2)] = c*( u1[i+shift] - u1[i-1+shift] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = c*( u1[i+1+shift] - u1[i+shift] )/dx;
   stn[abs(4*shift-4)] = c*( u1[i+2+shift] - u1[i+1+shift] )/dx;
   dfdx = -u[i]*weno_recon(stn); //reconstructed derivative
       
     }else if(i == imax-2){
       shift = 0;
              // WENO4
       stn[abs(4*shift-0)] = c*( u1[i-2+shift] - u1[i-3+shift] )/dx;
   stn[abs(4*shift-1)] = c*( u1[i-1+shift] - u1[i-2+shift] )/dx;
   stn[abs(4*shift-2)] = c*( u1[i+shift] - u1[i-1+shift] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = c*( u1[i+1+shift] - u1[i+shift] )/dx;
   stn[abs(4*shift-4)] = c*( u1[i+2+shift] - u1[i+1+shift] )/dx;
   dfdx = -u[i]*weno_recon(stn); //reconstructed derivative
       
     }else{
     
      if(sign(u[i]) >= 0 ) shift =0;
     if(sign(u[i]) < 0 ) shift = 1;

     stn[abs(4*shift-0)] = c*( u1[i-2+shift] - u1[i-3+shift] )/dx;
     stn[abs(4*shift-1)] = c*( u1[i-1+shift] - u1[i-2+shift] )/dx;
     stn[abs(4*shift-2)] = c*( u1[i+shift] - u1[i-1+shift] )/dx; 
     stn[abs(4*shift-3)] = c*( u1[i+1+shift] - u1[i+shift] )/dx;
     stn[abs(4*shift-4)] = c*( u1[i+2+shift] - u1[i+1+shift] )/dx;
     dfdx = -u1[i]*weno_recon(stn); //reconstructed derivative
     }
     
   uin = u[i];
   u_1 = u1[i];
   utemp = rk1(u_1, dt, dfdx);
   u2[i] = 0.75*uin + 0.25*utemp;

   } //second step

   for(int i=2; i<imax-1; ++i){
  /*Handling boundary condition*/
     if(i==1){
       shift = 1;
        // WENO4
       stn[abs(4*shift-0)] = c*( u2[i-2+shift] - u2[i-3+shift] )/dx;
   stn[abs(4*shift-1)] = c*( u2[i-1+shift] - u2[i-2+shift] )/dx;
   stn[abs(4*shift-2)] = c*( u2[i+shift] - u2[i-1+shift] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = c*( u2[i+1+shift] - u2[i+shift] )/dx;
   stn[abs(4*shift-4)] = c*( u2[i+2+shift] - u2[i+1+shift] )/dx;
   dfdx = -u[i]*weno_recon(stn); //reconstructed derivative
       
     }else if(i == imax-2){
       shift = 0;
              // WENO4
       stn[abs(4*shift-0)] = c*( u2[i-2+shift] - u2[i-3+shift] )/dx;
   stn[abs(4*shift-1)] = c*( u2[i-1+shift] - u2[i-2+shift] )/dx;
   stn[abs(4*shift-2)] = c*( u2[i+shift] - u2[i-1+shift] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = c*( u2[i+1+shift] - u2[i+shift] )/dx;
   stn[abs(4*shift-4)] = c*( u2[i+2+shift] - u2[i+1+shift] )/dx;
   dfdx = -u[i]*weno_recon(stn); //reconstructed derivative
       
     }else{

      if(sign(u[i]) >= 0 ) shift =0;
     if(sign(u[i]) < 0 ) shift = 1;

     stn[abs(4*shift-0)] = c*( u2[i-2+shift] - u2[i-3+shift] )/dx;
     stn[abs(4*shift-1)] = c*( u2[i-1+shift] - u2[i-2+shift] )/dx;
     stn[abs(4*shift-2)] = c*( u2[i+shift] - u2[i-1+shift] )/dx; 
     stn[abs(4*shift-3)] = c*( u2[i+1+shift] - u2[i+shift] )/dx;
     stn[abs(4*shift-4)] = c*( u2[i+2+shift] - u2[i+1+shift] )/dx;
     dfdx = -u2[i]*weno_recon(stn); //reconstructed derivative
     }
     
   uin = u[i];
   u_1 = u2[i];

   utemp = rk1(u_1, dt, dfdx);
   
   u1[i] = 0.33333f*uin + 0.66666f*utemp;
   
   } //third step
   
   swap(u,u1);
   
   output(iter, xg, u);
 
   iter +=1;
     
 }while(iter<itermax+1);

 //WENO output
 fp.open("cd-weno-ssprk-recon.csv");
 fp<<"x, fx\n";
 for(int i=0; i<imax; ++i){
   fp<<i*dx<<",\t"<<u[i]<<"\n";
 }
 fp.close();
 
}
