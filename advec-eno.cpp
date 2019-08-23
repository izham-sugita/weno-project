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

//------------------
float eno(float stn[5]) //calculate the right-hand side
{
  float fr;

  float select;

  float gamma1, gamma2, gamma3;
  
  gamma1 = 1.0/16.0;
  gamma2 = 5.0/8.0;
  gamma3 = 5.0/16.0;
    
  float p1, p2, p3;

  p1 = (3.0/8.0)*stn[0]+(-5.0/4.0)*stn[1]+(15.0/8.0)*stn[2];
  p2 = (-1.0/8.0)*stn[1]+(3.0/4.0)*stn[2]+(3.0/8.0)*stn[3];
  p3 = (3.0/8.0)*stn[2]+(3.0/4.0)*stn[3]+(-1.0/8.0)*stn[4];

  
  select = min( abs(p1), abs(p2) );
  select = min( select, abs(p3) );
  if( select==abs(p1) ) fr = p1;
  if( select==abs(p2) ) fr = p2;
  if( select==abs(p3) ) fr = p3;
  

  /*
  select = min( p1, p2 );
  select = min( select, p3 );
  fr = select;
  */

  
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
  for(int i=0; i<imax; ++i){
    u[i] = 1.0;
    //uex[i] = 0.0;
    
    //if(i*dx>=1.0f && i*dx<=2.0f){
    if(i*dx>=2.0f){      
      //u[i] = 1.0;
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

    for(int i=1; i<imax-1; ++i){
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
    //u1[i]= 0.0;
    //u2[i] = 0.0;

        u1[i]= uinit[i];
	u2[i]= uinit[i];
    u[i] = uinit[i];
  }

  float source;
  float utemp;
  float uin;
  float alpha;

  float stn[5];
  float fl, fr;
  float fleft, fright;
  float uleft, uright;
   
  iter =0;
 do{

   for(int i=3; i<imax-3; ++i){
   uin = u[i];
   alpha=0.5f;
/*------------------------------------------------------*/      
   stn[0] = u[i-2];
    stn[1] = u[i-1];
    stn[2] = u[i]; //interface [i,i+1]
    stn[3] = u[i+1];
    stn[4] = u[i+2];
    fr = eno(stn);
    
    /*for interface [i,i-1] */
      stn[0] = u[i-3];
    stn[1] = u[i-2];
    stn[2] = u[i-1]; //interface [i,i-1]
    stn[3] = u[i];
    stn[4] = u[i+1];
    fl = eno(stn);
/*--------------------------------------------------*/
    
   source = -c*(fr-fl)/dx;
   utemp = rkstep(uin,dt,alpha,source);
   u1[i] = utemp;
   } //first step

      for(int i=3; i<imax-3; ++i){
   uin = u[i];
   alpha=1.0f;
/*------------------------------------------------------*/      
   stn[0] = u1[i-2];
    stn[1] = u1[i-1];
    stn[2] = u1[i]; //interface [i,i+1]
    stn[3] = u1[i+1];
    stn[4] = u1[i+2];
    fr = eno(stn);
    
    /*for interface [i,i-1] */
      stn[0] = u1[i-3];
    stn[1] = u1[i-2];
    stn[2] = u1[i-1]; //interface [i,i-1]
    stn[3] = u1[i];
    stn[4] = u1[i+1];
    fl = eno(stn);
/*--------------------------------------------------*/
    
   source = -c*(fr-fl)/dx;
   utemp = rkstep(uin,dt,alpha,source);
   u2[i] = utemp;
   } //second step

   swap(u,u2);
   
     iter +=1;
     
 }while(iter<itermax+1);

 //WENO output
 fp.open("ENO.csv");
 fp<<"x, fx\n";
 for(int i=0; i<imax; ++i){
   fp<<i*dx<<",\t"<<u[i]<<"\n";
 }
 fp.close();
 
}
