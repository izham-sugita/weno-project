#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>
#include<sstream>
#include<iomanip>
#include<algorithm>

using namespace std;

/*function prototype*/
void output(int, vector<float>, vector<float> );
float rkstep(float, float, float,  float);
float cdweno2(float stn[3]);


void init(int imax, float dx,
	  vector<float> &u, vector<float> &xg)
{

  u.resize(imax);
  xg.resize(imax);
  
  for(int i=0; i<imax; ++i){
    u[i] = 0.0;
    xg[i] = i*dx;
    if(xg[i] >4.0 && xg[i]<6.0){
      u[i] = 1.0;
    }
  }
  
};

int main()
{

  vector<float> u, uinit, xg, u1, u2;
  int imax;
  float dx, dt;
  float c = 1.0;
  int i;
  int iter, itermax;
  
  cout<<"Enter imax "<<endl;
  cin>>imax;

  dx = 10.0/(imax-1);
  
  cout<<"Enter dt (dx ="<<dx<<") \n";
  cin>>dt;
  cout<<"Enter maximum iteration\n";
  cin>>itermax;
  
  uinit.resize(imax);
  u1.resize(imax);
  u2.resize(imax);

  //dx = 10.0/(imax-1);

  init(imax,dx,u,xg);

  uinit = u;  //vector container is COOL!
  u1 = u;
  u2 = u;

  /*
  i =0;
  do{
    cout<<uinit[i]<<", "<<u[i]<<endl;
    i+=1;
  }while(i<u.size());
  */

  float uin, dfdx, utemp;
  float stn[3];
  float eps = 1.0e-6;
  
  /*Main loop*/
  iter = 0;
  do{

    /*1st RK step*/
    for(i=1; i<imax-1; ++i){

      
      stn[0] = ( c*(u[i+1] - u[i])/dx );
      stn[1] = ( c*(u[i] - u[i-1])/dx  );
      
      //dfdx = -cdweno2(stn);
      
      
      
      float b1 = (u[i] - u[i-1])*(u[i] - u[i-1]);
      float b2 = (u[i+1] - u[i])*(u[i+1] - u[i]);
      float wtilde1 = 0.5/( (eps + b1)*(eps + b1) );
      float wtilde2 = 0.5/( (eps + b2)*(eps + b2) );
      float w1 = wtilde1/(wtilde1 + wtilde2);
      float w2 = wtilde2/(wtilde1 + wtilde2);

      //cout<<"w1="<<w1<<", w2="<<w2<<endl;
      
      dfdx = -(w1*stn[0]+w2*stn[1]) ;
	
      uin = u[i];
      utemp = rkstep(uin, dt, 0.5, dfdx);
      u1[i] = utemp;
    }

      /*2nd RK step*/
    for(i=1; i<imax-1; ++i){
            
      stn[0] = ( c*(u1[i+1] - u1[i])/dx );
      stn[1] = ( c*(u1[i] - u1[i-1])/dx );
      //dfdx = -cdweno2(stn);

       float b1 = (u1[i] - u1[i-1])*(u1[i] - u1[i-1]);
      float b2 = (u1[i+1] - u1[i])*(u1[i+1] - u1[i]);
      float wtilde1 = 0.5/( (eps + b1)*(eps + b1) );
      float wtilde2 = 0.5/( (eps + b2)*(eps + b2) );
      float w1 = wtilde1/(wtilde1 + wtilde2);
      float w2 = wtilde2/(wtilde1 + wtilde2);

      dfdx = -(w1*stn[0]+w2*stn[1]) ;
      uin = u[i];
      utemp = rkstep(uin, dt, 1.0, dfdx);
      u2[i] = utemp;
    }

    swap(u,u2);

    output(iter, xg, u);
    
    iter +=1;
  }while(iter < itermax);

  
  
}
