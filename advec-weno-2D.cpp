#include "basic.h"

using namespace std; // basic.h provides basic C++ lib. for my project

//WENO5 reconstruction procedure.
float weno_recon(float stn[5]);
/*sign function*/
int direction(float val);

float rk1(float u, float dt, float rhs);

void output2D(int steps, int imax, int jmax, 
              vector< vector<float> > x, 
              vector< vector<float> > y, 
              vector< vector<float> > z, 
              vector< vector<float> > fx);


int main()
{

  int imax =201, jmax = 201;
  int i, j;

  float dx = 1.0/(imax-1);
  float dy = 1.0/(jmax-1);
  float dt = 0.0001;
  int itermax = 1000;
  
  cout<<"dx="<<dx<<endl;
  cout<<"Enter dt"<<endl;
  cin>>dt;
  cout<<"Enter maximum iteration"<<endl;
  cin>>itermax;
  
  float xc = 0.5;
  float yc = 0.5;

  float omega = 2.0*pi; // rad/secs
  
  vector< vector<float> > u, v, f, f1, f2, fn;
  vector< vector<float> > xg, yg;
  
  u.resize(imax);
  v.resize(imax);
  f.resize(imax);
  f1.resize(imax);
  f2.resize(imax);
  fn.resize(imax);

  xg.resize(imax);
  yg.resize(imax);

  for(i =0; i<imax; ++i){
    u[i].resize(jmax);
    v[i].resize(jmax);
    f[i].resize(jmax);
    f1[i].resize(jmax);
    f2[i].resize(jmax);
    fn[i].resize(jmax);

    xg[i].resize(jmax);
    yg[i].resize(jmax);
  }

  for(i=0; i<imax; ++i){
    for(j=0; j<jmax; ++j){
      xg[i][j] = i*dx;
      yg[i][j] = j*dy;

      u[i][j] = -(yg[i][j] - yc)*omega;
      v[i][j] = (xg[i][j] - xc)*omega;

      if( xg[i][j] >=0.2 && xg[i][j] <=0.4){
	if(yg[i][j] >=0.2 && yg[i][j]<=0.4){
	  f[i][j] = 1.0;
	}
      }else{
	f[i][j] = 0.0;
      }
    }
  }

  /*Initiate value*/
  f1 = f;
  f2 = f;
  
  int iter = 0;
  
  float dfdx, dfdy;
  int shift;
  float stn[5];
  float source;
  float fin, ftemp;

  int istart=2;
  int iend = imax-2;
  int jstart=2;
  int jend = jmax-2;
  int steps=10;
  
  //cout<<direction(u[10][10])<<endl;
  
  do{

    /*First block, first step RK*/
    {
      //start of i-j loop
      for(i=istart; i<iend; ++i){
	for(j=jstart; j<jend; ++j){

	  //dfdx reconstruction
	  /*Handling boundary condition*/
     if(i==istart){
       shift = 1;
        // WENO4
       stn[abs(4*shift-0)] = ( f[i-2+shift][j] - f[i-3+shift][j] )/dx;
   stn[abs(4*shift-1)] = ( f[i-1+shift][j] - f[i-2+shift][j] )/dx;
   stn[abs(4*shift-2)] = ( f[i+shift][j] - f[i-1+shift][j] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = ( f[i+1+shift][j] - f[i+shift][j] )/dx;
   stn[abs(4*shift-4)] = ( f[i+2+shift][j] - f[i+1+shift][j] )/dx;
   dfdx = -u[i][j]*weno_recon(stn); //reconstructed derivative
       
     }else if(i == iend-1){
       shift = 0;
       // WENO4
       stn[abs(4*shift-0)] = ( f[i-2+shift][j] - f[i-3+shift][j] )/dx;
   stn[abs(4*shift-1)] = ( f[i-1+shift][j] - f[i-2+shift][j] )/dx;
   stn[abs(4*shift-2)] = ( f[i+shift][j] - f[i-1+shift][j] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = ( f[i+1+shift][j] - f[i+shift][j] )/dx;
   stn[abs(4*shift-4)] = ( f[i+2+shift][j] - f[i+1+shift][j] )/dx;
   dfdx = -u[i][j]*weno_recon(stn); //reconstructed derivative
       
     }else{
     
     if(direction(u[i][j] ) >= 0 ) shift =0;
     if(direction(u[i][j] ) < 0    ) shift = 1;
     
     stn[abs(4*shift-0)] = ( f[i-2+shift][j] - f[i-3+shift][j] )/dx;
     stn[abs(4*shift-1)] = ( f[i-1+shift][j] - f[i-2+shift][j] )/dx;
     stn[abs(4*shift-2)] = ( f[i+shift][j] - f[i-1+shift][j] )/dx; // space I_[i] 
     stn[abs(4*shift-3)] = ( f[i+1+shift][j] - f[i+shift][j] )/dx;
     stn[abs(4*shift-4)] = ( f[i+2+shift][j] - f[i+1+shift][j] )/dx;
     dfdx = -u[i][j]*weno_recon(stn); //reconstructed derivative
     }	  

     /*dfdy reconstruction*/
         if(j==jstart){
       shift = 1;
        // WENO4
       stn[abs(4*shift-0)] = ( f[i][j-2+shift] - f[i][j-3+shift] )/dy;
   stn[abs(4*shift-1)] = ( f[i][j-1+shift] - f[i][j-2+shift] )/dy;
   stn[abs(4*shift-2)] = ( f[i][j+shift] - f[i][j-1+shift] )/dy; // space I_[i] 
   stn[abs(4*shift-3)] = ( f[i][j+1+shift] - f[i][j+shift] )/dy;
   stn[abs(4*shift-4)] = ( f[i][j+2+shift] - f[i][j+1+shift] )/dy;
   dfdy = -v[i][j]*weno_recon(stn); //reconstructed derivative
       
     }else if(j == jend-1){
	   shift = 0;
     // WENO4
       stn[abs(4*shift-0)] = ( f[i][j-2+shift] - f[i][j-3+shift] )/dy;
   stn[abs(4*shift-1)] = ( f[i][j-1+shift] - f[i][j-2+shift] )/dy;
   stn[abs(4*shift-2)] = ( f[i][j+shift] - f[i][j-1+shift] )/dy; // space I_[i] 
   stn[abs(4*shift-3)] = ( f[i][j+1+shift] - f[i][j+shift] )/dy;
   stn[abs(4*shift-4)] = ( f[i][j+2+shift] - f[i][j+1+shift] )/dy;
   dfdy = -v[i][j]*weno_recon(stn); //reconstructed derivative
       
              
     }else{
     
     if(direction(v[i][j] ) >= 0 ) shift =0;
     if(direction(v[i][j] ) < 0    ) shift = 1;

            stn[abs(4*shift-0)] = ( f[i][j-2+shift] - f[i][j-3+shift] )/dy;
   stn[abs(4*shift-1)] = ( f[i][j-1+shift] - f[i][j-2+shift] )/dy;
   stn[abs(4*shift-2)] = ( f[i][j+shift] - f[i][j-1+shift] )/dy; // space I_[i] 
   stn[abs(4*shift-3)] = ( f[i][j+1+shift] - f[i][j+shift] )/dy;
   stn[abs(4*shift-4)] = ( f[i][j+2+shift] - f[i][j+1+shift] )/dy;
   dfdy = -v[i][j]*weno_recon(stn); //reconstructed derivative
          
     }	  

	 /*Time integration*/
	 source = dfdx + dfdy; //minus is included in reconstruction step
	 fin = f[i][j];
	 ftemp = rk1(fin, dt, source);
	 f1[i][j] = ftemp;
	 	 
	}
      } //end of i-j loop

    }

    /*Second block, second step RK*/
    {

           //start of i-j loop
      for(i=istart; i<iend; ++i){
	for(j=jstart; j<jend; ++j){

	  //dfdx reconstruction
	  /*Handling boundary condition*/
     if(i==istart){
       shift = 1;
        // WENO4
       stn[abs(4*shift-0)] = ( f1[i-2+shift][j] - f1[i-3+shift][j] )/dx;
   stn[abs(4*shift-1)] = ( f1[i-1+shift][j] - f1[i-2+shift][j] )/dx;
   stn[abs(4*shift-2)] = ( f1[i+shift][j] - f1[i-1+shift][j] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = ( f1[i+1+shift][j] - f1[i+shift][j] )/dx;
   stn[abs(4*shift-4)] = ( f1[i+2+shift][j] - f1[i+1+shift][j] )/dx;
   dfdx = -u[i][j]*weno_recon(stn); //reconstructed derivative
       
     }else if(i == iend-1){
       shift = 0;
        // WENO4
       stn[abs(4*shift-0)] = ( f1[i-2+shift][j] - f1[i-3+shift][j] )/dx;
   stn[abs(4*shift-1)] = ( f1[i-1+shift][j] - f1[i-2+shift][j] )/dx;
   stn[abs(4*shift-2)] = ( f1[i+shift][j] - f1[i-1+shift][j] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = ( f1[i+1+shift][j] - f1[i+shift][j] )/dx;
   stn[abs(4*shift-4)] = ( f1[i+2+shift][j] - f1[i+1+shift][j] )/dx;
   dfdx = -u[i][j]*weno_recon(stn); //reconstructed derivative
       
     }else{
     
     if(direction(u[i][j] ) >= 0 ) shift =0;
     if(direction(u[i][j] ) < 0    ) shift = 1;
     
     stn[abs(4*shift-0)] = ( f1[i-2+shift][j] - f1[i-3+shift][j] )/dx;
     stn[abs(4*shift-1)] = ( f1[i-1+shift][j] - f1[i-2+shift][j] )/dx;
     stn[abs(4*shift-2)] = ( f1[i+shift][j] - f1[i-1+shift][j] )/dx; // space I_[i] 
     stn[abs(4*shift-3)] = ( f1[i+1+shift][j] - f1[i+shift][j] )/dx;
     stn[abs(4*shift-4)] = ( f1[i+2+shift][j] - f1[i+1+shift][j] )/dx;
     dfdx = -u[i][j]*weno_recon(stn); //reconstructed derivative
     }	  

     /*dfdy reconstruction*/
         if(j==jstart){
       shift = 1;
        // WENO4
       stn[abs(4*shift-0)] = ( f1[i][j-2+shift] - f1[i][j-3+shift] )/dy;
   stn[abs(4*shift-1)] = ( f1[i][j-1+shift] - f1[i][j-2+shift] )/dy;
   stn[abs(4*shift-2)] = ( f1[i][j+shift] - f1[i][j-1+shift] )/dy; // space I_[i] 
   stn[abs(4*shift-3)] = ( f1[i][j+1+shift] - f1[i][j+shift] )/dy;
   stn[abs(4*shift-4)] = ( f1[i][j+2+shift] - f1[i][j+1+shift] )/dy;
   dfdy = -v[i][j]*weno_recon(stn); //reconstructed derivative
       
     }else if(j == jend-1){
	   shift = 0;
     // WENO4
       stn[abs(4*shift-0)] = ( f1[i][j-2+shift] - f1[i][j-3+shift] )/dy;
   stn[abs(4*shift-1)] = ( f1[i][j-1+shift] - f1[i][j-2+shift] )/dy;
   stn[abs(4*shift-2)] = ( f1[i][j+shift] - f1[i][j-1+shift] )/dy; // space I_[i] 
   stn[abs(4*shift-3)] = ( f1[i][j+1+shift] - f1[i][j+shift] )/dy;
   stn[abs(4*shift-4)] = ( f1[i][j+2+shift] - f1[i][j+1+shift] )/dy;
   dfdy = -v[i][j]*weno_recon(stn); //reconstructed derivative
       
              
     }else{
     
     if(direction(v[i][j] ) >= 0 ) shift =0;
     if(direction(v[i][j] ) < 0    ) shift = 1;

            stn[abs(4*shift-0)] = ( f1[i][j-2+shift] - f1[i][j-3+shift] )/dy;
   stn[abs(4*shift-1)] = ( f1[i][j-1+shift] - f1[i][j-2+shift] )/dy;
   stn[abs(4*shift-2)] = ( f1[i][j+shift] - f1[i][j-1+shift] )/dy; // space I_[i] 
   stn[abs(4*shift-3)] = ( f1[i][j+1+shift] - f1[i][j+shift] )/dy;
   stn[abs(4*shift-4)] = ( f1[i][j+2+shift] - f1[i][j+1+shift] )/dy;
   dfdy = -v[i][j]*weno_recon(stn); //reconstructed derivative
          
     }	  

	 /*Time integration*/
	 source = dfdx + dfdy; //minus is included in reconstruction step
	 fin = f1[i][j];
	 ftemp = rk1(fin, dt, source);

	 f2[i][j] = 0.75*f[i][j] + 0.25*ftemp;
	 
	}
      } //end of i-j loop
      
    }

    /*Third block(final), third step RK*/
    {
           //start of i-j loop
      for(i=istart; i<iend; ++i){
	for(j=jstart; j<jend; ++j){

	  //dfdx reconstruction
	  /*Handling boundary condition*/
     if(i==istart){
       shift = 1;
        // WENO4
       stn[abs(4*shift-0)] = ( f2[i-2+shift][j] - f2[i-3+shift][j] )/dx;
   stn[abs(4*shift-1)] = ( f2[i-1+shift][j] - f2[i-2+shift][j] )/dx;
   stn[abs(4*shift-2)] = ( f2[i+shift][j] - f2[i-1+shift][j] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = ( f2[i+1+shift][j] - f2[i+shift][j] )/dx;
   stn[abs(4*shift-4)] = ( f2[i+2+shift][j] - f2[i+1+shift][j] )/dx;
   dfdx = -u[i][j]*weno_recon(stn); //reconstructed derivative
       
     }else if(i == iend-1){
       shift = 0;
        // WENO4
       stn[abs(4*shift-0)] = ( f2[i-2+shift][j] - f2[i-3+shift][j] )/dx;
   stn[abs(4*shift-1)] = ( f2[i-1+shift][j] - f2[i-2+shift][j] )/dx;
   stn[abs(4*shift-2)] = ( f2[i+shift][j] - f2[i-1+shift][j] )/dx; // space I_[i] 
   stn[abs(4*shift-3)] = ( f2[i+1+shift][j] - f2[i+shift][j] )/dx;
   stn[abs(4*shift-4)] = ( f2[i+2+shift][j] - f2[i+1+shift][j] )/dx;
   dfdx = -u[i][j]*weno_recon(stn); //reconstructed derivative
       
     }else{
     
     if(direction(u[i][j] ) >= 0 ) shift =0;
     if(direction(u[i][j] ) < 0    ) shift = 1;
     
     stn[abs(4*shift-0)] = ( f2[i-2+shift][j] - f2[i-3+shift][j] )/dx;
     stn[abs(4*shift-1)] = ( f2[i-1+shift][j] - f2[i-2+shift][j] )/dx;
     stn[abs(4*shift-2)] = ( f2[i+shift][j] - f2[i-1+shift][j] )/dx; // space I_[i] 
     stn[abs(4*shift-3)] = ( f2[i+1+shift][j] - f2[i+shift][j] )/dx;
     stn[abs(4*shift-4)] = ( f2[i+2+shift][j] - f2[i+1+shift][j] )/dx;
     dfdx = -u[i][j]*weno_recon(stn); //reconstructed derivative
     }	  

     /*dfdy reconstruction*/
         if(j==jstart){
       shift = 1;
        // WENO4
       stn[abs(4*shift-0)] = ( f2[i][j-2+shift] - f2[i][j-3+shift] )/dy;
   stn[abs(4*shift-1)] = ( f2[i][j-1+shift] - f2[i][j-2+shift] )/dy;
   stn[abs(4*shift-2)] = ( f2[i][j+shift] - f2[i][j-1+shift] )/dy; // space I_[i] 
   stn[abs(4*shift-3)] = ( f2[i][j+1+shift] - f2[i][j+shift] )/dy;
   stn[abs(4*shift-4)] = ( f2[i][j+2+shift] - f2[i][j+1+shift] )/dy;
   dfdy = -v[i][j]*weno_recon(stn); //reconstructed derivative
       
     }else if(j == jend-1){
	  shift = 0;
     // WENO4
       stn[abs(4*shift-0)] = ( f2[i][j-2+shift] - f2[i][j-3+shift] )/dy;
   stn[abs(4*shift-1)] = ( f2[i][j-1+shift] - f2[i][j-2+shift] )/dy;
   stn[abs(4*shift-2)] = ( f2[i][j+shift] - f2[i][j-1+shift] )/dy; // space I_[i] 
   stn[abs(4*shift-3)] = ( f2[i][j+1+shift] - f2[i][j+shift] )/dy;
   stn[abs(4*shift-4)] = ( f2[i][j+2+shift] - f2[i][j+1+shift] )/dy;
   dfdy = -v[i][j]*weno_recon(stn); //reconstructed derivative
       
              
     }else{
     
     if(direction(v[i][j] ) >= 0 ) shift =0;
     if(direction(v[i][j] ) < 0    ) shift = 1;

            stn[abs(4*shift-0)] = ( f2[i][j-2+shift] - f2[i][j-3+shift] )/dy;
   stn[abs(4*shift-1)] = ( f2[i][j-1+shift] - f2[i][j-2+shift] )/dy;
   stn[abs(4*shift-2)] = ( f2[i][j+shift] - f2[i][j-1+shift] )/dy; // space I_[i] 
   stn[abs(4*shift-3)] = ( f2[i][j+1+shift] - f2[i][j+shift] )/dy;
   stn[abs(4*shift-4)] = ( f2[i][j+2+shift] - f2[i][j+1+shift] )/dy;
   dfdy = -v[i][j]*weno_recon(stn); //reconstructed derivative
          
     }	  

	 /*Time integration*/
	 source = dfdx + dfdy; //minus is included in reconstruction step
	 fin = f2[i][j];
	 ftemp = rk1(fin, dt, source);

	 f1[i][j] = 0.3333333*f[i][j] + 0.6666666*ftemp;
	 
	}
      } //end of i-j loop
      
    }

    /*Update*/
    swap(f,f1);

    if(iter%steps == 0){
    output2D(iter, imax, jmax, 
              xg, 
              yg, 
              f, 
              f);
    }
    
    iter +=1;
  }while(iter<itermax);
  
  
  /*file output*/
  /*
  ofstream fp;
  fp.open("init-2D.csv", ios::out);
  fp<<"x,y,z,u,v,f\n";
  for(i=0; i<imax; ++i){
    for(j=0; j<jmax; ++j){
      fp<<xg[i][j]<<","
	<<yg[i][j]<<","
	<<0.0<<","
	<<u[i][j]<<","
	<<v[i][j]<<","
	<<f[i][j]<<"\n";
    }
  }
  fp.close();
  */
  
}
