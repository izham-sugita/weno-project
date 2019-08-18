#include<iostream>
#include<cmath>
#include<vector>
#include<fstream>

using namespace std;


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

  
 float wtilde1 = gamma1/(1.0e-6 + b1)*(1.0e-6 + b1);
 float wtilde2 = gamma2/(1.0e-6 + b2)*(1.0e-6 + b2);
 float wtilde3 = gamma3/(1.0e-6 + b3)*(1.0e-6 + b3);
  

 float w1 = wtilde1/(wtilde1+wtilde2+wtilde3);
 float w2 = wtilde2/(wtilde1+wtilde2+wtilde3);
 float w3 = wtilde3/(wtilde1+wtilde2+wtilde3);
    

  fr = w1*p1 + w2*p2 + w3*p3; // flux at xi+1/2
  
  //fr = gamma1*p1 + gamma2*p2 + gamma3*p3; // flux at xi+1/2
  //fr = p1; // flux at xi+1/2
  //fr = p2; // flux at xi+1/2 //p2 <-problem?
  //fr = p3; // flux at xi+1/2 //p3 <-problem?
 
  return fr;
}

int main()
{

  ofstream fp; //for general use
  
  float stn[5];
  int imax;
  cout<<"Enter imax: "<<endl;
  cin>>imax;  
  float pi = 4.0*atan(1.0);
  float dx = pi/(imax-1);

  float x[imax];
  float fx[imax];

  //cout<<"x, fx"<<endl;
  for(int i=0; i<imax; ++i)
    {
      x[i] = i*dx;
      
      //fx[i] = sin(x[i]); //smooth function

      fx[i] = 0.0;
      if(x[i]>=0.25*pi and x[i]<=0.75*pi) fx[i]=1.0;
    }

  vector<float> xp, fint, ref;
    
  //  cout<<"x, fx, sin(x), err"<<endl;
  //interpolated values
  //on the plus
  int i=2;
  do{

    /*
    stn[0] = fx[i];
    stn[1] = fx[i+1];
    stn[2] = fx[i+2];
    stn[3] = fx[i+3];
    stn[4] = fx[i+4];
float xg = x[i+2] + 0.5*dx;    

*/
    
    stn[0] = fx[i-2];
    stn[1] = fx[i-1];
    stn[2] = fx[i];
    stn[3] = fx[i+1];
    stn[4] = fx[i+2];
    float xg = x[i] + 0.5*dx;    

    float fr = flux(stn);

    float dum;
    
    //cout<<xg<<", "<<fr<<", "<<sin(xg)<<", "
    //	<<fr-sin(xg)<<endl;
    if(xg>=0.25*pi and xg<=0.75*pi){
      dum=1.0;
      ref.push_back(dum);
    }
    else{dum =0.0; ref.push_back(dum);}
    //cout<<xg<<"\t"<<fr<<"\t"<<dum<<endl;
    
    xp.push_back(xg);
    fint.push_back(fr);

    i +=1;
  }while(i <imax-2);

  //interpolation on the minus
  vector<float> xl, fint_L;
  //i = imax-1;
  i = imax-2;
  do{

    /*
    stn[0] = fx[i];
    stn[1] = fx[i-1];
    stn[2] = fx[i-2];// point i+1
    stn[3] = fx[i-3];
    stn[4] = fx[i-4];
    float xn = x[i-2] - 0.5*dx;
    */
    
    stn[0] = fx[i-1];
    stn[1] = fx[i];
    stn[2] = fx[i+1]; 
    stn[3] = fx[i+2];
    stn[4] = fx[i+3];
    
    float xn = x[i+1] - 0.5*dx;
    
    float fl = flux(stn);
    

    xl.push_back(xn);
    fint_L.push_back(fl);
    
    i -=1;
  }while(i>2);

  cout<<endl;

  fp.open("weno-test-minus.csv");
  fp <<"x, fint_L\n";
  for(i=xl.size()-1; i>=0; --i){
    //cout<<xl[i]<<",\t"<<fint_L[i]<<"\n";
    fp<<xl[i]<<",\t"<<fint_L[i]<<"\n";
  }
  fp.close();
  
  //ofstream fp;
  fp.open("weno-test.csv", ios::out);
  fp << "x, fint, analytical \n";
  i=0;
  do{

    fp<<xp[i]<<",\t"<<fint[i]<<",\t"<<ref[i]<<"\n";
    i +=1;
    
  }while(i<xp.size());
  fp.close();

  //the pairs [xp,fint] and [xl,fint_L]
  int xpsize,xlsize,fintsize,fintLsize;
  xpsize = xp.size();
  xlsize = xl.size();
  fintsize = fint.size();
  fintLsize = fint_L.size();

  cout<<endl;
  cout<<"xpsize: "<<xpsize<<endl;
  cout<<"xlsize: "<<xlsize<<endl;

  cout<<"fintsize: "<<fintsize<<endl;
  cout<<"fintLsize: "<<fintLsize<<endl;

  cout<<endl;
  //cout<<"xp, xl\n";
  for(i=0; i<xlsize; ++i){
    //cout<<xp[i]<<", "<<xl[i]<<"\n";
  }

  //average interpolation, how to sort (xl,fint_L)?
  //just reverse
  vector<float> xlcopy, fint_Lcopy;
  for(i=xlsize-1; i>=0; --i){
    float copy = xl[i];
    xlcopy.push_back(copy);
    copy = fint_L[i];
    fint_Lcopy.push_back(copy);
  }

  fp.open("weno-flux_ave.csv", ios::out);
  cout<<endl;
  cout<<"xp, xl, fint, fint_L, flux_ave\n";
  fp<<"xp,fint, fint_L, flux_ave\n";
  vector<float> fl_ave;
  i=0;
  do{

    float dum = xp[i];
    for(int ii=0; ii<xlsize; ++ii){
      if(dum == xl[ii]){
	float ave=0.5*(fint[i]+fint_L[ii]);
	fl_ave.push_back(ave);

	/*
	cout<<dum<<", "<<xl[ii]<<", "<<fint[i]
	    <<", "<<fint_L[ii]<<", "<< ave
	    <<"\n";
	*/
	
	fp<<xp[i]<<", "
	  <<fint[i]<<", "
	  <<fint_L[ii]<<", "
	  <<ave<<"\n";
      }
    }
    i +=1;
    
  }while(i<xlsize-1);
  fp.close();
  
}
