#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include<algorithm>

using namespace std;

float cdweno2(float stn[2]){
  float dfdx;
  float gamma1 = 0.5;
  float gamma2 = 0.5;

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

  //initiating
  ofstream fp;
  float stn[2];
  int imax;
  int i;
  
  cout<<"Enter imax"<<endl;
  cin>>imax;

  float dx = 10.0/(imax-1);

  vector<float> x, fx;
  x.resize(imax);
  fx.resize(imax);
  for(i=0; i<imax; ++i){
    x[i] = i*dx;
    if(x[i] >=2.0 && x[i]<=4.0){
      fx[i] = 1.0;
    }
    else{
      fx[i] = 0.0;
    }
  }


  //interpolation step
  vector<float> xp, fint;
  float flux;
  float xg;
  i=2;
  do{

    stn[0] = fx[i];
    stn[1] = fx[i+1];

    xg = x[i] + 0.5*dx;
    flux = cdweno2(stn);

    xp.push_back(xg);
    fint.push_back(flux);
    i +=1;
  }while(i<imax-2);

  fp.open("cdweno2-recon.csv",ios::out);
  fp <<"xp, fint\n";
  i=0;
  do{
    fp<<xp[i]<<","<<fint[i]<<"\n";
    i +=1;
  }while(i<xp.size());
  fp.close();
  
}
