#include<iostream> // for cout,cin
#include<cmath>
#include<fstream> //for ifstream, ofstream
#include<vector> // for vector
#include<algorithm> // min_element

using namespace std;

float flux2(float stn[5]) //calculate the right-hand side
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

float flux1(float stn[5]) //calculate the right-hand side
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
  double b1 = 0.3333*(4.0*stn[0]*stn[0]
	       -19.0*stn[0]*stn[1]
	       +25.0*stn[1]*stn[1]
	       +11.0*stn[0]*stn[2]
	       -31.0*stn[1]*stn[2]
	       +10.0*stn[2]*stn[2]);

  double b2 = 0.3333*(4.0*stn[1]*stn[1]
	       -13.0*stn[1]*stn[2]
	       +13.0*stn[2]*stn[2]
	       +5.0*stn[1]*stn[3]
	       -13.0*stn[2]*stn[3]
	       +4.0*stn[3]*stn[3]);

  double b3 = 0.3333*(10.0*stn[2]*stn[2]
	       -31.0*stn[2]*stn[3]
	       +25.0*stn[3]*stn[3]
	       +11.0*stn[2]*stn[4]
	       -19.0*stn[3]*stn[4]
	       +4.0*stn[4]*stn[4]);

  
  double wtilde1 = gamma1/((1.0e-6 + b1)*(1.0e-6 + b1));
  double wtilde2 = gamma2/((1.0e-6 + b2)*(1.0e-6 + b2));
  double wtilde3 = gamma3/((1.0e-6 + b3)*(1.0e-6 + b3));
  

 double w1 = wtilde1/(wtilde1+wtilde2+wtilde3);
 double w2 = wtilde2/(wtilde1+wtilde2+wtilde3);
 double w3 = wtilde3/(wtilde1+wtilde2+wtilde3);

 // double eps= 1.0e-6;
 // cout<<wtilde1<<" "<<wtilde2<<" "<<wtilde3<<" "<<eps <<endl;

  fr = w1*p1 + w2*p2 + w3*p3; // flux at xi+1/2
 
  return fr;
}

float eno(float stn[5]) //calculate the right-hand side
{
  float fr;
     
  float p1, p2, p3;
  float select;

  p1 = (3.0/8.0)*stn[0]+(-5.0/4.0)*stn[1]+(15.0/8.0)*stn[2];
  p2 = (-1.0/8.0)*stn[1]+(3.0/4.0)*stn[2]+(3.0/8.0)*stn[3];
  p3 = (3.0/8.0)*stn[2]+(3.0/4.0)*stn[3]+(-1.0/8.0)*stn[4];

  //fr = w1*p1 + w2*p2 + w3*p3; // flux at xi+1/2

  select = min(abs(p1),abs(p2));
  select = min(abs(p3),abs(select));

  if(select==abs(p1)) fr = p1;
  if(select==abs(p2)) fr = p2;
  if(select==abs(p3)) fr = p3;
  
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

  vector<float> x, fx;
  x.resize(imax);
  fx.resize(imax);
  for(i=0; i<imax; ++i){
    x[i] = i*dx;
    //     fx[i] = 0.0;
    //if(x[i]>=0.25*pi and x[i]<=0.75*pi) fx[i]=1.0; //function with jump

    fx[i] = 1.0;
    if(x[i]>=0.5*pi) fx[i]=0.0; //function with jump
    
  }
  fp.open("ref.csv",ios::out);
  fp<<"x, fx\n";
  i=0;
  do{
    fp<<x[i]<<", "<<fx[i]<<"\n";
    i +=1;
  }while(i<imax);
  fp.close();

  //forward sweep
  vector<float> xp, fint;
  float xg, fr, dum;
  i=2;
  do{

    stn[0] = fx[i-2];
    stn[1] = fx[i-1];
    stn[2] = fx[i];
    stn[3] = fx[i+1];
    stn[4] = fx[i+2];

    xg = x[i] + 0.5*dx;
    //fr = flux2(stn);
    //fr = flux1(stn);
    fr = eno(stn);

    xp.push_back(xg);
    fint.push_back(fr);
    
    i +=1;
  }while(i < imax-2);

  //backward sweep
  vector<float> xl, fl;
  i=imax-3;
  do{

    stn[0] = fx[i+2];
    stn[1] = fx[i+1];
    stn[2] = fx[i];
    stn[3] = fx[i-1];
    stn[4] = fx[i-2];

    xg = x[i] - 0.5*dx;
    //fr = flux2(stn);
    //fr = flux1(stn);
    fr = eno(stn);

    xl.push_back(xg);
    fl.push_back(fr);
    
    i -=1;
  }while(i > 2);

  fp.open("weno-classic-fr.csv", ios::out);
  fp <<"xl, fl \n";
  i=0;
  do{
    fp<<xl[i]<<", "<<fl[i]<<"\n";
    i +=1;
  }while(i<xl.size());
  fp.close();

  fp.open("weno-classic-fl.csv", ios::out);
  fp <<"xp, fint \n";
  i=0;
  do{
    fp<<xp[i]<<", "<<fint[i]<<"\n";
    i +=1;
  }while(i<xp.size());
  fp.close();


  //combining interpolation and taking average
  vector<float> xc, fc;
  float fleft, fright, ave;
  for(i=2; i<imax-2; ++i){
    stn[0] = fx[i-2];
    stn[1] = fx[i-1];
    stn[2] = fx[i]; //interface [i,i+1]
    stn[3] = fx[i+1];
    stn[4] = fx[i+2];

    xg = x[i] + 0.5*dx; // exact position

    //fleft = flux2(stn);
    //fleft = flux1(stn);
    fleft = eno(stn);
    
    stn[0] = fx[i+3];
    stn[1] = fx[i+2];
    stn[2] = fx[i+1]; // interface [i,i+1]
    stn[3] = fx[i];
    stn[4] = fx[i-1];

    //fright = flux2(stn);
    //fright = flux1(stn);
    fright = eno(stn);

    ave = 0.5*(fleft + fright);

    xc.push_back(xg);
    fc.push_back(ave);
    
  }

  cout<<"Total interface nodes: "<<xc.size()<<endl;
  fp.open("weno-ave.csv", ios::out);
  fp<<"xc, fc\n";
  i=0;
  do{
    fp<<xc[i]<<", "<<fc[i]<<"\n";
    i +=1;
  }while(i<xc.size());
  fp.close();  
  
}
