#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<iomanip>

using namespace std;

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

int main()
{

  string filename;
  int timestepmax =10;
  int iter=0;
  int imax =100;
  
  vector<float> x, fx;
  x.resize(imax);
  fx.resize(imax);

  for(int i=0; i<imax; ++i){
    x[i] = 1.0;
    fx[i] = 1.0;
  }
  
  do{

    output(iter, x, fx);

    iter +=1;
  }while(iter<timestepmax);
  
}
