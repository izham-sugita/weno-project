/*
This file contains utility functions such as file output,
simple functions etc.
The content will be compiled first and linked later.
 */
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<iomanip>

using namespace std;


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

/*sign function*/
template<typename T>int sign(T val){
  return (T(0) < val) - (val < T(0));
}

/*Runge-Kutta step function*/
float rkstep(float u, float dt,
	     float alpha, float source)
{
  u = u + alpha*dt*source;
  return u;
}
