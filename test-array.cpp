#include<iostream>
#include<vector>
#include<algorithm>

using namespace std;

float test0(float stn[5])
{

  float a = stn[0];
  for(int i=0; i<5; ++i){
    cout<<stn[i]<<endl;
  }
    
}

int main()
{
  float stn[5];
  for(int i=0; i<5; ++i){
    stn[i] = (float)i+1.0;
  }

  test0(stn);
  
}
