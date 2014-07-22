#include <math.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include<time.h>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstring>

using namespace std;

ifstream IN;
int NV;


int main(int argc,char *argv[]){
	string fileName=argv[1];
	vector <int> primoVertice;
	vector <int> secondoVertice;
	NV=0;
	IN.open(fileName.data(),ios::in);
	while(!IN.eof()){
		double temp1,temp2;
		IN>>temp1;
		if(temp1>NV)NV=temp1;
		IN>>temp2;
		if(temp2>NV)NV=temp2;
		
		primoVertice.push_back(temp1);
		secondoVertice.push_back(temp2);
	}
	cout<<"numero di vertici="<<NV+1<<endl;
}
