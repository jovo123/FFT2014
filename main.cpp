//Author Josephine Vo
//Data 10/8/2014
//Fast Fourier Transform implementation using main (not object oriented)
//compile with g++, run with ./a.out

#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

const double PI = 3.1415926536;

vector<double> fft(int N, vector<double> &Data);

unsigned int reverseBits(unsigned int n);

int main(){

	vector<double> test;
	double input;
	int length;
	
	cout << "Enter how many samples you have: ";
	cin >> length;
	length = 2*length;
	for ( int i = 0; i < length; i += 2){
		cout << "Enter a sample (real):";
		cin >> input;
		test.push_back(input);
		test.push_back(0.0);
	}
	

	test = fft(length, test);
	length = test.size();
	
	cout << "FFT frequency Data ( real, imaginary ): ";
	for ( int j = 0; j < length; j += 2){
		cout << "( " << test[j] << ", " << test[j+1] << " ) ";
	}
		
	cout << endl;
	test.clear();
	return 0;	
		 
};

vector<double> fft(int N, vector<double> &Data){

	//check that there are 2^m inputs
	int k = ceil(log2(N));
	for ( int aug = N; aug < pow(2,k); aug++){
		Data.push_back(0.0); 
		N++;
	}
	

	//binary inversion (note that the indexes 
    //start from 0 witch means that the
    //real part of the complex is on the even-indexes 
    //and the complex part is on the odd-indexes
	int j=0;
    for (int i=0;i<N/2;i+=2) {
        if (j > i) {
            //swap the real part
            swap(Data[j],Data[i]);
            //swap the complex part
            swap(Data[j+1],Data[i+1]);
            // checks if the changes occurs in the first half
            // and use the mirrored effect on the second half
            if((j/2)<(N/4)){
                //swap the real part
                swap(Data[(N-(i+2))],Data[(N-(j+2))]);
                //swap the complex part
                swap(Data[(N-(i+2))+1],Data[(N-(j+2))+1]);
            }
        }
        int m = N/2;
        while (m >= 2 && j >= m) {
            j -= m;
            m = m/2;
        }
        j += m;
    }
	
	cout << "swapped vector: ";
	for ( int j = 0; j < N; j += 2){
		cout << "( " << Data[j] << ", " << Data[j+1] << " ) ";
	}

	//Danielson-Lanczos; 
	int mmax = 2;
	//loop through the stages (while mmax = 2, 4, ... N/2 is the number of btfy blocks per stage)
	while(N > mmax){ 
		int istep = 2 * mmax; //step to index for finding odd part of each btfy
		//loop through butterfly blocks
		for( int m = 0; m < mmax; m+=2){ 
			//compute a butterfly in each block			
			for (int i = 0; i < N; i += istep){ 
				int j = i + mmax;	
				double theta = (-2*(PI)/mmax)*((i/2)%mmax);
				double Wr = cos(theta);
				double Wi = sin(theta);
				double tempr = Data[j]*Wr + Data[j+1]*Wi; //real 
				double tempi = Data[j+1]*Wr - Data[j]*Wi;//imaginary
				Data[j] = Data[i] - tempr;
				Data[j+1] = Data[i+1] - tempi;
				Data[i] = Data[i] + tempr;
				Data[i+1] = Data[i+1] + tempi;
			}
		}
		mmax = istep;
	}
	
	cout << endl;
	return Data;
	
};
