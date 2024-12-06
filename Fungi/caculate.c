#include <stdio.h>
#include <math.h>
void addme(double * a, double * b,int *n, double * result )
{
	int count = 0;
	for (int num1 = 0; num1 < *n; ++num1)
	{
		for (int num2=num1+1; num2 < *n; ++num2)
		{
			double c=0;
		for (int m=0;m<20 ;++m)
		{  
			int m1 = num1*20+m;
			int m2 = num2*20+m;
			c=c+(a[m1]-b[m2])*(a[m1]-b[m2]);
			
			
	
		}
			c = sqrt(c);

			result[count] = c;
						count = count+1;
		}
	
	};


};
