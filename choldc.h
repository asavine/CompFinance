#pragma once

//	Choleski decomposition, inspired by Numerical Recipes

#include "matrix.h"

template <class T, class U>
void choldc(const matrix<U>& in, matrix<T>& out)
{
	int n = in.rows();
	double sum;

	for(int i=0; i<n; ++i)
	{
		auto* ai = in[i];
		for(int j=i; j<n; ++j)
		{
			sum = ai[j];
			auto* aj = in[j];
			auto* pj = out[j];
			for(int k=i-1; k>=0; --k)
			{
				sum -= ai[k] * aj[k];
			}
			if(i == j)
			{      
                if(sum < - 1.0e-15)
                {
                    throw "choldc : matrix not positive definite";
					return false;
				}
                if(fabs(sum) <= 1.0e-15) sum = 0.0;
				pi[i] = sqrt(sum);
			}
			else
			{
				if(pi == 0.0)
                {
                    pj[i] = 0.0;
                }
                else
                {
                    pj[i] = sum/pi[i];
                }
			}
		}
	}
}