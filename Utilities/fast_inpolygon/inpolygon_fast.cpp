
#include <mex.h>
#include <iostream>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#ifndef round
#define round(x) (x>0 ? (double)(int)(x+0.5) : (double)(int)(x-0.5))
#endif

template<typename _Tp> class LessThanIdx
{
public:
	LessThanIdx( const _Tp* _arr ) : arr(_arr) {}
	bool operator()(int a, int b) const { return arr[a] < arr[b]; }
	const _Tp* arr;
};


double* sort_array(double* arr, int length, int** index)
{
	double* sorted=(double*)malloc(sizeof(double)*length);
	*index=(int*)malloc(sizeof(int)*length);

	for(int j = 0; j < length; j++ )
		(*index)[j] = j;

	// matlab compatible stable sorting: Otherwise, errors occur.
	std::stable_sort( (*index), (*index) + length-1, LessThanIdx<double>(arr) );

	for(int j = 0; j < length; j++ )
		sorted[j]=arr[(*index)[j]];

	return sorted;
}


void find_min_max(double* arr, int length, double* minval, double* maxval)
{
	register int i;
	*maxval = arr[0];
	*minval = *maxval;

	for(i=0; i<length;i++)
	{
		const double arrVal=arr[i];
		if (arrVal >  *maxval ) *maxval = arrVal;
		else if (arrVal <  *minval ) *minval = arrVal;
	}
}


void test_points_in_polygon(double* X, double* Y, int lengthPts, double* PolyX, double* PolyY, int nnode, bool* In, bool* On)
{
	int* ex, *ey;
	int i,k,n,numIn=0;
	const double TOL = 1.0e-12;
	double yMin=0, yMax=0, xMin=0, xMax=0, tol, dxyX, dxyY;
	double* Xt=X, *Yt=Y;
	double* pxt=PolyX, *pyt=PolyY;
	int* Ind;

	ex=(int*)malloc(sizeof(int)*nnode);
	ey=(int*)malloc(sizeof(int)*nnode);

	for (i=0;i<nnode-1;i++)
	{
		ex[i]=i;
		ey[i]=i+1;
	}
	ex[nnode-1]=nnode-1;
	ey[nnode-1]=0;
	
	n=lengthPts;

	find_min_max(X, n, &xMin, &xMax);
	find_min_max(Y, n, &yMin, &yMax);
	
	dxyX=xMax-xMin;
	dxyY=yMax-yMin;

	if (dxyX>dxyY)
	{
		double* temp=Xt;
		Xt=Yt;
		Yt=temp;
		temp=pxt;
		pxt=pyt;
		pyt=temp;
	}
	
	tol = TOL*min(dxyX,dxyY);
	
	double* ys= sort_array(Yt, n, &Ind);
	double *xs=(double*)malloc(sizeof(double)*n);
	for (int i=0;i <n;i++)
		xs[i]=Xt[Ind[i]];
	
	int *cn=(int*)malloc(sizeof(int)*n);
	int *on=(int*)malloc(sizeof(int)*n);
	memset(cn, 0, sizeof(char)*n);
	memset(on, 0, sizeof(char)*n);
	
	for (k=0;k<nnode;k++)
	{
		int j;
		const int n1=ex[k];
		const int n2=ey[k];

		double y1=pyt[n1];
		double y2=pyt[n2];
		double x1, x2, xmin, xmax;
		int start;
		
		if (y1<y2)
		{
			x1=pxt[n1];
			x2=pxt[n2];
		}
		else
		{
			double yt=y1;
			y1=y2;
			y2=yt;
			x1=pxt[n2];
			x2=pxt[n1];
		}

		if (x1>x2)
		{
			xmin = x2;
			xmax = x1;
		}
		else
		{
			xmin = x1;
			xmax = x2;
		}

		if (ys[0]>=y1)
			start=1;
		else if (ys[n-1]<y1)
			start=n+1;
		else
		{
			int lower = 1;
			int upper = n;
			for (j=1; j<=n; j++)
			{
				start = round(0.5*(lower+upper));
				if (ys[start]<y1)
					lower = start;
				else if (ys[start-1]<y1)
					break;
				else
					upper = start;
			}
		}

		for (j=(start-1); j<n; j++)
		{
			double yCur=ys[j];
			double xCur;

			if (yCur<=y2)
			{
				xCur=xs[j];
				if (xCur>=xmin)
				{
					if (xCur<=xmax)
					{
						on[j]= on[j] || (abs((y2-yCur)*(x1-xCur)-(y1-yCur)*(x2-xCur))<tol);

						if ( (yCur<y2) && ((y2-y1)*(xCur-x1)<(yCur-y1)*(x2-x1) ) )
							cn[j] = !cn[j];
					}
				}
				else if(yCur<y2)
				{
					cn[j]=!cn[j];
				}
			}
			else
				break;
		}
	}

	for (k=0; k<n; k++)
	{
		In[Ind[k]]=cn[k] | on[k];
		On[Ind[k]]=on[k];
	}

	free(ex);
	free(ey);
	free(xs);
	free(ys);
	free(Ind);
	free(cn);
	free(on);
}



void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double* X=(double*)mxGetPr(prhs[0]);
	double* Y=(double*)mxGetPr(prhs[1]);
	double* PolyX=(double*)mxGetPr(prhs[2]);
	double* PolyY=(double*)mxGetPr(prhs[3]);
	
	const mwSize* dims= mxGetDimensions(prhs[0]);
	const mwSize* dims2= mxGetDimensions(prhs[2]);

	plhs[0]=mxCreateNumericMatrix(dims[0],dims[1], mxLOGICAL_CLASS, mxREAL);
	plhs[1]=mxCreateNumericMatrix(dims[0],dims[1], mxLOGICAL_CLASS, mxREAL);

	bool* In=(bool*)mxGetPr(plhs[0]);
	bool* On=(bool*)mxGetPr(plhs[1]);

	test_points_in_polygon(X, Y, max(dims[0], dims[1]), PolyX, PolyY, max(dims2[0], dims2[1]), In, On);
}
