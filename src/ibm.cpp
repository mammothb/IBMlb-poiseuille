#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fluid.h"
#include "mesh.h"
#include "ibm.h"

using namespace std;

/**
 *	Phi_2 function to calculate Dirac delta function with support of two nodes
 *	@param double r, Distance between nodes (Lagrangian - Eulerian)
 *  @return d, distance weighting
 */
double phi_2(double r)
{
	double phi = 0.0;
	if (0.0 <= fabs(r) && fabs(r) <= 1.0) phi = 1.0 - fabs(r);
  if (1.0 <= fabs(r)) phi = 0.0;
  return phi;
}

/**
 *	Phi_3 function to calculate Dirac delta function with four support
 *	@param double r, Distance between nodes (Lagrangian - Eulerian)
 *  @return double d, distance weighting 700, 10100
 */
double phi_3(double r)
{
	double phi = 0.0;
	if (0.0 <= fabs(r) && fabs(r) <= 0.5) phi = (1.0 / 3.0 * (1.0 + sqrt(1.0 - 3.0 * r * r));
  if (0.5 <= fabs(r) && fabs(r) <= 1.5) phi = (1.0 / 6.0 * (5.0 - 3.0 + fabs(r) - sqrt(-2.0 + 6.0 * fabs(r) - 3.0 * r * r));
  if (1.5 <= fabs(r)) phi = 0.0;
  return phi;
}

/**
 *	Phi_4 function to calculate Dirac delta function with four support
 *	@param double r, Distance between nodes (Lagrangian - Eulerian)
 *  @return d, distance weighting
 */
double phi_4(double r)
{
	double phi = 0.0;
	if (0.0 <= fabs(r) && fabs(r) <= 1.0) phi = (0.125 * (3.0 - 2.0 * fabs(r) + sqrt(1.0 + 4.0 * fabs(r) - (4.0 * r * r))));
  if (1.0 <= fabs(r) && fabs(r) <= 2.0) phi = (0.125 * (5.0 - 2.0 * fabs(r) - sqrt(-7.0 + 12.0 * fabs(r) - (4.0 * r * r))));
  if(fabs(r) >= 2.0) phi = 0.0;
  return phi;
}

/**
 *	Dirac_2 function to calculate Dirac delta function with four support
 *	@param x, Vector distance between nodes (Lagrangian - Eulerian)
 *  @return d, distance weighting
 */
double dirac_2(double *x)
{
    double d = phi_2(x[0]) * phi_2(x[1]) * phi_2(x[2]);
    return d;
}

/**
 *	Dirac_3 function to calculate Dirac delta function with four support
 *	@param x, Vector distance between nodes (Lagrangian - Eulerian)
 *  @return d, distance weighting1
 */
double dirac_3(double *x)
{
    double d = phi_3(x[0]) * phi_3(x[1]) * phi_3(x[2]);
    return d;
}

/**
 *	Dirac_4 function to calculate Dirac delta function with four support
 *	@param x, Vector distance between nodes (Lagrangian - Eulerian)
 *  @return d, distance weighting1
 */
double dirac_4(double *x)
{
    double d = phi_4(x[0]) * phi_4(x[1]) * phi_4(x[2]);
    return d;
}

/**
 * TODO: Implement spread function
 */
void spread(fluid fluido, mesh membrana, int x, int y, int z)
{
	// Browse all nodes in the fluid
	double pos[3]={0.0,0.0,0.0}, distancia[3]={0.0,0.0,0.0}, delta, df[3]={0.0,0.0,0.0}, fNodo[3]={0.0,0.0,0.0};
	double a, A, b, B, c, C;
	int nodos = membrana.darNumeroNodos();
	for(int u=0;u<nodos;u++)
	{
		membrana.darPosNodo(u, pos);
		membrana.darFuerzaNodo(u, fNodo);
		a = pos[0]-3.0;
		A = pos[0]+3.0;
		b = pos[1]-3.0;
		B = pos[1]+3.0;
		c = pos[2]-3.0;
		C = pos[2]+3.0;

		for(int i = (int) a;i<A;i++)
			for(int j = (int) b;j<B;j++)
				for(int k=(int) c;k<C;k++)
				{
					int posj = j;
					distancia[0]=pos[0]-i;
					distancia[1]=pos[1]-j;
					distancia[2]=pos[2]-k;
					delta = dirac_4(distancia);
					df[0]=fNodo[0]*delta;
					df[1]=fNodo[1]*delta;
					df[2]=fNodo[2]*delta;
					if(j<1)
					{
						posj = j+y-1;
					}else if(j>=y)
					{
						posj = j-y;
					}
					fluido.addFuerza(i,posj,k,df);
				}
		df[0] = 0.0;
		df[1] = 0.0;
		df[2] = 0.0;
	}
}


/**
 * TODO: Implement the interpolation function
 */
void interpolation(fluid fluido, mesh membrana, int x, int y, int z)
{
	// Loop through all mesh nodes
	double pos[3], distancia[3], delta, a, A, b, B, c, C, ux=0.0, uy=0.0, uz=0.0;
	int nNodos = membrana.darNumeroNodos();
	for(int u=0;u<nNodos;u++)
	{
		membrana.darPosNodo(u, pos);
		a = pos[0]-3.0;
		A = pos[0]+3.0;
		b = pos[1]-3.0;
		B = pos[1]+3.0;
		c = pos[2]-3.0;
		C = pos[2]+3.0;

		for(int i = (int) a; i < A; i++)
			for(int j = (int) b; j < B; j++)
				for(int k= (int) c; k < C; k++)
				{
					int posj = j;
					distancia[0]=pos[0]-i;
					distancia[1]=pos[1]-j;
					distancia[2]=pos[2]-k;
					delta = dirac_4(distancia);
					if(j<1)
					{
						posj = j+y-1;
					}else if(j>=y)
					{
						posj = j-y;
					}
					ux+=delta*fluido.darVelocidad(i,posj,k,0);
					uy+=delta*fluido.darVelocidad(i,posj,k,1);
					uz+=delta*fluido.darVelocidad(i,posj,k,2);
				}
		membrana.setVelocidad(u,ux,uy,uz);
		ux=0.0;
		uy=0.0;
		uz=0.0;
	}
}
