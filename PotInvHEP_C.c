#include <stdio.h>
#include <math.h>

#define PI 3.14159265359

// Physical constants (atomic units where hbar = 1, m = 0.5)
double m = 1.0;
double hbar = 1.0;
double L = 10.0;        // Length of the box

#define xmin 0.0
#define xmax 0.4
#define xstep 0.001
#define xnum 500   // ((xmax-xmin)/xstep)
double x[xnum],fx[xnum],W[xnum];
double fx0=0.0;

#define numRes 8
double ERes[numRes];

FILE *arq;

// Potential function (inicial = 0)
double V(double x) {
    return 0.0;
}

// System of ODEs
void derivatives(double x, double fx, double W, double E, double *dfx_dx) {
    *dfx_dx = E+fx*fx-W;
}

// Runge-Kutta 4th Order Step
double rk4_step(double x, double fx, double W, double E, double h) {
    double k1_fx;
    double k2_fx;
    double k3_fx;
    double k4_fx;

    derivatives(x, fx, W, E, &k1_fx);

    derivatives(x + 0.5*h, fx + 0.5*h*k1_fx, W, E, &k2_fx);

    derivatives(x + 0.5*h, fx + 0.5*h*k2_fx, W, E, &k3_fx);

    derivatives(x + h, fx + h*k3_fx, W, E, &k4_fx);

    return fx += (h / 6.0) * (k1_fx + 2*k2_fx + 2*k3_fx + k4_fx);
}

// Finds the energy E where shoot(E) == 0 (Secant Method)
void find_Pot() {
	int j,ptr,ptrmax,iterations;
	double E;

	ptrmax=(xmax-xmin)/xstep;
	arq=fopen("arqpot.csv","w");
	for(ptr=0;ptr<ptrmax;ptr++){
		x[ptr]=xmin+ptr*xstep;
		W[ptr] = V(x[ptr]);
		fprintf(arq,"%f; %f; \n",x[ptr],W[ptr]);
	}
	printf("\n");
	for(j=0;j<numRes-1;j++) {
		E=ERes[numRes-1-1-j]-ERes[numRes-1];
		fprintf(arq,"Nivel=%f \n",E);
	    fx[0] = fx0;      // Boundary condition at x = 0

		for(ptr=1;ptr<ptrmax;ptr++){
			fx[ptr]=rk4_step(x[ptr], fx[ptr-1], W[ptr], E, xstep);
			W[ptr]=2.0*E+2.0*fx[ptr]*fx[ptr]-W[ptr];
			fprintf(arq,"%f; %f; %f; %f;\n",E,x[ptr],fx[ptr],W[ptr]);
		}
		fprintf(arq,"\n");
	}
	fprintf(arq,"*** Resultado final ***\n");
	for(ptr=1;ptr<ptrmax;ptr++){
		W[ptr] = W[ptr]+ERes[numRes-1];
		fprintf(arq,"%f; %f; \n",x[ptr],W[ptr]);
	}

	fclose(arq);
}

int main() {
    // Pesquisa do potencial 1D para os seguintes autovalores
	ERes[0]= 1.2;
	ERes[1]= 4.6;
	ERes[2]= 9.6;
	ERes[3]= 190.0;
	ERes[4]= 210.0;
	ERes[5]= 2550.0;
	ERes[6]= 3552.0;
	ERes[7]= 8360.0;
	find_Pot();
    return 0;
}