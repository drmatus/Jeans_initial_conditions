#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <limits.h>

float ran2(long *idum);
float gasdev(long *idum);
void odeint(float *, int, float, float, float, float, float, int *, int *, void (*), void (*));
void rkqs(float *, float *, int, float *, float, float, float *, float *, float*,  void (*));
void indexx(unsigned long, float *, unsigned long *);
void rank(unsigned long, unsigned long *, unsigned long *);

long idum;
float M_pl, R_pl, r_cut, rho_pl, D,mass_gas, rho_gas, Gamma;
float G, max_prob, r_cut_MAX, potxf;
float _X,_Y,_Z;
int kmax,kount;

float *xp, **yp, dxsav;


FILE *vel;

//For a given radius and mass, calculate the crossing time of
//the model.
float CrossTime(float R, float M){
    float gp, ff, tcr;
    gp = 2.0*2.2489E-15;
    ff = 32.*R/(3*M_PI);
    ff = pow(ff,3);
    tcr = sqrt(ff/(M*gp));
    return tcr;
}

void Print_Table(FILE *out, float *data, int n){
    int i;
    for(i=0; i<n; i++)
        fprintf(out, "% 5.10f ", data[i]);
    fprintf(out, "\n");
}

//The following funcions calculate the potenctial
//of the cilinder and the plummer sphere
float get_gravity_magnitude_at_radius(float R){ 
        float acc = 0.;
        if (R > 0.)
            acc = 4.*M_PI*G*rho_gas / Gamma * pow(D,2.) / R * (pow(1. + pow(R/D,2.),Gamma/2.) - 1.);
        return acc;
}

float potential_gas(float x, float y, float z){
    
    float sampling;
    float A0, A1, xmin, xmax, dx, xx, dA;
    float xl, xr, yl, yr;
    float r2 = x*x+y*y;
    float r = sqrt(r2);        
    int i;
/*        """
        sampling = 500

        xmin = 0.
        xmax = r.value_in(units.parsec)
        dx   = (xmax-xmin)/(sampling-1.)

        A = 0. | units.kms**2

        xx = xmin
        while xx < xmax:
            xl = xx | units.parsec
            xr = (xx+dx) | units.parsec
            yl = self.get_gravity_magnitude_at_radius(xl) 
            yr = self.get_gravity_magnitude_at_radius(xr) 
            dA = ((yl+yr)/2.)*(xr-xl)
            A += dA
            xx += dx
        """
*/
    sampling = 8.;

    A0 = 1.;// | units.kms**2
    A1 = 100.;// | units.kms**2
    //While the error is greater than some value, loop then return:
    while (pow(((A1-A0)/A0),2) > pow(0.004,2)){
        sampling *= 2;
        A0 = A1;

        xmin = 0.;
        xmax = r;//.value_in(units.parsec)
        dx   = (xmax-xmin)/(sampling-1.);
        A1 = 0.;// | units.kms**2
        xx = xmin;
        i=0;
        while (xx < xmax){
			//The potential is determined via integration
			//of the gravitational field
            xl = xx;// | units.parsec
            xr = (xx+dx);// | units.parsec
            yl = get_gravity_magnitude_at_radius(xl) ;
            yr = get_gravity_magnitude_at_radius(xr) ;
            dA = ((yl+yr)/2.)*(xr-xl);
            A1 += dA;
            xx += dx;
        }
    }
    return A1;
}

float Potencial_pl(float x, float y, float z){
    float r,r2;
    float ret;
    r2 = x*x+y*y+z*z;
    //if(r2 > r_cut*r_cut) return 0;
    r = sqrt(r2+R_pl*R_pl);
    ret = (G*M_pl)/r;
    return ret;
}

float Potencial (float x, float y, float z,float w){
    float a,b;
    a = Potencial_pl(x,y,z);
    b = potential_gas(x,y,z);
    float ret = a+w*b;
    return ret;
}

//Mass of the gas enclosed by the point (x,y)
float Mass_gas(float x, float y){
    float R2,a,b, RD;
    R2 = x*x+y*y;
    RD = R2/D;
    a = (1.+RD/D);
    b = pow(a,Gamma/2)-1;
    return b*mass_gas;
}

//Calculate the dencity from each
//component and the sum of both
float Den_Gas(float x, float y, float z){
    float R2, RD, a;
    R2 = x*x+y*y;
    RD = R2/(D*D);
    a = pow(1+RD,(Gamma-2)/2.);
    return rho_gas*a;
}

float Den_Plummer(float x, float y, float z){
    float r2,pot,a;
    r2 = x*x+y*y+z*z;
	if(r2 > r_cut*r_cut) return 0;
    a = 1.+r2/(R_pl*R_pl);
    pot = rho_pl*pow(a,-5./2.);
    return pot;
}

float Den(float x, float y, float z,float w){
    return Den_Plummer(x,y,z)+Den_Gas(x,y,z)*w;

}

//Calculate the gravitational force from
//each component and the sum of both
float Grav_pl(float x, float y, float z){
    float R2,RA,a,b;
    R2 = x*x+y*y+z*z;
    RA = R2+R_pl*R_pl;
    a = sqrt(R2)/RA;
    b = G*M_pl/sqrt(RA);
    return a*b;
}

float Grav_gas(float x, float y, float z){
    float R;
    R = sqrt(x*x+y*y);
    return 2*G*Mass_gas(x,y)/R;
}

float Grav_total(float x, float y, float z, float w,int i){
    float g1, g2;
    float tmp;
    tmp = x;
    if (i==1){
        tmp = y;
    }
    else if(i==2){
        tmp = z;
    }
    g1 = Grav_pl(x,y,z) * tmp/sqrt(x*x+y*y+z*z);
    g2 = Grav_gas(x,y,z) * tmp/sqrt(x*x+y*y);
    return g1+g2*w;
}

//Calculate  the probability that a star is at a
//distance r from the center of a plummer sphere
float P(float r){
    //return Den_Plummer(r,0,0);
    return r*r*Den_Plummer(r,0,0);
}

//Functions for each axis

//x axis:
void fx(float x, float *f, float *dxdy){
    float r2,R2,r,R;
    R2  = x*x+_Y*_Y;
    r2  = R2+_Z*_Z;
    r = sqrt(r2);
    R = sqrt(R2);
    *(dxdy+1) = Den(x,_Y,_Z,0)*(Grav_pl(x,_Y,_Z)*x/r);//+Grav_gas(x,_Y,_Z)*x/R);
}

//y axis:
void fy(float y, float *f, float *dxdy){
    float r2,R2,r,R;
    R2  = _X*_X+y*y;
    r2  = R2+_Z*_Z;
    r = sqrt(r2);
    R = sqrt(R2);
    *(dxdy+1) = Den(_X,y,_Z,0)*(Grav_pl(_X,y,_Z)*y/r);//+Grav_gas(_X,y,_Z)*y/R);

}
//z axis:
void fz(float z, float *f, float *dxdy){
    float r2,R2,r,R;
    R2  = _X*_X+_Y*_Y;
    r2  = R2+z*z;
    r = sqrt(r2);
    R = sqrt(R2);
    *(dxdy+1) = Den(_X,_Y,z,0)*Grav_pl(_X,_Y,z)*z/r;

}

//Generate n, equally spaced, bins fron xi to xf
float *Bins(float xi, float xf, int n){
    float *ret,dx;
    int i;
    ret = (float *)malloc(sizeof(float)*n);
    xi = log10(xi);
    xf = log10(xf);
    dx = (xf-xi)/n;
    for(i=0; i<n; i++)
        ret[i] = pow(10,xi+i*dx);
    return ret;
}

void Generate_Star_XYZ(float *x, float *y, float *z,float r_cut){
    float X1,X2,X3;
    float r,r2, z2;
    X2 = ran2(&idum);
    X3 = ran2(&idum);

    do{
    X1 = ran2(&idum);
    r2 = pow(X1,-2./3.)-1;
    r = pow(r2,-0.5);
    r2 = r*r;
    }while(r > r_cut);

    *z = (1-2*X2)*r;
    z2 = (*z) * (*z);
    *x = pow(r2-z2,0.5)*cos(2*M_PI*X3);
    *y = pow(r2-z2,0.5)*sin(2*M_PI*X3);

    *x = *x*R_pl;
    *y = *y*R_pl;
    *z = *z*R_pl;

}

void Generate_Star_Vel(float *vx, float *vy, float *vz, float *X){//, float pls){
    float *sig2,den,vesq;
    float potr,a;   
    int nok,nbad,i;
    float tmp[5];
   
    sig2 = (float *)malloc(sizeof(float)*3);
    for (i=0; i<3;sig2[i++]=0);

    
    _X = X[0]>0.0?X[0]:-X[0];
    _Y = X[1]>0.0?X[1]:-X[1];
    _Z = X[2]>0.0?X[2]:-X[2];
    

    /*
    _X=X[0];
    _Y=X[1];
    _Z=X[2];
    */

    kmax=0;

    a = Potencial(X[0],X[1],X[2],0);
    vesq = 2*fabs(potxf-a);
    //vesq = pls*pls;
    if (vesq == 0){
        printf("ERROR: vesq = 0.\n");
        exit(0);
    }
    //Calculate the value of the velocity dispersion
    odeint(sig2-1,1,_X,r_cut_MAX,1e-6, 0.0001,0,&nok, &nbad, fx,rkqs);
    odeint(sig2  ,1,_Y,r_cut_MAX,1e-6, 0.0001,0,&nok, &nbad, fy,rkqs);
    odeint(sig2+1,1,_Z,r_cut_MAX,1e-6, 0.0001,0,&nok, &nbad, fz,rkqs);

    den = Den(_X,_Y,_Z,0);
    for(i=0; i<3; sig2[i++]/=0.5*den);

    
    //Generate velocities until they're under the escape velocity
    do{
        *vx = sqrt(sig2[0])*gasdev(&idum);
        *vy = sqrt(sig2[1])*gasdev(&idum);
        *vz = sqrt(sig2[2])*gasdev(&idum);
    }while((*vx)*(*vx)+(*vy)*(*vy)+(*vz)*(*vz) > vesq);

    tmp[0] = sqrt(_X*_X+_Y*_Y+_Z*_Z);
    tmp[1] = sqrt(vesq);
    tmp[2] = sqrt(sig2[0]);
    tmp[3] = sqrt(sig2[1]);
    tmp[4] = sqrt(sig2[2]);

    Print_Table(vel,tmp,5);

    free(sig2);

}

float *AverageVel(float **dat, int N){
    int i;
    float *v;
    v = (float *)malloc(4*sizeof(float));
    for (i=0; i<3; v[i++]=0);
    for (i=0; i<N; i++){
        v[0]+= dat[i][5];
        v[1]+= dat[i][6];
        v[2]+= dat[i][7];
    }

    v[0] /= N;
    v[1] /= N;
    v[2] /= N;
    //printf("vx = %f\tvy = %f\tvz = %f\n",v[0],v[1],v[2]);
    return v;
}

void FixVel(float **dat, int N, float *v){
    int i;
    for (i=0; i<N; i++){
        dat[i][5] -= v[0];
        dat[i][6] -= v[1];
        dat[i][7] -= v[2];
    }
    free(v);
}

//WHY?!?!?!?!??!
//This function calculates the escape velocity using the particles
//inside an sphere of radius r^2 = x^2 + y^2 + z^2
/*
float *RealPot(float *x, float *y, float *z, int N){
    int i;
    float *ret, *r;
    unsigned long *in, *irank;
    FILE *vesc;
    float tmp[2];

    ret = (float *)malloc(N*sizeof(float));
    r   = (float *)malloc(N*sizeof(float));

    in      = (unsigned long *)malloc(N*sizeof(unsigned long));
    irank   = (unsigned long *)malloc(N*sizeof(unsigned long));

    //calculate the radial distance from the center to each particle
    for(i=0;i<N;i++)
        r[i] = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);

    //Count how many particles are inside that radius
    //indexx and rank routines taken from Numerical Recipes in C, 2nd edition.
    indexx(N,r-1,in-1);
    rank(N,in-1,irank-1);

    //for(i=0; i<N; i++) irank[i]-=1;

    //Now, the escape velocity for each particle is saved:
    for(i=0; i<N; i++){
        ret[i] = sqrt(2*0.5*irank[i]*G/r[i]);
    }

    vesc = fopen("vesc.log","w");
    for (i=0; i<N; i++){
        tmp[0] = r[i];
        tmp[1] = ret[i];
        Print_Table(vesc,tmp,2);
    }
    fclose(vesc);

    return ret;
}
*/
int main(){
    float data[8],x,y,z,m,tcr;
    float **todos, *X,*Y,*Z, *VX,*VY,*VZ,*vesc;
    int i,N;
    FILE *tabla;
    char nombre[80];
    float magico = sqrt(1.695);


    //Get the parameters:
    //Number of stars
    printf("Number of stars: ");
    scanf("%d", &N);
    //Plummer radius in parsecs
    printf("Plummer radius [pc]: ");
    scanf("%f", &R_pl);
    //Name of the model:
    fscanf(fopen("name","r"),"%s", nombre);
    strcat(nombre,".pos");
    //Seed for the random number generator:
    i = open("/dev/urandom",O_RDONLY);
    read(i,&idum,sizeof(long)-5);
    close(i);
    //idum = -41;
    if (idum > 0) idum *= -1;

    //Mass of one star:
    m = 0.5;
    //Plummer sphere parameters:
    //mass: Number of stars x 0.5 M_sun
    M_pl = m*N;
    r_cut = 10;

    rho_pl = 3*M_pl/(4.*M_PI*R_pl*R_pl*R_pl);

    //gas column parameters:
    r_cut_MAX = 8.5; //pc, the limit of the simulation area
    Gamma =0.225 ;
    D = 5e-6;
    mass_gas = 53.0690998914; //Via fitting the line-mass function
    rho_gas =mass_gas*Gamma/(2*M_PI*D*D); //rho0

    //Other parameters:
    
    G = 0.0043; //gravity constant in M_sun, pc and kms
    max_prob = P(sqrt(2./3.)*R_pl); //Peak of prob. distribution (Plummer)

    tabla = fopen(nombre,"w");
    printf("\n");

    potxf = Potencial(r_cut_MAX,0,0,0);
    printf("seed: %ld\nRK\n\n", idum);
    todos = (float **)malloc(sizeof(float *)*N);
    X = (float *)malloc(sizeof(float)*N);
    Y = (float *)malloc(sizeof(float)*N);
    Z = (float *)malloc(sizeof(float)*N);
    VX = (float *)malloc(sizeof(float)*N);
    VY = (float *)malloc(sizeof(float)*N);
    VZ = (float *)malloc(sizeof(float)*N);

    //Logfile for the escape velocities:
    vel = fopen("vel.log","w");
    if (!tabla)
        printf("No Tabla\n");
    
    tcr = CrossTime(R_pl,M_pl)/1000.;
    printf("Crossing Time: %f kyr\n",tcr);
    printf("Recomended timestep: %f kyr\n",tcr/20.);
    printf("Generating Positions\n");
    for (i=0; i<N; i++){
        todos[i] = (float *)malloc(sizeof(float)*8);
        Generate_Star_XYZ(X+i,Y+i,Z+i,r_cut);
        todos[i][0] = m;   //Star's mass
        todos[i][1] = 0.5; //Star's radius
        todos[i][2] = X[i];
        todos[i][3] = Y[i];
        todos[i][4] = Z[i];
        /*
        X[i] = x;
        Y[i] = y;  //Star's position
        Z[i] = z;
        */

        printf("%d/%d done.\r",i,N);
        fflush(stdout);
    }
    //vesc = RealPot(X,Y,Z,N);
    printf("%d/%d done.\n",N,N);
    printf("Generating Velocities\n");
    for(i=0; i<N; i++){
        Generate_Star_Vel(VX+i,VY+i,VZ+i,todos[i]+2);//,vesc[i]);
        //1.695 es un factor magico que hace
        //que todo funcione. Ni idea por que
        //Update: esta relacionado con el uso
        //de "Standard Units":
        //M=G=R_v=1
        todos[i][5] = VX[i]/magico;
        todos[i][6] = VY[i]/magico;
        todos[i][7] = VZ[i]/magico;
        
        printf("%d/%d done.\r",i,N);
        fflush(stdout);
    }
    printf("%d/%d done.\n",N,N);
    
    float *v;
    v = AverageVel(todos,N);
    FixVel(todos,N,v);
    v = AverageVel(todos,N);

    printf("All done.         \n");
    for (i=0; i<N; i++){
        Print_Table(tabla, todos[i], 8);
    }

    //Close the files!!
    fclose(tabla);
    fclose(vel);

    //Free the arrays!!!
    for (i=0; i<N; i++) free(todos[i]);
    free(todos);
    free(X);
    free(Y);
    free(Z);
    free(VX);
    free(VY);
    free(VZ);
    free(vesc);

}
