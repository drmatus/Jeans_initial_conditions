#include <math.h>
#include "nrutil.h"

#define TINY 1.0e-25 //A small number.
#define FREERETURN {free_vector(d,1,n);free_vector(c,1,n);return;}

/*Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y, and an error estimate dy. If P(x) is the polynomial of degree N-1 such that P(xa_i)=ya_i. i=1,...,n, then the returned value y=P(x)
 *  */

void polint (float xa[], float ya[], int n, float x, float *y, float *dy){
    int i,m,ns=1;
    float den,dif,dift,ho,hp,w;
    float *c,*d;
    dif=fabs(x-xa[1]);
    c=vector(1,n);
    d=vector(1,n);
    for (i=1;i<=n;i++) { //Here we find the index ns of the closest table entry,
        if ( (dift=fabs(x-xa[i])) < dif) {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i]; //and initialize the tableau of c’s and d’s
        d[i]=ya[i];
    }
    *y=ya[ns--];  //This is the initial approximation to y.
    for (m=1;m<n;m++) { //For each column of the tableau,
        for (i=1;i<=n-m;i++) { //we loop over the current c’s and d’s and update
            ho=xa[i]-x;        //them.
            hp=xa[i+m]-x;
            w=c[i+1]-d[i];
            if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
            //This error can occur only if two input xa’s are (to within roundoff) identical.
            d[i]=hp*den; //Here the c’s and d’s are updated.
            c[i]=ho*den;
        }
        *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
        /*After each column in the tableau is completed, we decide which correction, c or d,
         * we want to add to our accumulating value of y, i.e., which path to take through the
         * tableau—forking up or down. We do this in such a way as to take the most “straight
         * line” route through the tableau to its apex, updating ns accordingly to keep track of
         * where we are. This route keeps the partial approximations centered (insofar as possible)
         * on the target x. The last dy added is thus the error indication.
         */
      }
      free_vector(d,1,n);
      free_vector(c,1,n);
}

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[])
/*Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f (xi ), with
 * x1 < x2 < . . . < xN , and given values yp1 and ypn for the first derivative of the interpolating
 * function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
 * the second derivatives of the interpolating function at the tabulated points xi . If yp1 and/or
 * ypn are equal to 1 × 1030 or larger, the routine is signaled to set the corresponding boundary
 * condition for a natural spline, with zero second derivative on that boundary.*/
{
    int i,k;
    float p,qn,sig,un,*u;
    u=vector(1,n-1);
    if (yp1 > 0.99e30)     //The lower boundary condition is set either to be “natural”
        y2[1]=u[1]=0.0;
        else {                 //or else to have a specified first derivative.
            y2[1] = -0.5;
            u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
        }
        for (i=2;i<=n-1;i++) { //This is the decomposition loop of the tridiagonal al-
            sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]); //gorithm. y2 and u are used for tem-
            p=sig*y2[i-1]+2.0; //porary storage of the decomposed
            y2[i]=(sig-1.0)/p; //factors.
            u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
            u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
        }
        if (ypn > 0.99e30)     //The upper boundary condition is set either to be “natural”
            qn=un=0.0;
        else {                 //or else to have a specified first derivative.
            qn=0.5;
            un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
        }
        y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
        for (k=n-1;k>=1;k--)   //This is the backsubstitution loop of the tridiagonal algorithm.
            y2[k]=y2[k]*y2[k+1]+u[k];
        free_vector(u,1,n-1);
}

void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)
/*Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai ’s in order),
 * and given the array y2a[1..n], which is the output from spline above, and given a value of
 * x, this routine returns a cubic-spline interpolated value y.*/
{
    void nrerror(char error_text[]);
    int klo,khi,k;
    float h,b,a;
    klo=1;                    //We will find the right place in the table by means of
    khi=n;                    //bisection. This is optimal if sequential calls to this
    while (khi-klo > 1) {     //routine are at random values of x. If sequential calls
        k=(khi+klo) >> 1;     //are in order, and closely spaced, one would do better
        if (xa[k] > x) khi=k; //to store previous values of klo and khi and test if
        else klo=k;           //they remain appropriate on the next call.
    }   //klo and khi now bracket the input value of x.
    h=xa[khi]-xa[klo];
    if (h == 0.0) nrerror("Bad xa input to routine splint"); //The xa’s must be dis-
    a=(xa[khi]-x)/h;                                         //tinct.
    b=(x-xa[klo])/h; //Cubic spline polynomial is now evaluated.
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void ratint(float xa[], float ya[], int n, float x, float *y, float *dy)
/*Given arrays xa[1..n] and ya[1..n], and given a value of x, this routine returns a value of
 * y and an accuracy estimate dy. The value returned is that of the diagonal rational function,
 * evaluated at x, which passes through the n points (xai , yai ), i = 1...n.*/
{
    int m,i,ns=1;
    float w,t,hh,h,dd,*c,*d;
    c=vector(1,n);
    d=vector(1,n);
    hh=fabs(x-xa[1]);
    for (i=1;i<=n;i++) {
        h=fabs(x-xa[i]);
        if (h == 0.0) {
            *y=ya[i];
            *dy=0.0;
            FREERETURN
        } else if (h < hh) {
            ns=i;
            hh=h;
        }
        c[i]=ya[i];
        d[i]=ya[i]+TINY; //The TINY part is needed to prevent a rare zero-over-zero condition.
    }
    *y=ya[ns--];
    for (m=1;m<n;m++) {
        for (i=1;i<=n-m;i++) {
            w=c[i+1]-d[i];
            h=xa[i+m]-x; //h will never be zero, since this was tested in the initializing loop.
            t=(xa[i]-x)*d[i]/h;
            dd=t-c[i+1];
            if (dd == 0.0) nrerror("Error in routine ratint");
            /*This error condition indicates that the interpolating function has a pole at the
            requested value of x.*/
            dd=w/dd;
            d[i]=c[i+1]*dd;
            c[i]=t*dd;
        }
        *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
    FREERETURN
}
