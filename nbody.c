#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <string.h>

//                                                                            N-BODY SYSTEM SIMULATION 

//                                                                                 Elia Di Santo 

#define G 6.67e-11
#define R 10
#define M 100

// Funzione per generare un numero casuale tra 0 e 1
double rand01() {
    return (double) rand() / RAND_MAX;
}

// Funzione per generare un numero casuale uniforme tra a e b
double rand_uniform(double a, double b) {
    return a + (b - a) * rand01();
}

int main() {

  srand(time(NULL));

  FILE* fp;
  FILE* fc;
  FILE* fl;
  FILE* fd;
  FILE* fg;
  FILE* fv;
  FILE* fh;
  FILE* fq;
  FILE* ft;
  FILE* ff;
  FILE* fr;
  FILE* fn;
  FILE* fm;
  FILE* f14;
    
  int n = 1000; // Numero di stelle

  double* x;
  double* y;
  double* z;
  double* densities;
  double* contatore;
  double* rr;
  double* vx;
  double* vy;
  double* vz;
  double* vv;
  double* uu;
  double* rho;
  double* m;

  double rhor,ms;

  int scelta;

  printf("\n1 sferica unif\n2 sferica non unif\n3 sferoidale\n\n");
  scanf("%d",&scelta);

  x = (double*)calloc(n+1, sizeof(double));
  y = (double*)calloc(n+1, sizeof(double));
  z = (double*)calloc(n+1, sizeof(double));
  vx = (double*)calloc(n+1, sizeof(double));
  vy = (double*)calloc(n+1, sizeof(double));
  vz = (double*)calloc(n+1, sizeof(double));
  vv = (double*)calloc(n+1, sizeof(double));
  uu = (double*)calloc(n+1, sizeof(double));
  densities = (double*)calloc(n+1, sizeof(double));
  contatore = (double*)calloc(n+1, sizeof(double));
  rr = (double*)calloc(n+1, sizeof(double));
  rho = (double*)calloc(n+1, sizeof(double));
  m = (double*)calloc(n+1, sizeof(double));
    
  
  if((fp=fopen("uscita.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
  if((fc=fopen("uscita1.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
  if((fl=fopen("uscita2.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
   if((fd=fopen("uscita3.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
   if((fg=fopen("uscita4.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
   if((fv=fopen("uscita5.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
   if((fh=fopen("uscita6.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
    if((fq=fopen("uscita7.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
    if((ft=fopen("uscita8.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
   if((ff=fopen("uscita9.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
   if((fr=fopen("uscita10.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
   if((fn=fopen("uscita11.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
   if((fm=fopen("uscita12.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
   if((f14=fopen("uscita14.dat","w+"))==NULL){
    printf("Errore nell'apertura del file\n");
    exit(EXIT_FAILURE);
  }
   
   double r;
   
   int i,j,a;
   
   double ii, jj, kk, d, vii, vjj, vkk, vd, U, llll=0,theta,phi,xxx,yyy,zzz;
   if(scelta==1){
   for (i = 0; i < n; i++) {                                                 // i=0
     if(i!=0){
       do {
	 
	 
	 
	 ii = rand_uniform(-R, R);
	 
	 jj = rand_uniform(-R, R);
	 
	 kk = rand_uniform(-R, R);
	 
	 d = sqrt(ii*ii + jj*jj + kk*kk);
	 
	 //if(i<n/2){
	 //m[i] = 1;
	 //} else {
	 //m[i] = 0.4;
	 //}
	 
	 m[i] = rand_uniform(0,1);
	 
	 
       } while (d > R); // Ripeti finché il punto non si trova all'interno della sfera
       
       x[i] = ii;
       y[i] = jj;
       z[i] = kk;
       rr[i] = d;

     } else {
       x[i] =0;
       y[i] = 0;
       z[i] = 0;
       rr[i] = 0;
       m[i] = M;
     }
     if(i==0){
       fprintf(fn,"%lf\t%lf\t%lf\t%lf\t%lf\n",x[i],y[i],z[i],rr[i],m[i]);
     }
     if(i<n/2 && i!=0){
       fprintf(fc,"%lf\t%lf\t%lf\t%lf\t%lf\n", x[i], y[i], z[i], rr[i], m[i]);
     } else {
       fprintf(fg,"%lf\t%lf\t%lf\t%lf\t%lf\n", x[i], y[i], z[i], rr[i], m[i]);
     }
   }
   }
   
   /////////////////////////////////////////////// (x_0,y_0,z_0) //////////////////////////////////////////////////////
   
   if(scelta==2){
   for(i=0;i<n;i++){
     if(i!=0){
       do{
	 r = R*((double)rand())/(RAND_MAX);
	 theta = 2*M_PI*((double)rand())/(RAND_MAX);
	 phi = M_PI*((double)rand())/(RAND_MAX);
	 
	 xxx = r*sin(phi)*cos(theta);
	 yyy = r*sin(phi)*sin(theta);
	 zzz = r*cos(phi);
	 m[i] = rand_uniform(0,1);
	 d = sqrt(xxx*xxx + yyy*yyy + zzz*zzz);
       }while(d>R);
       x[i] = xxx;
       y[i] = yyy;
       z[i] = zzz;
       rr[i] = d;
     } else {
       x[i] = 0;
       y[i] = 0;
       z[i] = 0;
       rr[i] = 0;
     }
     rhor = 3*M_PI*m[i]/(rr[i]*rr[i]*rr[i]);
     if(i==0){
       fprintf(fn,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x[i],y[i],z[i],rr[i],m[i],rhor);
     }
     if(i<n/2 && i!=0){
       fprintf(ff,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x[i],y[i],z[i],rr[i],m[i],rhor);
     } else {
       fprintf(fr,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x[i], y[i], z[i], rr[i],m[i],rhor);
     }
   }
   ms=0;
   for(i=1;i<n;i++){
     ms+=m[i];
   }
   printf("\n%lf\n",ms);
   ms=0;
   }
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   if(scelta==3){
   FILE* f13;
   if((f13=fopen("uscita13.dat","w+"))==NULL){
     printf("Errore nell'apertura del file\n");
     exit(EXIT_FAILURE);
   }
   
   
   double xx,yy,zz,nu,mu;
   
   for (i = 0; i < n; i++) {
     do{
       if(i!=0){
	 nu = rand_uniform(-M_PI, M_PI);
	 mu = rand_uniform(0,0.2);
	 phi = rand_uniform(-M_PI,M_PI);
	 r = rand_uniform(-R,R);
	 xx = r*cosh(mu)*cos(nu)*cos(phi);
	 yy = r*cosh(mu)*cos(nu)*sin(phi);
	 zz = r*sinh(mu)*sin(nu);
	 m[i] = rand_uniform(0,1);
	 d = sqrt(xx*xx+yy*yy+zz*zz);
       } else {
	 xx=0;
	 yy=0;
	 zz=0;
	 fprintf(fn,"%lf\t%lf\t%lf\n",x[i],y[i],z[i]);
       }
     }while(d>R && zz>0.2);
     //printf("%d\n",i);
     
     x[i]=xx;
     y[i]=yy;
     z[i]=zz;
     rr[i]=d;
     
     fprintf(f13,"%lf\t%lf\t%lf\t%lf\t%lf\n", x[i], y[i], z[i],rr[i],m[i]);
   }
   
   }
   
   
   //////////////////////////////// CALCOLATORE DI POTENZIALE PER OGNI STELLA ////////////////////////////////////////////
   for(i=0;i<n;i++){
     for(j=0;j<n;j++){
       if(i!=j){
	 U += G*m[i]*m[j]/(fabs(rr[i]-rr[j]));
       }
     }
     uu[i] = sqrt(2*U);
     fprintf(fl,"%lf\t%lf\t%lf\n",rr[i],U,uu[i]);
     U=0;
   }
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   
   ///////////////////////////////////////////// SELETTORE DI v_0 ///////////////////////////////////////////////////////
   for(i=0;i<n;i++){
     if(i!=0){
       do{
	 vii = rand_uniform(-uu[i],uu[i]);
	 vjj = rand_uniform(-uu[i],uu[i]);
	 vkk = rand_uniform(-uu[i],uu[i]);
	 
	 vd = sqrt(vii*vii + vjj*vjj + vkk*vkk);
       } while(vd > uu[i]);
       vx[i] = vii;
       vy[i] = vjj;
       vz[i] = vkk;
       
       vv[i] = vd;

     } else {
       vx[i] = 0;
       vy[i] = 0;
       vz[i] = 0;
       
       vv[i] = 0;
       
     }

     fprintf(fd,"%lf\t%lf\t%lf\t%lf\n",vx[i],vy[i],vz[i],vv[i]);
     
   }
   
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   
   //////////////////////////////////////////// RK4 ////////////////////////////////////////////////////////////////////
   
   double dt=10,t=0,T=10000,aa=0,bb=0,cc=0,b=0,b1=0,b2=0,b3=0,rho_tot=0,q=0,Q,K,E,EE,p,tff,rt,rl1,rl2,rl3,rl4,counter1,counter2,nv,rhor1=0;
   double dvx=0,dvy=0,dvz=0,dx=0,dy=0,dz=0;
   int k,n_stelle=0,ppp=0,out=0,in=0;
   
   double Ax0,Bx0,C0,D0,Ax1,Bx1,C1,D1,Ax2,Bx2,C2,D2,Ax3,Bx3,C3,D3; 
   double Ay0,By0,Cy0,Dy0,Ay1,By1,Cy1,Dy1,Ay2,By2,Cy2,Dy2,Ay3,By3,Cy3,Dy3; 
   double Az0,Bz0,Cz0,Dz0,Az1,Bz1,Cz1,Dz1,Az2,Bz2,Cz2,Dz2,Az3,Bz3,Cz3,Dz3; 
   
   printf("%d\t%d\t%d\n\n\n",n/3,n/3,n/17);
   ms=0;
   for(i=0;i<n;i++){
     ms+=m[i];
   }
   p = 3*(ms+M)/(4*M_PI*R*R*R);
   tff = sqrt((3*M_PI)/(32*G*p));
   printf("\nT_ff = %lf\n\n",tff);
   ms=0;
   for(i=0;i<n;i++){
     for(j=0;j<n;j++){
       if(rr[i]<R && rr[j]<R){
	 if(i!=j){
	   
	   U += 0.5*G*m[i]*m[j]/(fabs(rr[i]-rr[j]));
	   
	 }
       }
     }
     if(rr[i]<R){
       K += m[i]*vv[i]*vv[i];
     }
   }
   
   fprintf(fq,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,Q,U,T,q,E);
   
   for(j=0;j<R;j++){
     contatore[j]=0;
   }
   
   do{
     
     for(i=0;i<n;i++){
       if(rr[i]<=R){    //
	 if(i!=0){
	   for(j=0;j<n;j++){
	     if(i!=j && rr[j]<=R){   //
	       
	       aa+=x[j]-x[i];
	       bb+=y[j]-y[i];
	       cc+=z[j]-z[i];
	       b+=rr[i]-rr[j];
	       ms+=m[j];
	       
	     }
	   }
	   
	   // Passo 0                                                       dvxx0[i] = Ax0 / dxx0[i] = Bx0 / vv0[i] = Cx0 / drr[i] = Dx0
	   
	   Ax0 = ms*G*(aa/(fabs(b*b*b)))*dt;
	   Bx0 = vx[i]*dt;
	   
	   Ay0 = ms*G*(bb/(fabs(b*b*b)))*dt;
	   By0 = vy[i]*dt;
	   
	   Az0 = ms*G*(cc/(fabs(b*b*b)))*dt;
	   Bz0 = vz[i]*dt;
	   
	   
	   C0 = vv[i];
	   D0 = C0*dt;
	   
	   // Passo 1
	   
	   b1=0;
	   
	   for(j=0;j<n;j++){
	     if(i!=j && rr[j]<=R){   //
	       
	       b1+=pow(fabs(rr[i]-rr[j]+D0*0.5),3);
	       
	     }
	   }
	   Ax1 = ms*G*((aa+Bx0*0.5*(n-1))/((b1)))*dt;
	   Bx1 = Bx0 + Ax0*0.5*dt;
	   
	   Ay1 = ms*G*((bb+By0*0.5*(n-1))/((b1)))*dt;
	   By1 = By0 + Ay0*0.5*dt;
	   
	   Az1 = ms*G*((cc+Bz0*0.5*(n-1))/((b1)))*dt;
	   Bz1 = Bz0 + Az0*0.5*dt;
	   
	   
	   C1 = sqrt((vx[i]+Ax0*0.5*dt)*(vx[i]+Ax0*0.5*dt) + (vy[i]+Ay0*0.5*dt)*(vy[i]+Ay0*0.5*dt) + (vz[i]+Az0*0.5*dt)*(vz[i]+Az0*0.5*dt));
	   D1 = C0*dt + D0*0.5*dt;
	   
	   // Passo 2
	   
	   b2=0;
	   
	   for(j=0;j<n;j++){
	     if(i!=j && rr[j]<=R){         //
	       
	       b2+=pow(fabs(rr[i]-rr[j]+D1*0.5),3);
	       
	     }
	   }
	   
	   Ax2 = ms*G*((aa+Bx1*0.5*(n-1))/((b2)))*dt;
	   Bx2 = Bx0 + Ax1*0.5*dt;
	   
	   Ay2 = ms*G*((bb+By1*0.5*(n-1))/((b2)))*dt;
	   By2 = By0 + Ay1*0.5*dt;
	   
	   Az2 = ms*G*((cc+Bz1*0.5*(n-1))/((b2)))*dt;
	   Bz2 = Bz0 + Az1*0.5*dt;
	   
	   C2 = sqrt((vx[i]+Ax1*0.5*dt)*(vx[i]+Ax1*0.5*dt) + (vy[i]+Ay1*0.5*dt)*(vy[i]+Ay1*0.5*dt) + (vz[i]+Az1*0.5*dt)*(vz[i]+Az1*0.5*dt));
	   D2 = C0*dt + D1*0.5*dt;   
	   
	   // Passo 3
	   
	   b3=0;
	   
	   for(j=0;j<n;j++){
	     if(i!=j && rr[j]<=R){
	       
	       b3+=pow(fabs(rr[i]-rr[j]+D2),3);
	       
	     }
	   }
	   
	   Ax3 = ms*G*((aa+Bx2*(n-1))/(b3))*dt;
	   Bx3 = Bx0 + Ax2*dt;
	   
	   Ay3 = ms*G*((bb+By2*(n-1))/(b3))*dt;
	   By3 = By0 + Ay2*dt;
	   
	   Az3 = ms*G*((cc+Bz2*(n-1))/(b3))*dt;
	   Bz3 = Bz0 + Az2*dt;
	   
	   
	   C3 = sqrt((vx[i]+Ax2*dt)*(vx[i]+Ax2*dt) + (vy[i]+Ay2*dt)*(vy[i]+Ay2*dt) + (vz[i]+Az2*dt)*(vz[i]+Az2*dt));
	   D3 = C0*dt + D2*dt;   
	   
	   
	   // Aggiornamento variabili
	   
	   
	   vx[i] += (Ax0+2*Ax1+2*Ax2+Ax3)/6;
	   vy[i] += (Ay0+2*Ay1+2*Ay2+Ay3)/6;
	   vz[i] += (Az0+2*Az1+2*Az2+Az3)/6;
	   
	   x[i] += (Bx0+2*Bx1+2*Bx2+Bx3)/6;
	   y[i] += (By0+2*By1+2*By2+By3)/6;
	   z[i] += (Bz0+2*Bz1+2*Bz2+Bz3)/6;
	   
	   rr[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
	   vv[i] = sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	 }
	 if(i==0 && t==T-dt){
	   fprintf(fn,"%lf\t%lf\t%lf\t%lf\t%lf\n",x[i],y[i],z[i],rr[i],t);
	 }
	 if(i<n/2 && i!=0 && t==T-dt){
	   rhor = 3*M_PI*m[i]/(rr[i]*rr[i]*rr[i]);
	   fprintf(fv,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x[i],y[i],z[i],rr[i],t,rhor);
	 }
	 if(i>=n/2 && t==T-dt){
	   rhor = 3*M_PI*m[i]/(rr[i]*rr[i]*rr[i]);
	   fprintf(fh,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x[i],y[i],z[i],rr[i],t,rhor);
	 }
	 if(i==1){
	   fprintf(fm,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x[i],y[i],z[i],rr[i],t,rhor);
	 }
	 /*if(t==T-dt){
	   fprintf(fg,"%lf\t%lf\t%lf\t%lf\t%lf\n",x[i],y[i],z[i],rr[i],t);
	   }*/
	 if(t==0){
	   
	   for(k=1;k<=R;k++){
	     
	     if(rr[i]<=k && rr[i] > k-1){
	       
	       n_stelle += 1; 
	       
	     }
	     if(rr[i] > R){
	       
	       ppp+=1;
	     }
	     //printf("%d\t%d\t%d\t%d\n",i,k,n_stelle,ppp);
	     contatore[i]=k;
	     n_stelle = 0;
	     ppp = 0;
	   }
	   
	   contatore[i] = k;    // Dice in che shell sta la stella
	   
	 }
	 
	 for(j=0;j<R;j++){
	   if(rr[i] < R-j+1 && rr[i] >= R-j){
	  contatore[j] += 1;                    // Raccoglie il numero di stelle in ogni shell
	   }
	 }
	 
	 
	 if(rr[i]>R){
	   out+=1;
	 } else {
	   in+=1;
	 }

	 
	 aa=0;
	 bb=0;
	 cc=0;
	 b=0;
	 ms=0;
       
     }
     }
     
     
     for(j=0;j<R;j++){
       if(j!=R-1){
	 rho[j] = -(4/3)*M_PI*(contatore[j])/(pow((R-j),3)) + (4/3)*M_PI*(contatore[j])/(pow((R-j-1),3));
       } else {
	 rho[j] = (4/3)*M_PI*(contatore[j])/(pow((R-j),3));
       }
       rho_tot += rho[j];
       //printf("%lf\n",rho[j]);
     }
     
     //fprintf(fq,"%lf\t%lf\n",t,rho_tot);
     for(j=0;j<R;j++){
       contatore[j] = 0;
     }
     printf("%lf\t%d\t%d\n",t,in,out);
     rho_tot=0;
     out=0;
     in=0;
     
     U=0;
     K=0;
     
     for(i=0;i<n;i++){
       
       for(j=0;j<n;j++){
	 if(rr[i]<R && rr[j]<R){
	   if(i!=j){
	     
	     U += 0.5*G*m[i]*m[j]/(fabs(rr[i]-rr[j]));
	     
	   }
	 }
       }
       if(rr[i]<R){
	 K += m[i]*vv[i]*vv[i];
       }
     }
     
     E = K + U;
     EE = K-U;
     Q = K/U;
     if(t!=0){
       q+=Q/t;
     }
     fprintf(fq,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,Q,U,T,q,E,EE);
     
     
     rt=0;
     rl1=0;
     rl2=0;
     rl3=0;
     rl4=0;
     nv=0;
     counter1=0;
     counter2=0;
     do{
       
       for(i=0;i<n;i++){
	 
	 if(rr[i]<=rt && i<n/2){
	   counter1 += 1;
	 }
	 if(rr[i]<=rt && i>=n/2){
	   counter2 += 1;
	 }
       }
       if(rl1==0){
	 if(counter1 >= n/10 && counter1 <= n/4){
	   rl1=rt;
	 }
       }
       if(rl2==0){
	 if(counter1 >= n/4){
	   rl2=rt;
	 }
       }
       if(rl3==0){
	 if(counter2 >= n/10 && counter2 <= n/4){
	   rl3=rt;
	 }
       }
       if(rl4==0){
	 if(counter2 >= n/4){
	   rl4=rt;
	 }
       }
       if(rl1!=0 && rl2!=0 && rl3!=0 && rl4!=0){
	 fprintf(ft,"%lf\t%lf\t%lf\t%lf\t%lf\n",t,rl1,rl2,rl3,rl4);
       }
       counter1 = 0;
       counter2 = 0;
       
       rt+=0.1;
       
     }while(rt<R);
     
     rhor=0;
     
     ms=0;
     for(i=0;i<n;i++){
       if(rr[i]<=R/5){
	 ms+=m[i];
       }
     }
     
     for(i=0;i<n;i++){
       if(rr[i]<=R/5){
	 rhor=3*ms/(4*M_PI*R*R*R/125);
       }
     }
     ms=0;
     for(i=0;i<n;i++){
       if(rr[i]<=R/2){
	 ms+=m[i];
       }
     }
     for(i=0;i<n;i++){
       if(rr[i]<=R/2){
	 rhor1=3*ms/(4*M_PI*R*R*R/8);
       }
     }
     fprintf(f14,"%lf\t%lf\t%lf\n",t,rhor,rhor1);
     rhor1=0;
     rhor=0;
     ms=0;
     
     t+=dt;
     
     
   } while(t<T);
   
}

