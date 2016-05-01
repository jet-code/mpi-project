#include "m6378.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// CPU real-time instruction counter for timing

// #define constants
// XL, XR, UL, UR, UOFX, FOFX

float f(float x) {
//  d^2/dx^2 u(x) = f(x) ;
  return(FOFX);
 }

float exact(float x) {
//  u(x);
  return(UOFX);
 }

float slope(float y2, float y1, float x2, float x1) {
// compute slope m of line y = mx + b
  return((y2 - y1) / (x2 - x1));
 }

float yintercept(float y, float m, float x) {
// compute y-intercept b of line y = mx + b
  return(y - m*x);
 }

main( int argc, char *argv[] ) {

  int  np, myid, id;
  int I, nits;
  float xl, xl_local, xr, xr_local, ul, ur;
  float *uold = NULL, *unew = NULL; 
  float *utmp = NULL, *utotal = NULL;
  float x, xstar, urr;
  float R, e[2], mb[2];
  float dx, dx2, dt, odx2;
  int i, n;
  // timing vars
  int msec;
  clock_t start, diff;
  start = clock();

  M6378_Init( &argc, argv );
  M6378_Comm_size( &np );
  M6378_Comm_rank( &myid );

  if( argc != 3 ) {
    if( myid == 0 )
      fprintf( stderr, "Usage: %s local-grid-size number-iters\n", argv[0] );
    M6378_Finalize();
    exit(1);
  }
  
  I    = atoi(argv[1]);
  nits = atoi(argv[2]);
  
  xl = XL;
  xr = XR;
  ul = UL;
  ur = UR;
  
  dx = (xr - xl) / (np*I);
  dx2 = dx*dx;
  odx2 = 1.0/(dx2);
  dt = 0.25*dx2;
  
  xl_local = xl + (myid*I)*dx;
  xr_local = xl + (myid*I + I)*dx;

  uold = (float*)malloc( (I+1)*sizeof(float) );
  unew = (float*)malloc( (I+1)*sizeof(float) );
  

  if( !uold || !unew ) {
     fprintf( stderr, "Memory allocation error on id = %d.\n", myid );
     if( uold ) free( uold );
     if( unew ) free( unew );
     M6378_Finalize();
     exit(1);
  }
  
  if (myid == 0) {
     
     // Left boundary value
     uold[0] = ul;
  
     for( i = 1; i < I; i++ ) {
         uold[i] = 0.0;
      }
  
     // Right boundary value
     uold[I] = 0.0;
  }
   
  else if (myid > 0 && myid < (np-1)) {
  
     // Left boundary value
     uold[0] = 0.0;
  
     for( i = 1; i < I; i++ ) {
         uold[i] = 0.0;
      }
  
     // Right boundary value
     uold[I] = 0.0;
  }
     
  else if (myid == (np-1)) {
     
     // Left boundary value
     uold[0] = 0.0;
  
     for( i = 1; i < I; i++ ) {
         uold[i] = 0.0;
      }
  
     // Right boundary value
     uold[I] = ur;
  }
  
  //
  
  for(n = 0; n < nits; n++) {
      
      
      x = xl_local;
      
      for(i = 1; i < I; i++) {
          
          unew[i] = uold[i] + dt*((uold[i+1] - 2*uold[i] + uold[i-1])*odx2 - f(x));
          x += dx;
      }
      
      utmp = uold;
      uold = unew;
      unew = utmp;
      
  }
  
  // Share interface values between adjacent processes
  
  if (myid > 0) {
     
      if(!M6378_Wait(myid-1))
         goto ErrorExit;
        
      if(!M6378_Isend(myid-1, &uold[1], sizeof(float)))
         goto ErrorExit;
  }
  
  if (myid < (np-1)) {
     
      if(!M6378_Recv(myid+1, &urr, (int*)0))
         goto ErrorExit;
      
      // Compute residual R = D^2/dx^2 u(x) - f(x) 
      //  at the right interface of each process.
      
      xstar = xr_local;
      
      R = ((urr - 2*uold[I] + uold[I-1]) / dx2) - f(xstar);
      
      // Use D^2/dx^2 e(x) = -R*delta_at_interface_point
      //  to compute error at interface point
      e[1] = R*(dx*(xr - xstar)*(xstar - xl)) / (xr - xl);
      
      // slope and y-intercept for line on the left side of x_star
      mb[0] = slope(e[1], ul, xstar, xl);
      mb[1] = yintercept(ul, mb[0], xl);
      
      // error at left interface
      e[0] = mb[0]*xl_local + mb[1];
      
      
      printf("\n l1, myid = %d, e0 = %f , e1 = %f \n", myid, e[0], e[1]);
      
      // broadcast error correction to other processes
      
      for(id = 0; id < myid; id++) {
          
          if(!M6378_Wait(id))
            goto ErrorExit;
        
          if(!M6378_Isend(id, &mb[0], 2*sizeof(float)))
            goto ErrorExit;
       } 
      
      // slope and y-intercept for line on the right side of x_star
      mb[0] = slope(ur, e[1], xr, xstar);
      mb[1] = yintercept(ur, mb[0], xr);
      
      for(id = myid+1; id < np; id++) {
          
          if(!M6378_Wait(id))
            goto ErrorExit;
        
          if(!M6378_Isend(id, &mb[0], 2*sizeof(float)))
            goto ErrorExit;
       }
       
  } else {
        // Set process np-1 error to zero.
        e[0] = 0.0;
        e[1] = 0.0;
       
  }
   
  // Receive error corrections from other processes
  for(id = 0; id < np-1; id++) {
       
      if(id == myid) 
         continue;
       
      if(!M6378_Recv(id, &mb[0], (int*)0))
         goto ErrorExit;
       
       
      e[0] += mb[0]*(xl_local) + mb[1];
      e[1] += mb[0]*(xr_local) + mb[1];
       
      printf("\n l2, myid = %d, e0 = %f , e1 = %f \n", myid, e[0], e[1]);
  }
   
  // Piecewise linear interpolation of error: e^n(x) = m*x + b
  
    mb[0] = slope(e[1], e[0], xr_local, xl_local);
    mb[1] = yintercept(e[1], mb[0], xr_local);
  
    // u(x) = u^n(x) + e^n(x)
    
    x = xl_local;
    
    for(i = 0; i <= I; i++) {
        
        uold[i] = uold[i] + (mb[0]*x + mb[1]);
        x += dx;
    }

  // Send sub-domain results to master process
  if (myid > 0) {
     
    if(!M6378_Wait(0))
       goto ErrorExit;
       
    if(!M6378_Isend(0, &uold[1], I*sizeof(float)))
       goto ErrorExit;
   }
  
  // Collect final results
  if (myid == 0) {
     
     float x = xl;
     utotal = (float*)malloc(((np*I)+1)*sizeof(float));
     
     for (i = 1; i <= I; i++) {
         
         utotal[i] = uold[i];
      }
     
     for( id = 1; id < np; id++ ) {
         
         if( !M6378_Recv(id, &utotal[(id*I)+1], (int*)0))
            goto ErrorExit;
        }
     
     for (i = 0; i <= np*I; i++) {
        
        printf(" u_total_%d(%f) = %f,  e_%d = %f\n",i,x,utotal[i],i,(exact(x)-utotal[i]));
        x = x + dx;
      }
     
     }


  free( uold ); free( unew ); 
  
  if (myid == 0) {
     free( utotal );
     // Compute time
     diff = clock() - start;
     msec = diff * 1000 / CLOCKS_PER_SEC;
     printf("\nCompute time: %d s %d ms\n\n", msec/1000, msec%1000);
   }
  M6378_Finalize();
  exit(0);
 
ErrorExit:
  free( uold ); free( unew );
  
  if (myid == 0) {
     free( utotal );
   }
  M6378_Finalize();
  exit(1);
}
