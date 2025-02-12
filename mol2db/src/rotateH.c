/*
c----------------------------------------------------------------------
c
c       Copyright (C) 2001 Binqing Wei and Brian K. Shoichet
c              Northwestern University Medical School
c                         All Rights Reserved.
c----------------------------------------------------------------------
c */
// rotateH(float* A,float* B,float* C,float* D,float ang)
// pass in: Cartesian coordinates of 4 points (A,B,C,D) and an angle (in degree)
// pass out: Cartesian coordinates of the 4 points, point D's coordinates is
// transformed by rotating about BC axis for the specified angle.
//  assuming: A, B, C are not colinear.
// based on algorithm used in Babel
// Binqing Wei , 08-23-00

#include <stdio.h>
#include <math.h>

void rotateH(float* A,float* B,float* C,float* D, float ang)
{
#define SQUARE(x) ((x)*(x))
#define PI 3.1415926535897932384626433f

  float v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z;
  float c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z;
  float c1mag,c2mag,radang,costheta,m[9];
  float x,y,z,mag,rotang,sn,cs,t,tx,ty,tz;

  //calculate the torsion angle

  v1x = A[0] - B[0]; v2x = B[0] - C[0];
  v1y = A[1] - B[1]; v2y = B[1] - C[1];
  v1z = A[2] - B[2]; v2z = B[2] - C[2];
  v3x = C[0] - D[0];
  v3y = C[1] - D[1];
  v3z = C[2] - D[2];

  c1x = v1y*v2z - v1z*v2y;   c2x = v2y*v3z - v2z*v3y;
  c1y = -v1x*v2z + v1z*v2x;  c2y = -v2x*v3z + v2z*v3x;
  c1z = v1x*v2y - v1y*v2x;   c2z = v2x*v3y - v2y*v3x;
  c3x = c1y*c2z - c1z*c2y;
  c3y = -c1x*c2z + c1z*c2x;
  c3z = c1x*c2y - c1y*c2x; 

  c1mag = SQUARE(c1x)+SQUARE(c1y)+SQUARE(c1z);
  c2mag = SQUARE(c2x)+SQUARE(c2y)+SQUARE(c2z);
  if (c1mag*c2mag < 0.01) costheta = 1.0; //avoid div by zero error
  else costheta = (c1x*c2x + c1y*c2y + c1z*c2z)/(sqrt(c1mag*c2mag));

  if (costheta < -0.999999) costheta = -0.999999f;
  if (costheta >  0.999999) costheta =  0.999999f;
                                  
  if ((v2x*c3x + v2y*c3y + v2z*c3z) > 0.0) radang = -acos(costheta);
  else                                     radang = acos(costheta);

        //
        // now we have the torsion angle (radang) - set up the rot matrix
        //

        //find the difference between current and requested
  rotang = ang*PI/180; 

  sn = sin(rotang); cs = cos(rotang);t = 1 - cs;
  //normalize the rotation vector
  mag = sqrt(SQUARE(v2x)+SQUARE(v2y)+SQUARE(v2z));
  x = v2x/mag; y = v2y/mag; z = v2z/mag;

  //set up the rotation matrix
  m[0]= t*x*x + cs;     m[1] = t*x*y + sn*z;  m[2] = t*x*z - sn*y;
  m[3] = t*x*y - sn*z;  m[4] = t*y*y + cs;    m[5] = t*y*z + sn*x;
  m[6] = t*x*z + sn*y;  m[7] = t*y*z - sn*x;  m[8] = t*z*z + cs;

  //
  //now the matrix is set - time to rotate the atoms
  //
  tx = B[0]; ty = B[1]; tz = B[2];
  D[0] -= tx; D[1] -= ty; D[2] -= tz;
  x = D[0]*m[0] + D[1]*m[1] + D[2]*m[2];
  y = D[0]*m[3] + D[1]*m[4] + D[2]*m[5];
  z = D[0]*m[6] + D[1]*m[7] + D[2]*m[8];
  D[0] = x; D[1] = y; D[2] = z;
  D[0] += tx; D[1] += ty; D[2] += tz;

  return;  
}
