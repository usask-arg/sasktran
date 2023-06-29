/*
   readeph0.c: Dummy version of 'readeph.c' for use with NOVAS-C
               Version 2.0

   U. S. Naval Observatory
   Astronomical Applications Dept.
   3450 Massachusetts Ave., NW
   Washington, DC  20392-5420

   The real version of 'readeph' provides a heliocentric state vector 
   (position and velocity) of a minor planet at a specified date.  
   This dummy version is for use with NOVAS-C when positions of the 
   minor planets are not of interest.
*/

/*
   Function prototype.
*/

   double *readeph (int mp, char *name, double jd, 

                    int *error );


/********readeph */

double *readeph (int mp, char *name, double jd,

                 int *error )
/*
------------------------------------------------------------------------

   PURPOSE:
      This is a dummy version of function 'readeph'.  It returns a
      pointer to an array of six doubles, all set to zero.

   REFERENCES:
      None.

   INPUT
   ARGUMENTS:
      mp (int)
         The number of the asteroid for which the position in desired.
      name (char*)
         The name of the asteroid.
      jd (double)
         The Julian date on which to find the postition and velocity.
    
   OUTPUT
   ARGUEMENTS:
      *error (int)
         Error code.

   RETURNED
   VALUE:
      (double *)
         A six element array giving first the position in AU and then
         the velocity of the asteroid in AU/day.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      None.

   VER./DATE/
   PROGRAMMER:
      V1.0/06-97/JAB (USNO/AA)
      V1.1/08-98/JAB (USNO/AA): Support new 'readeph' argument list.

   NOTES:
      None.

------------------------------------------------------------------------
*/
{
   static double pv[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

   *error = 0;
   return (pv);
}
