C **** NOTE NPN1MAXVAL = 100 for Win32 and operational builds. Works for micron size particles. 
C ****      NPN1MAXVAL = 280 For non operational x64 build. Allows much bigger particles to be handled
C ****      The variable is It is defined in the Intel Fortran Pre-processor options, see the project Fortran options
      PARAMETER (NPN1MAXVAL=280)
      PARAMETER (NPN1=NPN1MAXVAL,
     &           NPNG1=4*NPN1, NPNG2=2*NPNG1, NPN2=2*NPN1,  
     &           NPL=NPN2+1, NPN3=NPN1+1,  
     &           NPN4=NPN1-25, NPN5=2*NPN4, NPN6=NPN4+1, NPL1=NPN5+1)


