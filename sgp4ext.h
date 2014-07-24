#ifndef _sgp4ext_
#define _sgp4ext_
/*     ----------------------------------------------------------------
*
*                                 sgp4ext.h
*
*    this file contains extra routines needed for the main test program for sgp4.
*    these routines are derived from the astro libraries.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2007
*                              by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              20 apr 07  david vallado
*                           misc documentation updates
*    changes :
*              14 aug 06  david vallado
*                           original baseline
*       ----------------------------------------------------------------      */
#pragma once

#include <string.h>
#include <math.h>

#include "ast2body.h"//DEBUG: 23/07/2014 Alek Lidtke: included ast2body to have access to newtonu
#include "sgp4unit.h"

#define undefined 999999.1
#define infinite 999999.9
// ------------------------- function declarations -------------------------

//DEBUG: 23/07/2014 Alek Lidtke: commented some functions out as they are already defined elsewhere.
//double  sgn
//        (
//          double x
//        );
//
//double  mag
//        (
//          double x[3]
//        );
//
//void    cross
//        (
//          double vec1[3], double vec2[3], double outvec[3]
//        );
//
//double  dot
//        (
//          double x[3], double y[3]
//        );
//
//double  angle
//        (
//          double vec1[3],
//          double vec2[3]
//        );
//
//void    newtonnu
//        (
//          double ecc, double nu,
//          double& e0, double& m
//        );
//
//double  asinh
//        (
//          double xval
//        );

void    rv2coe
        (
          double r[3], double v[3], double mu,
          double& p, double& a, double& ecc, double& incl, double& omega, double& argp,
          double& nu, double& m, double& arglat, double& truelon, double& lonper
        );

//void    jday
//        (
//          int year, int mon, int day, int hr, int minute, double sec,
//          double& jd
//        );
//
//void    days2mdhms
//        (
//          int year, double days,
//          int& mon, int& day, int& hr, int& minute, double& sec
//        );
//
//void    invjday
//        (
//          double jd,
//          int& year, int& mon, int& day,
//          int& hr, int& minute, double& sec
//        );

#endif

