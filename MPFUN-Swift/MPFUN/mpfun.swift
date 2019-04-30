//
//  MPFUN.swift
//  MPFUN-Swift
//
//  Created by Mike Griebling on 28 Apr 2019.
//  Copyright © 2019 Computer Inspirations. All rights reserved.
//

import Foundation

class mpfun {
    
    //----------------------------------------------------------------------------
    
    //   Integer constants:
    
    //   Name     Default   Description
    //   mpndpw       14    Largest n such that 10^n <= mpbdx (see below).
    //   mpldb         6    Logical device number for output of error messages.
    //   mpnbt        48    Number of significant bits in one mantissa word.
    //   mpnpr         7    Limit on summation loop to prevent overflow of exact value.
    //                        See usage in mpmul and mpdiv of module MPFUNB.
    //   mpnstr     2048    Maximum length of certain input character strings.
    //                        See usage in mpinp of module MPFUNC.
    //   mpoutl       80    Length of output lines.  See usage in mpout of MPFUNC.
    //   mprknd        8    Kind parameter for double precision.
    
    //   Double precision constants:
    
    //   Name     Default   Description
    //   mpbdx      2^48    2^mpnbt, the radix for MP numbers.
    //   mpbx2      2^96    Square of radix.
    //   mpb13x     2^13    Constant for checking the 40-bit restriction.
    //                        See usage in mpdmc40, mpmuld40 and mpdivd40 of MPFUNB.
    //   mpdpw  Log10(2^48) DP approximation to number of digits per mantissa word.
    //   mpexpmx   2^31/48  Largest permissible exponent, corresponding to a maximum
    //                        binary exponent of 2^31, or, in other words, a maximum
    //                        MP value of 2^(2^31) or approximately 10^646456993.
    //   mprdx   2^(-48)    Reciprocal of radix.
    //   mprx2   2^(-96)    Reciprocal of square of radix.
    //   mpb24x     2^24    Square root of radix.
    //   mpr24x  2^(-24)    Reciprocal of square root of radix.
    
//    integer, public:: mpndpw, mpldb, mpnbt, mpnpr, mpnstr, mpoutl, mprknd
//    parameter (mpndpw = 14, mpldb = 6, mpnbt = 48, mpnpr = 7, mpnstr = 2048, &
//    mpoutl = 80, mprknd = 8)
//    real (mprknd), public:: mpbdx, mpbx2, mpdpw, mpexpmx, mprdx, mprx2, mpb24x, &
//    mpr24x, mpb13x
//    parameter (mpbdx = 2.d0 ** mpnbt, mpbx2 = mpbdx**2, mpdpw = 14.449439791871d0, &
//    mpexpmx = 2.d0**31 / 48.d0, mprdx = 0.5d0 ** mpnbt, mprx2 = mprdx**2, &
//    mpb24x = 2.d0 ** 24, mpr24x = 0.5d0 ** 24, mpb13x = 2.d0 ** 13)
    
 
                
}
