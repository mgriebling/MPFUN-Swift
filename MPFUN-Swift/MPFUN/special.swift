//
//  sfuncs.swift - Special functions
//  MPFUN-Swift
//
//  Created by Mike Griebling on 6 May 2019.
//  Copyright © 2019 Computer Inspirations. All rights reserved.
//

import Foundation

extension MPFUN {
    
    static func mpberner (_ nb1 : Int, _ nb2 : Int, berne : inout Array<MPRNumber>, _ mpnw : Int) {
        
        //  This returns the even Bernouli numbers B(2*k), from B(2) = 1/6 up to
        //  B(2*nb2).  The array berne must be dimensioned as shown below.
        
        var ia, na, mpnw1  : Int
        // var berne(0:nb1+5,nb2),
        var t1 = MPRNumber(repeating: 0, count: mpnw+7)
        var t2 = t1; var t3 = t1; var t4 = t1; var t5 = t1
        
        //  End of declaration
        
        if mpnw < 4 || berne[0][1] < Double(mpnw + 4) || berne[0][nb2] < Double(mpnw + 4) {
            print("*** MPBERN: uninitialized or inadequately sized arrays")
            mpabrt (62)
        }
        
        mpnw1 = mpnw + 1
        t1[0] = Double(mpnw + 7)
        t2[0] = Double(mpnw + 7)
        t3[0] = Double(mpnw + 7)
        t4[0] = Double(mpnw + 7)
        t5[0] = Double(mpnw + 7)
        mpmuld (mppicon, 2.0, &t1, mpnw1)
        mpmul (t1, t1, &t2, mpnw1)
        mpdmc (-2.0, 0, &t1, mpnw1)
        
        for k in 1...nb2 {
            mpmuld (t1, Double(2*k - 1), &t3, mpnw1)
            mpmuld (t3, Double(2*k), &t4, mpnw1)
            mpdiv (t4, t2, &t1, mpnw1)
            t1[2] = -t1[2]
            mpdmc (2.0 * Double(k), 0, &t3, mpnw1)
            mpzetar (t3, t4, mpnw1)
            mpmul (t1, t4, &t5, mpnw1)
            mproun (&t5, mpnw)
            
            //   The next few lines (to !+) are necessary, rather than a simple to
            //   mpeq, to avoid a Fortran rank-mismatch error.
            
            //  mpeq (t5, berne(0,k), mpnw)
            
            ia = sign (1.0, t5[2])
            na = min (Int (abs (t5[2])), mpnw)
            berne[1][k] = Double(mpnw)
            berne[2][k] = Double(sign (Double(na), Double(ia)))
            
            for i in 2...na+2 {
                berne[i+1][k] = t5[i+1]
            }
            
            berne[na+4][k] = 0.0
            berne[na+5][k] = 0.0
        }
        
    } // mpberner
    
    static func mpbesseljr (_ anu : MPRNumber, _ t : MPRNumber, _ z : inout MPRNumber, _ mpnw : Int) {
        
        //   This evaluates the function BesselJ (ANU, T).  ANU must be nonnegative and
        //   not greater than 10^6 (this limit can be adjusted below).  To compensate
        //   for an unsually large amount of internal cancelation in these formulas, all
        //   computations are performed to 3*mpnw/2 words precision.
        
        //   In the parameter statement below:
        //     itrmx = limit of number of iterations in series; default = 100000.
        //     dasy = factor used to decide if asymptic series is used; default = 25.
        //     anumx = upper limit of anu argument; default = 1000.
        
        var mpnw1, ndp, nu, n1 : Int
        var t0 = MPRNumber(repeating:0, count: 3*mpnw/2+6)
        var t1 = t0; var t2 = t0; var t3 = t0
        var t4 = t0; var t5 = t0; var t6 = t0
        let itrmx = 100000; let dasy = 25.0; let anumx = 1.0e6
        
        // End of declaration
        
        if mpnw < 4 || anu[0] < Double(mpnw+4) || anu[0] < abs(anu[2]) + 4 || t[0] < Double(mpnw+4) || t[0] < abs(t[2]) + 4 || z[0] < Double(mpnw+6) {
            print ("*** MPBESSELJR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        if (anu[2] < 0.0 || anu[3] > 0.0 || (anu[3] == 0.0 && anu[4] > anumx)) {
            print ("*** MPBESSELJR: First argument must be >= 0 and <= \(anumx)")
            mpabrt (65)
        }
        
        mpnw1 = 3 * mpnw / 2
        t0[0] = Double(mpnw1 + 6)
        t1[0] = Double(mpnw1 + 6)
        t2[0] = Double(mpnw1 + 6)
        t3[0] = Double(mpnw1 + 6)
        t4[0] = Double(mpnw1 + 6)
        t5[0] = Double(mpnw1 + 6)
        t6[0] = Double(mpnw1 + 6)
        
        //   Select either the direct or the asymptotic series.
        
        if t[3] < 0.0 || t[3] == 0.0 && t[4] < dasy * Double(mpnw - 2) {
            t2[1] = Double(mpnw1)
            t2[2] = 1.0
            t2[3] = 0.0
            t2[4] = 1.0
            t2[5] = 0.0
            t2[6] = 0.0
            mpadd (anu, t2, &t0, mpnw1)
            mpgammar (t0, t1, mpnw1)
            mpdiv (t2, t1, &t3, mpnw1)
            mpeq (t3, &t1, mpnw1)
            mpeq (t1, &t0, mpnw1)
            mpmul (t, t, &t3, mpnw1)
            mpmuld (t3, 0.25, &t2, mpnw1)
            
            var flag = false
            for i in 1...itrmx {
                mpmul (t1, t2, &t3, mpnw1)
                mpdivd (t3, Double(i), &t4, mpnw1)
                mpdmc (Double(i), 0, &t5, mpnw1)
                mpadd (anu, t5, &t6, mpnw1)
                mpdiv (t4, t6, &t1, mpnw1)
                t1[2] = -t1[2]
                mpadd (t0, t1, &t3, mpnw1)
                mpeq (t3, &t0, mpnw1)
                if t1[2] == 0.0 || t1[3] < t0[3] - Double(mpnw1) { flag = true; break; /* goto 100 */ }
            }
            
            if !flag {
                print ("*** MPBESSELJR: loop overflow 1")
                mpabrt (66)
            }
            
            // 100 continue
            
            mpmuld (t, 0.5, &t1, mpnw1)
            mppower (t1, anu, t2, mpnw1)
            mpmul (t0, t2, &t3, mpnw1)
            mpeq (t3, &t0, mpnw1)
        } else {
            t0[1] = Double(mpnw1)
            t0[2] = 1.0
            t0[3] = 0.0
            t0[4] = 1.0
            t0[5] = 0.0
            t0[6] = 0.0
            t1[1] = Double(mpnw1)
            t1[2] = 0.0
            t1[3] = 0.0
            t1[4] = 0.0
            t2[1] = Double(mpnw1)
            t2[2] = 1.0
            t2[3] = 0.0
            t2[4] = 1.0
            t2[5] = 0.0
            t2[6] = 0.0
            mpmul (anu, anu, &t3, mpnw1)
            mpmuld (t3, 4.0, &t5, mpnw1)
            
            var flag = false
            for i in 1...itrmx {
                mpdmc (Double(2*i - 1), 0, &t4, mpnw1)
                mpmul (t4, t4, &t6, mpnw1)
                mpsub (t5, t6, &t4, mpnw1)
                mpmul (t2, t4, &t3, mpnw1)
                mpdivd (t3, 8.0 * Double(i), &t4, mpnw1)
                mpdiv (t4, t, &t2, mpnw1)
                if i % 2 == 0 {
                    mpeq (t2, &t3, mpnw1)
                    if i % 4 == 2 { t3[2] = -t3[2] }
                    mpadd (t0, t3, &t4, mpnw1)
                    mpeq (t4, &t0, mpnw1)
                } else {
                    mpeq (t2, &t3, mpnw1)
                    if i % 4 == 3 { t3[2] = -t3[2] }
                    mpadd (t1, t3, &t4, mpnw1)
                    mpeq (t4, &t1, mpnw1)
                }
                if t2[2] == 0.0 || t2[3] < t0[3] - Double(mpnw1) && t2[3] < t1[3] - Double(mpnw1) {
                    flag = true
                    break
                    // goto 110
                }
            }
            
            if !flag {
                print ("*** MPBESSELJR: loop overflow 2")
                mpabrt (66)
            }
            // 110 continue
            
            mpeq (mppicon,  &t2, mpnw1)
            mpmul (t2, anu, &t4, mpnw1)
            mpmuld (t4, 0.5, &t3, mpnw1)
            mpsub (t, t3, &t4, mpnw1)
            mpmuld (t2, 0.25, &t3, mpnw1)
            mpsub (t4, t3, &t5, mpnw1)
            mpcssnr (t5, t3, &t4, mpnw1)
            mpmul (t3, t0, &t5, mpnw1)
            mpmul (t4, t1, &t6, mpnw1)
            mpsub (t5, t6, &t3, mpnw1)
            t4[1] = Double(mpnw1)
            t4[2] = 1.0
            t4[3] = 0.0
            t4[4] = 2.0
            t4[5] = 0.0
            t4[6] = 0.0
            mpmul (t2, t, &t5, mpnw1)
            mpdiv (t4, t5, &t6, mpnw1)
            mpsqrt (t6, &t4, mpnw1)
            mpmul (t4, t3, &t0, mpnw1)
        }
        
        mproun (&t0, mpnw)
        mpeq (t0, &z, mpnw)
        
    } //  mpbesseljr

    static func mpgammar (_ t : MPRNumber, _ z: inout MPRNumber, _ mpnw : Int) {
        
        //   This evaluates the gamma function, using an algorithm of R. W. Potter.
        //   The argument t must not exceed 10^8 in size (this limit is set below),
        //   must not be zero, and if negative must not be integer.
        
        //   In the parameter statement below:
        //     itrmx = limit of number of iterations in series; default = 100000.
        //     con1 = 1/2 * log (10) to DP accuracy.
        //     dmax = maximum size of input argument.
        
        var j, k, mpnw1, ndp, neps, nt, n1, n2, n3 : Int
        var alpha, d1, d2, d3 : Double
        let al2 = 0.69314718055994530942; let dmax = 1.0e8; let itrmx = 100000
        var f1 = MPRNumber(repeating:0, count:9)
        var sum1 = MPRNumber(repeating:0, count:mpnw+7)
        var sum2 = sum1; var tn = sum1; var t1 = sum1; var t2 = sum1; var t3 = sum1
        var t4 = sum1; var t5 = sum1; var t6 = sum1
        
        // End of declaration
        
        if mpnw < 4 || t[0] < Double(mpnw + 4) || t[0] < abs (t[2]) + 4 || z[0] < Double(mpnw + 6) {
            print ("*** MPGAMMAR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        if (t[2] == 0 || t[3] > 0 || (t[3] == 0 && t[4] > dmax) || (t[2] < 0.0 && t[3] == 0.0 && abs(t[2]) == 1.0)) {
            print("*** MPGAMMAR: input argument must have absolute value <= \(dmax)",
                "must not be zero, and if negative must not be an integer.")
            mpabrt (65)
        }
        
        mpnw1 = mpnw + 1
        f1[0] = 9.0
        f1[1] = Double(mpnw1)
        f1[2] = 1.0
        f1[3] = 0.0
        f1[4] = 1.0
        f1[5] = 0.0
        f1[6] = 0.0
        sum1[0] = Double(mpnw + 7)
        sum2[0] = Double(mpnw + 7)
        tn[0] = Double(mpnw + 7)
        t1[0] = Double(mpnw + 7)
        t2[0] = Double(mpnw + 7)
        t3[0] = Double(mpnw + 7)
        t4[0] = Double(mpnw + 7)
        t5[0] = Double(mpnw + 7)
        t6[0] = Double(mpnw + 7)
        
        //   Find the integer and fractional parts of t.
        
        mpinfr (t, &t2, &t3, mpnw1)
        
        if t3[2] == 0.0 {
            
            //   If t is a positive integer, then apply the usual factorial recursion.
            
            mpmdc (t2, &d2, &n2, mpnw1)
            nt = Int(d2 * pow(2, Double(n2)))
            mpeq (f1, &t1, mpnw1)
            
            for i in 2...nt - 1 {
                mpmuld (t1, Double(i), &t2, mpnw1)
                mpeq (t2, &t1, mpnw1)
            }
            
            mproun (&t1, mpnw)
            mpeq (t1, &z, mpnw)
            return // goto 120
        } else if t[2] > 0.0 {
            
            //   Apply the identity Gamma[t+1] = t * Gamma[t] to reduce the input argument
            //   to the unit interval.
            
            mpmdc (t2, &d2, &n2, mpnw1)
            nt = Int(d2 * pow(2, Double(n2)))
            mpeq (f1, &t1, mpnw1)
            mpeq (t3, &tn, mpnw1)
            
            for i in 1...nt {
                mpdmc (Double(i), 0, &t4, mpnw1)
                mpsub (t, t4, &t5, mpnw1)
                mpmul (t1, t5, &t6, mpnw1)
                mpeq (t6, &t1, mpnw1)
            }
        } else {
            
            //   Apply the gamma identity to reduce a negative argument to the unit interval.
            
            mpsub (f1, t, &t4, mpnw1)
            mpinfr (t4, &t3, &t5, mpnw1)
            mpmdc (t3, &d3, &n3, mpnw1)
            nt = Int(d3 * pow(2, Double(n3)))
            
            mpeq (f1, &t1, mpnw1)
            mpsub (f1, t5, &t2, mpnw1)
            mpeq (t2, &tn, mpnw1)
            
            for i in 0...nt - 1 {
                //    t1 = t1 / (t + Double(i))
                mpdmc (Double(i), 0, &t4, mpnw1)
                mpadd (t, t4, &t5, mpnw1)
                mpdiv (t1, t5, &t6, mpnw1)
                mpeq (t6, &t1, mpnw1)
            }
        }
        
        //   Calculate alpha = bits of precision * log(2) / 2, then take the nearest integer
        //   value, so that d2 = 0.25 * alpha^2 can be calculated exactly in DP.
        
        alpha = aint (0.5 * Double(mpnbt) * al2 * Double(mpnw1 + 1))
        d2 = 0.25 * alpha*alpha
        
        mpeq (tn, &t2, mpnw1)
        mpdiv (f1, t2, &t3, mpnw1)
        mpeq (t3, &sum1, mpnw1)
        
        //   Evaluate the series with t.
        
        var flag = false
        for j in 1...itrmx {
            mpdmc (Double(j), 0, &t6, mpnw1)
            mpadd (t2, t6, &t4, mpnw1)
            mpmuld (t4, Double(j), &t5, mpnw1)
            mpdiv (t3, t5, &t6, mpnw1)
            mpmuld (t6, d2, &t3, mpnw1)
            mpadd (sum1, t3, &t4, mpnw1)
            mpeq (t4, &sum1, mpnw1)
            if t3[2] == 0.0 || t3[3] < sum1[3] - Double(mpnw1) { break; flag = true /* goto 100 */ }
        }
        
        if !flag {
            print ("*** MPGAMMAR: iteration limit execeeded \(itrmx)")
            mpabrt (67)
        }
        
        // 100 continue
        
        mpeq (tn, &t2, mpnw1)
        t2[2] = -t2[2]
        mpdiv (f1, t2, &t3, mpnw1)
        mpeq (t3, &sum2, mpnw1)
        
        //   Evaluate the same series with -t.
        
        flag = false
        for j in 1...itrmx {
            mpdmc (Double(j), 0, &t6, mpnw1)
            mpadd (t2, t6, &t4, mpnw1)
            mpmuld (t4, Double(j), &t5, mpnw1)
            mpdiv (t3, t5, &t6, mpnw1)
            mpmuld (t6, d2, &t3, mpnw1)
            mpadd (sum2, t3, &t4, mpnw1)
            mpeq (t4, &sum2, mpnw1)
            if t3[2] == 0.0 || t3[3] < sum2[3] - Double(mpnw1) {  break; flag = true /* goto 110 */ }
        }
        
        if !flag {
            print ("*** MPGAMMAR: iteration limit execeeded \(itrmx)")
            mpabrt (67)
        }
        
        // 110 continue
        
        //   Compute sqrt (mppic * sum1 / (tn * sin (mppic * tn) * sum2))
        //   and (alpha/2)^tn terms.
        
        mpeq (mppicon, &t2, mpnw1)
        mpmul (t2, tn, &t3, mpnw1)
        mpcssnr (t3, t4, t5, mpnw1)
        mpmul (t5, sum2, &t6, mpnw1)
        mpmul (tn, t6, &t5, mpnw1)
        mpmul (t2, sum1, &t3, mpnw1)
        mpdiv (t3, t5, &t6, mpnw1)
        t6[2] = -t6[2]
        mpsqrt (t6, &t2, mpnw1)
        
        mpdmc (0.5 * alpha, 0, &t3, mpnw1)
        mplog (t3, t4, mpnw1)
        mpmul (tn, t4, &t5, mpnw1)
        mpexp (t5, t6, mpnw1)
        mpmul (t2, t6, &t3, mpnw1)
        
        mpmul (t1, t3, &t4, mpnw1)
        
        //   Round to mpnw words precision.
        
        mproun (&t4, mpnw)
        mpeq (t4, &z, mpnw)
        
        // 120 continue
        
    } // mpgammar

    static func mpincgammar (_ s : MPRNumber, _ z : MPRNumber, _ g : inout MPRNumber, _ mpnw : Int) {
        
        //  This returns the incomplete gamma function, using a combination of formula
        //  8.7.3 of the DLMF (for modest-sized z) and formula 8.11.2 (for large z).
        
        var mpnw1, n : Int
        let dmax = 40.0; let itrmax = 1000000
        var f1 = MPRNumber(repeating:0, count:9)
        var t0 = MPRNumber(repeating:0, count:mpnw+7)
        var t1 = t0; var t2 = t0; var t3 = t0
        var t4 = t0; var t5 = t0; var t6 = t0
        
        // End of declaration
        
        if mpnw < 4 || s[0] < Double(mpnw+4) || s[0] < abs(s[2]) + 4 || z[0] < Double(mpnw+4) || z[0] < abs(z[2]) + 4 || g[0] < Double(mpnw+6) {
            print ("*** MPINCGAMMAR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        mpnw1 = mpnw + 1
        t0[0] = Double(mpnw + 7)
        t1[0] = Double(mpnw + 7)
        t2[0] = Double(mpnw + 7)
        t3[0] = Double(mpnw + 7)
        t4[0] = Double(mpnw + 7)
        t5[0] = Double(mpnw + 7)
        t6[0] = Double(mpnw + 7)
        f1[0] = 9.0
        f1[1] = Double(mpnw1)
        f1[2] = 1.0
        f1[3] = 0.0
        f1[4] = 1.0
        f1[5] = 0.0
        f1[6] = 0.0
        
        // if (abs (z) < dmax * mpnw) {
        
        if z[3] < 0.0 || (z[3] == 0.0 && z[4] < dmax * Double(mpnw)) {
            
            //  t1 = gamma (s)
            
            mpgammar (s, &t1, mpnw1)
            
            //  t2 = 1.0 / (s * t1)
            
            mpmul (s, t1, &t3, mpnw1)
            mpdiv (f1, t3, &t2, mpnw1)
            
            //   t0 = t2
            
            mpeq (t2, &t0, mpnw1)
            
            var flag = false
            for k in 1...itrmax {
                
                //    t2 = t2 * z / (s + Double(k))
                
                mpmul (t2, z, &t5, mpnw1)
                mpdmc (Double(k), 0, &t3, mpnw1)
                mpadd (s, t3, &t4, mpnw1)
                mpdiv (t5, t4, &t2, mpnw1)
                
                //    t0 = t0 + t2
                
                mpadd (t0, t2, &t3, mpnw1)
                mpeq (t3, &t0, mpnw1)
                
                if t2[2] == 0.0 || t2[3] < t0[3] - Double(mpnw)  { flag = true; break /* goto 100 */ }
            }
            
            if !flag {
                print ("*** MPINCGAMMAR: iteration limit exceeded: \(itrmax)")
                mpabrt (101)
            }
            
            // 100 continue
            
            //   gammainc = t1 * (1.0 - z ** s / exp (z) * t0)
            
            mppower (z, s, &t2, mpnw1)
            mpexp (z, t3, mpnw1)
            mpdiv (t2, t3, &t4, mpnw1)
            mpmul (t4, t0, &t5, mpnw1)
            mpsub (f1, t5, &t2, mpnw1)
            mpmul (t1, t2, &t3, mpnw1)
            mpeq (t3, &t1, mpnw1)
        } else {
            //  t0 = mpreal (1.0, mpnw)
            
            t0[2] = 1.0
            t0[3] = 0.0
            t0[4] = 1.0
            t0[5] = 0.0
            t0[6] = 0.0
            
            //  t1 = mpreal (1.0, mpnw)
            
            t1[2] = 1.0
            t1[3] = 0.0
            t1[4] = 1.0
            t1[5] = 0.0
            t1[6] = 0.0
            
            var flag = false
            for k in 1...itrmax {
                //    t1 = t1 * (s - Double(k)) / z
                
                mpdmc (Double(k), 0, &t2, mpnw1)
                mpsub (s, t2, &t3, mpnw1)
                mpmul (t1, t3, &t4, mpnw1)
                mpdiv (t4, z, &t1, mpnw1)
                
                //    t0 = t0 + t1
                
                mpadd (t0, t1, &t2, mpnw1)
                mpeq (t2, &t0, mpnw1)
                
                if t1[2] == 0.0 || t1[3] < t0[3] - Double(mpnw) { flag = true; break /* goto 110 */ }
            }
            
            if !flag {
                print ("*** MPINCGAMMAR: iteration limit exceeded: \(itrmax)")
                mpabrt (101)
            }
            
            // 110 continue
            
            //  gammainc = z ** (s - 1.0) / exp (z) * t0
            
            mpsub (s, f1, &t2, mpnw1)
            mppower (z, t2, &t3, mpnw1)
            mpexp (z, t4, mpnw1)
            mpdiv (t3, t4, &t2, mpnw1)
            mpmul (t2, t0, &t1, mpnw1)
        }
        
        // 200 continue
        
        mproun (&t1, mpnw)
        mpeq (t1, &g, mpnw)
        
    } // mpincgammar

    static func mpzetar (_ ss : MPRNumber, _ zz : inout MPRNumber, _ mpnw : Int) {
        
        //   This returns the zeta function at positive real argument SS using an algorithm
        //   due to Peter Borwein.
        
        var iss, j, mpnw1, n, nwds : Int
        var d1 : Double
        let itrmax = 1000000; let dfrac = 16.0; let dlogb = 33.27106466687737
        var t1 = MPRNumber(repeating: 0, count:mpnw+7); var t2 = t1; var t3 = t1; var t4 = t1
        var t5 = t1; var tn = t1; var tt = t1; var s = t1
        var f1 = MPRNumber(repeating: 0, count:9)
        var sgn : Double
        
        //  End of declaration
        
        if mpnw < 4 || ss[0] < Double(mpnw+4) || ss[0] < abs(ss[2]) + 4 || zz[0] < Double(mpnw+6) {
            print ("*** MPZETAR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   Check if argument is 1 -- undefined.
        
        if ss[2] == 1.0 && ss[3] == 0.0 && ss[4] == 1.0 {
            print ("*** MPZETAR: argument is 1")
            mpabrt (63)
        }
        
        mpnw1 = mpnw + 1
        s[0] = Double(mpnw + 7)
        t1[0] = Double(mpnw + 7)
        t2[0] = Double(mpnw + 7)
        t3[0] = Double(mpnw + 7)
        t4[0] = Double(mpnw + 7)
        t5[0] = Double(mpnw + 7)
        tn[0] = Double(mpnw + 7)
        tt[0] = Double(mpnw + 7)
        
        //   Set f1 = 1.
        
        f1[0] = 9.0
        f1[1] = Double(mpnw1)
        f1[2] = 1.0
        f1[3] = 0.0
        f1[4] = 1.0
        f1[5] = 0.0
        f1[6] = 0.0
        
        //   Check if argument is zero.  If so, the result is -1/2.
        
        if (ss[2] == 0.0) {
            mpdmc (-0.5, 0, &t1, mpnw1)
            mproun (&t1, mpnw)
            mpeq (t1, &zz, mpnw)
            return
            //goto 200
        }
        
        //   Check if argument is negative.
        
        if (ss[2] < 0.0) {
            
            //   Check if argument is a negative even integer.  If so, the result is zero.
            
            mpmuld (ss, 0.5, &t1, mpnw1)
            mpinfr (t1, &t2, &t3, mpnw1)
            if (t3[2] == 0.0) {
                t1[1] = Double(mpnw1)
                t1[2] = 0.0
                t1[3] = 0.0
                t1[4] = 0.0
                mproun (&t1, mpnw)
                mpeq (t1, &zz, mpnw)
                return
                //goto 200
            }
            
            //   Otherwise compute zeta(1-ss), and later apply the reflection formula.
            
            mpsub (f1, ss, &tt, mpnw1)
        } else {
            mpeq (ss, &tt, mpnw1)
        }
        
        //  Check if argument is large enough that computing with definition is faster.
        
        // if (tt .gt. mpreald (dlogb * mpnw / log (32.0 * mpnw), mpnw)) {
        
        d1 = dlogb * Double(mpnw) / log (32.0 * Double(mpnw))
        if tt[2] >= 1.0 && (tt[3] > 1.0 || tt[3] == 0.0 && tt[4] > d1) {
            
            //  t1 = mpreal (1.0, mpnw)
            
            t1[1] = Double(mpnw1)
            t1[2] = 1.0
            t1[3] = 0.0
            t1[4] = 1.0
            t1[5] = 0.0
            t1[6] = 0.0
            
            for i in 2...itrmax {
                
                //    t2 = mpreal (Double(i), mpnw) ** tt
                
                mpdmc (Double(i), 0, &t4, mpnw1)
                mppower (t4, tt, t2, mpnw1)
                
                //    t3 = 1.0 / t2
                
                mpdiv (f1, t2, &t3, mpnw1)
                
                //    t1 = t1 + t3
                
                mpadd (t1, t3, &t2, mpnw1)
                mpeq (t2, &t1, mpnw1)
                
                if t3[2] == 0.0 || t3[3] < Double(-mpnw) {
                    mproun (&t1, mpnw)
                    mpeq (t1, &zz, mpnw)
                    return
                    // goto 200
                }
            }
            
            print ("*** MPZETAR: iteration limit exceeded \(itrmax)")
            mpabrt (101)
        }
        
        n = Int(dfrac * Double(mpnw))
        
        // tn = mpreal (2.0, mpnw) ** n
        
        tn[0] = Double(mpnw + 7)
        mpdmc (1.0, n, &tn, mpnw1)
        
        // t1 = - tn
        
        mpeq (tn, &t1, mpnw1)
        t1[2] = -t1[2]
        
        // t2 = mpreal (0.0, mpnw)
        
        t2[2] = 0.0
        t2[3] = 0.0
        t2[4] = 0.0
        
        // s = mpreal (0.0, mpnw)
        
        s[1] = Double(mpnw1)
        s[2] = 0.0
        s[3] = 0.0
        s[4] = 0.0
        
        sgn = 1.0
        
        for j in 0...2 * n - 1 {
            //  t3 = mpreal (Double(j + 1), mpnw) ** tt
            
            mpdmc (Double(j + 1), 0, &t4, mpnw1)
            mppower (t4, tt, t3, mpnw1)
            
            //  s = s + sgn * t1 / t3
            
            mpdiv (t1, t3, &t4, mpnw1)
            if (sgn < 0.0) { t4[2] = -t4[2] }
            mpadd (s, t4, &t5, mpnw1)
            mpeq (t5, &s, mpnw1)
            
            sgn = -sgn
            
            if (j < n - 1) {
                //    t2 = mpreal (0.0, mpnw)
                
                t2[2] = 0.0
                t2[3] = 0.0
                t2[4] = 0.0
                
            } else if (j == n - 1) {
                //    t2 = mpreal (1.0, mpnw)
                
                t2[2] = 1.0
                t2[3] = 0.0
                t2[4] = 1.0
                t2[5] = 0.0
                t2[6] = 0.0
                
            } else {
                //     t2 = t2 * Double(2 * n - j) / Double(j + 1 - n)
                
                mpmuld (t2, Double(2 * n - j), &t3, mpnw1)
                mpdivd (t3, Double(j + 1 - n), &t2, mpnw1)
                
            }
            //  t1 = t1 + t2
            
            mpadd (t1, t2, &t3, mpnw1)
            mpeq (t3, &t1, mpnw1)
        }
        
        // t1 = - s / (tn * (1.0 - mpreal (2.0, mpnw) ** (1.0 - tt)))
        
        mpsub (f1, tt, &t3, mpnw1)
        t2[2] = 1.0
        t2[3] = 0.0
        t2[4] = 2.0
        t2[5] = 0.0
        t2[6] = 0.0
        mppower (t2, t3, &t4, mpnw1)
        mpsub (f1, t4, &t2, mpnw1)
        mpmul (tn, t2, &t3, mpnw1)
        mpdiv (s, t3, &t1, mpnw1)
        t1[2] = -t1[2]
        
        //   If original argument was negative, apply the reflection formula.
        
        if (ss[2] < 0.0) {
            mpgammar (tt, &t3, mpnw1)
            mpmul (t1, t3, &t2, mpnw1)
            mpmul (mppicon, tt, &t1, mpnw1)
            mpmuld (t1, 0.5, &t3, mpnw1)
            mpcssnr (t3, t4, &t5, mpnw1)
            mpmul (t2, t4, &t1, mpnw1)
            mpmuld (mppicon, 2.0, &t2, mpnw1)
            mppower (t2, tt, &t3, mpnw1)
            mpdiv (t1, t3, &t2, mpnw1)
            mpmuld (t2, 2.0, &t1, mpnw1)
        }
        
        // 200 continue
        
        // zetapbr = t1
        
        mproun (&t1, mpnw)
        mpeq (t1, &zz, mpnw)
    } // mpzetar

    static func mpzetaemr (_ nb1: Int, _ nb2: Int, _ berne: [MPRNumber], s: MPRNumber, _ z: inout MPRNumber, _ mpnw: Int) {
        
        //  This evaluates the Riemann zeta function, using the combination of
        //  the definition formula (for large s), and an Euler-Maclaurin scheme
        //  (see formula 25.2.9 of the DLMF.  The array berne contains precomputed
        //  Bernoulli numbers.  Its dimensions must be as shown below.
        
        var i, ia, iss, k, mpnw1, na, nn, n1, n2 : Int
        var d1, d2 : Double
        let itrmax = 1000000; let dfrac = 8.5; let dlogb = 33.27106466687737
        var t1 = MPRNumber(repeating: 0, count:mpnw+7); var t2 = t1; var t3 = t1; var t4 = t1
        var t5 = t1; var t6 = t1; var tt = t1; var s = t1; var t0 = t1; var t7 = t1; var t8 = t1
        var t9 = t1
        var f1 = MPRNumber(repeating: 0, count:9)
        
        // End of declaration
        
        if mpnw < 4 || s[0] < Double(mpnw + 4) || s[0] < abs (s[2]) + 4 || z[0] < Double(mpnw + 6) {
            print ("*** MPZETAEMR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   Check if argument is 1 -- undefined.
        
        if (s[2] == 1.0 && s[3] == 0.0 && s[4] == 1.0) {
            print ("*** MPZETAEMR: argument is 1")
            mpabrt (63)
        }
        
        //   Check if berne array has been initialized.
        
        if Int(berne[0][1]) < mpnw + 4 || berne[0][1] < abs(berne[2][1]) + 4 ||
            Int(berne[0][nb2]) < mpnw + 4 || berne[0][nb2] < abs(berne[2][nb2]) + 4 ||
            nb2 < Int(mpndpw * mpnw) {
            print ("*** MPZETAEMR: Bernoulli coefficient array must be initialized",
                   "with at least \(Int (mpndpw * mpnw)) entries.")
            mpabrt (62)
        }
        
        i = 0
        k = 0
        mpnw1 = mpnw + 1
        t0[0] = Double(mpnw + 7)
        t1[0] = Double(mpnw + 7)
        t2[0] = Double(mpnw + 7)
        t3[0] = Double(mpnw + 7)
        t4[0] = Double(mpnw + 7)
        t5[0] = Double(mpnw + 7)
        t6[0] = Double(mpnw + 7)
        t7[0] = Double(mpnw + 7)
        t8[0] = Double(mpnw + 7)
        t9[0] = Double(mpnw + 7)
        tt[0] = Double(mpnw + 7)
        
        //   Set f1 = 1.
        
        f1[0] = 9.0
        f1[1] = Double(mpnw1)
        f1[2] = 1.0
        f1[3] = 0.0
        f1[4] = 1.0
        f1[5] = 0.0
        f1[6] = 0.0
        
        //   Check if argument is zero.  If so, result is - 1/2.
        
        if (s[2] == 0.0) {
            mpdmc (-0.5, 0, &t1, mpnw)
            mproun (&t1, mpnw)
            mpeq (t1, &z, mpnw)
            return
            // goto 200
        }
        
        //   Check if argument is negative.
        
        if s[2] < 0.0 {
            
            //   Check if argument is a negative even integer.  If so, the result is zero.
            
            mpmuld (s, 0.5, &t1, mpnw1)
            mpinfr (t1, &t2, &t3, mpnw1)
            if (t3[2] == 0.0) {
                t1[1] = Double(mpnw1)
                t1[2] = 0.0
                t1[3] = 0.0
                t1[4] = 0.0
                mproun (&t1, mpnw)
                mpeq (t1, &z, mpnw)
                return
                // goto 200
            }
            
            //   Otherwise compute zeta(1-s), and later apply the reflection formula.
            
            mpsub (f1, s, &tt, mpnw1)
        } else {
            mpeq (s, &tt, mpnw1)
        }
        
        //  Check if argument is large enough that computing with definition is faster.
        
        // if (tt .gt. mpreald (dlogb * mpnw / log (32.0 * mpnw), mpnw)) {
        
        d1 = dlogb * Double(mpnw) / log (32.0 * Double(mpnw))
        if (tt[3] > 1.0 || (tt[3] == 0.0 && tt[4] > d1)) {
            
            //  t1 = mpreal (1.0, mpnw)
            
            t1[1] = Double(mpnw1)
            t1[2] = 1.0
            t1[3] = 0.0
            t1[4] = 1.0
            t1[5] = 0.0
            t1[6] = 0.0
            
            for i in 2...itrmax {
                
                //    t2 = mpreal (Double(i), mpnw) ** tt
                
                mpdmc (Double(i), 0, &t4, mpnw1)
                mppower (t4, tt, t2, mpnw1)
                
                //    t3 = 1.0 / t2
                
                mpdiv (f1, t2, &t3, mpnw1)
                
                //    t1 = t1 + t3
                
                mpadd (t1, t3, &t2, mpnw1)
                mpeq (t2, &t1, mpnw1)
                
                if t3[2] == 0.0 || t3[3] < Double(-mpnw) {
                    mproun (&t1, mpnw)
                    mpeq (t1, &z, mpnw)
                    return
                    // goto 200
                }
            }
            
            print ("*** MPZETAEMR: iteration limit exceeded \(itrmax)")
            mpabrt (101)
        }
        
        // t0 = mpreal (1.0, mpnw)
        
        t0[1] = Double(mpnw1)
        t0[2] = 1.0
        t0[3] = 0.0
        t0[4] = 1.0
        t0[5] = 0.0
        t0[6] = 0.0
        
        nn = Int(dfrac) * mpnw1
        
        for k in 2...nn {
            //  t1 = mpreal (Double(k), mpnw) ** tt
            
            mpdmc (Double(k), 0, &t2, mpnw1)
            mppower (t2, tt, t1, mpnw1)
            
            //  t0 = t0 + 1.0 / t1
            
            mpdiv (f1, t1, &t2, mpnw1)
            mpadd (t0, t2, &t3, mpnw1)
            mpeq (t3, &t0, mpnw1)
        }
        
        // t0 = t0 + Double(nn) / (t1 * (tt - 1.0)) - 0.5d0 / t1
        
        mpdmc (Double(nn), 0, &t2, mpnw1)
        mpsub (tt, f1, &t3, mpnw1)
        mpmul (t1, t3, &t4, mpnw1)
        mpdiv (t2, t4, &t3, mpnw1)
        mpadd (t0, t3, &t2, mpnw1)
        mpdmc (0.5, 0, &t3, mpnw1)
        mpdiv (t3, t1, &t4, mpnw1)
        mpsub (t2, t4, &t0, mpnw1)
        
        // t3 = tt
        
        mpeq (tt, &t3, mpnw1)
        
        // t2 = t3 / (12.0 * Double(nn) * t1)
        
        mpmuld (t1, 12.0 * Double(nn), &t4, mpnw1)
        mpdiv (t3, t4, &t2, mpnw1)
        
        // t5 = Double(nn) * t1
        
        mpmuld (t1, Double(nn), &t5, mpnw1)
        
        // t9 = Double(nn) ** 2
        
        mpdmc (Double(nn), 0, &t6, mpnw1)
        mpmul (t6, t6, &t9, mpnw1)
        
        var flag = false
        for k in 2...min (nb2, itrmax) {
            //  t3 = t3 * (tt + Double(2*k - 2)) * (tt + Double(2*k - 3)) / &
            //         (Double(2 * k - 1) * Double(2 * k - 2))
            
            mpdmc (Double(2 * k - 2), 0, &t4, mpnw1)
            mpadd (tt, t4, &t6, mpnw1)
            mpdmc (Double(2 * k - 3), 0, &t7, mpnw1)
            mpadd (tt, t7, &t8, mpnw1)
            mpmul (t6, t8, &t7, mpnw1)
            mpmul (t3, t7, &t4, mpnw1)
            mpdmc (Double(2 * k - 1), 0, &t6, mpnw1)
            mpdmc (Double(2 * k - 2), 0, &t7, mpnw1)
            mpmul (t6, t7, &t8, mpnw1)
            mpdiv (t4, t8, &t3, mpnw1)
            
            //  t5 = t5 * t9
            
            mpmul (t5, t9, &t6, mpnw1)
            mpeq (t6, &t5, mpnw1)
            
            //  t7 = t3 * berne(k) / (Double(2 * k) * t5)
            
            //   The next few lines (to !+) are necessary, rather than a simple call
            //   to mpmul, to avoid a Fortran rank-mismatch error.
            
            //   mpmul (t3, berne(0,k), t4, mpnw1)
            
            ia = sign (1.0, berne[2][k])
            na = min (Int (abs (berne[2][k])), mpnw1)
            t8[1] = Double(mpnw1)
            t8[2] = Double(sign (Double(na), Double(ia)))
            
            for i in 2...na + 2 {
                t8[i+1] = berne[i+1][k]
            }
            
            t8[na+4] = 0.0
            t8[na+5] = 0.0
            mpmul (t3, t8, &t4, mpnw1)
            //+
            mpmuld (t5, Double(2 * k), &t6, mpnw1)
            mpdiv (t4, t6, &t7, mpnw1)
            
            //  t2 = t2 + t7
            
            mpadd (t2, t7, &t4, mpnw1)
            mpeq (t4, &t2, mpnw1)
            
            if t7[2] == 0.0 || t7[3] < t2[3] - Double(mpnw) {
                flag = true
                break
                // goto 110
            }
        }
        
        if !flag {
            print ("*** MPZETAEMR: iteration limit exceeded \(min (nb2, itrmax))")
            mpabrt (101)
        }
        
        //110 continue
        
        // zetaem = t0 + t2
        
        mpadd (t0, t2, &t1, mpnw1)
        
        //   If original argument was negative, apply the reflection formula.
        
        if (s[2] < 0.0) {
            mpgammar (tt, &t3, mpnw1)
            mpmul (t1, t3, &t2, mpnw1)
            mpmul (mppicon, tt, &t1, mpnw1)
            mpmuld (t1, 0.5, &t3, mpnw1)
            mpcssnr (t3, t4, t5, mpnw1)
            mpmul (t2, t4, &t1, mpnw1)
            mpmuld (mppicon, 2.0, &t2, mpnw1)
            mppower (t2, tt, t3, mpnw1)
            mpdiv (t1, t3, &t2, mpnw1)
            mpmuld (t2, 2.0, &t1, mpnw1)
        }
        
        // 200 continue
        
        mproun (&t1, mpnw)
        mpeq (t1, &z, mpnw)
        
    } // mpzetaemr
            
    
}
