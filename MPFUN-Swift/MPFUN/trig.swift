//
//  trig.swift - Trigonometric functions
//  MPFUN-Swift
//
//  Created by Mike Griebling on 6 May 2019.
//  Copyright © 2019 Computer Inspirations. All rights reserved.
//

import Foundation

extension MPFUN {
    
    static func mpagmr (_ a: MPRNumber, _ b: MPRNumber, _ c: inout MPRNumber, _ mpnw: Int) {
        
        //   This performs the arithmetic-geometric mean (AGM) iterations on A and B.
        //   The AGM algorithm is as follows: Set a_0 = a and b_0 = b, then iterate
        
        //    a_{k+1} = (a_k + b_k)/2
        //    b_{k+1} = sqrt (a_k * b_k)
        
        //   until convergence (i.e., until a_k = b_k to available precision).
        //   The result is returned in C.
        
        var mpnw1 : Int
        let itrmx = 50
        var s0 = MPRNumber(repeating: 0, count:mpnw+7)
        var s1 = s0; var s2 = s0; var s3 = s0
        
        if mpnw < 4 || a[0] < Double(mpnw + 4) || b[0] < abs (a[2]) + 4 || c[0] < Double(mpnw + 6) {
            print ("*** MPAGMR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        mpnw1 = mpnw + 1
        mpeq (a, &s1, mpnw1)
        mpeq (b, &s2, mpnw1)
        
        var flag = false
        for _ in 1...itrmx {
            mpadd (s1, s2, &s0, mpnw1)
            mpmuld (s0, 0.5, &s3, mpnw1)
            mpmul (s1, s2, &s0, mpnw1)
            mpsqrt (s0, &s2, mpnw1)
            mpeq (s3, &s1, mpnw1)
            mpsub (s1, s2, &s0, mpnw1)
            
            //   Check for convergence.
            
            if s0[2] == 0.0 || s0[3] < 1.0 - Double(mpnw1) {
                flag = true
                break // goto 100
            }
        }
        
        if !flag {
            print ("*** MPAGMR: Iteration limit exceeded.")
            mpabrt (5)
        }
        
        // 100 continue
        
        mproun (&s1, mpnw)
        mpeq (s1, &c, mpnw)
        
    } // mpagmr
    
    static func mpang (_ x: MPRNumber, _ y: MPRNumber, _ a: inout MPRNumber, _ mpnw: Int) {
        
        //   This computes the MPR angle A subtended by the MPR pair (X, Y) considered as
        //   a point in the x-y plane.  This is more useful than an arctan or arcsin
        //   routine, since it places the result correctly in the full circle, i.e.
        //   -Pi < A <= Pi.  Pi must be precomputed to at least MPNW words precision
        //   and the stored in the array in module MPMODA.
        
        //   The Taylor series for Sin converges much more slowly than that of Arcsin.
        //   Thus this routine does not employ Taylor series, but instead computes
        //   Arccos or Arcsin by solving Cos (a) = x or Sin (a) = y using one of the
        //   following Newton iterations, both of which converge to a:
        
        //     z_{k+1} = z_k - [x - Cos (z_k)] / Sin (z_k)
        //     z_{k+1} = z_k + [y - Sin (z_k)] / Cos (z_k)
        
        //   The first is selected if Abs (x) <= Abs (y); otherwise the second is used.
        //   These iterations are performed with a maximum precision level MPNW that
        //   is dynamically changed, approximately doubling with each iteration.
        
        //   If the precision level MPNW exceeds MPNWX words, this static func calls
        //   MPANGX instead.  By default, MPNWX = 100 (approx. 1450 digits).
        
        var iq, ix, iy,kk, mpnw1, mq, nx, ny, n1, n2 : Int
        var t1, t2, t3 : Double
        let cl2 = 1.4426950408889633 // let cpi = 3.141592653589793
        let mprxx = 1e-14; let mpnwx = 100; let nit = 3
        var s0 = MPRNumber(repeating: 0, count:mpnw+7)
        var s1 = s0; var s2 = s0; var s3 = s0; var s4 = s0; var s5 = s0
        
        // End of declaration
        
        if mpnw < 4 || x[0] < Double(mpnw+4) || x[0] < abs(x[2]) + 4 || y[0] < Double(mpnw+4) || y[0] < abs(y[2]) + 4 || a[0] < Double(mpnw + 6) {
            print ("*** MPANG: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ix = sign (1.0, x[2])
        nx = min (Int (abs (x[2])), mpnw)
        iy = sign (1.0, y[2])
        ny = min (Int (abs (y[2])), mpnw)
        mpnw1 = mpnw + 1
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        s4[0] = Double(mpnw + 7)
        s5[0] = Double(mpnw + 7)
        
        //   Check if both X and Y are zero.
        
        if (nx == 0 && ny == 0) {
            print ("*** MPANG: Both arguments are zero.")
            mpabrt (7)
        }
        
        //   Check if Pi has been precomputed.
        
        if mpnw1 > Int(mppicon[1]) {
            print ("*** MPANG: Pi must be precomputed to precision \(mpnw1) words.",
                "See documentation for details.")
            mpabrt (8)
        }
        
        //   If the precision level mpnw exceeds mpnwx words, mpangx.
        
        if (mpnw > mpnwx) {
            mpangx (x, y, &a, mpnw)
            return
        }
        
        //   Check if one of X or Y is zero.
        
        if (nx == 0) {
            mpeq (mppicon, &s0, mpnw1)
            if (iy > 0) {
                mpmuld (s0, 0.5, &a, mpnw)
            } else {
                mpmuld (s0, -0.5, &a, mpnw)
            }
            return
        } else if (ny == 0) {
            if (ix > 0) {
                a[1] = Double(mpnw)
                a[2] = 0.0
                a[3] = 0.0
            } else {
                mpeq (s0, &a, mpnw)
            }
            return
        }
        
        //   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.
        
        t1 = Double(mpnw1)
        mq = Int(cl2 * log (t1) + 1.0 - Double(mprxx))
        
        //   Normalize x and y so that x^2 + y^2 = 1.
        
        mpmul (x, x, &s0, mpnw1)
        mpmul (y, y, &s1, mpnw1)
        mpadd (s0, s1, &s2, mpnw1)
        mpsqrt (s2, &s3, mpnw1)
        mpdiv (x, s3, &s1, mpnw1)
        mpdiv (y, s3, &s2, mpnw1)
        
        //   Compute initial approximation of the angle.
        n1 = 0; n2 = 0; t2 = 0
        mpmdc (s1, &t1, &n1, mpnw1)
        mpmdc (s2, &t2, &n2, mpnw1)
        n1 = max (n1, -mpnbt)
        n2 = max (n2, -mpnbt)
        t1 = t1 * pow(2.0, Double(n1))
        t2 = t2 * pow(2.0, Double(n2))
        t3 = atan2 (t2, t1)
        mpdmc (t3, 0, &s5, mpnw1)
        
        //   The smaller of x or y will be used from now on to measure convergence.
        //   This selects the Newton iteration (of the two listed above) that has the
        //   largest denominator.
        
        if (abs (t1) <= abs (t2)) {
            kk = 1
            mpeq (s1, &s0, mpnw1)
        } else {
            kk = 2
            mpeq (s2, &s0, mpnw1)
        }
        
        mpnw1 = 4
        iq = 0
        
        //   Perform the Newton-Raphson iteration described above with a dynamically
        //   changing precision level MPNW (one greater than powers of two).
        
        for k in 1...mq {
            mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
            
            //100  continue
            while true {
                mpcssnr (s5, &s1, &s2, mpnw1)
                
                if (kk == 1) {
                    mpsub (s0, s1, &s3, mpnw1)
                    mpdiv (s3, s2, &s4, mpnw1)
                    mpsub (s5, s4, &s1, mpnw1)
                } else {
                    mpsub (s0, s2, &s3, mpnw1)
                    mpdiv (s3, s1, &s4, mpnw1)
                    mpadd (s5, s4, &s1, mpnw1)
                }
                mpeq (s1, &s5, mpnw1)
                
                if (k == mq - nit && iq == 0) {
                    iq = 1
                    // goto 100
                } else {
                    break
                }
            }
        }
        
        //   Restore original precision level.
        
        mproun (&s5, mpnw)
        mpeq (s5, &a, mpnw)
        
        // 120 continue
        
    } // mpang
    
    static func mpangx (_ x: MPRNumber, _ y: MPRNumber, _ a: inout MPRNumber, _ mpnw: Int) {
        
        //   This computes the MPR angle A subtended by the MPR pair (X, Y) considered as
        //   a point in the x-y plane.  This is more useful than an arctan or arcsin
        //   routine, since it places the result correctly in the full circle, i.e.
        //   -Pi < A <= Pi.  Pi and Log(2) must be precomputed to at least MPNW words
        //   precision and the stored in the array in module MPMODA.
        
        //   This routine simply calls mpclogx.  For modest precision, use mpang.
        
        var mp7 : Int
        var s0 = MPRNumber(repeating: 0, count: 2*mpnw+14); var s1 = s0
        
        // End of declaration
        
        if mpnw < 4 || x[0] < Double(mpnw+4) || x[0] < abs(x[2]) + 4 || y[0] < Double(mpnw+4) || y[0] < abs(y[2]) + 4 || a[0] < Double(mpnw+6) {
            print ("*** MPANGX: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        mp7 = mpnw + 7
        s0[0] = Double(mp7)
        s0[mp7] = Double(mp7)
        s1[0] = Double(mp7)
        s1[mp7] = Double(mp7)
        mpeq (x, &s0, mpnw)
        var tmp = s0
        mpeq (y, &tmp, mpnw)
        s0[mp7...] = tmp[0...]
        mpclogx (s0, &s1, mpnw)
        mpeq (MPRNumber(s1[mp7...]), &a, mpnw)
        
    } // mpangx

    static func mpcagm (_ a: MPRNumber, _ b: MPRNumber, _ c: inout MPRNumber, _ mpnw: Int) {
        
        //   This performs the arithmetic-geometric mean (AGM) iterations on A and B
        //   for MPC arguments A and B.
        //   The AGM algorithm is as follows: Set a_0 = a and b_0 = b, { iterate
        
        //    a_{k+1} = (a_k + b_k)/2
        //    b_{k+1} = sqrt (a_k * b_k)
        
        //   until convergence (i.e., until a_k = b_k to available precision).
        //   The result is returned in C.
        
        var la, lb, lc, mp7, mpnw1 : Int
        let itrmx = 50
        var s0 = MPRNumber(repeating:0, count:2*mpnw+14)
        var s1 = s0; var s2 = s0; var s3 = s0
        
        la = Int(a[0])
        lb = Int(b[0])
        lc = Int(c[0])
        if mpnw < 4 || a[0] < abs(a[2]) + 4 || a[la] < abs(a[la+2]) + 4 || b[0] < abs(b[2]) + 4 || b[lb] < abs(b[lb+2]) + 4 ||
            Int(c[0]) < mpnw + 6 || c[lc] < Double(mpnw + 6) {
            print ("*** MPCAGM: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        mp7 = mpnw + 7
        s0[0] = Double(mp7)
        s0[mp7] = Double(mp7)
        s1[0] = Double(mp7)
        s1[mp7] = Double(mp7)
        s2[0] = Double(mp7)
        s2[mp7] = Double(mp7)
        s3[0] = Double(mp7)
        s3[mp7] = Double(mp7)
        mpnw1 = mpnw + 1
        mpceq (a, &s1, mpnw1)
        mpceq (b, &s2, mpnw1)
        
        var flag = false
        for _ in 1...itrmx {
            mpcadd (s1, s2, &s0, mpnw1)
            mpmuld (s0, 0.5, &s3, mpnw1)
            var t = MPRNumber(s3[mp7...])
            mpmuld (MPRNumber(s0[mp7...]), 0.5, &t, mpnw1)
            s3[mp7...] = t[0...]
            mpcmul (s1, s2, &s0, mpnw1)
            mpcsqrt (s0, &s2, mpnw1)
            mpceq (s3, &s1, mpnw1)
            mpcsub (s1, s2, &s0, mpnw1)
            
            //   Check for convergence.
            
            if ((s0[2] == 0.0 || s0[3] < 1.0 - Double(mpnw1)) && (s0[mp7+2] == 0.0 || s0[mp7+3] < 1.0 - Double(mpnw1))) {
                flag = true
                break // goto 100
            }
        }
        
        if !flag {
        print ("*** MPCAGM: Iteration limit exceeded.")
        mpabrt (5)
        }
        
        // 100 continue
        
        mproun (&s1, mpnw)
        var t = MPRNumber(s1[mp7...])
        mproun (&t, mpnw); s1[mp7...] = t[0...]
        mpceq (s1, &c, mpnw)
    } // mpcagm

    static func mpcexp (_ a: MPRNumber, _ b: inout MPRNumber, _ mpnw : Int) {
        
        //   This computes Exp[A], for MPC A.
        
        //   The formula is:  E^a1 * (Cos[a2] + I * Sin[a2]), where a1 and a2 are
        //   the real and imaginary parts of A.
        
        //   If the precision level MPNW exceeds MPNWX words, this static func calls
        //   MPCEXP instead.  By default, MPNWX = 300 (about 4300 digits).
        
        var la, lb, mpnw1 : Int
        let mpnwx = 300
        var s0 = MPRNumber(repeating: 0, count:mpnw+7)
        var s1 = s0; var s2 = s0; var s3 = s0; var s4 = s0
        
        // End of declaration
        
        la = Int(a[0])
        lb = Int(b[0])
        if mpnw < 4 || a[0] < abs(a[2]) + 4 || a[la] < abs(a[la+2]) + 4 || Int(b[0]) < mpnw + 6 || Int(b[lb]) < mpnw + 6 {
            print ("*** MPCEXP: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   If the precision level mpnw exceeds mpnwx, mpcexpx.
        
        if (mpnw > mpnwx) {
            mpcexpx (a, &b, mpnw)
            return
        }
        
        mpnw1 = mpnw + 1
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        s4[0] = Double(mpnw + 7)
        
        mpexp (a, &s0, mpnw1)
        mpcssnr (MPRNumber(a[la...]), &s1, &s2, mpnw1)
        mpmul (s0, s1, &s3, mpnw1)
        mpmul (s0, s2, &s4, mpnw1)
        
        mproun (&s3, mpnw)
        mproun (&s4, mpnw)
        mpeq (s3, &b, mpnw)
        var t = MPRNumber(b[lb...])
        mpeq (s4, &t, mpnw)
        b[lb...] = t[0...]
        
        // 100 continue
        
    } // mpcexp
    
    static func mpcexpx (_ a: MPRNumber, _ b: inout MPRNumber, _ mpnw : Int) {
        
        //   This computes the exponential of the MPC number A and returns the MPC
        //   result in B.
        
        //   This routine employs the following Newton iteration, which converges to b:
        
        //     x_{k+1} = x_k + x_k * [a - Log (x_k)]
        
        //   These iterations are performed with a maximum precision level MPNW that
        //   is dynamically changed, approximately doubling with each iteration.
        //   For modest levels of precision, use mpcexp.
        
        var iq, la, lb, mpnw1, mp7, mq, nb, n0, n1 : Int
        var t0, t1, t2 : Double
        let cl2 = 1.4426950408889633; let nit = 3; let mprxx = 1e-14
        var s0 = MPRNumber(repeating: 0, count:2*mpnw+14)
        var s1 = s0; var s2 = s0; var s3 = s0; var s4 = s0
        var r1 = MPRNumber(repeating: 0, count:mpnw+7); var r2 = r1
        
        // End of declaration
        
        la = Int(a[0])
        lb = Int(b[0])
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || a[la] < abs (a[la+2]) + 4
            || Int(b[0]) < mpnw + 6 || Int(b[lb]) < mpnw + 6) {
            print ("*** MPCEXPX: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        // ia = sign (1.0, a[2])
        t1 = 0; n1 = 0
        nb = min (Int (abs (b[2])), mpnw)
        mpmdc (a, &t1, &n1, mpnw)
        
        //   Check for overflows and underflows.
        
        if n1 > 30 {
            if t1 > 0.0 {
                print ("*** MPCEXPX: Real part of argument is too large.")
                mpabrt (34)
            } else {
                b[1] = Double(mpnw)
                b[2] = 0.0
                b[3] = 0.0
                b[nb+1] = Double(mpnw)
                b[nb+2] = 0.0
                b[nb+3] = 0.0
                return
            }
        }
        
        t1 = t1 * pow(2.0, Double(n1))
        if abs (t1) > 1488522236.0 {
            if t1 > 0 {
                print ("*** MPCEXPX: Real part of argument is too large.")
                mpabrt (34)
            } else {
                b[1] = Double(mpnw)
                b[2] = 0.0
                b[3] = 0.0
                b[nb+1] = Double(mpnw)
                b[nb+2] = 0.0
                b[nb+3] = 0.0
                return
            }
        }
        
        //   Check if imaginary part is too large to compute meaningful cos/sin values.
        
        mpmdc (MPRNumber(a[la...]), &t1, &n1, mpnw)
        if (n1 >= mpnbt * (mpnw - 1)) {
            print ("*** MPCEXPX: imaginary part is too large to compute cos or sin.")
            mpabrt (28)
        }
        
        mpnw1 = mpnw + 1
        mp7 = mpnw + 7
        s0[0] = Double(mp7)
        s0[mp7] = Double(mp7)
        s1[0] = Double(mp7)
        s1[mp7] = Double(mp7)
        s2[0] = Double(mp7)
        s2[mp7] = Double(mp7)
        s3[0] = Double(mp7)
        s3[mp7] = Double(mp7)
        s4[0] = Double(mp7)
        s4[mp7] = Double(mp7)
        r1[0] = Double(mp7)
        r2[0] = Double(mp7)
        
        //   Check if Pi and Log(2) have been precomputed.
        
        if mpnw1 > Int(mplog2con[1]) {
            print ("*** MPCEXPX: Pi and Log(2) must be precomputed to precision \(mpnw1) words.",
                "See documentation for details.")
            mpabrt (53)
        }
        
        //   Reduce imaginary part to within -pi and pi.
        
        mpeq (a, &s4, mpnw1)
        mpmuld (mppicon, 2.0, &s0, mpnw1)
        mpdiv (MPRNumber(a[la...]), s0, &s1, mpnw1)
        mpnint (s1, &s2, mpnw1)
        mpmul (s2, s0, &s1, mpnw1)
        var t = MPRNumber(s4[mp7...])
        mpsub (MPRNumber(a[la...]), s1, &t, mpnw1)
        s4[mp7...] = t[0...]
        
        //   Check if imaginary part is -pi; if so correct to +pi.
        
        mpadd (MPRNumber(s4[mp7...]), mppicon, &s2, mpnw1)
        if s2[2] <= 0.0 || Int(s2[3]) < -mpnw {
            var t = MPRNumber(s4[mp7...])
            mpeq (mppicon, &t, mpnw1)
            s4[mp7...] = t[0...]
        }
        
        //   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.
        
        t2 = Double(mpnw1)
        mq = Int(cl2 * log (t2) + 1.0 - Double(mprxx))
        
        //   Compute initial approximation of Exp (A) (DP accuracy is OK).
        
        mpnw1 = 4
        
        // The following code (between here and iq = 0) is the equivalent of:
        //   mpcexp (s4, s3, mpnw1)
        
        mpdiv (s4, mplog2con, &s0, mpnw1)
        mpinfr (s0, &s1, &s2, mpnw1)
        mpmdc (s1, &t1, &n1, mpnw1)
        n1 = Int (t1 * pow(2.0, Double(n1)))
        t0 = 0; n0 = 0
        mpmdc (s2, &t0, &n0, mpnw1)
        n0 = min (max (n0, -100), 0)
        t0 = pow(2.0, (t0 * pow(2.0, Double(n0))))
        mpdmc (t0, n1, &s0, mpnw1)
        
        mpmdc (MPRNumber(s4[mp7...]), &t0, &n0, mpnw1)
        t0 = t0 * pow(2.0, Double(n0))
        t1 = cos (t0)
        t2 = sin (t0)
        if (abs (t1) < 1e-14) { t1 = 0.0 }
        if (abs (t2) < 1e-14) { t2 = 0.0 }
        mpdmc (t1, 0, &s1, mpnw1)
        t = MPRNumber(s1[mp7...])
        mpdmc (t2, 0, &t, mpnw1); s1[mp7...] = t[0...]
        mpmul (s0, s1, &s3, mpnw1)
        t = MPRNumber(s3[mp7...])
        mpmul (s0, MPRNumber(s1[mp7...]), &t, mpnw1); s3[mp7...] = t[0...]
        iq = 0
        
        //   Perform the Newton-Raphson iteration described above with a dynamically
        //   changing precision level MPNW (one greater than powers of two).
        
        for k in 0...mq {
            if (k > 1) { mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1 }
            
            // 100  continue
            
            while true {
                mpclogx (s3, &s0, mpnw1)
                
                //   Check if we need to add or subtract 2*pi to the output of imaginary part,
                //   in order to remain consistent with previous iterations.
                
                mpsub (MPRNumber(s4[mp7...]), MPRNumber(s0[mp7...]), &r1, 4)
                if (r1[2] != 0.0 && r1[3] >= -1.0) {
                    var t = MPRNumber(s0[mp7...])
                    if (r1[2] > 0.0) {
                        mpadd (MPRNumber(s0[mp7...]), mppicon, &r1, mpnw1)
                        mpadd (r1, mppicon, &t, mpnw1)
                    } else {
                        mpsub (MPRNumber(s0[mp7...]), mppicon, &r1, mpnw1)
                        mpsub (r1, mppicon, &t, mpnw1)
                    }
                    s0[mp7...] = t[0...]
                }
                
                mpcsub (s4, s0, &s1, mpnw1)
                mpcmul (s3, s1, &s2, mpnw1)
                mpcadd (s3, s2, &s1, mpnw1)
                mpceq (s1, &s3, mpnw1)
                if (k == mq - nit && iq == 0) {
                    iq = 1 // goto 100
                } else {
                    break
                }
            }
        }
        
        //   Restore original precision level.
        
        mproun (&s1, mpnw)
        t = MPRNumber(s1[mp7...])
        mproun (&t, mpnw); s1[mp7...] = t[0...]
        mpceq (s1, &b, mpnw)
        
//        130 continue
        
    } // mpcexpx

    static func mpclog (_ a: MPRNumber, _ b: inout MPRNumber, _ mpnw : Int) {
        
        //   This computes Log[A], for MPC A.
        
        //   The formula is:  1/2 * Log[r] + I * Theta, where r = a1^2 + a2^2,
        //   Theta is the angle corresponding to (a1, a2), and a1 and a2 are the
        //   real and imaginary parts of A.
        
        //   If the precision level MPNW exceeds MPNWX words, this static func calls
        //   MPCLOGX instead.  By default, MPNWX = 30 (about 400 digits).
        
        var la, lb, mpnw1 : Int
        let mpnwx = 30
        var s0 = MPRNumber(repeating:0, count:mpnw+7)
        var s1 = s0; var s2 = s0; var s3 = s0; var s4 = s0
        
        // End of declaration
        
        la = Int(a[0])
        lb = Int(b[0])
        if mpnw < 4 || a[0] < abs(a[2]) + 4 || a[la] < abs(a[la+2]) + 4 || Int(b[0]) < mpnw + 6 || Int(b[lb]) < mpnw + 6 {
            print ("*** MPCLOG: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   If precision level mpnw exceeds mpnwx words, mpclogx.
        
        if mpnw > mpnwx {
            mpclogx (a, &b, mpnw)
            return //goto 100
        }
        
        mpnw1 = mpnw + 1
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        s4[0] = Double(mpnw + 7)
        
        mpmul (a, a, &s0, mpnw1)
        mpmul (MPRNumber(a[la...]), MPRNumber(a[la...]), &s1, mpnw1)
        mpadd (s0, s1, &s2, mpnw1)
        mplog (s2, &s3, mpnw1)
        mpmuld (s3, 0.5, &s0, mpnw1)
        mpang (a, MPRNumber(a[la...]), &s1, mpnw1)
        
        mproun (&s0, mpnw)
        mproun (&s1, mpnw)
        mpeq (s0, &b, mpnw)
        var t = MPRNumber(b[lb...])
        mpeq (s1, &t, mpnw); b[lb...] = t[0...]
        
        // 100 continue
    } // mpclog
    
    static func mpclogx (_ a: MPRNumber, _ b: inout MPRNumber, _ mpnw : Int) {
        
        //   This computes the natural logarithm of the MP number A and returns the MP
        //   result in B.  Pi and Log(2) must be precomputed to at least MPNW words
        //   precision and the stored in the arrays in module MPMODA.
        
        //   This uses the following algorithm, which is due to Salamin and Brent.  If
        //   A is extremely close to 1, use a Taylor series.  Otherwise select n such
        //   that z = a * 2^n is at least 2^m, where m is the number of bits of desired
        //   precision in the result.  Then
        
        //   Log(x) = Pi / [2 AGM (1, 4/x)]
        
        //   For modest precision, or if A is close to 2, use mplog.
        
        var iss, la, lb, mpnw1, mp7, na1, na2, n1, n2 : Int
        var st, tol, t1, tn : Double
        let itrmax = 1000000; let rtol = pow(0.5, Double(7))
        var f1 = MPRNumber(repeating:0, count:19); var f4 = f1
        var s0 = MPRNumber(repeating:0, count:2*mpnw+14); var s1 = s0; var s2 = s0
        var s3 = s0; var s4 = s0
        
        // End of declaration
        
        la = Int(a[0])
        lb = Int(b[0])
        if mpnw < 4 || a[0] < abs(a[2]) + 4 || a[la] < abs(a[la+2]) + 4 || Int(b[0]) < mpnw+6 || Int(b[lb]) < mpnw+6 {
            print ("*** MPCLOGX: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //      let ia1 = sign (1.0, a[2])
        na1 = min (Int (abs (a[2])), mpnw)
        //      let ia2 = sign (1.0, a[la+2])
        na2 = min (Int (abs (a[la+2])), mpnw)
        
        if na1 == 0 && na2 == 0 {
            print ("*** MPCLOGX: Argument is zero.")
            mpabrt (50)
        }
        
        //   Check if input is exactly one.
        
        if a[2] == 1.0 && a[3] == 0.0 && a[4] == 1.0 && a[la+2] == 0.0 {
            b[1] = Double(mpnw)
            b[2] = 0.0
            b[3] = 0.0
            b[4] = 0.0
            b[lb+1] = Double(mpnw)
            b[lb+2] = 0.0
            b[lb+3] = 0.0
            b[lb+4] = 0.0
            return // goto 120
        }
        
        mpnw1 = mpnw + 1
        mp7 = mpnw + 7
        s0[0] = Double(mp7)
        s0[mp7] = Double(mp7)
        s1[0] = Double(mp7)
        s1[mp7] = Double(mp7)
        s2[0] = Double(mp7)
        s2[mp7] = Double(mp7)
        s3[0] = Double(mp7)
        s3[mp7] = Double(mp7)
        s4[0] = Double(mp7)
        s4[mp7] = Double(mp7)
        
        //   Check if Pi and Log(2) have been precomputed.
        
        if mpnw1 > Int(mplog2con[1]) {
            print ("*** MPCLOGX: Pi and Log(2) must be precomputed to precision \(mpnw1) words.",
                "See documentation for details.")
            mpabrt (53)
        }
        
        f1[0] = 9.0
        f1[1] = Double(mpnw1)
        f1[2] = 1.0
        f1[3] = 0.0
        f1[4] = 1.0
        f1[5] = 0.0
        f1[6] = 0.0
        f1[9] = 9.0
        f1[10] = Double(mpnw1)
        f1[11] = 0.0
        f1[12] = 0.0
        f1[13] = 0.0
        
        f4[0] = 9.0
        f4[1] = Double(mpnw1)
        f4[2] = 1.0
        f4[3] = 0.0
        f4[4] = 4.0
        f4[5] = 0.0
        f4[6] = 0.0
        f4[9] = 9.0
        f4[10] = Double(mpnw1)
        f4[11] = 0.0
        f4[12] = 0.0
        f4[13] = 0.0
        
        //   If the argument is sufficiently close to 1, employ a Taylor series.
        
        mpcsub (a, f1, &s0, mpnw1)
        
        if (s0[2] == 0.0 || s0[3] <= min(-2.0, -rtol * Double(mpnw1))) && (s0[mp7+2] == 0.0 || s0[mp7+3] <= min(-2.0, -rtol * Double(mpnw1))) {
            mpceq (s0, &s1, mpnw1)
            mpceq (s1, &s2, mpnw1)
            //       i1 = 1
            iss = 1
            if (s0[2] == 0.0) {
                tol = s0[mp7+3] - Double(mpnw1)
            } else if (s0[mp7+2] == 0.0) {
                tol = s0[3] - Double(mpnw1)
            } else {
                tol = max (s0[3], s0[mp7+3]) - Double(mpnw1)
            }
            
            for i1 in 2...itrmax {
                iss = -iss
                st = Double(iss * i1)
                mpcmul (s1, s2, &s3, mpnw1)
                mpceq (s3, &s2, mpnw1)
                mpdivd (s3, st, &s4, mpnw1)
                var t = MPRNumber(s4[mp7...])
                mpdivd (MPRNumber(s3[mp7...]), st, &t, mpnw1)
                s4[mp7...] = t[0...]
                mpcadd (s0, s4, &s3, mpnw1)
                mpceq (s3, &s0, mpnw1)
                if (s4[2] == 0.0 || s4[3] < tol) && (s4[mp7+2] == 0.0 || s4[mp7+3] < tol) {
                    //goto 110
                    //  Restore original precision level.
                    mproun (&s0, mpnw)
                    t = MPRNumber(s0[mp7...])
                    mproun (&t, mpnw); s0[mp7...] = t[0...]
                    mpceq (s0, &b, mpnw)
                    return
                }
            }
            
            print ("*** MPCLOGX: Iteration limit exceeded: \(itrmax)")
            mpabrt (54)
        }
        
        //   Multiply the input by a large power of two.
        n1 = 0; t1 = 0
        mpmdc (a, &t1, &n1, mpnw1)
        n2 = mpnbt * (mpnw1 / 2 + 2) - n1
        tn = Double(n2)
        mpdmc (1.0, n2, &s1, mpnw1)
        mpmul (a, s1, &s0, mpnw1)
        var t = MPRNumber(s0[mp7...])
        mpmul (MPRNumber(a[la...]), s1, &t, mpnw1)
        s0[mp7...] = t[0...]
        
        //   Perform AGM iterations.
        
        mpceq (f1, &s1, mpnw1)
        mpcdiv (f4, s0, &s2, mpnw1)
        mpcagm (s1, s2, &s3, mpnw1)
        
        //   Compute Pi / (2 * A), where A is the limit of the AGM iterations.
        
        mpmuld (s3, 2.0, &s0, mpnw1)
        t = MPRNumber(s0[mp7...])
        mpmuld (MPRNumber(s3[mp7...]), 2.0, &t, mpnw1)
        s0[mp7...] = t[0...]
        mpeq (mppicon, &s3, mpnw1)
        t = MPRNumber(s3[mp7...])
        mpdmc (0.0, 0, &t, mpnw1)
        s3[mp7...] = t[0...]
        mpcdiv (s3, s0, &s1, mpnw1)
        
        //   Subtract TN * Log(2).
        
        mpeq (mplog2con, &s3, mpnw1)
        mpmuld (s3, tn, &s2, mpnw1)
        mpsub (s1, s2, &s0, mpnw1)
        
        //   Check if imaginary part is -pi; if so correct to +pi.
        
        mpadd (MPRNumber(s1[mp7...]), mppicon, &s2, mpnw1)
        t = MPRNumber(s0[mp7...])
        if (s2[2] <= 0.0 || Int(s2[3]) < -mpnw) {
            mpeq (mppicon, &t, mpnw1)
        } else {
            mpeq (MPRNumber(s1[mp7...]), &t, mpnw1)
        }
        s0[mp7...] = t[0...]
        
        // 110 continue
        
        //  Restore original precision level.
        
        mproun (&s0, mpnw)
        t = MPRNumber(s0[mp7...])
        mproun (&t, mpnw); s0[mp7...] = t[0...]
        mpceq (s0, &b, mpnw)
        
        // 120 continue
        
    } // mpclogx
    
    static func mpcpowcc (_ a: MPRNumber, _ b: MPRNumber, _ c: inout MPRNumber,_ mpnw : Int) {
        
        //   This computes A^B, where A and B are MPC.
        
        var la, lb, lc : Int
        var s1 = MPRNumber(repeating:0, count:2*mpnw+12); var s2 = s1
        
        // End of declaration
        
        la = Int(a[0])
        lb = Int(b[0])
        lc = Int(c[0])
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || a[la] < abs(a[la+2]) + 4 || b[0] < abs (b[2]) + 4 || b[lb] < abs(b[lb+2]) + 4
            || Int(c[0]) < mpnw + 6 || Int(c[lc]) < mpnw + 6 {
            print ("*** MPCPOWCC: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let l3 = mpnw + 6
        let dl3 = Double(l3)
        s1[0] = dl3
        s1[l3] = dl3
        s2[0] = dl3
        s2[l3] = dl3
        mpclog (a, &s1, mpnw)
        mpcmul (s1, b, &s2, mpnw)
        mpcexp (s2, &c, mpnw)
        
    } // mpcpowcc
    
    static func mpcpowcr (_ a: MPRNumber, _ b: MPRNumber, _ c: inout MPRNumber, _ mpnw : Int) {
        
        //   This computes A^B, where A is MPC and B is MPR.
        
        var la, lc, l3 : Int
        var s1 = MPRNumber(repeating:0, count:2*mpnw+12); var s2 = s1
        
        // End of declaration
        
        la = Int(a[0])
//        lb = Int(b[0])
        lc = Int(c[0])
        if mpnw < 4 || a[0] < abs(a[2]) + 4 || a[la] < abs(a[la+2]) + 4 || b[0] < abs (b[2]) + 4 || Int(c[0]) < mpnw+6 || Int(c[lc]) < mpnw+6 {
            print ("*** MPCPOWCR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        l3 = mpnw + 6
        let dl3 = Double(l3)
        s1[0] = dl3
        s1[l3] = dl3
        s2[0] = dl3
        s2[l3] = dl3
        mpclog (a, &s1, mpnw)
        mpmul (b, s1, &s2, mpnw)
        var t = MPRNumber(s2[l3...])
        mpmul (b, MPRNumber(s1[l3...]), &t, mpnw); s2[l3...] = t[0...]
        mpcexp (s2, &c, mpnw)
        
    } //  mpcpowcr
    
    static func mpcpowrc (_ a: MPRNumber, _ b: MPRNumber, _ c: inout MPRNumber, _ mpnw : Int) {
        
        //   This computes A^B, where A is MPR and and B is MPC.
        
        var lb, lc, l3 : Int
        var s1 = MPRNumber(repeating:0, count:2*mpnw+12); var s2 = s1
        
        // End of declaration
        
        // la = Int(a[0])
        lb = Int(b[0])
        lc = Int(c[0])
        if mpnw < 4 || a[0] < abs(a[2]) + 4 || b[0] < abs (b[2]) + 4 || b[lb] < abs(b[lb+2]) + 4 || Int(c[0]) < mpnw + 6 || Int(c[lc]) < mpnw + 6 {
            print ("*** MPCPOWRC: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        l3 = mpnw + 6
        s1[0] = Double(l3)
        s2[0] = Double(l3)
        s2[l3] = Double(l3)
        mplog (a, &s1, mpnw)
        mpmul (s1, b, &s2, mpnw)
        var t = MPRNumber(s2[l3...])
        mpmul (s1, MPRNumber(b[lb...]), &t, mpnw); s2[l3...] = t[0...]
        mpcexp (s2, &c, mpnw)
    } // mpcpowrc
    
    static func mpcsshr (_ a: MPRNumber, _ x: inout MPRNumber, _ y: inout MPRNumber, _ mpnw : Int) {
        
        //   This computes the hyperbolic cosine and sine of the MPR number A and
        //   returns the two MPR results in X and Y, respectively.  If the argument
        //   is very close to zero, a Taylor series is used; otherwise this routine
        //   calls mpexp.
        
        var mpnw1, mpnw2 : Int
        let itrmx = 1000000  // let mpnwx = 700
        var f = MPRNumber(repeating:0, count:10)
        var s0 = MPRNumber(repeating:0, count:mpnw+7)
        var s1 = s0; var s2 = s0; var s3 = s0
        var t2 : Double
        
        // End of declaration
        
        if (mpnw < 4 || Int(a[0]) < mpnw+4 || a[0] < abs (a[2]) + 4 || Int(x[0]) < mpnw+6 || Int(y[0]) < mpnw+6) {
            print ("*** MPCSSHR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        mpnw1 = mpnw + 1
        f[0] = 9.0
        f[1] = Double(mpnw)
        f[2] = 1.0
        f[3] = 0.0
        f[4] = 1.0
        f[5] = 0.0
        f[6] = 0.0
        
        //   If argument is very small, compute the sinh using a Taylor series.
        //   This avoids accuracy loss that otherwise occurs by using exp.
        
        if s0[3] < -1.0 {
            mpeq (a, &s0, mpnw1)
            mpmul (s0, s0, &s2, mpnw1)
            mpnw2 =  mpnw1
            
            //   The working precision used to compute each term can be linearly reduced
            //   as the computation proceeds.
            var flag = false
            for j in 1...itrmx {
                t2 = Double((2 * j) * (2 * j + 1))
                mpmul (s2, s1, &s3, mpnw2)
                mpdivd (s3, t2, &s1, mpnw2)
                mpadd (s1, s0, &s3, mpnw1)
                mpeq (s3, &s0, mpnw1)
                
                //   Check for convergence of the series, and adjust working precision
                //   for the next term.
                
                if s1[2] == 0.0 || s1[3] < s0[3] - Double(mpnw1) { flag = true; break; /* goto 110 */ }
                mpnw2 = min (max (mpnw1 + Int (s1[3] - s0[3]) + 1, 4), mpnw1)
            }
            
            if !flag {
                print ("*** MPCSSHR: Iteration limit exceeded.")
                mpabrt (29)
            }
            
            //110 continue
            
            mpmul (s0, s0, &s2, mpnw1)
            mpadd (f, s2, &s3, mpnw1)
            mpsqrt (s3, &s1, mpnw1)
            mproun (&s1, mpnw)
            mpeq (s1, &x, mpnw)
            mproun (&s0, mpnw)
            mpeq (s0, &y, mpnw)
        } else {
            mpexp (a, &s0, mpnw1)
            mpdiv (f, s0, &s1, mpnw1)
            mpadd (s0, s1, &s2, mpnw1)
            mpmuld (s2, 0.5, &s3, mpnw1)
            mproun (&s3, mpnw)
            mpeq (s3, &x, mpnw)
            mpsub (s0, s1, &s2, mpnw1)
            mpmuld (s2, 0.5, &s3, mpnw1)
            mproun (&s3, mpnw)
            mpeq (s3, &y, mpnw)
        }
        
        //100 continue
        
    } // mpcsshr

    
    static func mpcssnr (_ a: MPRNumber, _ x: inout MPRNumber, _ y: inout MPRNumber, _ mpnw : Int) {
        
        //   This computes the cosine and sine of the MPR number A and returns the
        //   two MPR results in X and Y, respectively.  Pi must be precomputed to
        //   at least MPNW words precision and the stored in the array MPPICON in
        //   module MPMODA.
        
        //   This routine uses the conventional Taylor series for Sin (s):
        
        //   Sin (s) =  s - s^3 / 3// + s^5 / 5// - s^7 / 7// ...
        
        //   where the argument S has been reduced to (-pi, pi).  To further
        //   accelerate the series, the reduced argument is divided by 2^NQ, where NQ
        //   is computed as int (sqrt (0.5d0 * N)), where N is the precision in bits.
        //   After convergence of the series, the double-angle formulas for cos are
        //   applied NQ times.
        
        //   If the precision level MPNW exceeds MPNWX, this static func calls
        //   MPCSSNX instead.  By default, mpnwx = 100000 (approx. 1450000 digits).
        
        var iss, mpnw1, mpnw2 : Int
        var na, nq, n1 : Int
        let itrmx = 1000000; let mpnwx = 100000
        var t1, t2 : Double
        var f1 = MPRNumber(repeating:0, count:9); var f2 = f1
        var s0 = MPRNumber(repeating:0, count:mpnw+7); var s1 = s0; var s2 = s1; var s3 = s1
        var s4 = s1; var s5 = s1; var s6 = s1
        
        // End of declaration
        
        if mpnw < 4 || Int(a[0]) < mpnw+4 || a[0] < abs (a[2]) + 4 || Int(x[0]) < mpnw+6 || Int(y[0]) < mpnw+6 {
            print ("*** MPCSSNR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   If the precision level mpnw exceeds mpnwx, mpcssx.
        
        if mpnw > mpnwx {
            mpcssnx (a, &x, &y, mpnw)
            return // goto 120
        }
        
        // ia = sign (1.0, a[2])
        na = min (Int (abs (a[2])), mpnw)
        if na == 0 {
            x[1] = Double(mpnw)
            x[2] = 1.0
            x[3] = 0.0
            x[4] = 1.0
            y[1] = Double(mpnw)
            y[2] = 0.0
            y[3] = 0.0
            return // goto 120
        }
        
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        s4[0] = Double(mpnw + 7)
        s5[0] = Double(mpnw + 7)
        s6[0] = Double(mpnw + 7)
        mpnw1 = mpnw + 1
        
        //   Set f1 = 1 and f2 = 1/2.
        
        f1[0] = 9.0
        f1[1] = Double(mpnw)
        f1[2] = 1.0
        f1[3] = 0.0
        f1[4] = 1.0
        f1[5] = 0.0
        f1[6] = 0.0
        f2[0] = 9.0
        f2[1] = Double(mpnw)
        f2[2] = 1.0
        f2[3] = -1.0
        f2[4] = 0.5 * Double(mpbdx)
        f2[5] = 0.0
        f2[6] = 0.0
        
        //   Check if Pi and Sqrt(2)/2 have been precomputed in data statements to the
        //   requested precision.
        
        if mpnw1 > Int(mppicon[1]) {
            print ("*** MPCSSNR: Pi and Sqrt(2)/2 must be precomputed to precision \(mpnw1) words).",
                "See documentation for details.")
            mpabrt (27)
        }
        
        //   Check if argument is too large to compute meaningful cos/sin values.
        t1 = 0; n1 = 0
        mpmdc(a, &t1, &n1, mpnw)
        if n1 >= mpnbt * (mpnw - 1) {
            print ("*** MPCSSNR: argument is too large to compute cos or sin.")
            mpabrt (28)
        }
        
        //   Reduce to between - Pi and Pi.
        
        mpmuld (mppicon, 2.0, &s0, mpnw1)
        mpdiv (a, s0, &s1, mpnw1)
        mpnint (s1, &s2, mpnw1)
        mpmul (s0, s2, &s4, mpnw1)
        mpsub (a, s4, &s3, mpnw1)
        
        //   Check if reduced argument is zero.  If so { cos = 1 and sin = 0.
        
        if (s3[2] == 0.0) {
            s0[1] = Double(mpnw1)
            s0[2] = 1.0
            s0[3] = 0.0
            s0[4] = 1.0
            s0[5] = 0.0
            s0[6] = 0.0
            s1[1] = Double(mpnw1)
            s1[2] = 0.0
            s1[3] = 0.0
            // goto 115
            //   Restore original precision level.
            
            mproun (&s0, mpnw)
            mproun (&s1, mpnw)
            mpeq (s0, &x, mpnw)
            mpeq (s1, &y, mpnw)
            return
        }
        
        //   Determine nq to scale reduced argument, then divide by 2^nq.
        //   If reduced argument is very close to zero, then nq = 0.
        
        if s3[3] >= -1.0 {
            nq = Int (sqrt (0.5 * Double(mpnw1 * mpnbt)))
        } else {
            nq = 0
        }
        
        // write (6, *) "nq =", nq
        
        mpdivd (s3, pow(2.0, Double(nq)), &s0, mpnw1)
        mpeq (s0, &s1, mpnw1)
        
        //   Compute the sin of the reduced argument of s1 using a Taylor series.
        
        mpmul (s0, s0, &s2, mpnw1)
        mpnw2 =  mpnw1
        iss = Int(s0[2])
        
        //   The working precision used to compute each term can be linearly reduced
        //   as the computation proceeds.
        var flag = false
        for i1 in 1...itrmx {
            t2 = Double(-(2 * i1) * (2 * i1 + 1))
            mpmul (s2, s1, &s3, mpnw2)
            mpdivd (s3, t2, &s1, mpnw2)
            mpadd (s1, s0, &s3, mpnw1)
            mpeq (s3, &s0, mpnw1)
            
            //   Check for convergence of the series, and adjust working precision
            //   for the next term.
            
            if (s1[2] == 0.0 || s1[3] < s0[3] - Double(mpnw1)) { flag = true; break; /* goto 110 */ }
            mpnw2 = min (max (mpnw1 + Int (s1[3] - s0[3]) + 1, 4), mpnw1)
        }
        
        if !flag {
            print ("*** MPCSSNR: Iteration limit exceeded.")
            mpabrt (29)
        }
        
        // 110 continue
        
        if nq > 0 {
            
            //   Apply the formula cos(2*x) = 2*cos^2(x) - 1 NQ times to produce
            //   the cosine of the reduced argument, except that the first iteration is
            //   cos(2*x) = 1 - 2*sin^2(x), since we have computed sin(x) above.
            //   Note that these calculations are performed as 2 * (cos^2(x) - 1/2) and
            //   2 * (1/2 - sin^2(x)), respectively, to avoid loss of precision.
            
            mpmul (s0, s0, &s4, mpnw1)
            mpsub (f2, s4, &s5, mpnw1)
            mpmuld (s5, 2.0, &s0, mpnw1)
            
            for _ in 2...nq {
                mpmul (s0, s0, &s4, mpnw1)
                mpsub (s4, f2, &s5, mpnw1)
                mpmuld (s5, 2.0, &s0, mpnw1)
            }
            
            //   Compute sin of result and correct sign.
            
            mpmul (s0, s0, &s4, mpnw1)
            mpsub (f1, s4, &s5, mpnw1)
            mpsqrt (s5, &s1, mpnw1)
            if iss < 1 { s1[2] = -s1[2] }
        } else {
            
            //   In case nq = 0, compute cos of result.
            
            mpeq (s0, &s1, mpnw1)
            mpmul (s0, s0, &s4, mpnw1)
            mpsub (f1, s4, &s5, mpnw1)
            mpsqrt (s5, &s0, mpnw1)
        }
        
        // 115 continue
        
        //   Restore original precision level.
        
        mproun (&s0, mpnw)
        mproun (&s1, mpnw)
        mpeq (s0, &x, mpnw)
        mpeq (s1, &y, mpnw)
        
        // 120 continue
        
    } //mpcssnr
    
    static func mpcssnx (_ a: MPRNumber, _ x: inout MPRNumber, _ y: inout MPRNumber, _ mpnw : Int) {
        
        //   This computes the cosine and sine of the MPR number A and returns the
        //   two MPR results in X and Y, respectively.  Pi and Log(2) must be precomputed to at
        //   least MPNW words precision and the stored in the array in module MPMODA.
        
        //   This routine merely calls mpcexp.  For modest levels of precision, use mpcssn.
        
        var mp7 : Int
        var f = MPRNumber(repeating:0, count:9)
        var s0 = MPRNumber(repeating:0, count:mpnw+7); var s1 = s0
        
        // End of declaration
        
        if (mpnw < 4 || Int(a[0]) < mpnw + 4 || a[0] < abs (a[2]) + 4 || Int(x[0]) < mpnw + 6 || Int(y[0]) < mpnw + 6) {
            print ("*** MPCSSNX: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        mp7 = mpnw + 7
        s0[0] = Double(mp7)
        s0[mp7] = Double(mp7)
        s1[0] = Double(mp7)
        s1[mp7] = Double(mp7)
        f[0] = 9.0
        f[1] = Double(mpnw)
        f[2] = 0.0
        f[3] = 0.0
        mpeq (f, &s0, mpnw)
        var t = MPRNumber(s0[mp7...])
        mpeq (a, &t, mpnw)
        s0[mp7...] = t[0...]
        mpcexpx (s0, &s1, mpnw)
        mpeq (s1, &x, mpnw)
        mpeq (MPRNumber(s1[mp7...]), &y, mpnw)
        
    } // mpcssnx
    
    static func ceiling (_ x : Double) -> Int {
        let xa = abs(x)
        let xi = Int(xa)  // Int() aka floor()
        let delta = xa - Double(xi)
        if x < 0 {
            if delta > 0 { return -(xi + 1) }
            return -xi
        } else {
            if delta > 0 { return xi + 1 }
            return xi
        }
    }
    
    static func mpegammaq (_ egamma : inout MPRNumber, _ mpnw: Int) {
        
        //   This computes Euler"s gamma to available precision (MPNW mantissa words).
        //   The algorithm is the following, which is an improvement to a scheme due to
        //   Sweeney (see https://www.davidhbailey.com/dhbpapers/const.pdf):
        
        //   Select N such that 1/(2^N * Exp(2^N)) < desired epsilon. Then compute
        //   Gamma = 2^N/Exp(2^N) * (Sum_{m >= 0} 2^(m*N)/(m+1)// * H(m+1)) - N * Log(2),
        //   where H(m) = 1 + 1/2 + ... + 1/m.
        
        var mpnw1, neps, nn : Int
        let itrmx = 1000000
        var f = MPRNumber(repeating:0, count:9)
        var s0 = MPRNumber(repeating:0, count:mpnw+7); var s1 = s0; var s2 = s1; var s3 = s1
        var s4 = s1; var s5 = s1; var s6 = s1; var s7 = s1
        
        // End of declaration.
        
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        s4[0] = Double(mpnw + 7)
        s5[0] = Double(mpnw + 7)
        s6[0] = Double(mpnw + 7)
        s7[0] = Double(mpnw + 7)
        mpnw1 = mpnw + 1
        
        //   Check if Log(2) has been precomputed.
        
        if mpnw1 > Int(mplog2con[1]) {
            print ("*** MPEGAMMA: Log(2] must be precomputed to precision \(mpnw1) words.",
                "See documentation for details.")
            mpabrt (35)
        }
        
        //   Compute eps and nn based on precision level.
        
        neps = -mpnw1 - 1
        nn = ceiling (log (Double(mpnw1 * mpnbt + mpnbt) * log (2.0)) / log (2.0))
        
        //   Initialize s0 through s4 to 1.
        
        s0[1] = Double(mpnw)
        s0[2] = 1.0
        s0[3] = 0.0
        s0[4] = 1.0
        s0[5] = 0.0
        s0[6] = 0.0
        
        s1[1] = Double(mpnw)
        s1[2] = 1.0
        s1[3] = 0.0
        s1[4] = 1.0
        s1[5] = 0.0
        s1[6] = 0.0
        
        s2[1] = Double(mpnw)
        s2[2] = 1.0
        s2[3] = 0.0
        s2[4] = 1.0
        s2[5] = 0.0
        s2[6] = 0.0
        
        s3[1] = Double(mpnw)
        s3[2] = 1.0
        s3[3] = 0.0
        s3[4] = 1.0
        s3[5] = 0.0
        s3[6] = 0.0
        
        s4[1] = Double(mpnw)
        s4[2] = 1.0
        s4[3] = 0.0
        s4[4] = 1.0
        s4[5] = 0.0
        s4[6] = 0.0
        
        s7[1] = Double(mpnw)
        s7[2] = 1.0
        s7[3] = 0.0
        s7[4] = 2.0
        s7[5] = 0.0
        s7[6] = 0.0
        
        //   Set s7 = 2^nn.
        
        mpdmc (1.0, nn, &s7, mpnw1)
        
        //  Set f = 1.
        
        f[0] = 9.0
        f[1] = Double(mpnw1)
        f[2] = 1.0
        f[3] = 0.0
        f[4] = 1.0
        f[5] = 0.0
        f[6] = 0.0
        
        var flag = false
        for m in 1...itrmx {
            mpmul (s7, s0, &s5, mpnw1)
            mpeq (s5, &s0, mpnw1)
            mpdmc (Double(m + 1), 0, &s5, mpnw1)
            mpdiv (f, s5, &s6, mpnw1)
            mpadd (s1, s6, &s5, mpnw1)
            mpeq (s5, &s1, mpnw1)
            mpmuld (s2, Double(m+1), &s5, mpnw1)
            mpeq (s5, &s2, mpnw1)
            mpmul (s0, s1, &s5, mpnw1)
            mpdiv (s5, s2, &s3, mpnw1)
            mpadd (s3, s4, &s5, mpnw1)
            mpeq (s5, &s4, mpnw1)
            if Int(s3[3] - s4[3]) < neps { flag = true; break /* goto 100 */ }
        }
        
        if !flag {
            print ("*** MPEGAMMA: Loop end error.")
            mpabrt (36)
        }
        
        // 100 continue
        
        mpexp (s7, &s5, mpnw1)
        mpdiv (s7, s5, &s6, mpnw1)
        mpmul (s6, s4, &s5, mpnw1)
        mpmuld (mplog2con, Double(nn), &s6, mpnw1)
        mpsub (s5, s6, &s0, mpnw1)
        
        //   Restore original precision level.
        
        mproun (&s0, mpnw)
        mpeq (s0, &egamma, mpnw)
        
    } // mpegammaq
    
    static func mpexp (_ a: MPRNumber, _ b: inout MPRNumber, _ mpnw : Int) {
        
    }
    
    //   This computes the exponential function of the MPR number A and returns
    //   the MPR result in B.  Log(2) must be precomputed to at least MPNW words
    //   precision and the stored in the array MPLOG2CON in module MPMODA.
    
    //   This routine uses a modification of the Taylor series for Exp (t):
    
    //   Exp (t) =  (1 + r + r^2 / 2// + r^3 / 3// + r^4 / 4// ...) ^ q * 2 ^ n
    
    //   where the argument T has been reduced to within the closest factor of Log(2).
    //   To further accelerate the series, the reduced argument is divided by 2^NQ.
    //   After convergence of the series, the result is squared NQ times.  NQ = 12
    //   by default.
    
    //   If the precision level MPNW exceeds MPNWX words, this static func calls
    //   MPEXPX instead.  By default, MPNWX = 700 (approx. 10100 digits).
    
//    implicit none
//    integer i, ia, itrmx, j, mpnw, mpnwx, mpnw1, mpnw2, na, nq, nz, n1, n2
//    real (mprknd) t1, t2
//    parameter (itrmx = 1000000, mpnwx = 700)
//    real (mprknd) a(0:), b(0:), f(0:8), s0(0:mpnw+6), &
//    s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), s4(0:mpnw+6)
//    
//    // End of declaration
//    
//    if (mpnw < 4 || a(0] < mpnw + 4 || a(0] < abs (a(2]) + 4 || &
//    b[0] < mpnw + 6) {
//    write (mpldb, 1]
//    1 format ("*** MPEXP: uninitialized or inadequately sized arrays"]
//    mpabrt (99]
//    }
//    
//    ia = sign (1.0, a(2])
//    na = min (int (abs (a(2])), mpnw)
//    mpmdc (a, t1, n1, mpnw)
//    
//    //   Check for overflows and underflows.
//    
//    if (n1 > 30) {
//    if (t1 > 0.0) {
//    write (mpldb, 2)
//    2   format ("*** MPEXP: Argument is too large.")
//    mpabrt (34)
//    } else {
//    b[1] = mpnw
//    b[2] = 0.0
//    b[3] = 0.0
//    goto 130
//    }
//    }
//    
//    t1 = t1 * 2.0 ** n1
//    if (abs (t1) > 1488522236.0) {
//    if (t1 > 0) {
//    write (mpldb, 2)
//    mpabrt (34)
//    } else {
//    b[1] = mpnw
//    b[2] = 0.0
//    b[3] = 0.0
//    goto 130
//    }
//    }
//    
//    //   If the precision level mpnw exceeds mpnwx words, mpexpx.
//    
//    if (mpnw > mpnwx) {
//    mpexpx (a, b, mpnw)
//    goto 130
//    }
//    
//    s0[0] = mpnw + 7
//    s1[0] = mpnw + 7
//    s2[0] = mpnw + 7
//    s3[0] = mpnw + 7
//    s4[0] = mpnw + 7
//    mpnw1 = mpnw + 1
//    
//    //   Set f1 = 1.
//    
//    f[0] = 9.0
//    f[1] = mpnw1
//    f[2] = 1.0
//    f[3] = 0.0
//    f[4] = 1.0
//    f[5] = 0.0
//    f[6] = 0.0
//    
//    //   Check if Log(2) has been precomputed.
//    
//    if (mpnw1 > mplog2con(1)) {
//    write (mpldb, 3) mpnw1
//    3 format ("*** MPLOG: Log(2) must be precomputed to precision",i9," words."/ &
//    "See documentation for details.")
//    mpabrt (35)
//    }
//    
//    //   Compute the reduced argument A" = A - Log(2) * Nint [A / Log(2)].  Save
//    //   NZ = Nint [A / Log(2)] for correcting the exponent of the final result.
//    
//    mpdiv (a, mplog2con, s0, mpnw1)
//    mpnint (s0, s1, mpnw1)
//    mpmdc (s1, t1, n1, mpnw1)
//    nz = nint (t1 * 2.0 ** n1)
//    mpmul (mplog2con, s1, s2, mpnw1)
//    mpsub (a, s2, s0, mpnw1)
//    
//    //   Check if the reduced argument is zero.
//    
//    if (s0[2] == 0.0) {
//    s0[1] = mpnw1
//    s0[2] = 1.0
//    s0[3] = 0.0
//    s0[4] = 1.0
//    s0[5] = 0.0
//    s0[6] = 0.0
//    goto 120
//    }
//    
//    //   Divide the reduced argument by 2 ^ NQ.
//    
//    nq = max (nint (dble (mpnw * mpnbt)** 0.4d0), 1)
//    mpdivd (s0, 2.0**nq, s1, mpnw1)
//    
//    //   Compute Exp using the usual Taylor series.
//    
//    mpeq (f, s2, mpnw1)
//    mpeq (f, s3, mpnw1)
//    mpnw2 =  mpnw1
//    
//    //   The working precision used to compute each term can be linearly reduced
//    //   as the computation proceeds.
//    
//    do j = 1, itrmx
//    t2 = dble (j)
//    mpmul (s2, s1, s0, mpnw2)
//    mpdivd (s0, t2, s2, mpnw2)
//    mpadd (s3, s2, s0, mpnw1)
//    mpeq (s0, s3, mpnw1)
//    
//    //   Check for convergence of the series, and adjust working precision
//    //   for the next term.
//    
//    if (s2[2] == 0.0 || s2[3] < s0[3] - mpnw1) goto 100
//    mpnw2 = min (max (mpnw1 + int (s2[3] - s0[3]) + 1, 4), mpnw1)
//    }
//    
//    write (mpldb, 4)
//    4 format ("*** MPEXP: Iteration limit exceeded.")
//    mpabrt (36)
//    
//    100 continue
//    
//    //   Raise to the (2 ^ NQ)-th power.
//    
//    do j = 1, nq
//    mpmul (s0, s0, s1, mpnw1)
//    mpeq (s1, s0, mpnw1)
//    }
//    
//    //   Multiply by 2 ^ NZ.
//    
//    120 continue
//    
//    mpdmc (1.0, nz, s2, mpnw1)
//    mpmul (s0, s2, s1, mpnw1)
//    
//    //   Restore original precision level.
//    
//    mproun (s1, mpnw)
//    mpeq (s1, b, mpnw)
//    
//    130 continue
//    
//    return
//    end static func mpexp
    
    static func mpexpx ( _ a: MPRNumber, _ b: inout MPRNumber, _ mpnw : Int) {
        
    }
    
    //   This computes the exponential of the MPR number A and returns the MPR
    //   result in B.
    
    //   This routine employs the following Newton iteration, which converges to b:
    
    //     x_{k+1} = x_k + x_k * [a - Log (x_k)]
    
    //   These iterations are performed with a maximum precision level MPNW that
    //   is dynamically changed, approximately doubling with each iteration.
    //   For modest levels of precision, use mpexp.
    
//    implicit none
//    integer i, ia, iq, k, mpnw, mpnw1, mq, na, nit, n1, n2
//    real (mprknd) cl2, t1, t2, mprxx
//    parameter (cl2 = 1.4426950408889633d0, nit = 3, mprxx = 1d-14)
//    real (mprknd) a(0:), b(0:), s0(0:mpnw+6), s1(0:mpnw+6), &
//    s2(0:mpnw+6), s3(0:mpnw+6)
//    
//    // End of declaration
//    
//    if (mpnw < 4 || a(0] < mpnw + 4 || a(0] < abs (a(2]) + 4 || &
//    b[0] < mpnw + 6) {
//    write (mpldb, 1)
//    1 format ("*** MPEXPX: uninitialized or inadequately sized arrays")
//    mpabrt (99)
//    }
//    
//    ia = sign (1.0, a(2])
//    na = min (int (abs (a(2])), mpnw)
//    mpmdc (a, t1, n1, mpnw)
//    
//    //   Check for overflows and underflows.
//    
//    if (n1 > 30) {
//    if (t1 > 0.0) {
//    write (mpldb, 2)
//    2   format ("*** MPEXPX: Argument is too large.")
//    mpabrt (34)
//    } else {
//    b[1) = mpnw
//    b[2) = 0.0
//    b[3) = 0.0
//    goto 130
//    }
//    }
//    
//    t1 = t1 * 2.0 ** n1
//    if (abs (t1) > 1488522236.0) {
//    if (t1 > 0) {
//    write (mpldb, 2)
//    mpabrt (34)
//    } else {
//    b[1) = mpnw
//    b[2) = 0.0
//    b[3) = 0.0
//    goto 130
//    }
//    }
//    
//    s0[0) = mpnw + 7
//    s1[0) = mpnw + 7
//    s2[0) = mpnw + 7
//    s3[0) = mpnw + 7
//    mpnw1 = mpnw + 1
//    
//    //   Check if Pi and Log(2) have been precomputed.
//    
//    if (mpnw1 > mplog2con(1)) {
//    write (mpldb, 3) mpnw1
//    3 format ("*** MPEXPX: Pi and Log(2) must be precomputed to precision",i9," words."/ &
//    "See documentation for details.")
//    mpabrt (53)
//    }
//    
//    //   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.
//    
//    t2 = mpnw1
//    mq = cl2 * log (t2) + 1.0 - mprxx
//    
//    //   Compute initial approximation of Exp (A) (DP accuracy is OK).
//    
//    mpnw1 = 4
//    
//    // The following code (between here and iq = 0) is the equivalent of:
//    //   mpexp (a, s3, mpnw1)
//    
//    mpdiv (a, mplog2con, s0, mpnw1)
//    mpinfr (s0, s1, s2, mpnw1)
//    mpmdc (s1, t1, n1, mpnw1)
//    n1 = int (t1 * 2.0**n1)
//    mpmdc (s2, t2, n2, mpnw1)
//    n2 = min (max (n2, -100), 0)
//    t2 = 2.0 ** (t2 * 2.0**n2)
//    mpdmc (t2, n1, s3, mpnw1)
//    
//    iq = 0
//    
//    //   Perform the Newton-Raphson iteration described above with a dynamically
//    //   changing precision level MPNW (one greater than powers of two).
//    
//    do k = 0, mq
//    if (k > 1) mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
//    
//    100  continue
//    
//    mplogx (s3, s0, mpnw1)
//    mpsub (a, s0, s1, mpnw1)
//    mpmul (s3, s1, s2, mpnw1)
//    mpadd (s3, s2, s1, mpnw1)
//    mpeq (s1, s3, mpnw1)
//    if (k == mq - nit && iq == 0) {
//    iq = 1
//    goto 100
//    }
//    }
//    
//    //   Restore original precision level.
//    
//    mproun (s1, mpnw)
//    mpeq (s1, b, mpnw)
//    
//    130 continue
//    
//    return
//    end static func mpexpx
    
    static func mpinitran ( _ mpnw: Int) {
        
    }
    
    //   This routine computes pi, log(2) sqrt(2)/2, and stores this data in the
    //   proper arrays in module MPFUNA.  MPNW is the largest precision level
    //   (in words) that will be subsequently required for this run at the user level.
    
//    implicit none
//    integer mpnw, nwds, nwds6
//    
//    //   Add three words to mpnw, since many of the routines in this module
//    //   increase the working precision level by one word upon entry.
//    
//    nwds = mpnw + 3
//    
//    //  Compute pi, log(2) and sqrt(2)/2.
//    
//    nwds6 = nwds + 6
//    mplog2con(0) = nwds6
//    mplog2con(1) = 0
//    mppicon(0) = nwds6
//    mppicon(1) = 0
//    
//    mppiq (mppicon, nwds)
//    mplog2q (mppicon, mplog2con, nwds)
//    
//    return
//    end static func mpinitran
    
    static func mplog ( _ a: MPRNumber, _ b: inout MPRNumber, _ mpnw : Int) {
        
    }
    
    //   This computes the natural logarithm of the MPR number A and returns the MPR
    //   result in B.
    
    //   The Taylor series for Log converges much more slowly than that of Exp.
    //   Thus this routine does not employ Taylor series (except if the argument
    //   is extremely close to 1), but instead computes logarithms by solving
    //   Exp (b) = a using the following Newton iteration:
    
    //     x_{k+1} = x_k + [a - Exp (x_k)] / Exp (x_k)
    
    //   These iterations are performed with a maximum precision level MPNW that
    //   is dynamically changed, approximately doubling with each iteration.
    
    //   If the precision level MPNW exceeds MPNWX words, this static func calls
    //   MPLOGX instead.  By default, MPNWX = 30 (approx. 430 digits).
    
//    implicit none
//    integer i, ia, iq, is, itrmax, i1, k, mpnw, mpnwx, mpnw1, mq, na, nit, n1
//    real (mprknd) alt, cl2, rtol, st, tol, t1, t2, mprxx
//    parameter (alt = 0.693147180559945309d0, cl2 = 1.4426950408889633d0, &
//    rtol = 0.5d0**7, itrmax = 1000000, nit = 3, mprxx = 1d-14, mpnwx = 30)
//    real (mprknd) a(0:), b(0:), f1(0:8), s0(0:mpnw+6), &
//    s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), s4(0:mpnw+6)
//    
//    // End of declaration
//    
//    if (mpnw < 4 || a(0) < mpnw + 4 || a(0) < abs (a(2)) + 4 || &
//    b[0) < mpnw + 6) {
//    write (mpldb, 1)
//    1 format ("*** MPLOG: uninitialized or inadequately sized arrays")
//    mpabrt (99)
//    }
//    
//    ia = sign (1.0, a(2))
//    na = min (int (abs (a(2))), mpnw)
//    
//    if (ia .lt. 0 || na == 0) {
//    write (mpldb, 2)
//    2 format ("*** MPLOG: Argument is less than or equal to zero.")
//    mpabrt (50)
//    }
//    
//    //   Check if input is exactly one.
//    
//    if (a(2) == 1.0 && a(3) == 0.0 && a(4) == 1.0) {
//    b[1) = mpnw
//    b[2) = 0.0
//    b[3) = 0.0
//    goto 130
//    }
//    
//    //  If the precision level is more than mpnwx, then mplogx.
//    
//    if (mpnw > mpnwx) {
//    mplogx (a, b, mpnw)
//    goto 130
//    }
//    
//    mpnw1 = mpnw + 1
//    s0[0) = mpnw + 7
//    s1[0) = mpnw + 7
//    s2[0) = mpnw + 7
//    s3[0) = mpnw + 7
//    s4[0) = mpnw + 7
//    
//    f1[0) = 9.0
//    f1[1) = mpnw1
//    f1[2) = 1.0
//    f1[3) = 0.0
//    f1[4) = 1.0
//    f1[5) = 0.0
//    f1[6) = 0.0
//    
//    //   If the argument is sufficiently close to 1, employ a Taylor series.
//    
//    mpsub (a, f1, s0, mpnw1)
//    
//    if (s0[2) == 0.0 || s0[3) <= min (-2.0, - rtol * mpnw1)) {
//    mpeq (s0, s1, mpnw1)
//    mpeq (s1, s2, mpnw1)
//    i1 = 1
//    is = 1
//    tol = s0[3) - mpnw1
//    
//    do i1 = 2, itrmax
//    is = - is
//    st = is * i1
//    mpmul (s1, s2, s3, mpnw1)
//    mpeq (s3, s2, mpnw1)
//    mpdivd (s3, st, s4, mpnw1)
//    mpadd (s0, s4, s3, mpnw1)
//    mpeq (s3, s0, mpnw1)
//    if (s4[2) == 0.0 || s4[3) < tol) goto 120
//    }
//    
//    write (mpldb, 3) itrmax
//    3 format ("*** MPLOG: Iteration limit exceeded",i10)
//    mpabrt (54)
//    }
//    
//    //   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.
//    
//    t2 = mpnw
//    mq = cl2 * log (t2) + 1.0 - mprxx
//    
//    //   Compute initial approximation of Log (A).
//    
//    mpmdc (a, t1, n1, mpnw)
//    t1 = log (t1) + n1 * alt
//    mpdmc (t1, 0, s3, mpnw)
//    mpnw1 = 4
//    iq = 0
//    
//    //   Perform the Newton-Raphson iteration described above with a dynamically
//    //   changing precision level MPNW (one greater than powers of two).
//    
//    do k = 0, mq
//    if (k > 1) mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
//    
//    110  continue
//    
//    mpexp (s3, s0, mpnw1)
//    mpsub (a, s0, s1, mpnw1)
//    mpdiv (s1, s0, s2, mpnw1)
//    mpadd (s3, s2, s1, mpnw1)
//    mpeq (s1, s3, mpnw1)
//    if (k == mq - nit && iq == 0) {
//    iq = 1
//    goto 110
//    }
//    }
//    
//    //   Restore original precision level.
//    
//    120 continue
//    
//    mproun (s3, mpnw)
//    mpeq (s3, b, mpnw)
//    
//    130 continue
//    
//    return
//    end static func mplog
    
    static func mplogx ( _ a: MPRNumber, _ b: inout MPRNumber, _ mpnw : Int) {
        
    }
    
    //   This computes the natural logarithm of the MP number A and returns the MP
    //   result in B.  Pi and Log(2) must be precomputed to at least MPNW words
    //   precision and the stored in the arrays in module MPMODA.
    
    //   This uses the following algorithm, which is due to Salamin and Brent.  If
    //   A is extremely close to 1, use a Taylor series.  Otherwise select n such
    //   that z = a * 2^n is at least 2^m, where m is the number of bits of desired
    //   precision in the result.  Then
    
    //   Log(x) = Pi / [2 AGM (1, 4/x)]
    
    //   For modest precision, or if A is close to 2, use mplog.
    
//    implicit none
//    integer i, ia, is, itrmax, i1, i2, k, mpnw, mpnw1, na, n1, n2
//    real (mprknd) st, rtol, tol, t1, t2, tn
//    parameter (itrmax = 1000000, rtol = 0.5d0**7)
//    real (mprknd) a(0:), b(0:), f1(0:8), f4(0:8), s0[0:mpnw+6), &
//    s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), s4(0:mpnw+6)
//    
//    // End of declaration
//    
//    if (mpnw < 4 || a(0) < mpnw + 4 || a(0) < abs (a(2)) + 4 || &
//    b[0) < mpnw + 6) {
//    write (mpldb, 1)
//    1 format ("*** MPLOGX: uninitialized or inadequately sized arrays")
//    mpabrt (99)
//    }
//    
//    ia = sign (1.0, a(2))
//    na = min (int (abs (a(2))), mpnw)
//    
//    if (ia .lt. 0 || na == 0) {
//    write (mpldb, 2)
//    2 format ("*** MPLOGX: Argument is less than or equal to zero.")
//    mpabrt (50)
//    }
//    
//    //   Check if input is exactly one.
//    
//    if (a(2) == 1.0 && a(3) == 0.0 && a(4) == 1.0) {
//    b[1) = mpnw
//    b[2) = 0.0
//    b[3) = 0.0
//    b[4) = 0.0
//    goto 120
//    }
//    
//    mpnw1 = mpnw + 1
//    s0[0) = mpnw + 7
//    s1[0) = mpnw + 7
//    s2[0) = mpnw + 7
//    s3[0) = mpnw + 7
//    s4[0) = mpnw + 7
//    
//    //   Check if Pi and Log(2) have been precomputed.
//    
//    if (mpnw1 > mplog2con(1)) {
//    write (mpldb, 3) mpnw1
//    3 format ("*** MPLOGX: Pi and Log(2) must be precomputed to precision",i9," words."/ &
//    "See documentation for details.")
//    mpabrt (53)
//    }
//    
//    f1[0) = 9.0
//    f1[1) = mpnw1
//    f1[2) = 1.0
//    f1[3) = 0.0
//    f1[4) = 1.0
//    f1[5) = 0.0
//    f1[6) = 0.0
//    
//    f4[0) = 9.0
//    f4[1) = mpnw1
//    f4[2) = 1.0
//    f4[3) = 0.0
//    f4[4) = 4.0
//    f4[5) = 0.0
//    f4[6) = 0.0
//    
//    //   If the argument is sufficiently close to 1, employ a Taylor series.
//    
//    mpsub (a, f1, s0, mpnw1)
//    
//    if (s0[2) == 0.0 || s0[3) <= min (-2.0, - rtol * mpnw1)) {
//    mpeq (s0, s1, mpnw1)
//    mpeq (s1, s2, mpnw1)
//    i1 = 1
//    is = 1
//    tol = s0[3) - mpnw1
//    
//    do i1 = 2, itrmax
//    is = - is
//    st = is * i1
//    mpmul (s1, s2, s3, mpnw1)
//    mpeq (s3, s2, mpnw1)
//    mpdivd (s3, st, s4, mpnw1)
//    mpadd (s0, s4, s3, mpnw1)
//    mpeq (s3, s0, mpnw1)
//    if (s4[2) == 0.0 || s4[3) < tol) goto 110
//    }
//    
//    write (mpldb, 4) itrmax
//    4 format ("*** MPLOGX: Iteration limit exceeded",i10)
//    mpabrt (54)
//    }
//    
//    //   Multiply the input by a large power of two.
//    
//    mpmdc (a, t1, n1, mpnw1)
//    n2 = mpnbt * (mpnw1 / 2 + 2) - n1
//    tn = n2
//    mpdmc (1.0, n2, s1, mpnw1)
//    mpmul (a, s1, s0, mpnw1)
//    
//    //   Perform AGM iterations.
//    
//    mpeq (f1, s1, mpnw1)
//    mpdiv (f4, s0, s2, mpnw1)
//    mpagmr (s1, s2, s3, mpnw1)
//    
//    //   Compute Pi / (2 * A), where A is the limit of the AGM iterations.
//    
//    mpmuld (s3, 2.0, s0, mpnw1)
//    mpeq (mppicon, s3, mpnw1)
//    mpdiv (s3, s0, s1, mpnw1)
//    
//    //   Subtract TN * Log(2).
//    
//    mpeq (mplog2con, s3, mpnw1)
//    mpmuld (s3, tn, s2, mpnw1)
//    mpsub (s1, s2, s0, mpnw1)
//    
//    110 continue
//    
//    //  Restore original precision level.
//    
//    mproun (s0, mpnw)
//    mpeq (s0, b, mpnw)
//    
//    120 continue
//    
//    return
//    end static func mplogx
    
    static func mplog2q (_ pi: MPRNumber, _ alog2: inout MPRNumber, _ mpnw : Int) {
        
    }
    
    //   This computes log(2) to mpnw words precision, using an algorithm due to Salamin
    //   and Brent:  Select n > 2^m, where m is the number of bits of desired precision
    //   precision in the result.  Then
    
    //   Log(2) = Pi / [2 AGM (1, 4/x)]
    
    //   Where AGM (a, b) denotes the arithmetic-geometric mean:  Set a_0 = a and
    //   b_0 = b, then iterate
    //    a_{k+1} = (a_k + b_k)/2
    //    b_{k+1} = sqrt (a_k * b_k)
    //   until convergence (i.e., until a_k = b_k to available precision).
    
//    implicit none
//    integer i, mpnw, mpnw1, n, n1, n48
//    real (mprknd) cpi, st, t1, t2, tn
//    parameter (cpi = 3.141592653589793238d0)
//    real (mprknd) alog2(0:), f1(0:8), f4(0:8), pi(0:), &
//    s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), s4(0:mpnw+6)
//    
//    // End of declaration
//    
//    if (mpnw < 4 || pi(0) < mpnw + 4 || pi(0) < abs (pi(2)) + 4 || &
//    alog2(0) < mpnw + 6) {
//    write (mpldb, 1)
//    1 format ("*** MPLOG2Q: uninitialized or inadequately sized arrays")
//    mpabrt (99)
//    }
//    
//    //   Define sections of the scratch array.
//    
//    s0[0) = mpnw + 7
//    s1[0) = mpnw + 7
//    s2[0) = mpnw + 7
//    s3[0) = mpnw + 7
//    s4[0) = mpnw + 7
//    mpnw1 = mpnw + 1
//    
//    //   Unless precision is very high, just copy log2 from table.
//    
//    if (mpnw1 <= mplog2con(1)) {
//    mpeq (mplog2con, alog2, mpnw)
//    goto 100
//    }
//    
//    //   Check if Pi has been precomputed.
//    
//    mpmdc (pi, t1, n1, mpnw)
//    if (n1 /= 1 || abs (t1 * 2.0**n1 - cpi) > mprdx || abs (pi(2)) < mpnw) {
//    write (mpldb, 2) mpnw
//    2 format ("*** MPLOG2Q: Pi must be precomputed to precision",i9," words."/ &
//    "See documentation for details.")
//    mpabrt (53)
//    }
//    
//    //   Define sections of the scratch array.
//    
//    s0[0) = mpnw + 7
//    s1[0) = mpnw + 7
//    s2[0) = mpnw + 7
//    s3[0) = mpnw + 7
//    s4[0) = mpnw + 7
//    mpnw1 = mpnw + 1
//    
//    //   Set f1 = 1.
//    
//    f1[0) = 9.0
//    f1[1) = mpnw1
//    f1[2) = 1.0
//    f1[3) = 0.0
//    f1[4) = 1.0
//    f1[5) = 0.0
//    f1[6) = 0.0
//    
//    //   Set f4 = 4.
//    
//    f4[0) = 9.0
//    f4[1) = mpnw1
//    f4[2) = 1.0
//    f4[3) = 0.0
//    f4[4) = 4.0
//    f4[5) = 0.0
//    f4[6) = 0.0
//    
//    //   Set s4 to 2^(n/2), where n is the number of bits desired. n48 = n/48.
//    //   Note that this value can be directly set in the first few words of s4,
//    //   avoiding explicit exponentiation.
//    
//    n = mpnbt * (mpnw1 / 2 + 2)
//    n48 = n / mpnbt
//    s4[1) = mpnw1
//    s4[2) = 1.0
//    s4[3) = n48
//    s4[4) = 1.0
//    s4[5) = 0.0
//    s4[6) = 0.0
//    
//    //   Perform AGM iterations.
//    
//    mpeq (f1, s1, mpnw1)
//    mpdiv (f4, s4, s2, mpnw1)
//    mpagmr (s1, s2, s3, mpnw1)
//    
//    //   Set Log(2) = Pi / (2 * N * S3), where S3 is the limit of the AGM iterations.
//    
//    mpmuld (s3, 2.0 * n, s1, mpnw1)
//    mpdiv (pi, s1, s2, mpnw1)
//    mproun (s2, mpnw)
//    mpeq (s2, alog2, mpnw)
//    
//    100 continue
//    
//    return
//    end static func mplog2q
    
    static func mppiq (_ pi: inout MPRNumber, _ mpnw: Int) {
        
    }
    
    //   This computes Pi to available precision (MPNW mantissa words).
    //   The algorithm that is used for computing Pi, which is due to Salamin
    //   and Brent, is as follows:
    
    //   Set  A_0 = 1,  B_0 = 1/Sqrt(2)  and  D_0 = Sqrt(2) - 1/2.
    
    //   Then from k = 1 iterate the following operations:
    
    //   A_k = 0.5 * (A_{k-1} + B_{k-1})
    //   B_k = Sqrt (A_{k-1} * B_{k-1})
    //   D_k = D_{k-1} - 2^k * (A_k - B_k) ^ 2
    
    //   Then  P_k = (A_k + B_k) ^ 2 / D_k  converges quadratically to Pi.
    //   In other words, each iteration approximately doubles the number of correct
    //   digits, providing all iterations are done with the maximum precision.
    //   The constant cl2 (below) = 1 / log(2) (DP approximation).
    
//    implicit none
//    integer i, k, mpnw, mpnw1, mq
//    real (mprknd) pi(0:), s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), &
//    s3(0:mpnw+6), s4(0:mpnw+6), f(0:8)
//    real (mprknd) cl2, t1
//    parameter (cl2 = 1.4426950408889633d0)
//    
//    // End of declaration
//    
//    if (mpnw < 4 || pi(0) < mpnw + 6) {
//    write (mpldb, 1)
//    1 format ("*** MPPIQ: uninitialized or inadequately sized arrays")
//    mpabrt (99)
//    }
//    
//    s0[0) = mpnw + 7
//    s1[0) = mpnw + 7
//    s2[0) = mpnw + 7
//    s3[0) = mpnw + 7
//    s4[0) = mpnw + 7
//    mpnw1 = mpnw + 1
//    
//    //   Unless precision is very high, just copy pi from table.
//    
//    if (mpnw1 <= mppicon(1)) {
//    mpeq (mppicon, pi, mpnw)
//    goto 100
//    }
//    
//    //   Determine the number of iterations required for the given precision level.
//    //   This formula is good only for this Pi algorithm.
//    
//    t1 = mpnw1 * log10 (mpbdx)
//    mq = cl2 * (log (t1) - 1.0) + 1.0
//    
//    //   Initialize as above.
//    
//    s0[1) = mpnw
//    s0[2) = 1.0
//    s0[3) = 0.0
//    s0[4) = 1.0
//    s0[5) = 0.0
//    s0[6) = 0.0
//    f[0) = 9.0
//    f[1) = mpnw1
//    f[2) = 1.0
//    f[3) = 0.0
//    f[4) = 2.0
//    f[5) = 0.0
//    f[6) = 0.0
//    mpsqrt (f, s2, mpnw1)
//    mpmuld (s2, 0.5d0, s1, mpnw1)
//    f[3) = -1.0
//    f[4) = 0.5d0 * mpbdx
//    mpsub (s2, f, s4, mpnw1)
//    
//    //   Perform iterations as described above.
//    
//    do k = 1, mq
//    mpadd (s0, s1, s2, mpnw1)
//    mpmul (s0, s1, s3, mpnw1)
//    mpsqrt (s3, s1, mpnw1)
//    mpmuld (s2, 0.5d0, s0, mpnw1)
//    mpsub (s0, s1, s2, mpnw1)
//    mpmul (s2, s2, s3, mpnw1)
//    t1 = 2.0 ** k
//    mpmuld (s3, t1, s2, mpnw1)
//    mpsub (s4, s2, s3, mpnw1)
//    mpeq (s3, s4, mpnw1)
//    }
//    
//    //   Complete computation.
//    
//    mpadd (s0, s1, s2, mpnw1)
//    mpmul (s2, s2, s2, mpnw1)
//    mpdiv (s2, s4, s2, mpnw1)
//    mpeq (s2, s0, mpnw1)
//    
//    //   Restore original precision level.
//    
//    mproun (s0, mpnw)
//    mpeq (s0, pi, mpnw)
//    
//    100 continue
//    
//    return
//    end static func mppiq
    
    static func mppower (_ a: MPRNumber, _ b: MPRNumber, _ c: inout MPRNumber, _ mpnw : Int) {
        
    }
    
    //   This computes C = A ^ B, where A, B and C are MPR.  It first checks if
    //   B is the quotient of two integers up to 10^7 in size, in which case it
    //   calls MPNPWR and MPNRTR.  Otherwise it calls MPLOG and MPEXP.
    
//    implicit none
//    integer i, mpnw, mpnw1, n1
//    real (mprknd) a1, a2, a3, a4, a5, a6, q1, t0, t1, t2, t3, mprxx, mprxx2
//    parameter (mprxx = 1.d-14, mprxx2 = 5.d-10)
//    real (mprknd) a(0:), b(0:), c(0:), s0(0:mpnw+6), &
//    s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6)
//    
//    //  End of declaration
//    
//    if (mpnw < 4 || a(0) < mpnw + 4 || b[0) < abs (a(2)) + 4 || &
//    c[0) < mpnw + 6) {
//    write (mpldb, 1)
//    1 format ("*** MPPOWER: uninitialized or inadequately sized arrays")
//    mpabrt (99)
//    }
//    
//    //   Check if A <= 0 (error), or A = 1 or B = 0 or B = 1.
//    
//    if (a(2) <= 0.0) {
//    write (mpldb, 2)
//    2 format ("*** MPPOWER: A^B, where A is less than zero.")
//    mpabrt (61)
//    } else if ((a(2) == 1.0 && a(3) == 0.0 && a(4) == 1.0) &
//    || b[2) == 0.0) {
//    c[1) = mpnw
//    c[2) = 1.0
//    c[3) = 0.0
//    c[4) = 1.0
//    c[5) = 0.0
//    c[6) = 0.0
//    goto 200
//    } else {if (b[2) == 1.0 && b[3) == 0.0 && b[4) == 1.0) {
//    mpeq (a, c, mpnw)
//    goto 200
//    }
//    
//    s0[0) = mpnw + 6
//    s1[0) = mpnw + 6
//    s2[0) = mpnw + 6
//    s3[0) = mpnw + 6
//    
//    //   Check if B is rational using the extended Euclidean algorithm in DP.
//    
//    mpmdc (b, t1, n1, mpnw)
//    
//    if (n1 >= -mpnbt && n1 <= mpnbt) {
//    t0 = abs (t1 * 2.0**n1)
//    t1 = max (t0, 1.0)
//    t2 = min (t0, 1.0)
//    a1 = 1.0
//    a2 = 0.0
//    a3 = 0.0
//    a4 = 1.0
//    
//    do i = 1, 20
//    q1 = aint (t1 / t2)
//    a5 = a1 - q1 * a3
//    a6 = a2 - q1 * a4
//    t3 = t2
//    t2 = t1 - q1 * t2
//    t1 = t3
//    a1 = a3
//    a2 = a4
//    a3 = a5
//    a4 = a6
//    if (t2 < mprxx2) goto 100
//    }
//    
//    goto 110
//    }
//    
//    100 continue
//    
//    a3 = abs (a3)
//    a4 = abs (a4)
//    
//    //  If b = a3/a4 or a4/a3 (except for sign) or then mpnpwr and mpnrtr.
//    
//    if (abs (t0 - a3 / a4) / t0 < mprxx) {
//    a3 = sign (a3, b[2))
//    mpdmc (a3, 0, s0, mpnw)
//    mpdmc (a4, 0, s1, mpnw)
//    mpdiv (s0, s1, s2, mpnw)
//    mpsub (b, s2, s0, mpnw)
//    if (s0[2) == 0.0 || s0[3) < b[3) + 1 - mpnw) {
//    mpnpwr (a, int (a3), s0, mpnw)
//    mpnrtr (s0, int (a4), c, mpnw)
//    goto 200
//    }
//    } else if (abs (t0 - a4 / a3) / t0 < mprxx) {
//    a4 = sign (a4, b[2))
//    mpdmc (a4, 0, s0, mpnw)
//    mpdmc (a3, 0, s1, mpnw)
//    mpdiv (s0, s1, s2, mpnw)
//    mpsub (b, s2, s0, mpnw)
//    if (s0[2) == 0.0 || s0[3) < b[3) + 1 - mpnw) {
//    mpnpwr (a, int (a4), s0, mpnw)
//    mpnrtr (s0, int (a3), c, mpnw)
//    goto 200
//    }
//    }
//    
//    110 continue
//    
//    //  Call mplog and mpexp.
//    
//    mplog (a, s0, mpnw)
//    mpmul (s0, b, s1, mpnw)
//    mpexp (s1, c, mpnw)
//    
//    200 continue
//    
//    return
//    end static func mppower
    
}
