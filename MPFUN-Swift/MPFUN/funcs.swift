//
//  MPFUN.swift
//  MPFUN-Swift
//
//  Created by Mike Griebling on 28 Apr 2019.
//  Copyright © 2019 Computer Inspirations. All rights reserved.
//

import Foundation

extension MPFUN {
    
    static public func sign (_ a: Double, _ b: Double) -> Int {
        let a = Int(abs(a))
        if b < 0 { return -a }
        return a
    }

    static public func aint (_ a: Double) -> Double {
        return Double(Int(a))
    }
    
    public typealias MPReal = Array<Double>
    public typealias MPComplex = Array<Complex64>
    // Note: Complex numbers extend the MPReal array and are positioned at
    // the a[0]th offset.
    
    static func mpabrt (_ ier : Int) {
        //  This routine terminates execution.  Users may wish to replace the
        //  default STOP with a call to a system routine that provides a traceback.
        assertionFailure("*** MPABRT: Execution terminated, error code = \(ier)")
    }
    
    static func mpcabs (_ a : MPReal, _ b : inout MPReal, _ mpnw : Int) {
        
        //   This routine returns the absolute value of the MPC argument A (the
        //   result is of type MPR).
        
        var la, mpnw1 : Int
        var s0 = MPReal(repeating: 0, count: mpnw+7)
        var s1 = s0; var s2 = s0
        
        // End of declaration
        
        la = Int(a[0])
        // lb = Int(b[0])
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || a[la] < abs (a[la+2]) + 4 || b[0] < Double(mpnw + 6) {
            print ("*** MPCABS: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        mpnw1 = mpnw + 1
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        mpmul (a, a, &s0, mpnw1)
        mpmul (MPReal(a[la...]), MPReal(a[la...]), &s1, mpnw1)
        mpadd (s0, s1, &s2, mpnw1)
        mpsqrt (s2, &s0, mpnw1)
        mproun (&s0, mpnw)
        mpeq (s0, &b, mpnw)
    } // mpcabs
    
    static func mpceq (_ a : MPReal, _ b : inout MPReal, _ mpnw : Int) {
        
        //   Sets the MPC number B equal to A.
        
        var ia, la, lb, na : Int
        
        // End of declaration
        
        la = Int(a[0])
        lb = Int(b[0])
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || a[la] < abs (a[la+2]) + 4 || Int(b[0]) < mpnw + 6 || Int(b[lb]) < mpnw + 6) {
            print ("*** MPCEQ: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ia = sign (1.0, a[2])
        na = min (Int (abs (a[2])), mpnw)
        if (na == 0)  {
            b[1] = Double(mpnw)
            b[2] = 0.0
            b[3] = 0.0
            // goto 110
        } else {
            b[1] = Double(mpnw)
            b[2] = Double(sign (Double(na), Double(ia)))
            
            for i in 2...na + 2 {
                b[i+1] = a[i+1]
            }
            
            b[na+4] = 0.0
            b[na+5] = 0.0
        }
        
        // 110 continue
        
        ia = sign (1.0, a[la+2])
        na = min (Int (abs (a[la+2])), mpnw)
        if (na == 0)  {
            b[lb+1] = Double(mpnw)
            b[lb+2] = 0.0
            b[lb+3] = 0.0
            return // goto 120
        }
        b[lb+1] = Double(mpnw)
        b[lb+2] = Double(sign (Double(na), Double(ia)))
        
        for i in 2...na + 2 {
            b[i+lb+1] = a[i+la+1]
        }
        
        b[na+lb+4] = 0.0
        b[na+lb+5] = 0.0
        
        // 120 continue
    } // mpceq
    
    static func mpcadd (_ a : MPReal, _ b : MPReal, _ c : inout MPReal, _ mpnw : Int) {
        
        //   This routine adds the MPC numbers A and B.
        
        var la, lb, lc : Int
        
        // End of declaration
        
        la = Int(a[0])
        lb = Int(b[0])
        lc = Int(c[0])
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || a[la] < abs (a[la+2]) + 4 ||
            b[0] < abs (b[2]) + 4 || b[lb] < abs (b[lb+2]) + 4 ||
            Int(c[0]) < mpnw + 6 || Int(c[lc]) < mpnw + 6) {
            print ("*** MPCADD: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        mpadd (a, b, &c, mpnw)
        var t = MPReal(c[lc...])
        mpadd (MPReal(a[la...]), MPReal(b[lb...]), &t, mpnw)
        c[lc...] = t[0...]
    } //  mpcadd
    
    static func mpcmul (_ a : MPReal, _ b : MPReal, _ c : inout MPReal, _ mpnw : Int) {
        //   This routine multiplies the MPC numbers A and B.
        
        var la, lb, lc, mpnw1 : Int
        var s0 = MPReal(repeating: 0, count: mpnw+7)
        var s1 = s0; var s2 = s0; var s3 = s0
        
        // End of declaration
        
        la = Int(a[0])
        lb = Int(b[0])
        lc = Int(c[0])
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || a[la] < abs (a[la+2]) + 4
            || b[0] < abs (b[2]) + 4 || b[lb] < abs (b[lb+2]) + 4 ||
            Int(c[0]) < mpnw + 6 || Int(c[lc]) < mpnw + 6) {
            print ("*** MPCMUL: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        mpnw1 = mpnw + 1
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        
        mpmul (a, b, &s0, mpnw1)
        mpmul (MPReal(a[la...]), MPReal(b[lb...]), &s1, mpnw1)
        mpsub (s0, s1, &s2, mpnw1)
        mpmul (a, MPReal(b[lb...]), &s0, mpnw1)
        mpmul (MPReal(a[la...]), b, &s1, mpnw1)
        mpadd (s0, s1, &s3, mpnw1)
        
        mproun (&s2, mpnw)
        mproun (&s3, mpnw)
        mpeq (s2, &c, mpnw)
        var t = MPReal(c[lc...])
        mpeq (s3, &t, mpnw)
        c[lc...] = t[0...]
    } // mpcmul
    
    static func mpcdiv (_ a : MPReal, _ b : MPReal, _ c : inout MPReal, _ mpnw : Int) {
        
        //   This routine divides the MPC numbers A and B.
        
        var la, lb, lc, mpnw1 : Int
        var s0 = MPReal(repeating: 0, count: mpnw+7)
        var s1 = s0; var s2 = s0; var s3 = s0; var s4 = s0
        
        // End of declaration
        
        la = Int(a[0])
        lb = Int(b[0])
        lc = Int(c[0])
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || a[la] < abs (a[la+2]) + 4
            || b[0] < abs (b[2]) + 4 || b[lb] < abs (b[lb+2]) + 4 ||
            Int(c[0]) < mpnw + 6 || Int(c[lc]) < mpnw + 6) {
            print ("*** MPCDIV: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        mpnw1 = mpnw + 1
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        s4[0] = Double(mpnw + 7)
        
        mpmul (a, b, &s0, mpnw1)
        mpmul (MPReal(a[la...]), MPReal(b[lb...]), &s1, mpnw1)
        mpadd (s0, s1, &s2, mpnw1)
        mpmul (a, MPReal(b[lb...]), &s0, mpnw1)
        s0[2] = -s0[2]
        mpmul (MPReal(a[la...]), b, &s1, mpnw1)
        mpadd (s0, s1, &s3, mpnw1)
        
        mpmul (b, b, &s0, mpnw1)
        mpmul (MPReal(b[lb...]), MPReal(b[lb...]), &s1, mpnw1)
        mpadd (s0, s1, &s4, mpnw1)
        mpdiv (s2, s4, &s0, mpnw1)
        mpdiv (s3, s4, &s1, mpnw1)
        
        
        mproun (&s0, mpnw)
        mproun (&s1, mpnw)
        mpeq (s0, &c, mpnw)
        var t = MPReal(c[lc...])
        mpeq (s1, &t, mpnw)
        c[lc...] = t[0...]
    } //  mpcdiv
    
    static func mpcsqrt (_ a : MPReal, _ b : inout MPReal, _ mpnw : Int) {
        
        //   This routine returns the square root of the MPC argument A.
        //   The formula is:
        
        //   1/Sqrt[2] * (Sqrt[r + a1] + I * a2 / Sqrt[r + a1])  if a1 >= 0, or
        //   1/Sqrt[2] * (|a2| / Sqrt[r - a1] + I * Sgn[a2] * Sqrt[r - a1]) if a1 < 0,
        
        //   where r = Sqrt[a1^2 + a2^2], and a1 and a2 are the real and imaginary
        //   parts of A.
        
        var lb, la, mpnw1 : Int
        var s0 = MPReal(repeating: 0, count: mpnw+7)
        var s1 = s0; var s2 = s0; var s3 = s0; var s4 = s0
        
        // End of declaration
        
        la = Int(a[0])
        lb = Int(b[0])
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || a[la] < abs (a[la+2]) + 4 || Int(b[0]) < mpnw + 6 || Int(b[lb]) < mpnw + 6) {
            print ("*** MPCSQRT: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        mpnw1 = mpnw + 1
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        s4[0] = Double(mpnw + 7)
        
        mpmul (a, a, &s0, mpnw1)
        mpmul (MPReal(a[la...]), MPReal(a[la...]), &s1, mpnw1)
        mpadd (s0, s1, &s2, mpnw1)
        mpsqrt (s2, &s0, mpnw1)
        
        if (a[2] >= 0.0) {
            mpadd (s0, a, &s1, mpnw1)
            mpsqrt (s1, &s0, mpnw1)
            mpdiv (MPReal(a[la...]), s0, &s1, mpnw1)
        } else {
            mpsub (s0, a, &s2, mpnw1)
            mpsqrt (s2, &s1, mpnw1)
            mpdiv (MPReal(a[la...]), s1, &s0, mpnw1)
            s0[2] = abs (s0[2])
            s1[2] = Double(sign (s1[2], a[la+2]))
        }
        
        mpdmc (0.5, 0, &s3, mpnw1)
        mpsqrt (s3, &s2, mpnw1)
        mpmul (s0, s2, &s3, mpnw1)
        mpmul (s1, s2, &s4, mpnw1)
        
        mproun (&s3, mpnw)
        mproun (&s4, mpnw)
        mpeq (s3, &b, mpnw)
        var t = MPReal(b[lb...])
        mpeq (s4, &t, mpnw)
        b[lb...] = t[0...]
    } // mpcsqrt

    
    static func mpcsub (_ a : MPReal, _ b : MPReal, _ c : inout MPReal, _ mpnw : Int) {
        
        //   This routine subtracts the MPC numbers A and B.
        
        let la = Int(a[0])
        let lb = Int(b[0])
        let lc = Int(c[0])
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || a[la] < abs (a[la+2]) + 4 || b[0] < abs (b[2]) + 4 || b[lb] < abs (b[lb+2]) + 4 ||
            Int(c[0]) < mpnw + 6 || Int(c[lc]) < mpnw + 6) {
            print ("*** MPCSUB: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        mpsub (a, b, &c, mpnw)
        var t = MPReal(c[lc...])
        mpsub (MPReal(a[la...]), MPReal(b[lb...]), &t, mpnw)
        c[lc...] = t[0...]
    } // mpcsub
    
  
    static func mpcpr (_ a : MPReal, _ b : MPReal, _ ic : inout Int, _ mpnw : Int) {
        
        //   This routine compares the MPR numbers A and B and returns in IC the value
        //   -1, 0, or 1 depending on whether A < B, A = B, or A > B.
        //   Note that the first and second words do NOT need to be the same for the
        //   result to be "equal".
        
        var s0 = MPReal(repeating: 0, count:mpnw+6)
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < abs (b[2]) + 4 {
            print ("*** MPCPR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        s0[0] = Double(mpnw + 6)
        mpsub (a, b, &s0, mpnw)
        if s0[2].isSignMinus {
            ic = -1
        } else if s0[2].isZero {
            ic = 0
        } else {
            ic = 1
        }
    } // mpcpr
    
    
    static func mpdmc (_ a : Double, _ n : Int, _ b: inout MPReal, _ mpnw: Int) {
        
        //   This routine converts the DP number A * 2^N to MPR form in B.
        
        //   NOTE however that if A = 0.1D0, for example, { B will NOT be the true
        //   multiprecision equivalent of 1/10, since 0.1d0 is not an exact binary value.
        
        //   Examples of exact binary values (good): 123456789.0, 0.25d0, -5.3125d0.
        //   Examples of inexact binary values (bad): 0.1d0, 1234567.8d0, -3333.3d0.
        
        var n1, n2 : Int
        var aa : Double
        
        // End of declaration
        
        if mpnw < 4 || b[0] < Double(mpnw + 6) {
            print("*** MPDMC: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   Check for zero.
        
        if a == 0 {
            b[1] = Double(mpnw)
            b[2] = 0
            b[3] = 0
            return
        }
        n1 = n / mpnbt
        n2 = n - mpnbt * n1
        aa = abs (a) * pow(2.0, Double(n2))
        
        //   Reduce AA to within 1 and MPBDX.
        
        if aa >= mpbdx {
            
            for k in 1...100 {
                aa = mprdx * aa
                if aa < mpbdx {
                    n1 = n1 + k
                    break //goto 120
                }
            }
            
        } else if aa < 1 {
            
            for k in 1...100 {
                aa = mpbdx * aa
                if aa >= 1 {
                    n1 = n1 - k
                    break //goto 120
                }
            }
            
        }
        
        //   Store successive sections of AA into B.
        
        // 120  continue
        
        b[3] = Double(n1)
        b[4] = aint (aa)
        aa = mpbdx * (aa - b[3+1])
        b[5] = aint (aa)
        aa = mpbdx * (aa - b[4+1])
        b[6] = aint (aa)
        b[7] = 0
        b[8] = 0
        
        var i = 6
        while i >= 3 {
            if b[i+1] != 0 {
                break //goto 140
            }
            i -= 1
        }
        
        // 140  continue
        
        b[1] = Double(mpnw)
        aa = Double(i - 2)
        b[2] = Double(sign (aa, a))
        
        //150 continue
        //return
    } // mpdmc
    
    static func mpdmc40 (_ a : Double, _ n : Int, _ b: inout MPReal, _ mpnw: Int) {
        
        //   This routine converts the DP number A * 2^N to MPR form in B.  In contrast
        //   to mpdmc, this routine only allows 40 significant bits (approximately
        //   12 significant decimal digits) in A.  If more nonzero bits are present,
        //   an error is flagged.
        
        //   Examples of exact binary values (good): 123456789.0, 0.25d0, -5.3125d0.
        //   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.
        
        var t1, t2 : Double
        
        // End of declaration
        
        if mpnw < 4 || b[0] < Double(mpnw + 6) {
            print("*** MPDMC40: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   This convoluted-looking code tests whether A has more than 40
        //   significant bits (actually whether the trailing 13 bits are zero).
        
        t1 = mpb13x * abs (a)
        t2 = abs (abs (a) + t1) - abs (t1)
        if t2 == abs (a) {
            mpdmc (a, n, &b, mpnw)
        } else {
            print("*** MPDMC40: DP value has more than 40 significant bits:",
                  "and thus very likely represents an unintended loss of accuracy.",
                  "Fix the issue, or else use functions mpprodd, mpquotd, mpreald or mpcmplxdc.",
                  "See documentation for details.")
            mpabrt (82)
        }
    } // mpdmc40
    
    
    static func mpeq (_ a : MPReal, _ b : inout MPReal, _ mpnw : Int) {
        
        //   Sets the MPR number B equal to the MPR number A.
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < Double(mpnw + 6) {
            print ("*** MPEQ: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let ia = sign (1, a[2])
        let na = min (Int (abs (a[2])), mpnw)
        if na == 0  {
            b[1] = Double(mpnw)
            b[2] = 0
            b[3] = 0
            return
        }
        b[1] = Double(mpnw)
        b[2] = Double(sign (Double(na), Double(ia)))
        
        for i in 2...na + 2 {
            b[i+1] = a[i+1]
        }
        
        b[na+4] = 0
        b[na+5] = 0
        
        // 110 continue
    } // mpeq
    
    static func mpinfr (_ a: MPReal, _ b: inout MPReal, _ c: inout MPReal, _ mpnw: Int) {
        
        //   Sets B to the integer part of the MPR number A and sets C equal to the
        //   fractional part of A.  Note this is NOT the quite same as the greatest
        //   integer function as often defined in some mathematical books and papers.
        //   Examples:  If A = 1.95, { B = 1., C = 0.95.
        //     If A = -3.25, { B = -3., C = -0.25.
        
        var ia, ma, na, nb, nc : Int
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < Double(mpnw + 6) || c[0] < Double(mpnw + 6) {
            print ("*** MPINFR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   Check if  A  is zero.
        
        ia = sign (1, a[2])
        na = min (Int (abs (a[2])), mpnw)
        ma = Int(a[3])
        if na == 0  {
            b[1] = Double(mpnw)
            b[2] = 0
            b[3] = 0
            c[1] = Double(mpnw)
            c[2] = 0
            c[3] = 0
            return
            //goto 120
        }
        
        if ma >= mpnw - 1 {
            print("*** MPINFR: Argument is too large.")
            mpabrt (40)
        }
        
        //   Place integer part in  B.
        
        nb = min (max (ma + 1, 0), na)
        if nb == 0 {
            b[1] = Double(mpnw)
            b[2] = 0.0
            b[3] = 0.0
        } else {
            b[1] = Double(mpnw)
            b[2] = Double(sign (Double(nb), Double(ia)))
            b[3] = Double(ma)
            b[nb+4] = 0.0
            b[nb+5] = 0.0
            
            for i in 3...nb + 2 {
                b[i+1] = a[i+1]
            }
        }
        
        //   Place fractional part in C.
        
        nc = na - nb
        if nc <= 0 {
            c[1] = Double(mpnw)
            c[2] = 0.0
            c[3] = 0.0
        } else {
            c[1] = Double(mpnw)
            c[2] = Double(sign (Double(nc), Double(ia)))
            c[3] = Double(ma - nb)
            c[nc+4] = 0.0
            c[nc+5] = 0.0
            
            for i in 3...nc + 2 {
                c[i+1] = a[i+nb+1]
            }
        }
        
        //   Fix up results.  B may have trailing zeros and C may have leading zeros.
        
        mproun (&b, mpnw)
        mproun (&c, mpnw)
        
        //120  continue
        // return
    } // mpinfr
    
    static func mpadd (_ a : MPReal, _ b : MPReal, _ c : inout MPReal, _ mpnw : Int) {
        
        //   This routine adds MPR numbers A and B to yield C.
        
        var ia, ib, ish, ixa, ixb, ixd, m1, m2, m3, m4, m5, na, nb, nd, nsh : Int
        var d = MPReal(repeating: 0, count: mpnw+7)
        var db : Double
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs(a[2]) + 4 || b[0] < abs(b[2]) + 4 || c[0] < Double(mpnw+6) {
            print ("*** MPADD: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ia = sign (1, a[2])
        ib = sign (1, b[2])
        na = min (Int (abs (a[2])), mpnw)
        nb = min (Int (abs (b[2])), mpnw)
        
        //   Check for zero inputs.
        
        if na == 0 {
            
            //   A is zero -- the result is B.
            
            c[1] = Double(mpnw)
            c[2] = Double(sign (Double(nb), Double(ib)))
            
            for i in 2...nb+2 {
                c[i+1] = b[i+1]
            }
            
            return
        } else if nb == 0 {
            
            //   B is zero -- the result is A.
            
            c[1] = Double(mpnw)
            c[2] = Double(sign (Double(na), Double(ia)))
            
            for i in 2...na+2 {
                c[i+1] = a[i+1]
            }
            
            return
        }
        
        if ia == ib {
            db = 1
        } else {
            db = -1
        }
        ixa = Int(a[3])
        ixb = Int(b[3])
        ish = ixa - ixb
        
        if ish >= 0 {
            
            //   A has greater exponent than B, so B must be shifted to the right.
            
            m1 = min (na, ish)
            m2 = min (na, nb + ish)
            m3 = na
            m4 = min (max (na, ish), mpnw + 1)
            m5 = min (max (na, nb + ish), mpnw + 1)
            
            if m1 >= 1 {
                for i in 1...m1 {
                    d[i+3] = a[i+3]
                }
            }
            
            if m2 >= m1+1 {
                for i in m1+1...m2 {
                    d[i+3] = a[i+3] + db * b[i+2-ish+1]
                }
            }
            
            if m3 >= m2+1 {
                for i in m2+1...m3 {
                    d[i+3] = a[i+3]
                }
            }
            
            if m2 >= m3+1 {
                for i in m3+1...m4 {
                    d[i+3] = 0
                }
            }
            
            if m5 >= m4+1 {
                for i in m4+1...m5 {
                    d[i+3] = db * b[i+2-ish+1]
                }
            }
            
            nd = m5
            ixd = ixa
            d[nd+4] = 0
            d[nd+5] = 0
        } else {
            
            //   B has greater exponent than A, so A must be shifted to the right.
            
            nsh = -ish
            m1 = min (nb, nsh)
            m2 = min (nb, na + nsh)
            m3 = nb
            m4 = min (max (nb, nsh), mpnw + 1)
            m5 = min (max (nb, na + nsh), mpnw + 1)
            
            if m1 >= 1 {
                for i in 1...m1 {
                    d[i+3] = db * b[i+3]
                }
            }
            
            if m2 >= m1+1 {
                for i in m1+1...m2 {
                    d[i+3] = a[i+2-nsh+1] + db * b[i+3]
                }
            }
            
            if m3 >= m2+1 {
                for i in m2+1...m3 {
                    d[i+3] = db * b[i+3]
                }
            }
            
            if m2 >= m3+1 {
                for i in m3+1...m4 {
                    d[i+3] = 0
                }
            }
            
            if m5 >= m4+1 {
                for i in m4+1...m5 {
                    d[i+3] = a[i+2-nsh+1]
                }
            }
            
            nd = m5
            ixd = ixb
            d[nd+4] = 0
            d[nd+5] = 0
        }
        
        //   Call mpnorm to fix up result and store in c.
        
        d[0] = Double(mpnw + 7)
        d[1] = Double(mpnw)
        d[2] = Double(sign (Double(nd), Double(ia)))
        d[3] = Double(ixd)
        
        mpnorm (d, &c, mpnw)
    }
    
    static func mpmul (_ a : MPReal, _ b : MPReal, _ c : inout MPReal, _ mpnw : Int) {
        
        //   This routine multiplies MPR numbers A and B to yield C.
        
        //   This routine returns up to MPNW mantissa words of the product.  If the
        //   complete double-long product of A and B is desired (for example in large
        //   integer applications), { MPNW must be at least as large as the sum of
        //   the mantissa lengths of A and B.  In other words, if the precision levels
        //   of A and B are both 64 words, { MPNW must be at least 128 words to
        //   produce the complete double-long product in C.
        
        var i1, i2, j3, ia, ib, na, nb, nc, n2 : Int
        let mpnwx = 200
        var a1, a2, c1, c2, dc, d2, t1, t2 : Double
        var d = MPReal(repeating: 0, count: mpnw+6)
        var b1 = d; var b2 = d
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < abs (b[2]) + 4 || c[0] < Double(mpnw + 6) {
            print ("*** MPMUL: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ia = sign (1.0, a[2])
        ib = sign (1.0, b[2])
        na = min (Int (abs (a[2])), mpnw)
        nb = min (Int (abs (b[2])), mpnw)
        nc = min (na + nb, mpnw)
        
        if na == 0 || nb == 0 {
            
            //   One of the inputs is zero -- result is zero.
            
            c[1] = Double(mpnw)
            c[2] = 0.0
            c[3] = 0.0
            return
        }
        
        if na == 1 && a[4] == 1.0 {
            
            //   A is 1 or -1 -- result is B or -B.
            
            c[1] = Double(mpnw)
            c[2] = Double(sign (Double(nb), Double(ia * ib)))
            c[3] = a[3] + b[3]
            
            for i in 3...nb + 2 {
                c[i+1] = b[i+1]
            }
            return
            
        } else if nb == 1 && b[4] == 1.0 {
            
            //   B is 1 or -1 -- result is A or -A.
            
            c[1] = Double(mpnw)
            c[2] = Double(sign (Double(na), Double(ia * ib)))
            c[3] = a[3] + b[3]
            
            for i in 3...na + 2 {
                c[i+1] = a[i+1]
            }
            return
        }
        
        if na > mpnwx && nb > mpnwx {
            
            //   Precision levels of both arguments are higher than mpnwx, so call mpmulx.
            
            mpmulx (a, b, &c, mpnw)
            return
        }
        
        d2 = a[3] + b[3]
        d[0] = Double(mpnw + 6)
        d[1] = Double(mpnw)
        
        // probably not needed since we already zeroed d
        for i in 1...nc + 4 {
            d[i+1] = 0.0
        }
        
        for i in 0...nb + 4 {
            b1[i] = mpb24x * aint (mpr24x * b[i+1])
            b2[i] = b[i+1] - b1[i]
        }
        
        //   Perform ordinary long multiplication algorithm.  Accumulate at most MPNW+4
        //   mantissa words of the product.
        
        for j in 3...na + 2 {
            a1 = mpb24x * aint (mpr24x * a[j+1])
            a2 = a[j+1] - a1
            j3 = j - 3
            n2 = min (nb + 2, mpnw + 4 - j3)
            
            for i in 3...n2 {
                dc = a1 * b2[i] + a2 * b1[i]
                c1 = mpbdx * aint (mprdx * dc)
                c2 = dc - c1
                d[i+j3] = d[i+j3] + mprdx * (a1 * b1[i] + c1)
                d[i+j3+1] = d[i+j3+1] + a2 * b2[i] + c2
            }
            
            //   Release carries periodically to avoid overflowing the exact integer
            //   capacity of double precision floating point words in D.
            
            if (j - 2) % mpnpr == 0 {
                i1 = max (3, j - mpnpr)
                i2 = n2 + j3
                
                for i in i1...i2 {
                    t1 = d[i+1]
                    t2 = Double(Int(mprdx * t1))
                    d[i+1] = t1 - mpbdx * t2
                    d[i] += t2
                }
            }
        }
        
        //   If D(3) is nonzero, shift the result one cell right.
        
        if d[3] != 0.0 {
            d2 = d2 + 1.0
            
            for i in stride(from:nc + 4, through:3, by:-1) {
                d[i+1] = d[i]
            }
        }
        d[2] = Double(sign (Double(nc), Double(ia * ib)))
        d[3] = d2
        
        //   Fix up result, since some words may be negative or exceed MPBDX.
        
        mpnorm (d, &c, mpnw)
        
        //    200 continue
        
        //    return
    } // mpmul
    
    static func mpmuld (_ a : MPReal, _ b : Double, _ c : inout MPReal, _ mpnw : Int) {
        
        //   This routine multiplies the MPR number A by the DP number B to yield C.
        
        //   Note, however, that if B = 0.1D0 for example (or any other value that
        //   is not either a whole number or exact binary fraction), { C will NOT
        //   be the true multiprecision product A * B.  This is because the double
        //   precision value 0.1d0 is only a 15-digit approximation of 1/10, and thus
        //   the product C will only be good to 15 digits or so.
        
        //   Examples of exact binary values (good): 35.0, 395.5, 0.125, -5.3125.
        //   Examples of inexact binary values (bad):  0.1, 0.8, -2.95, 33.3.
        
        var ia, ib, na, n1 : Int
        var a1, a2, bb, b1, b2, c1, c2, dc : Double
        var d = MPReal(repeating: 0, count: mpnw+6)
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || c[0] < Double(mpnw + 6) {
            print("*** MPMULD: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   Check for zero inputs.
        
        ia = sign (1.0, a[2])
        na = min (Int(abs (a[2])), mpnw)
        ib = sign (1.0, b)
        if na == 0 || b == 0.0 {
            c[1] = Double(mpnw)
            c[2] = 0.0
            c[3] = 0.0
            return
        }
        bb = abs (b)
        n1 = 0
        
        //   Reduce BB to within 1 and MPBDX.
        
        if bb >= mpbdx {
            for k in 1...100 {
                bb = mprdx * bb
                if bb < mpbdx {
                    n1 = n1 + k
                    break // goto 120
                }
            }
        } else if bb < 1.0 {
            for k in 1...100 {
                bb = mpbdx * bb
                if bb >= 1.0 {
                    n1 = n1 - k
                    break // goto 120
                }
            }
        }
        
        //   If B cannot be represented exactly in a single mantissa word, use MPMUL.
        
        // 120  continue
        
        if bb != aint (bb) {
            bb = Double(sign (bb, b))
            d[0] = Double(mpnw + 6)
            d[1] = Double(mpnw)
            mpdmc (bb, n1 * mpnbt, &d, mpnw)
            mpmul (a, d, &c, mpnw)
            return
        }
        
        //   Perform short multiply operation.
        
        b1 = mpb24x * aint (mpr24x * bb)
        b2 = bb - b1
        d[0] = Double(mpnw + 6)
        d[1] = Double(mpnw)
        d[2] = Double(sign (Double(na + 1), Double(ia * ib)))
        
        for i in 2...min (na + 5, mpnw + 4) {
            d[i+1] = 0.0
        }
        
        for i in 3...na + 2 {
            a1 = mpb24x * aint (mpr24x * a[i+1])
            a2 = a[i+1] - a1
            dc = a1 * b2 + a2 * b1
            c1 = mpbdx * aint (mprdx * dc)
            c2 = dc - c1
            d[i] = d[i] + mprdx * (a1 * b1 + c1)
            d[i+1] = a2 * b2 + c2
        }
        
        //   If d(3) is nonzero, shift all words right by one.
        
        if d[3] > 0.0 {
            for i in stride(from:na + 3, through:3, by:-1) {
                d[i+1] = d[i]
            }
            
            d[3] = a[3] + Double(n1 + 1)
        } else {
            d[3] = a[3] + Double(n1)
        }
        
        //   Fix up the result.
        
        mpnorm (d, &c, mpnw)
        
        d[3] = a[3] + Double(n1)
        
        //    140 continue
        
        //    return
    } // mpmuld
    
    static func mpmuld40 (_ a : MPReal, _ b : Double, _ c : inout MPReal, _ mpnw : Int) {
        
        //   This routine multiples the MP number A by the DP number B to yield C.
        //   In contrast to mpmuld, this routine only allows 40 significant bits
        //   (approximately 12 significant decimal digits) in B.  If more nonzero bits
        //   are present, an error is flagged.
        
        //   Examples of exact binary values (good): 123456789.0, 0.25d0, -5.3125d0.
        //   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.
        
        var t1, t2 : Double
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || c[0] < Double(mpnw + 6) {
            print("*** MPMULD40: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   This convoluted-looking code tests whether B has more than 40
        //   significant bits (actually whether the trailing 13 bits are zero).
        
        t1 = mpb13x * abs (b)
        t2 = abs (abs (b) + t1) - abs (t1)
        if t2 == abs (b) {
            mpmuld (a, b, &c, mpnw)
        } else {
           print ("*** MPMULD40: DP value has more than 40 significant bits: \(b)",
                "and thus very likely represents an unintended loss of accuracy.",
                    "Fix the issue, or else use functions mpprodd, mpquotd, mpreald or mpcmplxdc.",
                "See documentation for details.")
            mpabrt (83)
        }
        
        return
    } // mpmuld40
    
    static func mpdiv (_ a: MPReal, _ b: MPReal, _ c: inout MPReal, _ mpnw: Int) {
        
        //   This divides the MPR number A by the MP number B to yield C.
        
        var i2, i3, ia, ib, ij, iss, j, j3, md, na, nb, nc : Int
        let mpnwx = 200
        var a1, a2, b1, b2, c1, c2, dc, rb, t0, t1, t2 : Double
        var d = MPReal(repeating: 0, count: mpnw+7)
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < abs (b[2]) + 4 || c[0] < Double(mpnw + 6) {
            print("*** MPDIV: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ia = sign (1.0, a[2])
        ib = sign (1.0, b[2])
        na = min (Int (abs (a[2])), mpnw)
        nb = min (Int (abs (b[2])), mpnw)
        
        //   Check if dividend is zero.
        
        if na == 0 {
            c[1] = Double(mpnw)
            c[2] = 0.0
            c[3] = 0.0
            c[4] = 0.0
            return
            //goto 190
        }
        
        //   Check if divisor is zero.
        
        if nb == 0 {
            print ("*** MPDIV: Divisor is zero.")
            mpabrt (31)
        }
        
        if na > mpnwx && nb > mpnwx {
            
            //   Precision is very high, so call mpdivx.
            
            mpdivx (a, b, &c, mpnw)
            return
            // goto 200
        }
        
        //   Compute double precision approximation to trial divisor and its reciprocal.
        
        t0 = mpbdx * b[4]
        if nb >= 2 { t0 = t0 + b[5] }
        if nb >= 3 { t0 = t0 + mprdx * b[6] }
        rb = 1.0 / t0
        
        md = min (na + nb, mpnw)
        d[0] = Double(mpnw + 6)
        d[1] = Double(mpnw)
        d[2] = 0.0
        
        for i in 2...na + 1 {
            d[i+1] = a[i+2]
        }
        
        for i in na + 2...md + 4 {
            d[i+1] = 0.0
        }
        
        //   Perform ordinary long division algorithm.
        
        for j in 2...mpnw + 3 {
            
            //   Compute trial dividend and trial quotient term.
            
            t1 = mpbdx * mpbdx * d[j] + mpbdx * d[j+1] + d[j+2]
            if j <= mpnw + 2 { t1 = t1 + mprdx * d[j+3] }
            t0 = aint (rb * t1)
            
            //   Split trial quotient term into high and low order halves, 24-bits each.
            
            a1 = mpb24x * aint (mpr24x * t0)
            a2 = t0 - a1
            j3 = j - 3
            i2 = min (nb, mpnw + 2 - j3) + 2
            ij = i2 + j3
            
            //   Multiply trial quotient term by each term of divisor, saving high and low-
            //   order parts into appropriate terms of the running dividend.
            
            for i in 3...i2 {
                i3 = i + j3
                b1 = mpb24x * aint (mpr24x * b[i+1])
                b2 = b[i+1] - b1
                dc = a1 * b2 + a2 * b1
                c1 = mpbdx * aint (mprdx * dc)
                c2 = dc - c1
                d[i3] = d[i3] - mprdx * (a1 * b1 + c1)
                d[i3+1] = d[i3+1] - a2 * b2 - c2
            }
            
            //   Release carries periodically to avoid overflowing the exact integer
            //   capacity of double precision floating point words in D.
            
            if j - 1 % mpnpr == 0 {
                for i in j + 1...ij {
                    t1 = d[i]
                    t2 = aint (mprdx * t1)
                    d[i] = t1 - mpbdx * t2
                    d[i-1] = d[i-1] + t2
                }
            }
            
            d[j+1] = d[j+1] + mpbdx * d[j]
            d[j] = t0
            
            //  Continue computing quotient terms past the end of the input dividend, until
            //  trial running dividend is zero.
            
            if j >= na + 2 {
                if ij <= mpnw + 1 { d[ij+4] = 0.0 }
            }
        }
        
        //   Final bookkeeping.
        
        j = mpnw + 3
        
        // 170 continue
        
        d[j+1] = 0.0
        if d[2] == 0.0 {
            iss = 1
        } else {
            iss = 2
        }
        nc = min (j - 1, mpnw)
        d[nc+4] = 0.0
        d[nc+5] = 0.0
        
        for i in stride(from: j + 1, through: 3, by: -1) {
            d[i+1] = d[i-iss+1]
        }
        
        d[0] = Double(mpnw + 6)
        d[1] = Double(mpnw)
        d[2] = Double(sign (Double(nc), Double(ia * ib)))
        d[3] = a[3] - b[3] + Double(iss - 2)
        c[1] = Double(mpnw)
        
        //   Call mpnorm to fix up any remaining bugs and perform rounding.
        
        mpnorm (d, &c, mpnw)
        
        //190 continue
        //200 continue
        //return
    } // mpdiv
    
    static func mpdivd (_ a: MPReal, _ b: Double, _ c: inout MPReal, _ mpnw: Int) {
        
        //   This routine divides the MPR number A by the DP number B to yield C.
        
        //   NOTE however that if A = 0.1D0, for example, { C will NOT be the true
        //   multiprecision equivalent of the quotient, since 0.1d0 is not an exact
        //   binary value.
        
        //   Examples of exact binary values (good): 123456789.0, 0.25d0, -5.3125d0.
        //   Examples of inexact binary values (bad): 0.1d0, 1234567.8d0, -3333.3d0.
        
        var iss, ia, ib, j, md, na, nc, n1 : Int
        var a1, a2, bb, b1, b2, c1, c2, dc, rb, t0, t1 : Double
        var d = MPReal(repeating: 0, count: mpnw+6)
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || c[0] < Double(mpnw + 6) {
            print("*** MPDIVD: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ia = sign (1.0, a[2])
        na = min (Int (abs (a[2])), mpnw)
        ib = sign (1.0, b)
        
        //   Check if dividend is zero.
        
        if na == 0 {
            c[1] = Double(mpnw)
            c[2] = 0.0
            c[3] = 0.0
            return
        }
        
        //   Check if divisor is zero.
        
        if b == 0.0 {
            print ("*** MPDIVD: Divisor is zero.")
            mpabrt (32)
        }
        
        n1 = 0
        bb = abs (b)
        
        //   Reduce BB to within 1 and MPBDX.
        
        if bb >= mpbdx {
            
            for k in 1...100 {
                bb = mprdx * bb
                if bb < mpbdx {
                    n1 = n1 + k
                    break // goto 120
                }
            }
            
        } else if bb < 1.0 {
            
            for k in 1...100 {
                bb = mpbdx * bb
                if bb >= 1.0 {
                    n1 = n1 - k
                    break // goto 120
                }
            }
            
        }
        
        // 120 continue
        
        //   If B cannot be represented exactly in a single mantissa word, use MPDIV.
        
        if bb != aint (bb) {
            bb = Double(sign (bb, b))
            d[0] = Double(mpnw + 6)
            d[1] = Double(mpnw)
            mpdmc (bb, n1 * mpnbt, &d, mpnw)
            mpdiv (a, d, &c, mpnw)
            return
            // goto 190
        }
        
        //   Compute double precision approximation to trial divisor and its reciprocal.
        
        t0 = mpbdx * bb
        rb = 1.0 / t0
        b1 = mpb24x * aint (mpr24x * bb)
        b2 = bb - b1
        
        md = min (na + 1, mpnw)
        d[0] = Double(mpnw + 6)
        d[1] = Double(mpnw)
        d[2] = 0.0
        
        for i in 2...na + 1 {
            d[i+1] = a[i+2]
        }
        
        for i in na + 2...md + 4 {
            d[i+1] = 0.0
        }
        
        //   Perform ordinary short division algorithm.
        
        for j in 2...mpnw + 3 {
            
            //   Compute trial dividend and trial quotient term.
            
            t1 = mpbdx*mpbdx * d[j] + mpbdx * d[j+1] + d[j+2]
            if j <= mpnw + 2 { t1 = t1 + mprdx * d[j+3] }
            t0 = aint (rb * t1)
            
            //   Split trial quotient term into high and low order halves, 24-bits each.
            
            a1 = mpb24x * aint (mpr24x * t0)
            a2 = t0 - a1
            
            //   Multiply trial quotient term by each term of divisor, saving high and low-
            //   order parts into appropriate terms of the running dividend.
            
            dc = a1 * b2 + a2 * b1
            c1 = mpbdx * aint (mprdx * dc)
            c2 = dc - c1
            d[j] = d[j] - mprdx * (a1 * b1 + c1)
            d[j+1] = d[j+1] - a2 * b2 - c2
            d[j+1] = d[j+1] + mpbdx * d[j]
            d[j] = t0
            
            //  Continue computing quotient terms past the end of the input dividend, until
            //  trial running dividend is zero.
            
            if j >= na + 2 {
                if j <= mpnw + 1 { d[j+4] = 0.0 }
            }
        }
        
        //   Final bookkeeping.
        
        j = mpnw + 3
        
        // 170 continue
        
        d[j+1] = 0.0
        if d[2] == 0.0 {
            iss = 1
        } else {
            iss = 2
        }
        nc = min (j - 1, mpnw)
        d[nc+4] = 0.0
        d[nc+5] = 0.0
        
        for i in stride(from: j+1, through: 3, by: -1) {
            d[i+1] = d[i-iss+1]
        }
        
        d[0] = Double(mpnw + 6)
        d[1] = Double(mpnw)
        d[2] = Double(sign (Double(nc), Double(ia * ib)))
        d[3] = a[3] - Double(n1 + iss - 2)
        c[1] = Double(mpnw)
        
        //   Call mpnorm to fix up any remaining bugs and perform rounding.
        
        mpnorm (d, &c, mpnw)
        
        //190 continue
        //return
    } // mpdivd
    
    static func mpdivd40 (_ a: MPReal, _ b: Double, _ c: inout MPReal, _ mpnw: Int) {
        
        //   This routine divides the MPR number A by the DP number B to yield C.
        //   In contrast to mpdivd, this routine only allows 40 significant bits
        //   (approximately 12 significant decimal digits) in B.  If more nonzero bits
        //   are present, an error is flagged.
        
        //   Examples of exact binary values (good): 123456789.0, 0.25d0, -5.3125d0.
        //   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.
        
        var t1, t2 : Double
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || c[0] < Double(mpnw + 6) {
            print("*** MPDIVD40: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   This convoluted-looking code tests whether B has more than 40
        //   significant bits (actually whether the trailing 13 bits are zero).
        
        t1 = mpb13x * abs (b)
        t2 = abs (abs (b) + t1) - abs (t1)
        if t2 == abs (b) {
            mpdivd (a, b, &c, mpnw)
        } else {
            print ("*** MPDIVD40: DP value has more than 40 significant bits:",
                   "and thus very likely represents an unintended loss of accuracy.",
                   "Fix the issue, or else use functions mpprodd, mpquotd, mpreald or mpcmplxdc.",
                   "See documentation for details.")
            mpabrt (81)
        }
    } // mpdivd40
    
    static func mpmdc (_ a: MPReal, _ b: inout Double, _ n: inout Int, _ mpnw: Int) {
        
        //   This returns a DP approximation the MPR number A in the form B * 2^n.
        
        var n, na : Int
        var aa : Double
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 {
            print ("*** MPMDC: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        if a[2] == 0  {
            b = 0
            n = 0
            return
        }
        
        na = Int(abs (a[2]))
        aa = a[4]
        if na >= 2 { aa = aa + mprdx * a[5] }
        if na >= 3 { aa = aa + mprx2 * a[6] }
        if na >= 4 { aa = aa + mprdx * mprx2 * a[7] }
        
        n = mpnbt * Int(a[3])
        b = Double(sign (aa, Double(a[2])))
        
        //   Reduce b to within 1 and 2.
        
        na = Int(log (abs (b)) / log (2.0) + mprdx)
        b = b / pow(2.0, Double(na))
        n = n + na
        if abs (b) < 1 {
            b = 2 * b
            n = n - 1
        } else if abs (b) > 2 {
            b = 0.5 * b
            n = n + 1
        }
        
        //   100  continue
        //   return
    } // mpmdc
    
    static func mpnint (_ a: MPReal, _ b: inout MPReal, _ mpnw: Int) {
        
        //   This sets B to the nearest integer to the MPR number A.
        //   Examples:  If A = 1.49, B = 1.; if A = 3.5, B = 4; if A = -2.5, B = -3.
        
        var ia, ma, na : Int
        var s0 = MPReal(repeating: 0, count: mpnw+6)
        var s1 = s0
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < Double(mpnw + 6) {
            print ("*** MPNINT: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ia = sign (1.0, a[2])
        na = min (Int (abs (a[2])), mpnw)
        ma = Int(a[3])
        if na == 0 {
            
            //   A is zero -- result is zero.
            
            b[1] = Double(mpnw)
            b[2] = 0.0
            b[3] = 0.0
            return
        }
        
        if ma >= mpnw {
            
            //   A cannot be represented exactly as an integer.
            
            print ("*** MPNINT: Argument is too large.")
            mpabrt (56)
        }
        
        //   Add or subtract 1/2 from the input, depending on its sign, {
        //   return the greatest integer.
        
        s0[0] = Double(mpnw + 6)
        s1[0] = Double(mpnw + 6)
        
        mpdmc (0.5, 0, &s0, mpnw)
        if ia == 1 {
            mpadd (a, s0, &s1, mpnw)
        } else {
            mpsub (a, s0, &s1, mpnw)
        }
        mpinfr (s1, &b, &s0, mpnw)
        
        //110 continue
        //return
    } // mpnint
    
    static func mpnorm (_ d: MPReal, _ a: inout MPReal, _ mpnw : Int) {
        
        //   This converts the MP number in array D to the standard normalized form
        //   in A.
        
        //   MPNORM assumes that two extra mantissa words are input at the end of D.
        //   This reduces precision loss when it is necessary to shift the result to
        //   the left. All words up to index A(2) + 5 in A *must* have data, even if 0.
        
        var ia, na, n4 : Int
        var a2, t1, t2, t3 : Double
        var d = d  // create a local mutable copy
        
        // End of declaration
        
        if mpnw < 4 || d[0] < abs(d[2]) + 4 || a[0] < Double(mpnw + 6) {
            print ("*** MPNORM: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ia = sign (1, d[2])
        na = min (Int(abs(d[2])), mpnw)
        if na == 0  {
            a[1] = Double(mpnw)
            a[2] = 0
            a[3] = 0
            return
        }
        n4 = na + 4
        a2 = d[3]
        d[3] = 0
        
        repeat {
            t1 = 0
            
            for i in stride(from: n4, through: 3, by: -1) {
                t3 = t1 + d[i+1]
                t2 = mprdx * t3
                t1 = Double(Int(t2))
                if t2 < 0 && t1 != t2 { t1 = t1 - 1 }
                d[i+1] = t3 - t1 * mpbdx
            }
            
            d[3] = d[3] + t1
            
            if d[3] < 0 {
                
                //   D(2) is negative -- negate all words and re-normalize.
                
                ia = -ia
                d[4] = d[4] + mpbdx * d[3]
                d[3] = 0
                
                for i in 2...n4 {
                    d[i+1] = -d[i+1]
                }
            } else if d[3] > 0 {
                
                //   The fixup loops above "spilled" a nonzero number into D(2).  Shift the
                //   entire number right one cell.  The exponent and length of the result
                //   are increased by one.
                
                for i in stride(from: n4, through: 3, by: -1) {
                    a[i+1] = d[i]
                }
                
                na = min (na + 1, mpnw)
                a2 = a2 + 1
                break
            } else {
                for i in 3...n4 {
                    a[i+1] = d[i+1]
                }
                break
            }
        } while true
        
        //   Perform rounding and truncation.
        
        a[1] = Double(mpnw)
        a[2] = Double(sign (Double(na), Double(ia)))
        a[3] = a2
        
        mproun (&a, mpnw)
    } // mpnorm
    
    
    static func mpnpwr (_ a : MPReal, _ n : Int, _ b : inout MPReal, _ mpnw : Int) {
        
        //   This computes the N-th power of the MPR number A and returns the result
        //   in B.  When N is zero, 1 is returned.  When N is negative, the reciprocal
        //   of A ^ |N| is returned.
        
        var i, j, k, kk, kn, k0, k1, k2, mn, mpnw1, na, nn, nws : Int
        var t1 : Double
        let cl2 = 1.4426950408889633e0; let mprxx = 1e-14
        var s0 = MPReal(repeating: 0, count: mpnw+7)
        var s1 = s0
        var s2 = s0
        
        func reciprocal() {
            if n < 0 {
                mpdmc (1, 0, &s1, mpnw1)
                mpdiv (s1, s2, &s0, mpnw1)
                mpeq (s0, &s2, mpnw1)
            }
            
            //   Restore original precision level.
            
            mproun (&s2, mpnw)
            mpeq (s2, &b, mpnw)
        }
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < Double(mpnw + 6) {
            print ("*** MPNPWR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        na = min(Int(abs (a[2])), mpnw)
        if na == 0 {
            if n >= 0 {
                b[1] = Double(mpnw)
                b[2] = 0
                b[3] = 0
                return
            } else {
                print("*** MPNPWR: Argument is zero and N is negative or zero.")
                mpabrt (57)
            }
        }
        
        mpnw1 = mpnw + 1
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        
        nn = abs (n)
        if nn == 0 {
            mpdmc (1, 0, &b, mpnw)
            return
        } else if nn == 1 {
            mpeq (a, &s2, mpnw1)
            reciprocal()
            return
        } else if nn == 2 {
            mpmul (a, a, &s2, mpnw1)
            reciprocal()
            return
        }
        
        //   Determine the least integer MN such that 2 ^ MN > NN.
        
        t1 = Double(nn)
        mn = Int(cl2 * log (t1) + 1 + mprxx)
        mpdmc (1, 0, &s2, mpnw1)
        mpeq (a, &s0, mpnw1)
        kn = nn
        
        //   Compute B ^ N using the binary rule for exponentiation.
        
        for j in 1...mn {
            kk = kn / 2
            if kn != 2 * kk {
                mpmul (s2, s0, &s1, mpnw1)
                mpeq (s1, &s2, mpnw1)
            }
            kn = kk
            if j < mn {
                mpmul (s0, s0, &s1, mpnw1)
                mpeq (s1, &s0, mpnw1)
            }
        }
        
        //   Compute reciprocal if N is negative.
        
        reciprocal()
    } // mpnpwr
    
    static func mpnrtr (_ a : MPReal, _ n : Int, _ b : inout MPReal, _ mpnw : Int) {
        
        //   This computes the N-th root of the MPR number A and returns result in B.
        //   N must be at least one and must not exceed 2 ^ 30.
        
        //   This subroutine employs the following Newton-Raphson iteration, which
        //   converges to A ^ (-1/N):
        
        //    X_{k+1} = X_k + (X_k / N) * (1 - A * X_k^N)
        
        //   The reciprocal of the final approximation to A ^ (-1/N) is the N-th root.
        //   These iterations are performed with a maximum precision level MPNW that
        //   is dynamically changed, approximately doubling with each iteration.
        
        //   When N is large and A is very near one, the following binomial series is
        //   employed instead of the Newton scheme:
        
        //   (1 + x)^(1/N)  =  1  +  x / N  +  x^2 * (1 - N) / (2// N^2)  +  ...
        
        //   See the comment about the parameter NIT in MPDIVX.
        
        var ia, iq, k, mpnw1, mq, na, n1, n2, n3 : Int
        var t1, t2, tn : Double
        let alt = 0.693147180559945309; let cl2 = 1.4426950408889633
        let nit = 3; let n30 = Int(pow(2, 30.0)); let mprxx = 1e-14
        var f1 = MPReal(repeating: 0, count: 8)
        var s0 = MPReal(repeating: 0, count: mpnw+6)
        var s1 = s0; var s2 = s0; var s3 = s0
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs(a[2])+4 || b[0] < Double(mpnw+6) {
            print ("*** MPNRTR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ia = sign (1, a[2])
        na = min (Int (abs (a[2])), mpnw)
        
        if na == 0 {
            b[2] = 0
            b[3] = 0
            return
        }
        if ia < 0 {
            print ("*** MPNRTR: Argument is negative.")
            mpabrt (59)
        }
        
        if n <= 0 || n > n30 {
            print ("*** MPNRTR: Improper value of N")
            mpabrt (60)
        }
        
        //   If N = 1 or 2, call MPEQ or MPSQRT instead.
        
        if n == 1 {
            mpeq (a, &b, mpnw)
            return
        } else if n == 2 {
            mpsqrt (a, &b, mpnw)
            return
        }
        
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        mpnw1 = mpnw + 1
        
        //   Set f1 = 1.
        f1 = [9, Double(mpnw1), 1, 0, 1, 0, 0]
        
        //   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.
        
        t1 = Double(mpnw)
        mq = Int(cl2 * log (t1) + 1 - mprxx)
        
        //   Check how close A is to 1.
        
        mpsub (a, f1, &s0, mpnw1)
        if s0[2] == 0 {
            mpeq (f1, &b, mpnw)
            return
        }
        n1 = 0
        mpmdc (s0, &t1, &n1, mpnw1)
        n2 = Int(cl2 * log (abs (t1)))
        t1 = t1 * pow(0.5, Double(n2))
        n1 = n1 + n2
        
        if n1 <= -30 {
            t2 = Double(n)
            n2 = Int(cl2 * log (t2) + 1 + mprxx)
            n3 = -mpnbt * mpnw1 / n1
            if Double(n3) < 1.25 * Double(n2) {
                
                //   A is so close to 1 that it is cheaper to use the binomial series.
                
                mpdivd (s0, t2, &s1, mpnw1)
                mpadd (f1, s1, &s2, mpnw1)
                k = 0
                
                repeat {
                    //                100 continue
                    k = k + 1
                    t1 = Double(1 - k * n)
                    t2 = Double((k + 1) * n)
                    mpmuld (s1, t1, &s3, mpnw1)
                    mpdivd (s3, t2, &s1, mpnw1)
                    mpmul (s0, s1, &s3, mpnw1)
                    mpeq (s3, &s1, mpnw1)
                    mpadd (s1, s2, &s3, mpnw1)
                    mpeq (s3, &s2, mpnw1)
                    if s1[2] != 0 && s1[3] >= Double(-mpnw1) {
                        // goto 100
                    } else {
                        mpeq (s2, &s1, mpnw1)
                        mproun (&s1, mpnw)
                        mpeq (s1, &b, mpnw)
                        return
                    }
                } while true
            }
        }
        
        //   Compute the initial approximation of A ^ (-1/N).
        
        tn = Double(n)
        mpmdc (a, &t1, &n1, mpnw1)
        n2 = -n1 / Int(tn)
        t2 = exp (-1.0 / tn * (log (t1) + (Double(n1) + tn * Double(n2)) * alt))
        mpdmc (t2, n2, &s2, mpnw1)
        mpnw1 = 5
        iq = 0
        
        //   Perform the Newton-Raphson iteration described above with a dynamically
        //   changing precision level MPNW (one greater than powers of two).
        
        for k in 1...mq {
            if k > 2 { mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1 }
            
            repeat {
 //               110  continue
                mpnpwr (s2, n, &s0, mpnw1)
                mpmul (a, s0, &s1, mpnw1)
                mpsub (f1, s1, &s0, mpnw1)
                mpmul (s2, s0, &s1, mpnw1)
                mpdivd (s1, tn, &s0, mpnw1)
                mpadd (s2, s0, &s1, mpnw1)
                mpeq (s1, &s2, mpnw1)
                if k == mq - nit && iq == 0 {
                    iq = 1
//                    goto 110
                } else {
                    break
                }
            } while true
        }
        
        //   Take the reciprocal to give final result.
        
        mpdiv (f1, s2, &s1, mpnw)
        
        //   Restore original precision level.
//        130 continue
        mproun (&s1, mpnw)
        mpeq (s1, &b, mpnw)
    } //  mpnrtr
    
    static func mpoutw (_ anam : String, _ a : MPReal, _ mpnw : Int) {
        
        //   This outputs the words of A up to the end of the active mantissa.
        //   This is for internal debugging only; it should not be called by user.
        var na : Int
        
        na = min (Int (abs (a[2])), mpnw)
        print(anam)
        for i in 0...na+5 {
            print(a[i])
        }
//        write (iu, "(4f18.0)") (a[i], i = 0, na + 5)
    } //  mpoutw
    
    static func mproun (_ a : inout MPReal, _ mpnw : Int) {
        
        //   This performs rounding and truncation of the MPR number A.  It is called
        //   by MPNORM, and also by other subroutines when the precision level is
        //   modified.  It is not intended to be directly called by the user.
        //   The parameter MPEXPMX is the absolute value of the largest exponent word
        //   allowed for MP numbers (see system parameters at start of this module).
        
        var ia, k, na, n4 : Int
        var a2 : Double
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || a[0] < Double(mpnw + 6) {
            print("*** MPROUN: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   Check for initial zeroes.
        
        a2 = a[3]
        a[3] = 0
        ia = sign (1, a[2])
        na = min (Int(abs (a[2])), mpnw)
        n4 = na + 4
        
        if a[4] == 0 {
            
            //   Find the first nonzero word and shift the entire number left.  The length
            //   of the result is reduced by the length of the shift.
            var flag = false
            var i = 4
            while i <= n4 {
                if a[i+1] != 0  {
                    flag = true
                    break // goto 110
                }
                i += 1
            }
            
            if !flag {
                a[2] = 0
                a[3] = 0
                return
            }
            
            // 110 continue
            
            k = i - 3
            
            for i in 3...n4-k {
                a[i+1] = a[i+k+1]
            }
            
            a2 = a2 - Double(k)
            na = na - max (k - 2, 0)
            if k == 2 { a[na+4] = 0 }
        }
        
        //   Perform rounding.
        
        if na == mpnw {
            if a[na+4] >= 0.5 * mpbdx { a[na+3] += 1 }
            
            //   Release carries as far as necessary due to rounding.
            var flag = false
            for i in stride(from: na+2, through: 3, by: -1) { // na+2, 3, -1 {
                if a[i+1] < mpbdx { flag = true; break /* goto 140 */ }
                a[i+1] = a[i+1] - mpbdx
                a[i] = a[i] + 1
            }
            
            //   Release of carries due to rounding continued all the way to the start --
            //   i.e. number was entirely 9"s.
            if !flag {
                a[4] = a[3]
                na = 1
                a2 = a2 + 1
            }
        }
        
        // 140 continue
        
        if a[na+3] == 0 {
            
            //   At least the last mantissa word is zero.  Find the last nonzero word
            //   and adjust the length of the result accordingly.
            var flag = false
            var i = na + 2
            while i >= 3 {
                if a[i+1] != 0  {
                    flag = true; break //  goto 160
                }
                i -= 1
            }
            
            if !flag {
                a[2] = 0
                a[3] = 0
                return
            }
            
            // 160  continue
            
            na = i - 2
        }
        
        //   Check for overflow and underflow.
        
        if a2 < -mpexpmx {
            print ("*** MPROUN: Exponent underflow.")
            mpabrt (68)
        } else if a2 > mpexpmx {
            print ("*** MPROUN: Exponent overflow.")
            mpabrt (69)
        }
        
        //   Check for zero.
        
        if a[4] == 0 {
            a[1] = Double(mpnw)
            a[2] = 0
            a[3] = 0
        } else {
            a[1] = Double(mpnw)
            a[2] = Double(sign (Double(na), Double(ia)))
            a[3] = a2
            a[na+4] = 0
            a[na+5] = 0
        }
        
        //            170  continue
        //
        //            return
    } //  mproun
    
    static func mpsqrt (_ a : MPReal, _ b : inout MPReal, _ mpnw : Int) {
        
        //   This computes the square root of the MPR number A and returns the result in B.
        
        //   This subroutine employs the following Newton-Raphson iteration, which
        //   converges to 1 / Sqrt(A):
        
        //    X_{k+1} = X_k + 0.5 * (1 - X_k^2 * A) * X_k
        
        //   where the multiplication () * X_k is performed with only half of the
        //   normal level of precision.  These iterations are performed with a
        //   working precision level MPNW that is dynamically changed, approximately
        //   doubling with each iteration (except that at iteration NIT before the final
        //   iteration, the iteration is repeated without doubling the precision, in order to
        //   enhance accuracy) .  The final iteration is performed as follows
        //   (this is due to A. Karp):
        
        //    Sqrt(A) = (A * X_n) + 0.5 * [A - (A * X_n)^2] * X_n  (approx.)
        
        //   where the multiplications A * X_n and [] * X_n are performed with only
        //   half of the final level of precision.
        
        var ia, iq, mpnw1, mq, n, na, nw1, nw2, n2 : Int
        var t1, t2 : Double
        let cl2 = 1.4426950408889633; let mprxx = 1e-14;  let nit = 3
        var s0 = MPReal(repeating: 0, count: mpnw+7)
        var s1 = s0; var s2 = s0; var s3 = s0
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < Double(mpnw + 6) {
            print ("*** MPSQRT: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ia = sign (1, a[2])
        na = min (Int (abs (a[2])), mpnw)
        
        if na == 0 {
            b[1] = Double(mpnw)
            b[2] = 0
            b[3] = 0
            return
        }
        if ia < 0 {
            print("*** MPSQRT: Argument is negative.")
            mpabrt (70)
            return
        }
        
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        
        //   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.
        
        t1 = Double(mpnw)
        mq = Int(cl2 * log (t1) + Double(1 - mprxx))
        
        //   Compute the initial approximation of 1 / Sqrt(A).
        n = 0
        mpmdc (a, &t1, &n, mpnw)
        n2 = -n / 2
        t2 = sqrt (t1 * pow(2, Double(n + 2 * n2)))
        t1 = 1 / t2
        mpdmc (t1, n2, &s2, mpnw)
        mpdmc (1, 0, &s3, mpnw)
        
        mpnw1 = 5
        iq = 0
        nw1 = mpnw1
        nw2 = mpnw1
        
        //   Perform the Newton-Raphson iteration described above with a dynamically
        //   changing precision level MPNW (one greater than powers of two).
        
        for k in 1...mq - 1 {
            if k > 2 {
                nw1 = mpnw1
                mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
                nw2 = mpnw1
            }
            
            // 100  continue
            repeat {
                mpmul (s2, s2, &s0, nw2)
                mpmul (a, s0, &s1, nw2)
                mpsub (s3, s1, &s0, nw2)
                mpmul (s2, s0, &s1, nw1)
                mpmuld (s1, 0.5, &s0, nw1)
                mpadd (s2, s0, &s1, nw2)
                mpeq (s1, &s2, nw2)
                
                if k == mq - nit && iq == 0 {
                    iq = 1
                    // goto 100
                } else {
                    break
                }
            } while true
        }
        
        //   Perform last iteration using Karp"s trick.
        
        nw1 = mpnw1
        mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
        nw2 = mpnw1
        
        mpmul (a, s2, &s0, nw1)
        mpmul (s0, s0, &s1, nw2)
        mpsub (a, s1, &s3, nw2)
        mpmul (s3, s2, &s1, nw1)
        mpmuld (s1, 0.5, &s3, nw1)
        mpadd (s0, s3, &s2, nw2)
        
        //   Restore original precision level.
        
        mproun (&s2, mpnw)
        mpeq (s2, &b, mpnw)
        
//        120 continue
//        return
    } //  mpsqrt
    
    static func mpsub  (_ a : MPReal, _ b : MPReal, _ c : inout MPReal, _ mpnw : Int)  {
        
        //   This routine subtracts MPR numbers A and B to yield C.
        
        var nb : Int
        var s = MPReal(repeating: 0, count: mpnw+6)
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < abs (b[2]) + 4 || c[0] < Double(mpnw + 6) {
            print ("*** MPSUB: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        nb = min (abs (Int (b[2])), mpnw)
        s[0] = Double(mpnw + 6)
        s[1] = Double(mpnw)
        if b[2] == 0 {
            s[2] = 0
        } else if b[2] > 0 {
            s[2] = Double(-nb)
        } else {
            s[2] = Double(nb)
        }
        
        for i in 3...nb + 5 {
            s[i] = b[i]
        }
        
        mpadd (a, s, &c, mpnw)
    } //  mpsub
    
    // ***  The following are the extra-high precision routines:
    
    static func mpdivx (_ a : MPReal, _ b : MPReal, _ c : inout MPReal, _ mpnw : Int)  {
        
        //   This divides A by B and returns the result in C.
        
        //   This subroutine employs the following Newton-Raphson iteration, which
        //   converges to 1 / B:
        
        //    X_{k+1} = X_k + (1 - X_k * B) * X_k
        
        //   where the multiplication () * X_k is performed with only half of the
        //   normal level of precision.  These iterations are performed with a
        //   working precision level MPNW that is dynamically changed, approximately
        //   doubling with each iteration (except that at iteration NIT before the
        //   final iteration, the iteration is repeated without doubling the
        //   precision, in order to enhance accuracy).  The final iteration is
        //   performed as follows (this is due to A. Karp):
        
        //    A / B = (A * X_n) + [A - (A * X_n) * B] * X_n  (approx.)
        
        //   where the multiplications A * X_n and [] * X_n are performed with only
        //   half of the final level of precision.
        
        var iq, mpnw1, mq, n, na, nb, nw1, nw2 : Int
        var t1, t2 : Double
        let cl2 = 1.4426950408889633; let mprxx = 1e-14; let nit = 3
        var s0 = MPReal(repeating: 0, count: mpnw+6)
        var s1 = s0; var s2 = s0; var s3 = s0
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < abs (b[2]) + 4 || c[0] < Double(mpnw + 6) {
            print("*** MPDIVX: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //ia = sign (1, a[2])
        //ib = sign (1, b[2])
        na = min (Int (abs (a[2])), mpnw)
        nb = min (Int (abs (b[2])), mpnw)
        
        if na == 0 {
            // Note: Original code assigned to b[] which would do nothing
            c[1] = Double(mpnw)
            c[2] = 0
            c[3] = 0
            return
        }
        if nb == 0 {
            print ("*** MPDIVX: Divisor is zero.")
            mpabrt (33)
        }
        
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        s3[0] = Double(mpnw + 7)
        
        //   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.
        
        t1 = Double(mpnw)
        mq = Int(cl2 * log (t1) + 1 - mprxx)
        
        //   Compute the initial approximation of 1 / B.
        n = 0
        mpmdc (b, &t1, &n, mpnw)
        t2 = 1 / t1
        mpdmc (t2, -n, &s2, mpnw)
        mpdmc (1, 0, &s3, mpnw)
        
        mpnw1 = mpnw + 1
        mpnw1 = 5
        iq = 0
        nw1 = mpnw1
        nw2 = mpnw1
        
        //   Perform the Newton-Raphson iteration described above with a dynamically
        //   changing precision level MPNW (one greater than powers of two).
        
        for k in 1...mq - 1 {
            if k > 2 {
                nw1 = mpnw1
                mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
                nw2 = mpnw1
            }
            
            // 100  continue
            repeat {
                mpmul (b, s2, &s1, nw2)
                mpsub (s3, s1, &s0, nw2)
                mpmul (s2, s0, &s1, nw1)
                mpadd (s2, s1, &s0, nw2)
                mpeq (s0, &s2, nw2)
                if k == mq - nit && iq == 0 {
                    iq = 1
                    // goto 100
                } else {
                    break
                }
            } while true
        }
        
        //   Perform last iteration using Karp"s trick.
        
        nw1 = mpnw1
        mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
        nw2 = mpnw1
        
        mpmul (a, s2, &s0, nw1)
        mpmul (s0, b, &s1, nw2)
        mpsub (a, s1, &s3, nw2)
        mpmul (s3, s2, &s1, nw1)
        mpadd (s0, s1, &s2, nw2)
        
        //   Restore original precision level.
        
        mproun (&s2, mpnw)
        mpeq (s2, &c, mpnw)
        
//        120 continue
//        return
    } //  mpdivx
    
    static func mpfftcr (_ iss : Int, _ m : Int, _ n : Int, _ nsq : Int, _ x : inout MPComplex, _ y : inout MPReal) {
        
        //   This performs an N-point complex-to-real FFT, where N = 2^M.  X is the
        //   double complex input array, and Y is the double precision output array.
        //   The array X is used as a scratch array in MPFFT1, and so is overwritten.
        //   X and Y must be dimensioned as shown below.  IS is the sign of the FFT.
        //   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
        //   This routine is not intended to be called directly by the user.
        
        var ku, mx, n2, n4 : Int
        var dc1 = MPComplex(repeating:Complex64(), count:n/2)
        var ai, a1, a2, x1, x2 : Complex64

        mx = Int(mpuu1[1].re)

        //   Check if input parameters are invalid.

        if (iss != 1 && iss != -1) || m < 3 || m > mx {
            print ("*** MPFFTCR: Either the UU arrays have not been initialized or else one of the input parameters is invalid")
            mpabrt (677)
        }

        //n1 = Int(pow(2, Double(m / 2)))
        n2 = n / 2
        n4 = n / 4
        ai = Complex64(0, 1)

        //   Construct the input to MPFFT1.

        dc1[1] = 0.5 * Complex64((x[1] + x[n2+1]).re, (x[1] - x[n2+1]).re)
        if (iss == 1) {
            dc1[n4+1] = x[n4+1].conj
        } else {
            dc1[n4+1] = x[n4+1]
        }
        ku = n2

        if (iss != 1) {
            for k in 2...n4 {
                x1 = x[k]
                x2 = x[n2+2-k].conj
                a1 = x1 + x2
                a2 = ai * mpuu1[k+ku] * (x1 - x2)
                dc1[k] = 0.5 * (a1 + a2)
                dc1[n2+2-k] = 0.5 * (a1 - a2).conj
            }
        } else {
            for k in 2...n4 {
                x1 = x[k]
                x2 = x[n2+2-k].conj
                a1 = x1 + x2
                a2 = ai * mpuu1[k+ku].conj * (x1 - x2)
                dc1[k] = 0.5 * (a1 + a2)
                dc1[n2+2-k] = 0.5 * (a1 - a2).conj
            }
        }

        //   Perform a normal N/2-point FFT on DC1.
        // Not sure how to deal with this?
        assert(false, "Fix mpfftcr!")
       // mpfft1 (iss, m - 1, n1, n2 / n1, &dc1, &x)

        //   Copy DC1 to Y such that DC1(k) = Y(2k-1) + i Y(2k).

//        for k in 1...n / 2 {
//            y[2*k-1] = dc1[k].re
//            y[2*k] = dc1[k].im
//        }
    } //  mpfftcr
    
    static func mpfftrc (_ iss : Int, _ m : Int, _ n : Int, _ nsq : Int, _ x: MPReal, _ y: inout MPComplex) {
        
        //   This performs an N-point real-to-complex FFT, where N = 2^M.  X is the
        //   dobule precision input array, and Y is the double complex output array.
        //   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
        //   This routine is not intended to be called directly by the user.
        
        var mx, n2 : Int
        var dc1 = MPComplex(repeating:Complex64(), count:n/2)
        // var ai, a1, a2, z1, z2 : Complex64
        
        mx = Int(mpuu1[1].re)
        
        //   Check if input parameters are invalid.
        
        if (iss != 1 && iss != -1) || m < 3 || m > mx {
            print ("*** MPFFTRC: either the UU arrays have not been initialized",
                   "or else one of the input parameters is invalid \(iss) \(m) \(mx)")
            mpabrt (677)
        }
        
        //n1 = Int(pow(2, Double(m / 2)))
        n2 = n / 2
        // n21 = n2 + 1
//        n4 = n / 4
//        let ai = Complex64(0, -1)
        
        //   Copy X to DC1 such that DC1(k) = X(2k-1) + i X(2k).
        
        for k in 1...n2 {
            dc1[k] = Complex64(x[2*k-1], x[2*k])
        }
        
        //   Perform a normal N/2-point FFT on DC1.
        
        // Not sure how to deal with this?
        assert(false, "Fix mpfftrc!")
        // mpfft1 (iss, m - 1, n1, n2 / n1, &dc1, &y)
        
        //   Reconstruct the FFT of X.
        
//        y[1] = Complex64(2 * (dc1[1].re + dc1[1].im), 0)
//        if (iss == 1) {
//            y[n4+1] = 2 * dc1[n4+1]
//        } else {
//            y[n4+1] = 2 * dc1[n4+1].conj
//        }
//        y[n2+1] = Complex64(2 * (dc1[1].re - dc1[1].im), 0)
//        ku = n2
//
//        if (iss == 1) {
//            for k in 2...n4 {
//                z1 = dc1[k]
//                z2 = dc1[n2+2-k].conj
//                a1 = z1 + z2
//                a2 = ai * mpuu1[k+ku] * (z1 - z2)
//                y[k] = a1 + a2
//                y[n2+2-k] = (a1 - a2).conj
//            }
//        } else {
//            for k in 2...n4 {
//                z1 = dc1[k]
//                z2 = dc1[n2+2-k].conj
//                a1 = z1 + z2
//                a2 = ai * mpuu1[k+ku].conj * (z1 - z2)
//                y[k] = a1 + a2
//                y[n2+2-k] = (a1 - a2).conj
//            }
//        }
        
    } //  mpfftrc
    

    static func mpfft1 (_ iss : Int, _ m : Int, _ n1 : Int, _ n2 : Int, _ x: inout [MPComplex], _ y: inout [MPComplex]) {
        
        //   This routine performs a complex-to-complex FFT.  IS is the sign of the
        //   transform, N = 2^M is the size of the transform.  N1 = 2^M1 and N2 = 2^M2,
        //   where M1 and M2 are defined as below.  X is the input and output array,
        //   and Y is a scratch array.  X must have at N, and Y at least N + N1*MPNSP1,
        //   double complex cells.  The arrays MPUU1 and MPUU2 must have been
        //   initialized by calling MPINIFFT.  This routine is not intended to be called
        //   directly by the user.
        
        //   This employs the two-pass variant of the "four-step" FFT.  See the
        //   article by David H. Bailey in J. of Supercomputing, March 1990, p. 23-35.
        var iu, j2, ku, m1, m2, nr1, nr2 : Int
        let zero = MPComplex(repeating: Complex64(), count:mpnrow+mpnsp1+1)
        var z1 = [MPComplex](repeating: zero, count:n1+1)
        var z2 = z1
        
      //  let n = Int(pow(2.0, Double(m)))
        m1 = (m + 1) / 2
        m2 = m - m1
        nr1 = min (n1, mpnrow)
        nr2 = min (n2, mpnrow)
        ku = Int(mpuu2[m].re)
        
        for i in stride(from: 0, through: n1-1, by: nr1) { // 0, n1 - 1, nr1
            
            //   Copy NR1 rows of X (treated as a N1 x N2 complex array) into Z1.
            
            for j in 1...n2 {
                for k in 1...nr1 {
                    z1[k][j] = x[i+k][j]
                }
            }
            
            //   Perform NR1 FFTs, each of length N2.
            mpfft2 (iss, nr1, m2, n2, &z1, &z2)
            
            //   Multiply the resulting NR1 x N2 complex block by roots of unity and
            //   store transposed into the appropriate section of Y.
            
            iu = i + ku - n1 - 1
            if iss == 1 {
                for j in 1...n2 {
                    for k in 1...nr1 {
                        y[j][i+k] = mpuu2[iu+k+j*n1] * z1[k][j]
                    }
                }
            } else {
                for j in 1...n2 {
                    for k in 1...nr1 {
                        y[j][i+k] = mpuu2[iu+k+j*n1].conj * z1[k][j]
                    }
                }
            }
        }
        
        for i in stride(from:0, through:n2 - 1, by:nr2) {
            
            //   Copy NR2 rows of the Y array into Z2.
            
            for j in 1...n1 {
                for k in 1...nr2 {
                    z2[k][j] = y[i+k][j]
                }
            }
            
            //   Perform NR2 FFTs, each of length N1.
            
            mpfft2 (iss, nr2, m1, n1, &z2, &z1)
            
            //   Copy NR2 x N1 complex block back into X array.  It"s a little more
            //   complicated if M is odd.
            
            if m.isMultiple(of: 2) {
                for j in 1...n1 {
                    for k in 1...nr2 {
                        x[i+k][j] = z2[k][j]
                    }
                }
            } else {
                for j in 1...n1 / 2 {
                    j2 = 2 * j - 1
                    
                    for k in 1...nr2 {
                        x[i+k][j] = z2[k][j2]
                        x[i+k+n2][j] = z2[k][j2+1]
                    }
                }
            }
        }
    } //  mpfft1
    
    static func mpfft2 (_ iss : Int, _ ns : Int, _ m : Int, _ n : Int, _ x: inout [MPComplex], _ y: inout [MPComplex]) {
        
        //   This performs NS simultaneous N-point complex-to-complex FFTs, where
        //   N = 2^M.  X is the input and output array, and Y is a scratch array.
        //   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
        //   This routine is not intended to be called directly by the user.
        
        //   Perform the second variant of the Stockham FFT.
        var flag = false
        for l in stride(from: 1, through:m, by:2) {
            mpfft3 (iss, l, ns, m, n, &x, &y)
            if (l == m) { flag = true; break /* goto 100 */ }
            mpfft3 (iss, l + 1, ns, m, n, &y, &x)
        }
        
        if !flag { return }
            
            //   Copy Y to X.
            
        // 100 continue
        
        for j in 1...n {
            for i in 1...ns {
                x[i][j] = y[i][j]
            }
        }
        
        // 110 continue
        
    } //  mpfft2
    
    static func mpfft3 (_ iss : Int, _ l : Int, _ ns : Int, _ m : Int, _ n : Int, _ x: inout [MPComplex], _ y: inout [MPComplex]) {
        
        //   This performs the L-th iteration of the second variant of the Stockham FFT
        //   on the NS vectors in X.  X is input/output, and Y is a scratch array.
        //   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
        //   This routine is not intended to be called directly by the user.
        
        var i11, i12, i21, i22, li, lj, lk, ku, n1 : Int
        var u1, x1, x2 : Complex64
        
        //   Set initial parameters.
        
        n1 = n / 2
        lk = Int(pow(2, Double(l - 1)))
        li = Int(pow(2, Double(m - l)))
        lj = 2 * lk
        ku = li + 1
        
        for i in 0...li - 1 {
            i11 = i * lk + 1
            i12 = i11 + n1
            i21 = i * lj + 1
            i22 = i21 + lk
            if (iss == 1) {
                u1 = mpuu1[i+ku]
            } else {
                u1 = mpuu1[i+ku].conj
            }
            
            for k in 0...lk - 1 {
                for j in 1...ns {
                    x1 = x[j][i11+k]
                    x2 = x[j][i12+k]
                    y[j][i21+k] = x1 + x2
                    y[j][i22+k] = u1 * (x1 - x2)
                }
            }
        }
    } //  mpfft3
    
    static func mpinifft (_ mpnw : Int) {
        
        //   This computes the root of unity arrays UU1 and UU2, which are required by
        //   the FFT routines, and places this data in the proper arrays defined in
        //   module MPFUNA.  MPNW is the largest precision level (in words) that will be
        //   subsequently used for this run.
        
        var iu, ku, ln, m, mm, mm1, mm2, mq, nn, nn1, nn2, nq, nwds : Int
        var d1 : Double
        let cl2 = 1.4426950408889633; let mprxx = 1e-14
        var pi, t1, ti, tpn : Double
        
        //  Determine sizes for FFT arrays.  Three words are added to mpnw, since many
        //  routines in MPFUND in particular increase the working precision upon entry.
        
        nwds = mpnw + 3
        d1 = 1.5 * Double((nwds + 1))
        m = Int(cl2 * log (d1) + 1 - mprxx)
        mq = m + 2
        nq = Int(pow(2, Double(mq)))
        
        if (mq + nq > mplfftx) {
            print ("*** MPINIFFT: Insufficient space for arrays mpuu1 and mpuu2.",
                "At least \(mq + nq) double complex cells must be allocated for each of",
                "these arrays in module mpfuna. See documentation for details.")
            mpabrt (91)
        }
        
        mpuu1[1] = Complex64(Double(mq),0)
        ku = 2
        ln = 1
        pi = acos (-1)
        
        for _ in 1...mq {
            t1 = pi / Double(ln)
            
            for i in 0...ln - 1 {
                ti = Double(i) * t1
                mpuu1[i+ku] = Complex64(cos (ti), sin (ti))
            }
            
            ku = ku + ln
            ln = 2 * ln
        }
        
        // write (6, 2) ku - 1
        // 2 format ("MPINIFFT: Size of table mpuu1 =",i10)
        
        ku = mq + 1
        mpuu2[1] = Complex64(Double(mq),0)
        
        for k in 2...mq {
            mpuu2[k] = Complex64(0, 0)
        }
        
        for k in 2...mq - 1 {
            mpuu2[k] = Complex64(Double(ku),0)
            mm = k
            nn = Int(pow(2, Double(mm)))
            mm1 = (mm + 1) / 2
            mm2 = mm - mm1
            nn1 = Int(pow(2, Double(mm1)))
            nn2 = Int(pow(2, Double(mm2)))
            tpn = 2 * pi / Double(nn)
            
            for j in 0...nn2 - 1 {
                for i in 0...nn1 - 1 {
                    iu = ku + i + j * nn1
                    t1 = tpn * Double(i * j)
                    mpuu2[iu] = Complex64(cos (t1), sin (t1))
                }
            }
            
            ku = ku + nn
        }
        
        // write (6, 3) ku - 1
        // 3 format ("MPINIFFT: Size of table mpuu2 =",i10)
        
    } //  mpinifft
   
 
    static func mplconv (_ iq: Int, _ n: Int, _ nsq: Int, _ a: MPReal, _ b: MPReal, _ c: inout MPReal) {
        
        //   This computes the linear convolution of A and B, returning the result
        //   in C.  If IQ is 1, { it is presumed B = A; if IQ = 2, { A //= B.
        //   NSQ is a spacing parameter, which should be set to more than sqrt (3*n).
        
        var m1, m2, n1, n2, n4, nm : Int
        var c0, an, t1, t2 : Double
        let cl2 = 1.4426950408889633; let mprxx = 1e-14; let mpffterrmx = 0.375
        var d1 = MPReal(repeating:0, count:6*n+3); var d2 = d1; var d3 = d1
        var dc1 = MPComplex(repeating: Complex64(), count:3*n+nsq*mpnsp1+4); var dc2 = dc1
        
        t1 = Double(n)
        m1 = Int(cl2 * log (t1) + Double(1 - mprxx))
        n1 = Int(pow(2, Double(m1)))
        m2 = m1 + 1
        n2 = 2 * n1
        n4 = 2 * n2
        nm = min (2 * n, n2)
        
        if (abs (iq) == 1) {
            
            //   Compute the square of a -- only one forward FFT is needed.
            
            for i in 1...n {
                d1[i] = a[i]
            }
            
            for i in n + 1...n2 {
                d1[i] = 0
            }
            
            //   Perform a forward real-to-complex FFT on the vector in a.
            
            mpfftrc (1, m2, n2, nsq, d1, &dc1)
            
            //   Square the resulting complex vector.
            
            for i in 1...n1 + 1 {
                dc1[i] = dc1[i] * dc1[i]
            }
        } else {
            
            //   Compute the product of a and b -- two forward FFTs are needed.
            for i in 1...n {
                d1[i] = a[i]
                d2[i] = b[i]
            }
            
            for i in n + 1...n2 {
                d1[i] = 0
                d2[i] = 0
            }
            
            //   Perform forward real-to-complex FFTs on the vectors in a and b.
            
            mpfftrc (1, m2, n2, nsq, d1, &dc1)
            mpfftrc (1, m2, n2, nsq, d2, &dc2)
            
            //   Multiply the resulting complex vectors.
            
            for i in 1...n1 + 1 {
                dc1[i] = dc1[i] * dc2[i]
            }
        }
        
        //   Perform an inverse complex-to-real FFT on the resulting data.
        
        mpfftcr (-1, m2, n2, nsq, &dc1, &d3)
        
        //   Divide by N4.
        
        an = Double(1 / n4)
        c0 = 0
        
        for i in 1...nm {
            t1 = an * d3[i]
            //  t2 = anint (t1)
            if (d3[i] >= 0) {
                t2 = aint (t1 + 0.5)
            } else {
                t2 = aint (t1 - 0.5)
            }
            c[i] = t2
            c0 = max (c0, abs (Double ((t2 - t1))))
        }
        
        mpffterr = max (c0, mpffterr)
        if (c0 > mpffterrmx) {
            print ("*** MPLCONV: excessive rounding error = \(c0)")
            mpabrt (55)
        }
    } //  mplconv
    
    
    static func mpmulx (_ a: MPReal, _ b: MPReal, _ c: inout MPReal, _ mpnw: Int) {
        
        //   This routine multiplies MP numbers A and B to yield the MP product C,
        //   using a FFT-convolution technique.  Before calling MPMULX, the arrays
        //   UU1 and UU2 must be initialized by calling MPINIFFT.  For modest levels
        //   of precision, use MPMUL.
        
        var ia, ib, na, nb, nc, nn, nx : Int
        let mprxx = 1e-14
        var d = MPReal(repeating: 0, count:mpnw+8)
        var t0, t1, t2, t3, t4, t5 : Double
        let r16 = pow(0.5, Double(16)); let t16 = Double(1<<16); let r32 = pow(0.5, Double(32))
        let t32 = Double(1<<32); let r48 = pow(0.5, Double(48)); let t48 = Double(1<<48)
        var d1 = MPReal(repeating: 0, count:3*mpnw+17); var d2 = d1
        var d3 = MPReal(repeating: 0, count:6*mpnw+33)
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < abs (b[2]) + 4 || Int(c[0]) < mpnw + 6 {
            print ("*** MPMULX: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ia = sign (1, a[2])
        ib = sign (1, b[2])
        na = min (Int (abs (a[2])), mpnw)
        nb = min (Int (abs (b[2])), mpnw)
        nc = min (na + nb, mpnw)
        nn = 3 * max (na, nb)
        nx = Int(sqrt (Double(3 * nn)) + mprxx)
        
        //   Divide each word into three 16-bit chunks.
        
        for i in 0...na - 1 {
            t1 = a[i+4]
            d1[3*i] = aint (r32 * t1)
            t1 = t1 - t32 * d1[3*i]
            d1[3*i+1] = aint (r16 * t1)
            d1[3*i+2] = t1 - t16 * d1[3*i+1]
        }
        
        for i in 3 * na...nn - 1 {
            d1[i] = 0
        }
        
        //   If A is the same array as B, { there is no need to deal with B.
        
        //        if (loc (a) == loc (b)) {
        //            mplconv (1, nn, nx, d1, d2, &d3)
        //        } else {
        for i in 0...nb - 1 {
            t1 = b[i+4]
            d2[3*i] = aint (r32 * t1)
            t1 = t1 - t32 * d2[3*i]
            d2[3*i+1] = aint (r16 * t1)
            d2[3*i+2] = t1 - t16 * d2[3*i+1]
        }
        
        for i in 3 * nb...nn - 1 {
            d2[i] = 0
        }
        
        mplconv (2, nn, nx, d1, d2, &d3)
        //       }
        
        //   Release carries.
        
        for i in 0...3 * nc + 13 {
            d1[i] = 0
        }
        
        for i in 0...min (3 * nc + 9, 2 * nn - 1) {
            t0 = d3[i]
            t1 = aint (r48 * t0)
            t2 = t0 - t48 * t1
            t3 = aint (r32 * t2)
            t4 = t2 - t32 * t3
            t5 = aint (r16 * t4)
            d1[i+3] = t4 - t16 * t5
            d1[i+2] = d1[i+2] + t5
            d1[i+1] = d1[i+1] + t3
            d1[i] = d1[i] + t1
        }
        
        //  Recombine words, with proper offset.
        
        d[0] = 0
        d[1] = 0
        d[2] = 0
        d[3] = 0
        
        for i in 0...nc + 3 {
            d[i+4] = t32 * d1[3*i+2] + t16 * d1[3*i+3] + d1[3*i+4]
        }
        
        d[0] = Double(mpnw + 6)
        d[1] = Double(mpnw)
        d[2] = Double(sign (Double(nc), Double(ia * ib)))
        d[3] = a[3] + b[3] + 1
        
        //   Fix up the result.
        
        d1[0] = Double(mpnw + 6)
        mpnorm (d, &c, mpnw)
        
        // 190 continue
    } //  mpmulx
    
    //>   These two subroutines are for real*16 support.  See "Uncomment" below
    //>   for differences.
    
    static func mpmqc (_ a: MPReal, _ b: inout Double, _ n: inout Int, _ mpnw: Int) {
        
        //   This returns a quad (real*16) approximation the MPR number A in the form B * 2^n.
        
        var na : Int
        
        //>  Uncomment this line if real*16 is supported.
        // parameter (knd = kind (0.q0))
        //>  Otherwise uncomment this line.
        var aa, b : Double
        
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 {
            print ("*** MPMQC: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        if a[2] == 0  {
            b = 0
            n = 0
            return // goto 100
        }
        
        na = Int(abs (a[2]))
        aa = a[4]
        if (na >= 2) { aa = aa + mprdx * a[5] }
        if (na >= 3) { aa = aa + mprx2 * a[6] }
        if (na >= 4) { aa = aa + mprdx * mprx2 * a[7]}
        
        n = mpnbt * Int(a[3])
        b = Double(sign (aa, a[2]))
        
        //   Reduce b to within 1 and 2.
        
        na = Int(log(abs (b)) / log(2) + mprdx)
        b = b / Double(1<<na)
        n = n + na
        if (abs (b) < 1) {
            b = 2 * b
            n = n - 1
        } else if (abs (b) > 2) {
            b = 0.5 * b
            n = n + 1
        }
        
        // 100  continue
    } //  mpmqc
    
    static func mpqmc (_ a: Double, _ n: Int, _ b: inout MPReal, _ mpnw: Int) {
        
        //   This routine converts the quad (real*16) number A * 2^N to MPR form in B.
        
        //   NOTE however that if A = 0.1q0, for example, { B will NOT be the true
        //   multiprecision equivalent of 1/10, since 0.1q0 is not an exact binary value.
        
        //   Examples of exact binary values (good): 123456789.0, 0.25d0, -5.3125d0.
        //   Examples of inexact binary values (bad): 0.1d0, 1234567.8d0, -3333.3d0.
        
        var n1, n2 : Int
        
        //>  Uncomment this line if real*16 is supported.
        // parameter (knd = kind (0.q0))
        //>  Otherwise uncomment this line.
        
        var aa : Double
        
        // End of declaration
        
        if mpnw < 4 || Int(b[0]) < mpnw + 6 {
            print ("*** MPQMC: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        //   Check for zero.
        
        if a == 0 {
            b[1] = Double(mpnw)
            b[2] = 0
            b[3] = 0
            return // goto 150
        }
        n1 = n / mpnbt
        n2 = n - mpnbt * n1
        aa = abs (a) * pow(2.0, Double(n2))
        
        //   Reduce AA to within 1 and MPBDX.
        
        if (aa >= mpbdx) {
            
            for k in 1...350 {
                aa = mprdx * aa
                if (aa < mpbdx) {
                    n1 = n1 + k
                    break // goto 120
                }
            }
            
        } else if (aa < 1) {
            
            for k in 1...350 {
                aa = mpbdx * aa
                if (aa >= 1) {
                    n1 = n1 - k
                    break // goto 120
                }
            }
            
        }
        
        //   Store successive sections of AA into B.
        
        // 120  continue
        
        b[3] = Double(n1)
        b[4] = aint (aa)
        aa = mpbdx * (aa - b[3+1])
        b[5] = aint (aa)
        aa = mpbdx * (aa - b[4+1])
        b[6] = aint (aa)
        aa = mpbdx * (aa - b[5+1])
        b[7] = aint (aa)
        b[8] = 0
        b[9] = 0
        
        var i = 7
        while i >= 3 {
            if b[i+1] != 0  { break /* goto 140 */ }
            i -= 1
        }
        
        // 140  continue
        
        b[1] = Double(mpnw)
        aa = Double(i - 2)
        b[2] = Double(sign (aa, a))
        
        //150 continue
        // return
    } //  mpqmc
    
    
}
