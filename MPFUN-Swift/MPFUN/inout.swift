//*****************************************************************************

//  MPFUN-Fort: A thread-safe arbitrary precision computation package
//  Binary-decimal, decimal-binary and I/O functions (module MPFUNC)

//  Revision date:  9 Nov 2018

//  AUTHOR:
//     David H. Bailey
//     Lawrence Berkeley National Lab (retired) and University of California, Davis
//     Email: dhbailey@lbl.gov

//  COPYRIGHT AND DISCLAIMER:
//    All software in this package (c) 2017 David H. Bailey.
//    By downloading or using this software you agree to the copyright, disclaimer
//    and license agreement in the accompanying file DISCLAIMER.txt.

//  PURPOSE OF PACKAGE:
//    This package permits one to perform floating-point computations (real and
//    complex) to arbitrarily high numeric precision, by making only relatively
//    minor changes to existing Fortran-90 programs.  All basic arithmetic
//    operations and transcendental functions are supported, together with several
//    special functions.

//    In addition to fast execution times, one key feature of this package is a
//    100% THREAD-SAFE design, which means that user-level applications can be
//    easily converted for parallel execution, say using a threaded parallel
//    environment such as OpenMP.  There are NO global shared variables (except
//    static compile-time data), and NO initialization is necessary unless
//    extremely high precision (> 19,500 digits) is required.

//  DOCUMENTATION:
//    A detailed description of this package, and instructions for compiling
//    and testing this program on various specific systems are included in the
//    README file accompanying this package, and, in more detail, in the
//    following technical paper:

//    David H. Bailey, "MPFUN2015: A thread-safe arbitrary precision package,"
//    available at http://www.davidhbailey.com/dhbpapers/mpfun2015.pdf.

//  DESCRIPTION OF THIS MODULE (MPFUNC):
//    This module contains subroutines for binary-decimal and decimal-binary
//    conversion, together with low-level E-format and F-format conversion, and
//    basic input and output.

//
//  inout.swift
//  MPFUN-Swift
//
//  Created by Mike Griebling on 1 May 2019.
//  Copyright © 2019 Computer Inspirations. All rights reserved.
//

import Foundation

extension MPFUN {
    
    static func mpctomp (_ a: String, _ b: inout MPRNumber, _ mpnw: Int) {
        
        //  Converts the string A into the MPR number B.
        //  Restrictions: (a) no embedded blanks; (b) a leading digit (possibly
        //  zero) must be present; and (c) a period must be present.  An exponent
        //  (with "d" or "e") may optionally follow the numeric value.
        
        var i1, i2, iexp, isgn, iexpsgn, ix, j, kde, kend : Int
        var kper, lexp : Int
        var lnum, lnum1, lnum2, mpnw1, n1, n2 : Int
        var t1, t2 : Double
        var ai, ksgn, kexpsgn : Character
        var ca : String
        var kexpst = ""
        var knumst = ""
        var knumst2 = ""
        let lexpmx = 9

        let d10w = pow(10.0, Double(mpndpw))
        var f = MPRNumber(repeating: 0, count: 9)
        var s0 = MPRNumber(repeating: 0, count: mpnw+7)
        var s1 = s0; var s2 = s0
        var a = a.trimmingCharacters(in: CharacterSet.whitespaces)   // mutable version of input string
        
        func abort(_ error: Int) {
            print ("*** MPCTOMP: Syntax error in input string; code = \(error)",
                "Restrictions: (a) no embedded blanks; (b) a leading digit (possibly",
                "zero) must be present; and (c) a period must be present.  An exponent",
                "(with 'd' or 'e') may optionally follow the numeric value.")
        }
        
        // End of declaration
        
        if mpnw < 4 || b[0] < Double(mpnw + 6) {
            print("*** MPCTOMP: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        s2[0] = Double(mpnw + 7)
        f[0] = 9.0
        f[1] = Double(mpnw)
        
//        for i in 2...8 {
//            f[i] = 0
//        }
        
        mpnw1 = mpnw + 1
        kde = 0
        kend = a.count
        kper = 0
        ksgn = " "
        kexpsgn = "+"
        
        //   Locate:
        //     kstart = index of first nonblank character.
        //     kend = index of last nonblank character.
        
//        for i in 1...n {
//            if (a[i] != " ")  { goto 100 }
//        }
//
        //   Input is completely blank.
        if a.isEmpty {
            abort(4)
            mpabrt (41)
        }
        
//        100 continue
//
//        kstart = i
//
//        for i in n, kstart, -1 {
//            if (a[i] /= " ") goto 110
//        }
//
//        i = kstart
//
//        110 continue
//
//        kend = i
        
        //   Scan input for:
        //     kde = index of "d" or "e".
        //     kexpend = index of end of exponent.
        //     kexpst = index of start of exponent.
        //     kespsgn = index of sign of exponent.
        //     knumend1 = index of end of numeric part prior to period.
        //     knumend2 = index of end of numeric part after period.
        //     knumst1 = index of start of numeric part prior to period.
        //     knumst2 = index of start of numeric part after period.
        //     kper = index of period.
        //     ksgn = index of sign of number.
        
        for (i, ch) in a.enumerated() {
            if ch.isWhitespace {
                abort(2)
                mpabrt (41)
            } else if ch == "+" || ch == "-" {
                if i == 0 {
                    ksgn = ch
                } else if kde > 0 && kexpst.isEmpty && i < kend {
                    kexpsgn = ch
                } else {
                    abort(3)
                    mpabrt (41)
                }
            } else if ch == "e" || ch == "E" || ch == "d" || ch == "D" {
                if kde == 0 && kper > 0 && i < kend {
                    kde = i
                } else {
                    abort(4)
                    mpabrt (41)
                }
            } else if ch == "." {
                if kper == 0 && kde == 0 && !knumst.isEmpty { // knumst1 > 0 && knumst2 == 0 {
                    kper = i
                } else {
                    abort(5)
                    mpabrt (41)
                }
            } else if ch.isNumber {
                if kexpst.isEmpty {
                    if kper == 0 {
                        knumst.append(ch)
                    } else {
                        knumst2.append(ch)
                    }
                } else {
                    kexpst.append(ch)
                }
                if i == kend {
                    if kde != 0 && kexpst.isEmpty {
                        abort(6)
                        mpabrt (41)
                    }
                }
            } else {
                abort(7)
                mpabrt (41)
            }
        }
        
        // write (6, *) "kde, kend, kexpend, kexpst =", kde, kend, kexpend, kexpst
        // write (6, *) "kexpsgn, numend1, knumend2, knumst1 =", kexpsgn, knumend1, knumend2, knumst1
        // write (6, *) "knumst2, kper, ksgn, kstart =", knumst2, kper, ksgn, kstart
        
        //   Decode exponent.
        
        if !kexpst.isEmpty {
            lexp = kexpst.count
            if (lexp > lexpmx) {
                print ("*** MPCTOMP: exponent string is too long.")
                mpabrt (85)
            }
            
//            for i in 1...lexp {
//                ca[i:i] = a[i+kexpst-1]
//            }
            
            iexp = Int(mpdigin (kexpst, lexp))
            if kexpsgn == "-" { iexp = -iexp }
        } else {
            iexp = 0
        }
        
        //   Determine sign of number.
        
        if ksgn == " " {
            isgn = 1
        } else if ksgn == "+" {
            isgn = 1
        } else if ksgn == "-" {
            isgn = -1
        }
        
        //   Determine lengths of two sections of number.
        
        lnum1 = knumst.count
        lnum2 = knumst2.count
        lnum = lnum1 + lnum2
        knumst += knumst2       // combine numerical strings
        
        // write (6, *) "iexp, lnum1, lnum2 =", iexp, lnum1, lnum2
        
        //   Determine the number of chunks of digits and the left-over.
        
        n1 = lnum / mpndpw
        n2 = lnum % mpndpw
        
        //   Construct first (left-over) portion, right-justified in CA.
    
//        ix = knumst1 - 1
        ca = String(knumst.prefix(n2))
        knumst = String(knumst.dropFirst(n2))
//        for i in 1...n2 {
//            ix = ix + 1
//            ca[i+mpndpw-n2:i+mpndpw-n2] = a[ix]
//        }
        
        t1 = mpdigin (ca, mpndpw)
        if t1 > 0 {
            f[2] = 1.0
            f[3] = 0.0
            f[4] = t1
        } else {
            f[2] = 0.0
            f[3] = 0.0
            f[4] = 0.0
        }
        mpeq (f, &s0, mpnw1)
        
        //   Process remaining chunks of digits.
        
        for _ in 0..<n1 {
//            for i in 1...mpndpw {
//                ix = ix + 1
//                if (ix == kper) { ix = ix + 1 }
//                ca[i:i] = a[ix]
//            }
            ca = String(knumst.prefix(mpndpw))
            knumst = String(knumst.dropFirst(mpndpw))
            
            t1 = mpdigin (ca, mpndpw)
            if (t1 > 0) {
                f[2] = 1.0
                f[3] = 0.0
                f[4] = t1
            } else {
                f[2] = 0.0
                f[3] = 0.0
                f[4] = 0.0
            }
            
            mpmuld (s0, d10w, &s1, mpnw1)
            mpadd (s1, f, &s0, mpnw1)
        }
        
        //  Correct exponent.
        
        iexp = iexp - lnum2
        f[2] = 1.0
        f[3] = 0.0
        f[4] = 10.0
        mpnpwr (f, iexp, &s1, mpnw1)
        mpmul (s0, s1, &s2, mpnw1)
        if ksgn == "-" { s2[2] = -s2[2] }
        
        //   Restore original precision and exit.
        
        mproun (&s2, mpnw)
        mpeq (s2, &b, mpnw)
        
        // write (6, *) "mpctomp: output ="
        // mpout (6, 420, 400, b, mpnw)
        
    } // mpctomp
    
    static func mpdigin (_ ca : String, _ n : Int) -> Double {
        
        //   This converts the string CA of nonblank length N to double precision.
        //   CA may only be modest length and may only contain digits.  Blanks are ignored.
        //   This is intended for internal use only.
        var k : Int
        
        var d1 = 0.0
        for ch in ca {
            if !ch.isWhitespace {
                k = ch.wholeNumberValue ?? -1
                if k < 0 {
                    print ("*** MPDIGIN: non-digit in character string = \(ch)")
                    mpabrt (86)
                } else if k <= 9 {
                    d1 = 10.0 * d1 + Double(k)
                }
            }
        }
        
        return d1
    } //  mpdigin

    static func mpdigout (_ a : Double, _ n : Int) -> String {
        
        //   This converts the double precision input A to a character(32) string of
        //   nonblank length N.  A must be a whole number, and N must be sufficient
        //   to hold it.  This is intended for internal use only.
        
        var d1, d2 : Double
        var ca : String
        let digits = "0123456789"
        var k : Int
        
        // End of declaration
        
        ca = " "
        d1 = abs (a)
        
        for _ in 1...n {
            d2 = aint (d1 / 10.0)
            k = Int(d1 - 10.0 * d2)
            d1 = d2
            ca = String(digits[digits.index(digits.startIndex, offsetBy: k)]) + ca
        }
        
        return ca
    } // mpdigout

    static func mpeformat (_ a : MPRNumber, _ nb : Int, _ nd : Int, _ b : inout String, _ mpnw : Int) {
        
        //   Converts the MPR number A into character form in the character(1) array B.
        //   NB (input) is the length of the output string, and ND (input) is the
        //   number of digits after the decimal point.  The format is analogous to
        //   Fortran E format.  The result is left-justified among the NB cells of B.
        //   The condition NB >= ND + 10 must hold or an error message will result.
        //   NB cells must be available in array B.
        
        var ia, ixp, i1, i2, mpnw1, na, nexp, nl : Int
        var ca, b2 : String
        let digits = "0123456789"
        var aa, an, t1, d10w : Double
        var f = MPRNumber(repeating: 0, count: 9)
        var s0 = MPRNumber(repeating: 0, count: mpnw+7)
        var s1 = s0
        
        // End of declaration
        
        if mpnw < 4 || a[0] < abs (a[2]) + 4 || nb < (nd + 10) {
            print ("*** MPEFORMAT: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        d10w = pow(10, Double(mpndpw))
        ia = sign (1.0, a[2])
        na = min (Int (abs (a[2])), mpnw)
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        mpnw1 = mpnw + 1
        
        //   Set f = 10.
        
        f[0] = 9.0
        f[1] = Double(mpnw1)
        f[2] = 1.0
        f[3] = 0.0
        f[4] = 10.0
        f[5] = 0.0
        f[6] = 0.0
        
        //   Determine power of ten for exponent, and scale input to within 1 and 10.
        
        if (na > 0) {
            aa = a[4]
            if (na >= 2) { aa = aa + mprdx * a[5] }
            if (na >= 3) { aa = aa + mprx2 * a[6] }
            t1 = log10 (2.0) * Double(mpnbt) * a[3] + log10 (aa)
            
            if (t1 >= 0.0) {
                nexp = Int(t1)
            } else {
                nexp = Int(t1 - 1.0)
            }
            
            if (nexp == 0) {
                mpeq (a, &s1, mpnw1)
            } else if (nexp > 0) {
                mpnpwr (f, nexp, &s0, mpnw1)
                mpdiv (a, s0, &s1, mpnw1)
            } else if (nexp < 0) {
                mpnpwr (f, -nexp, &s0, mpnw1)
                mpmul (a, s0, &s1, mpnw1)
            }
            
            //   If we didn"t quite get it exactly right, multiply or divide by 10 to fix.
            
            // 100 continue
            repeat {
                if (s1[3] < 0.0) {
                    nexp = nexp - 1
                    mpmuld (s1, 10.0, &s0, mpnw1)
                    mpeq (s0, &s1, mpnw1)
                    // goto 100
                } else if (s1[4] >= 10.0) {
                    nexp = nexp + 1
                    mpdivd (s1, 10.0, &s0, mpnw1)
                    mpeq (s0, &s1, mpnw1)
                    // goto 100
                } else {
                    break
                }
            } while true
            
            s1[2] = abs (s1[2])
        } else {
            nexp = 0
            mpeq (a, &s1, mpnw1)
        }
        
        //   Insert sign and first digit.
        
        b2 = ""
        if ia == -1 {
            b2.append("-")
        }
        if na > 0 {
            an = s1[4]
        } else {
            an = 0.0
        }
        ca = mpdigout (an, 1)
        b2.append(ca.first!)
        b2.append(".")
        ixp = b2.count-1
        
        //   Set f = an.
        
        f[0] = 9.0
        f[1] = Double(mpnw1)
        f[2] = 1.0
        f[3] = 0.0
        f[4] = an
        f[5] = 0.0
        f[6] = 0.0
        mpsub (s1, f, &s0, mpnw1)
        mpmuld (s0, d10w, &s1, mpnw1)
        
        //   Calculate the number of remaining chunks.
        
        nl = nd / mpndpw + 1
        
        //   Insert the digits of the remaining words.
        
        for _ in 1...nl {
            if (s1[2] != 0.0 && s1[3] == 0.0) {
                an = s1[4]
                f[2] = 1.0
                f[3] = 0.0
                f[4] = an
            } else {
                f[2] = 0.0
                f[3] = 0.0
                f[4] = 0.0
                an = 0.0
            }
            
            ca = mpdigout (an, mpndpw)
            b2 += ca
//            for i in 1...mpndpw {
//                ix = ix + 1
//                if ix > nb + 50 {
//                    print ("MPEFORMAT: Insufficient space in B2 array.")
//                    mpabrt (84)
//                }
//                b2[ix] = ca[i]
//            }
            
            mpsub (s1, f, &s0, mpnw1)
            mpmuld (s0, d10w, &s1, mpnw1)
        }
        
        //   Round the result.
        
        if b2.count >= nd + 1 {
            i1 = b2[b2.index(b2.startIndex, offsetBy: nd+1)].wholeNumberValue ?? 0
            if i1 >= 5 {
                
                //   Perform rounding, beginning at the last digit (position ND).  If the rounded
                //   digit is 9, set to 0, then repeat at position one digit to left.  Continue
                //   rounding if necessary until the decimal point is reached.
                var flag = false
                for i in stride(from:nd-1, through:0, by:-1) {
                    // manipulating characters in a Swift string is ... awkward
                    let ix = b2.index(b2.startIndex, offsetBy: i)
                    let r = Range(uncheckedBounds: (lower: ix, upper: ix))
                    i2 = b2[ix].wholeNumberValue ?? 0  //index (digits, b2[i]) - 1
                    if i2 <= 8 {
                        b2.replaceSubrange(r, with: String(digits[digits.index(digits.startIndex, offsetBy: i2+1)]))
//                        b2.remove(at: ix)
//                        b2.insert(digits[digits.index(digits.startIndex, offsetBy: i2+1)], at: ix)
//                        b2[i] = digits[i2+2]
                        flag = true; break // goto 180
                    } else {
                        b2.replaceSubrange(r, with: "0")
                    }
                }
                
                //   We have rounded up all digits to the right of the decimal point.  If the
                //   digit to the left of the decimal point is a 9, { set that digit to 1
                //   and increase the exponent by one; otherwise increase that digit by one.
                if !flag {
                    let ix = b2.index(b2.startIndex, offsetBy: ixp-1)
                    let r = Range(uncheckedBounds: (lower: ix, upper: ix))
                    if b2[ix] == "9" {
                        b2.replaceSubrange(r, with: "1")
                        nexp = nexp + 1
                    } else {
                        i1 = b2[ix].wholeNumberValue ?? 0
                        let ix = digits.index(digits.startIndex, offsetBy: i1+1)
                        b2.replaceSubrange(r, with: String(digits[ix]))
                    }
                }
            }
        }
        
        // 180 continue
        
        //   Done with mantissa.  Insert exponent.
        
        // ix = nd + 1
        b2.append("e")
        if nexp < 0 {
            b2.append("-")
        }
        ca = mpdigout (Double(abs (nexp)), 10)
        
        var gk = 9
        for k in 0..<10 {
            let ik = ca.index(ca.startIndex, offsetBy: k)
            if ca[ik] != "0"  { gk = k; break /* goto 190 */ }
        }
        
//        gk = 10
        
        // 190 continue
        
        for i in gk..<10 {
            let ik = ca.index(ca.startIndex, offsetBy: i)
            b2.append(ca[ik])
        }
        
        for _ in 1...nb {
            b2.append(" ")
        }
        
        //   Copy entire b2 array to B.
        b = b2
//        for i in 1...nb {
//            b[i] = b2[i]
//        }
    } // mpeformat

    
    static func mpfformat (_ a : MPRNumber, _ nb : Int, _ nd : Int, _ b : inout String, _ mpnw : Int) {
        
        //   Converts the MPR number A into character form in the character(1) array B.
        //   NB (input) is the length of the output string, and ND (input) is the
        //   number of digits after the decimal point.  The format is analogous to
        //   Fortran F format; the result is right-justified among the NB cells of B.
        //   The condition NB >= ND + 10 must hold or an error message will result.
        //   However, if it is found during execution that there is not sufficient space,
        //   to hold all digits, the entire output field will be filled with asterisks.
        //   NB cells of type character(1) must be available in B.
        
        var ixp, j, mpnw1, nb2, nexp : Int
        var b2, ca : String
        var t1 : Double
        var f = MPRNumber(repeating: 0, count: 8)
        var s0 = MPRNumber(repeating: 0, count: mpnw+6)
        var s1 = s0
        
        // End of declaration
        
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || nb < nd + 10) {
            print ("*** MPFFORMAT: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        // ia = sign (1.0, a[2])
        // na = min (Int (abs (a[2])), mpnw)
        s0[0] = Double(mpnw + 7)
        s1[0] = Double(mpnw + 7)
        mpnw1 = mpnw + 1
        
        //   Set f = 10.
        
        f[0] = 9.0
        f[1] = Double(mpnw1)
        f[2] = 1.0
        f[3] = 0.0
        f[4] = 10.0
        f[5] = 0.0
        f[6] = 0.0
        
        nb2 = nb + 20
        b2 = ""
        mpeformat (a, nb2, nb, &b2, mpnw+1)
        
        //   Trim off trailing blanks.
        b2 = b2.trimmingCharacters(in: .whitespaces)
        //        for i in stride(from:nb2, to:1, by:-1) {
        //            if b2[i] != " " { goto 90 }
        //        }
        //
        //        90 continue
        
        nb2 = b2.count
        
        //   Look for the "e" in B2.
        var k = b2.startIndex
        if let index = b2.firstIndex(of: "e") {
            k = index
        } else {
            //        for k in 1...nb2 {
            //            if b2[k] == "e" { goto 100 }
            //        }
            print ("*** MPFFORMAT: Syntax error in output of mpeformat")
            mpabrt (84)
        }
        
        // 100 continue
        
        //   Check the sign of the exponent.
        k = b2.index(k, offsetBy: 1) //   k = k + 1
        if b2[k] == "-" {
            ixp = -1
            k = b2.index(k, offsetBy: 1) // k = k + 1
        } else {
            ixp = 1
        }
        j = 0
        ca = " "
        
        //   Copy the exponent into CA.
        
        while b2.index(k, offsetBy: 1) < b2.endIndex {
            k = b2.index(k, offsetBy: 1)
            j = j + 1
            if (j <= 16) { ca.append(b2[k]) }
        }
        
        t1 = mpdigin (ca, j)
        
        //   Check if there is enough space in the ouput array for all digits.
        b = ""
        if Int(t1) + nd + 3 > nb {
            b = b.padding(toLength: nb, withPad: "*", startingAt: 0)
//            for i in 1...nb {
//                b[i] = "*"
//            }
            
            // goto 200
        } else {
            nexp = ixp * Int(t1)
            
            //   Insert the sign of the number, if any.
            
            if b2.first! == "-" {
                b.append(b2.removeFirst())
            }
            
            if nexp == 0 {
                
                //   Exponent is zero.  Copy first digit, period and ND more digits.
                
                for _ in 1...nd + 2 {
                    b.append(b2.removeFirst())
                }
                
                //  goto 200
            } else if nexp > 0 {
                
                //   Exponent is positive.  Copy first digit, skip the period, then copy
                //   nexp digits.
                
                b.append(b2.removeFirst())
                b2.removeFirst()            // skip decimal point
                
                for _ in 1...nexp {
                    b.append(b2.removeFirst())
                }
                
                //   Insert the period.
                
                b.append(".")
                
                //   Copy nd more digits.
                
                for _ in 1...nd {
                    b.append(b2.removeFirst())
                }
                
                //goto 200
            } else {
                
                //   Exponent is negative.  Insert a zero, then a period, then nexp - 1
                //   zeroes, then the first digit, then the remaining digits up to ND total
                //   fractional digits.
                
                b.append("0")
                b.append(".")
                
                for _ in 1...nexp - 1 {
                    b.append("0")
                }
                
                b.append(b2.removeFirst())
                b2.removeFirst()            // skip decimal point
                
                for _ in nexp...nd {
                    b.append(b2.removeFirst())
                }
            }
        }
        
        //200 continue
        
        //   Right-justify in field.
        
        let ki = nb - b.count
        
//        for i in 1...i1 {
//            b[nb-i+1] = b[nb-i-ki+1]
//        }
//
        for _ in 1...ki {
            b.insert(" ", at: b.startIndex)
        }
    } // mpfformat

    static func mpinp (_ iu : InputStream, _ a : inout MPRNumber, _ mpnw : Int) {
        
        //   This routine reads the MPR number A from InputStream IU.  The digits of A
        //   may span more than one line, provided that a "\" appears at the end of
        //   a line to be continued (any characters after the "\" on the same line
        //   are ignored).  Individual input lines may not exceed 2048 characters in
        //   length, although this limit can be changed in the system parameters
        //   (parameter mpnstr) in module MPFUNA.  Embedded blanks are allowed anywhere.
        //   An exponent with "e" or "d" may optionally follow the numeric value.
        
        //   A scratch array below (CHR1) holds character data for input to mpctomp.
        //   It is dimensioned MPNW * (MPNDPW + 1) + 1000 (see below).  If more nonblank
        //   input characters than this are input, they are ignored.
        
        var i1, lnc1, lncx, ln1 : Int
        var line1, chr1 : String
        let validc = " 0123456789+-.dDeE"
        
        func get() -> Character {
            // returns 0 for errors and end-of-file
            // Note: Only works for 8-bit ASCII
            var buffer = Array<UInt8>(repeating: 0, count: 1)
            let _ = iu.read(&buffer, maxLength: 1)
            return Character(UnicodeScalar(buffer[0]))
        }
        
        func getLine() -> String {
            var line = ""
            if !iu.hasBytesAvailable {
                print ("*** MPINP: End-of-file encountered.")
                mpabrt (72)
                return ""
            } else {
                // find non-blank character
                var ch : Character
                repeat {
                    ch = get()
                    if !ch.isWhitespace { line.append(ch) }
                } while (!ch.isNewline || line.isEmpty) && iu.hasBytesAvailable
                return line
            }
        }
        
        // End of declaration
        
        if mpnw < 4 || a[0] < Double(mpnw + 6) {
            print ("*** MPINP: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        lnc1 = 0
        lncx = mpnw * (mpndpw + 1) + 1000
        
        //100 continue

//        lab1: while true {
//                read (iu, "(a)", end = 200) line1
//
//                //   Find the last nonblank character.
//
//                for i in mpnstr, 1, -1 {
//                    if (line1[i] != " ") { break lab1 /* goto 110 */ }
//                }
//
//                //   Input line is blank -- ignore.
//
//                // goto 100
//        }
        
        //110 continue
        
        // ln1 = i
        
        //   Scan input line, looking for valid characters.
        chr1 = ""
        line1 = getLine()
        for ch in line1 {
//            if ch == "\\" { line1 = getLine() /* goto 100 */ }
            let x = validc.firstIndex(of: ch)
            if x == nil && !ch.isWhitespace {
                print ("*** MPINP: Invalid input character = '\(ch)'")
                mpabrt (87)
            } else if !ch.isWhitespace {
                chr1.append(ch)
            }
        }
        
        mpctomp (chr1, &a, mpnw)
        
        //300 return
    } // mpinp

    static func mpout (_ iu : OutputStream, _ ln : Int, _ nd : Int, _ a : MPRNumber, _ mpnw : Int) {
        
        //   This routine writes MPR number A to logical unit IU in E(LN,ND) format.
        //   This is output on MPOUTL characters per line.  The value of MPOUTL is set
        //   in the system parameters at the start of module MPFUNA.

        var chr1 = ""
        
        // End of declaration
        
        mpeformat (a, ln, nd, &chr1, mpnw)
        
        //        write (cform1, 1) mpoutl
        //        1 format ("(",i8,"a1)")
        //        write (cform2, 2) mpoutl
        //        2 format ("(",i8,"a1,"\")")
        // convert the numeric string to ascii bytes
        if let d = chr1.data(using: .ascii) {
            // output the data string as bytes
            let bytes = Array<UInt8>(d)
//            if ln <= mpoutl {
            iu.write(bytes, maxLength: ln)
//               write (iu, fmt = cform1) (chr1(i), i = 1, ln)
//            } else if ln % mpoutl == 0 {
//                ln1 = mpoutl * (ln / mpoutl) - mpoutl
//                write (iu, fmt = cform2) (chr1(i), i = 1, ln1)
//                write (iu, fmt = cform1) (chr1(i), i = ln1 + 1, ln)
//            } else {
//                ln1 = mpoutl * (ln / mpoutl)
//                write (iu, fmt = cform2) (chr1(i), i = 1, ln1)
//                write (iu, fmt = cform1) (chr1(i), i = ln1 + 1, ln)
//            }
        } else {
            print ("*** MPOUT: Invalid ASCII character in string = '\(chr1)'")
            mpabrt (87)
        }
    } // mpout

    
    
}
