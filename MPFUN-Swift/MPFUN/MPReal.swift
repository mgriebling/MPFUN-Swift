//
//  MPReal.swift
//  MPFUN-Swift
//
//  Created by Mike Griebling on 4 May 2019.
//  Copyright © 2019 Computer Inspirations. All rights reserved.
//

import Foundation

public struct MPReal {
    
    var number : MPFUN.MPRNumber
    
    // default precision level in digits
    public static var digitPrecision = 1200 {
        didSet {
            mpwds = digitPrecision / MPFUN.mpndpw + 2
            mpwds6 = mpwds + 6
        }
    }
    
    // dependent variables that are altered from the precision
    static var mpwds = digitPrecision / MPFUN.mpndpw + 2
    static var mpwds6 = mpwds + 6
    
    init(_ numberWords : Int = mpwds6) {
        //   Set initial number = 0
        number = MPFUN.MPRNumber(arrayLiteral: Double(numberWords), 7, 1, 0, 0, 0, 0)
        number.reserveCapacity(numberWords)
    }
    
    
}
