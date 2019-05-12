//
//  MPReal.swift
//  MPFUN-Swift
//
//  Created by Mike Griebling on 4 May 2019.
//  Copyright © 2019 Computer Inspirations. All rights reserved.
//

import Foundation

public struct MPReal : Codable {
    
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
        number = MPFUN.MPRNumber(repeating: 0, count: numberWords)
        number[0] = Double(numberWords)
        number[1] = 7
        number[2] = 1
    }
    
    init(_ s : String, _ numberWords : Int = mpwds6) {
        self.init()
        MPFUN.mpctomp(s, &number, numberWords)
    }
    
}


extension MPReal : CustomStringConvertible {
    
    /// Support for conversion to string
    public var description: String {
        var s = ""
        MPFUN.mpeformat(number, 50, 40, &s, MPReal.mpwds6)
        return s
    }
    
}

extension MPReal : CustomDebugStringConvertible {
    
    /// Support for conversion to debug string
    public var debugDescription: String {
        var s = ""
        for i in 0..<number.count {
            s += "\(number[i]) "
        }
        return s
    }
    
}

extension MPReal : ExpressibleByStringLiteral {
    
    // allows things like "let a : MPReal = "1.234"
    public typealias ExtendedGraphemeClusterLiteralType = StringLiteralType
    public typealias UnicodeScalarLiteralType = Character
    public typealias FloatLiteral = Double
    public init (stringLiteral s: String) { self.init(s) }
    public init (extendedGraphemeClusterLiteral s: ExtendedGraphemeClusterLiteralType) { self.init(stringLiteral:s) }
    public init (unicodeScalarLiteral s: UnicodeScalarLiteralType) { self.init(stringLiteral:"\(s)") }
    
}

extension MPReal : Comparable {
    
    public static func < (lhs: MPReal, rhs: MPReal) -> Bool {
        var code = 0
        MPFUN.mpcpr(lhs.number, rhs.number, &code, MPReal.mpwds6)
        return code == -1
    }
    
}

extension MPReal : Hashable { }





