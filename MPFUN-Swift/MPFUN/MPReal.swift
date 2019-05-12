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
    
    public init(_ numberWords : Int = 0) {
        //   Set initial number = 0
        let words = numberWords == 0 ? MPReal.mpwds6 : numberWords
        number = MPFUN.MPRNumber(repeating: 0, count: words)
        number[0] = Double(words)
        number[1] = 7
        number[2] = 1
    }
    
    init(_ s : String, _ numberWords : Int = 0) {
        self.init(numberWords)
        MPFUN.mpctomp(s, &number, numberWords)
    }
    
    init(_ i: Int, _ numberWords : Int = 0) {
        //   Set f = i
        self.init(numberWords)
        number[0] = 9.0
 //       number[1] = Double(numberWords)
        number[2] = 1.0
        number[3] = 0.0
        number[4] = Double(i)
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

extension MPReal : AdditiveArithmetic {
    
    public static func - (lhs: MPReal, rhs: MPReal) -> MPReal {
        var result = MPReal()
        MPFUN.mpsub(lhs.number, rhs.number, &result.number, MPReal.mpwds6)
        return result
    }
    
    public static func + (lhs: MPReal, rhs: MPReal) -> MPReal {
        var result = MPReal()
        MPFUN.mpadd(lhs.number, rhs.number, &result.number, MPReal.mpwds6)
        return result
    }
    
    public static var zero: MPReal { return MPReal() }
    public static func += (lhs: inout MPReal, rhs: MPReal) { lhs = lhs + rhs }
    public static func -= (lhs: inout MPReal, rhs: MPReal) { lhs = lhs - rhs  }
    
}

extension MPReal : ExpressibleByIntegerLiteral {
    
    public init(integerLiteral value: Int) { self.init(value) }
    public typealias IntegerLiteralType = Int
    
}

extension MPReal : Numeric {
    
    public init?<T>(exactly source: T) where T : BinaryInteger {
        // we'll cheat for now
        self.init(source.description)
    }
    
    public var magnitude: MPReal {
        var mag = MPReal()
        MPFUN.mpcabs(number, &mag.number, MPReal.mpwds6)
        return mag
    }
    
    public static func * (lhs: MPReal, rhs: MPReal) -> MPReal {
        var result = MPReal()
        MPFUN.mpmul(lhs.number, rhs.number, &result.number, MPReal.mpwds6)
        return result
    }
    
    public static func *= (lhs: inout MPReal, rhs: MPReal) { lhs = lhs * rhs }
    public typealias Magnitude = MPReal
    
}

extension MPReal : SignedNumeric { }

extension MPReal : Strideable {
    
    public func distance(to other: MPReal) -> Double {
        return 0  // TBD
    }
    
    public func advanced(by n: Double) -> MPReal {
        return self // TBD
    }
    
    
    public typealias Stride = Double
    
}

extension MPReal : FloatingPoint {
    
    public mutating func round(_ rule: FloatingPointRoundingRule) {
        switch rule {
        case .awayFromZero: break
        case .down: break
        case .toNearestOrAwayFromZero: break
        case .toNearestOrEven: break
        case .towardZero: break
        case .up: break
        @unknown default: assert(false, "MPReal unknown rounding rule")
        }
    }
    
    public init(sign: FloatingPointSign, exponent: Int, significand: MPReal) {
        self.init()
    }
    
    public init(signOf: MPReal, magnitudeOf: MPReal) {
        let sign = signOf.sign
        let mag = magnitudeOf.magnitude
        self.init(sign: sign, exponent: 0, significand: mag)
    }
    
    public init<Source>(_ value: Source) where Source : BinaryInteger {
        self.init(exactly: value)!
    }
    
    public static var radix: Int {
        return 2
    }
    
    public static var nan: MPReal {
        return MPReal()
    }
    
    public static var signalingNaN: MPReal {
        return MPReal()
    }
    
    public static var infinity: MPReal {
       return MPReal()
    }
    
    public static var greatestFiniteMagnitude: MPReal {
        return MPReal()
    }
    
    public static var pi: MPReal {
       return MPReal()
    }
    
    public var ulp: MPReal {
        return MPReal()
    }
    
    public static var leastNormalMagnitude: MPReal {
        return MPReal()
    }
    
    public static var leastNonzeroMagnitude: MPReal {
        return MPReal()
    }
    
    public var sign: FloatingPointSign {
        let ia = MPFUN.sign(1, number[2])
        return ia < 0 ? .minus : .plus
    }
    
    public var exponent: Int {
        return Int(number[3])
    }
    
    public var significand: MPReal {
        return self
    }
    
    public static func / (lhs: MPReal, rhs: MPReal) -> MPReal {
        return MPReal()
    }
    
    public static func /= (lhs: inout MPReal, rhs: MPReal) {
    }
    
    public mutating func formRemainder(dividingBy other: MPReal) {
    }
    
    public mutating func formTruncatingRemainder(dividingBy other: MPReal) {
    }
    
    public mutating func formSquareRoot() {
        MPFUN.mpsqrt(self.number, &self.number, MPReal.mpwds6)
    }
    
    public mutating func addProduct(_ lhs: MPReal, _ rhs: MPReal) {
        self += lhs * rhs
    }
    
    public var nextUp: MPReal {
        return self
    }
    
    public func isEqual(to other: MPReal) -> Bool {
        var code = 0
        MPFUN.mpcpr(number, other.number, &code, MPReal.mpwds6)
        return code == 0
    }
    
    public func isLess(than other: MPReal) -> Bool {
        return self < other
    }
    
    public func isLessThanOrEqualTo(_ other: MPReal) -> Bool {
        return isEqual(to: other) || isLess(than: other)
    }
    
    public func isTotallyOrdered(belowOrEqualTo other: MPReal) -> Bool {
        return true
    }
    
    public var isNormal: Bool {
        return true  // TBD
    }
    
    public var isFinite: Bool {
        return false // TBD
    }
    
    public var isZero: Bool {
        return false
    }
    
    public var isSubnormal: Bool {
        return false
    }
    
    public var isInfinite: Bool {
        return false // TBD
    }
    
    public var isNaN: Bool {
        return false
    }
    
    public var isSignalingNaN: Bool {
        return false
    }
    
    public var isCanonical: Bool {
        return true
    }
    
    public typealias Exponent = Int
    
}



