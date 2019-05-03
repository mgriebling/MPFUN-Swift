//
//  BytePacker.swift
//  ProgramLTCE
//
//  Created by Michael Griebling on 28Jun2016.
//  Copyright © 2016 Solinst Canada. All rights reserved.
//

import Foundation

public typealias Byte = UInt8

//********************************************************************************
/**
     Defines the data endianness which can be either little- or big-endian.  The
     *nativeByteOrder* is this computer's native byte ordering.
     
     - Author:  Michael Griebling
     - Date:   	28 June 2016
 
 ******************************************************************************** */
public enum ByteOrder {
    case bigEndian
    case littleEndian
    
    /// Machine specific byte order
    static public let nativeByteOrder: ByteOrder = (Int(CFByteOrderGetCurrent()) == Int(CFByteOrderLittleEndian.rawValue)) ? .littleEndian : .bigEndian
}

//********************************************************************************
/**
     Native-endian raw byte-based initialization and extraction.  The *init* creates
     the associated class while the *bytes* attribute returns the bytes making up
     the class.  
 
     Note: It is possible to conform to this protocol for more complex types such
     as *String* (see String PackableType conformance).
     
     - Author:  Michael Griebling
     - Date:   	28 June 2016
 
 ******************************************************************************** */
public protocol PackableType {
    
    init(_: [Byte])
    init()                            // to allow zero-based initialization
    var bytes : [Byte] { get }
    var bytesReversed : Self { get }  // allow base types to opt out of endian reversals
    
}

extension Int : PackableType {}
extension Int8 : PackableType {}
extension Int16 : PackableType {}
extension Int32 : PackableType {}
extension Int64 : PackableType {}

extension UInt : PackableType {}
extension UInt8 : PackableType {}
extension UInt16 : PackableType {}
extension UInt32 : PackableType {}
extension UInt64 : PackableType {}
extension Double : PackableType {}
extension Float : PackableType {}
extension Array : PackableType {}

//********************************************************************************
/**
     This is an example of how a more complex type can conform to the *PackableType*
     protocol. The initialization uses the *String*'s byte-based initialization 
     while the *bytes* attribute returns the bytes from the *String*'s *Data* encoding.
     The *bytesReversed* attribute opts out of byte reversal due to endian
     changes by returning *self*.
     
     - Author:  Michael Griebling
     - Date:   	28 June 2016
 
 ******************************************************************************** */
extension String : PackableType {
    
    public init(_ bytes : [Byte], encoding: String.Encoding = .ascii) {
        self = String(bytes: bytes, encoding: encoding) ?? ""
    }
    public func getBytes(encoding: String.Encoding = .ascii) -> [Byte] {
        return (self.data(using: encoding) ?? Data()).bytes
    }
    public var bytes : [Byte]         { return getBytes() }
    public var bytesReversed : String { return self }  // no endian change
    
}

//********************************************************************************
/**
     The base implementation of the *PackableType* which applies to any types
     without specific overriding behaviour.
     
     - Author:  Michael Griebling
     - Date:   	28 June 2016
 
 ******************************************************************************** */
public extension PackableType {
    
    /// Initializes the PackableType self with the byte array data in *bytes*.
    ///
    /// - Parameter bytes: used to initialize self in a native-endian format.
    init(_ bytes: [Byte]) {
        self = bytes.withUnsafeBufferPointer({ $0.baseAddress!.withMemoryRebound(to: Self.self, capacity: 1) { $0.pointee } })
    }
    
    var bytesReversed : Self {
        return Self(self.bytes.reversed())
    }
    
    /// Returns the native-endian byte representation of self.
    var bytes: [Byte] {
        var value = self
        let valueArray = withUnsafePointer(to: &value) {
            Array(UnsafeBufferPointer(start: $0.withMemoryRebound(to: Byte.self, capacity: 1){$0}, count: MemoryLayout<Self>.size))
        }
        return valueArray
    }
   
}

//********************************************************************************
/**
     A static class encapsulation of the *PackableType* byte packing and unpacking.
     This class helps to make explicit whenever a *pack* or *unpack* operation
     requires more inputs such as which byte ordering to use.
 
     Note the difference between pack/unpack of a base type and an array of that
     base type.  Trying to pack/unpack an array will not work correctly without
     the use of this class.
     
     - Author:  Michael Griebling
     - Date:   	28 June 2016
 
 ******************************************************************************** */
public class BytePacker {
    
    public class func pack<T:PackableType> (_ value : T, byteOrder : ByteOrder = .nativeByteOrder) -> [Byte] {
        return byteOrder == .littleEndian ? value.bytes : value.bytesReversed.bytes
    }
    
    public class func pack<T:PackableType> (_ values : [T], byteOrder : ByteOrder = .nativeByteOrder) -> [Byte] {
        var valueArray = [Byte]()
        valueArray.reserveCapacity(values.count*MemoryLayout<T>.size)
        for item in values {
            let bytes = pack(item, byteOrder: byteOrder)  // with byte swapping
            valueArray.append(contentsOf: bytes)
        }
        return valueArray
    }
    
    
    public class func unpack<T:PackableType> (_ bytes: [Byte], byteOrder : ByteOrder = .nativeByteOrder) -> T {
        let x = T(bytes)
        return byteOrder == .littleEndian ? x : x.bytesReversed
    }
    
    public class func unpack<T:PackableType> (array : [Byte], byteOrder : ByteOrder = .nativeByteOrder) -> [T] {
        let typeSize = MemoryLayout<T>.size
        let pointCount = array.count/typeSize
        var result : [T] = []
        result.reserveCapacity(pointCount)
        for ptr in (0..<pointCount).map({ UnsafePointer<Byte>(array) + typeSize*$0 }) {
            let x = ptr.withMemoryRebound(to: T.self, capacity: 1, { $0.pointee })
            let y = byteOrder == .littleEndian ? x : x.bytesReversed
            result.append(y)
        }
        return result
    }
    
}

public extension Data {
    
    var bytes : [Byte] { return [Byte](self) }
    
}



