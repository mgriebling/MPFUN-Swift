
import Foundation

//********************************************************************************
//
// This source is Copyright (c) 2017 by Solinst Canada.  All rights reserved.
//
//********************************************************************************
/**
    This module implements a stream-based input/output facility.  The Stream provides a standard
    set of functions that can read/write various data types from/to an underlying physical medium.
  	This *MStream* module forms the basis of all the memory, UART, and SPI input/output functions.
  	Input/output of all standard data types has been abstracted to this *MStream* implementation
  	which can be reused for all data input/output, no matter what the medium that is connected to
  	the *MStream*.  The *endian* attribute controls whether little- or big-endian reads/writes are
    used so that portability to other processor architectures is simplified.  
 
    Two types of Stream input/output are supported:
 
    - Binary read/write functions where passed arguments are written/read as binary bytes
  	to/from the underlying stream with no conversions performed.  The endian writing settings
  	will be used for writing intrinsic types like *Int*, *Int16*, *Double*, etc.  Note:
  	No endian conversions are possible with the generic byte-based *write* and *read* so use
  	these operations with extreme caution and as a last resort since the resulting stream contents
  	will not be portable across computers with different native endian architectures.
  	The default endian modes are to read/write with big-endian byte ordering.
 
    - Text-based read/write functions where passed arguments are converted to/from
 	text strings and are written/read to/from the underlying stream as a series of character
 	strings.  These conversions will be most useful when reading/writing from/to a serial
 	terminal.  For optimal performance in an embedded system, these functions do not use
 	the C runtime *sprintf* or *printf* conversion functions.
 
 	Be very careful when mixing the binary- and ASCII-stream functions since the two are not
 	interchangeable.  That is if a stream is written as bytes it must be read as bytes and not as
 	an ASCII-stream.
 
  - Author:   Michael Griebling
  - Date:     2 March 2017
 
 ******************************************************************************** */

public class MStream {
    
    public enum errorType { case ok, invalid }
    
    private var lastError : errorType = .ok
    
    public var input : InputStream!   { didSet { input?.open() } }
    public var output : OutputStream! { didSet { output?.open() } }
    public var endian : ByteOrder = .nativeByteOrder
    public var encoding : String.Encoding = .ascii
    
    //********************************************************************************
    /**
         Returns *true* if the input stream has more data for reading.
         
         - Author:  Michael Griebling
         - Date:    11 May 2017
     
     ******************************************************************************** */
    public var hasBytesAvailable : Bool { return input?.hasBytesAvailable ?? false }
    
    //********************************************************************************
    /**
         Returns *true* if the output stream has more space for writing.
         
         - Author:  Michael Griebling
         - Date:    11 May 2017
     
     ******************************************************************************** */
    public var hasSpaceAvailable : Bool { return output?.hasSpaceAvailable ?? false }
    
    //********************************************************************************
    /**
         Initialize the stream as a writer.
         
         - Author:  Michael Griebling
         - Date:    11 May 2017
     
     ******************************************************************************** */
    public init (output : OutputStream?) {
        self.input = nil
        self.output = output
        open()
    }
    
    //********************************************************************************
    /**
         Initialize the stream as a reader.
         
         - Author:  Michael Griebling
         - Date:    11 May 2017
     
     ******************************************************************************** */
    public init (input : InputStream) {
        self.input = input
        self.output = nil
        open()
    }
    
    //********************************************************************************
    /**
         Close any open streams when deallocating the stream.
         
         - Author:  Michael Griebling
         - Date:    11 May 2017
     
     ******************************************************************************** */
    deinit { close() }
    
    //********************************************************************************
    /**
         Return the active stream error.
     
       - Author:   Michael Griebling
       - Date:    6 May 2010
     
     ******************************************************************************** */
    public var error : errorType { return lastError }
    
    //********************************************************************************
    /**
         Clear any error for stream.
     
       - Author:   Michael Griebling
       - Date:    6 May 2010
     
     ******************************************************************************** */
    public func clearError() { lastError = .ok }
    
    //********************************************************************************
    /**
         Close any open streams and go to sleep.
         
         - Author:  Michael Griebling
         - Date:    11 May 2017
     
     ******************************************************************************** */
    public func close() { output?.close(); input?.close() }
    
    //********************************************************************************
    /**
         Open any available streams.
         
         - Author:  Michael Griebling
         - Date:    11 May 2017
     
     ******************************************************************************** */
    public func open()  {
        output?.open();
        let err = output?.streamError
        if err != nil { print("Open error : \(err!)") }
        input?.open()
    }
    
    //********************************************************************************
    /**
         Output a generic argument *arg* with *size* bytes to the stream.  The byte-ordering for
        this function is unspecified since it relies on the underlying hardware read/write functions.
     
       - Author:   Michael Griebling
       - Date:    6 May 2010
     
     ******************************************************************************** */
    public func write(_ arg : [Byte]) {
        guard arg.count > 0 else { return }  // streams fails with a zero-byte write
        assert(output != nil, "Undefined output stream")
        let result = output.hasSpaceAvailable
        assert(result, "No space left in output stream")
        
        let size = output.write(arg, maxLength: arg.count)
        if size != arg.count { lastError = .invalid }
    }
    
    //********************************************************************************
    /**
         Output a generic *PackableType* *T* *arg* to the stream.  The endianness 
         of the write is controlled by the object's *endian* attribute.
         
         - Author:   Michael Griebling
         - Date:    6 May 2010
     
     ******************************************************************************** */
    public func write<T:PackableType> (_ arg : T) {
        write(BytePacker.pack(arg, byteOrder: endian))
    }
    
    //********************************************************************************
    /**
         Output a generic *PackableType* array of *T* *arg* to the stream.  The endianness
         of the write is controlled by the object's *endian* attribute.
     
         - Author:  Michael Griebling
         - Date:    11 May 2017
     
     ******************************************************************************** */
    public func write<T:PackableType> (_ arg : [T]) {
        write(BytePacker.pack(arg, byteOrder: endian))
    }

    private static let CBIT = 128
    private static let CMASK = CBIT-1
    
    //********************************************************************************
    /**
         Output an integer *lint* to the stream in compressed and portable format.
     
           - Author:   Michael Griebling
           - Date:    6 May 2010
     
     ******************************************************************************** */
    public func writeNum(_ lint : Int) {
        var lint = lint
        while lint < -64 || lint > 63 {
            write(UInt8((lint & MStream.CMASK) + MStream.CBIT))
            lint >>= 7
        }
        write(UInt8(lint & MStream.CMASK))
    } /* WriteNum */
    
    //********************************************************************************
    /**
        Output the variable-length string *str* to the stream using *encoding*.  If
        *noCount* is true, only the string is output without a leading byte count.
     
       - Author:  Michael Griebling
       - Date:    6 May 2010
     
     ******************************************************************************** */
    public func writeString(_ str : String, noCount: Bool = true) {
        let s = str.getBytes(encoding: encoding)
        if !noCount { writeNum(s.count) }
        write(s)
    }
    
    //********************************************************************************
    /**
         Returns an array of *size* bytes from the stream.  The byte-ordering for 
         this function is unspecified since it relies on the underlying hardware 
         read/write functions of the stream.
     
       - Author:   Michael Griebling
       - Date:    6 May 2010
     
     ******************************************************************************** */
    public func read(_ size : Int) -> [UInt8] {
        assert(input != nil, "Undefined input stream")
        assert(input.hasBytesAvailable, "No data left in input stream")
        
        var buffer = [UInt8](repeating: 0, count: size)
        let status = input.read(&buffer, maxLength: size)
        if status >= 0 && status < size { buffer = Array(buffer.dropLast(size - status)) }
        if status == -1 { lastError = .invalid }
        return buffer
    }
    
    //********************************************************************************
    /**
         Read a variable length string from the stream with *encoding*. Returns
         an empty string if the number of bytes does not match the selected encoding.
         This function can read strings written by *write(String)*. If the *size*
         is zero, it is assumed the size is part of the data in the input stream
         as written by *writeString* with *noCount* set *false*.
     
       - Author:   Michael Griebling
       - Date:    6 May 2010
     
     ******************************************************************************** */
    public func readString(size : Int = 0) -> String {
        let length = size == 0 ? readNum() : size
        let bytes = read(length)
        return String(bytes, encoding: encoding)
    }

    //********************************************************************************
    /**
         Read a generic *PackableType* type *T* *arg* from the stream.  An optional
         *endian* argument supports either little- or big-Endian reads.  By default
         native endian reads are performed.
         
         - Author:   Michael Griebling
         - Date:    6 May 2010
     
     ******************************************************************************** */
    public func read<T:PackableType>() -> T {
        let size = MemoryLayout<T>.size
        let s = read(size)
        return BytePacker.unpack(s, byteOrder: endian)
    }
    
    //********************************************************************************
    /**
         Return the uncompressed integer that was read from the stream.
     
       - Author:   Michael Griebling
       - Date:    6 May 2010
     
     ******************************************************************************** */
    public func readNum () -> Int {
        var s : Int
        var ch : Int
        var n : Int
        
        func readByte() -> Int { let x : Byte = read(); return Int(x) }
        
        s=0; n=0; ch=readByte()
        while ch >= MStream.CBIT {
            n += Int(ch-MStream.CBIT) << s; s += 7
            ch=readByte()
        }
        
        /* output the expanded number */
        return (n + (((ch & 63) - ((ch >> 5) << 5)) << s))
    }
    
    
    //********************************************************************************
    /**
         Return the decimal number corresponding to the hexadecimal input character 
         *ch*.
     
       - Author:   Michael Griebling
       - Date:    6 May 2010
     
     ******************************************************************************** */
    private func toHex (_ ch : UnicodeScalar) -> Int {
        if let x = Int(String(ch), radix: 16) { return x }
        return 0
    }
    
    //********************************************************************************
    /**
         Return the cardinal number that was read as a hexadecimal number string 
         with *len* characters from the stream.
     
       - Author:   Michael Griebling
       - Date:    6 May 2010
     
     ******************************************************************************** */
    public func readHexString(len : Int) -> Int {
        let str : String = read()
        if str.count != len { return 0 }
        return str.unicodeScalars.reduce(0) { (result, ch) in result<<4 | toHex(ch) } // convert string to hex number
    } /* end readHex() */
    
    //********************************************************************************
    /**
         	Output an integer string representation of *int* to the stream using at least *size* characters.
     		The string will be padded with spaces if *size* is greater than the resultant string.
     		A *size* of 0 will produce a minimum length string for the integer *int*.
     
       - Author:   Michael Griebling
       - Date:    6 May 2010
     
     ******************************************************************************** */
    public func writeString(_ int : Int, size : Int) {
    	// convert integer to a string
        let size = max(size, 0)  // eliminate negatives
        let str = String(format: size == 0 ? "%d" : "%-\(size)d", int)
        write(str)
    }
    
    //********************************************************************************
    /**
         Output a cardinal string representation of *uint* to the stream using at least *size* characters.
         The string will be padded with spaces if *size* is greater than the resultant string.
     	 A *size* of 0 (default) will produce a minimum length string for the cardinal *uint*.
     
       - Author:   Michael Griebling
       - Date:    6 May 2010
       - Bug:     10 Mar 2015 - MG - Corrected output for large unsigned numbers.
     
     ******************************************************************************** */
    public func writeString(_ uint : UInt, size : Int = 0) {
        let size = max(size, 0)  // eliminate negatives
        let str = String(format: size == 0 ? "%d" : "%-\(size)d", uint)
        write(str)
    }
    
    //********************************************************************************
    /**
         Output an unsigned integer *hex* to the stream as a hexadecimal string with 
         *size* characters.  The maximum hex string size is sixteen characters.
     
       - Author:   Michael Griebling
       - Date:    6 May 2010
     
     ******************************************************************************** */
    public func writeHexString(_ hex : UInt, digits: Int = MemoryLayout<UInt>.size * 2) {
        let maxSize = MemoryLayout<UInt>.size * 2    // two hex digits per byte
        var size = min(digits, maxSize)
        if size <= 0 { size = maxSize }

        let buf = String(format: "%0\(size)X", hex)
        write(buf)
    }
    
    //********************************************************************************
    /**
        	Output a date/time string representation of *dtime* using a formatting 
            string *format* to the stream. For example, fstr = "dd/MM/yyyy HH:mm:ss" 
            produces an output of "DD/MM/YYYY HH:MM:SS".
     
       - Author:   Michael Griebling
       - Date:    6 May 2010
     
     ******************************************************************************** */
    public func writeString(_ date: Date, usingFormat format: String) {
        let f = DateFormatter()
        f.timeZone = TimeZone(secondsFromGMT: 0)
        f.dateFormat = format
     	write(f.string(from: date))
     }
    
}


//********************************************************************************
/**
     A memory-based writer that writes all its output to a memory buffer.  The
     written memory can be accessed using the *data* property.  If it is
     necessary to use overwrites, please use the *RandomWriterStream* instead.
     Use the inherited *MStream* object's methods to perform writes.
     
     - Author:   Michael Griebling
     - Date:    11 May 2017
 
 ******************************************************************************** */
public class MemoryWriterStream : MStream {
    
//    /*******************************************************************************
//      Deprecated methods - please use RandomWriteStream instead
//     *******************************************************************************/
//    public func overwrite (_ arg : [Byte], at offset : Int) {
//        var buffer = data
//        buffer.replaceSubrange(offset..<offset+arg.count, with: arg)
//        
//        close(); open()
//        _ = write(buffer)
//    }
//
//    /*******************************************************************************
//     Deprecated methods - please use RandomWriteStream instead
//     *******************************************************************************/
//    public func overwrite (_ arg : [Byte], inRange offsetRange : CountableRange<Int> ) {
//        var buffer = data
//        buffer.replaceSubrange(offsetRange, with: arg)
//        close(); open()
//        _ = write(buffer)
//    }
//    
//    /*******************************************************************************
//     Deprecated methods - please use RandomWriteStream instead
//     *******************************************************************************/
//    public func removeBytes(from offset : Int) {
//        var buffer = data
//        let size = data.count
//        buffer.removeSubrange(offset..<size)
//        close()
//        open()
//        _ = write(buffer)
//    }
// 
//    /*******************************************************************************
//     Deprecated methods - please use RandomWriteStream instead
//     *******************************************************************************/
//    public func searchSequence(from sequence : [Byte] ) -> Int {
//        var index = 0
//        var buffer = data
//        
//        while index < (buffer.count - sequence.count) {
//            let subarray = buffer[index ..< (index + sequence.count)]
//            
//            if subarray.elementsEqual(sequence) {
//                return index
//            }
//            
//            index += 1
//        }
//        return buffer.count
//    }
    
    // Other methods
    public var data : [UInt8] {
        if let buffer = output?.property(forKey: Stream.PropertyKey.dataWrittenToMemoryStreamKey) as? Data {
            return [UInt8](buffer)
        }
        return []
    }
    
    public init() {
        super.init(output: OutputStream.toMemory())
    }
    
}

//********************************************************************************
/**
     A file-/memory-based writer stream that allows positionable writes.  
     For example, the position can be set to 10 and the following stream writes 
     will start at that position in the file.  If a position is set that is larger 
     than the current file length, the position is set to the maximum in the 
     current file.  Positions are based on byte indices.  The current file contents 
     are available using the *data* property.
 
     If no file *url* is passed during object creation, an internal memory stream
     of the file will be created which is accessible using *data*; otherwise, an
     an external file will be created at the *url* after the *close*.  The
     writer is cleared after a successful *close*.
 
     It is possible to insert new bytes into the stream at the current *position*
     by setting *insertBytes* to the number of bytes that are to be inserted.  For
     example, if replacing original contents like "123" at *position* with a new 
     string "456789", set *insertBytes* to the difference in lengths or 3 and write 
     the string "456789" to *position*.  Three extra bytes are inserted and the 
     remaining three characters overwrite the existing file contents.  
     **Note: Any stored file positions occuring after the insertion point are 
     invalidated by the insertion and need to be adjusted accordingly.**
     
     - Author:  Michael Griebling
     - Date:    30 June 2017
 
 ******************************************************************************** */
public class RandomWriterStream : MStream {
    
    private var file = [Byte]()
    private var url : URL?
    
    /// If non-zero, this many bytes are inserted into the writer stream.  If the
    /// stream is at the end of the file, *insertBytes* has no effect but will be
    /// adjusted until it is zero.
    public var insertBytes : Int = 0
    
    public var position : Int = 0 {
        didSet {
            // clamp position to legal limits
            position = min(max(0, position), file.count)
        }
    }
    
    public var data : [Byte] { return file }
    public override var hasSpaceAvailable: Bool { return true }
    
    public override func write(_ arg: [Byte]) {
        let len = arg.count
        if position == file.count {
            // append to the end of the file
            file += arg; position += len
            insertBytes = max(0, insertBytes - len)
        } else if position+len <= file.count {
            if insertBytes >= len {
                // insert len bytes
                file.insert(contentsOf: arg, at: position)
                insertBytes -= len
            } else if insertBytes > 0 {
                // insert *insertBytes* length of arg
                file.insert(contentsOf: arg[0..<insertBytes], at: position)
                
                // overwrite the rest of the arg
                file.replaceSubrange(position+insertBytes..<position+len, with: arg[insertBytes..<len])
                insertBytes = 0
            } else {
                // overwrite *arg* to the file at *position*
                file.replaceSubrange(position..<position+len, with: arg)
            }
            position += len
        } else {
            // overwrite with file extension
            let overlap = file.count - position
            write(Array(arg[0..<overlap]))
            write(Array(arg[overlap..<len]))
        }
    }
    
    /// Using a Boyer-Moore search algorithm (see *Array+Find* array extension in *SolinstUtils*)
    public func find (pattern: [Byte]) -> Int {
        if let pos = file.find(pattern: pattern) {
            return pos
        }
        return file.count
    }
    
    public override func close() {
        if let url = url, let output = FileWriterStream(url: url, append: false) {
            output.open()
            output.write(file)
            output.close()
            
            // clear memory contents
            file = []
            position = 0
            self.url = nil
        }
    }
    
    public init(url : URL? = nil) {
        super.init(output: nil)
        self.url = url
    }
    
}

//********************************************************************************
/**
     A writer stream that is connected to the standard output stream.  All writes
    appear on the debug console. Use the inherited *MStream* object's methods 
    to perform writes.
 
    This class provides an example of overriding the *write* and *hasSpaceAvailable*
    methods to extend the *MStream* base class.  When using a method override, pass
    nil to the *init* method for the *OutputStream*.  Normally, the *OutputStream*
    provides these two methods.
 
     - Author:   Michael Griebling
     - Date:    11 May 2017
 
 ******************************************************************************** */
public class StandardWriterStream : MStream {
    
    public override var hasSpaceAvailable: Bool {
        return true
    }
    
    public override func write(_ arg: [Byte]) {
        // try to convert buffer to a string
        let s = String(arg, encoding: encoding)
        print(s, terminator: "")
    }
    
    public init() {
        super.init(output: nil)
    }
    
}

//********************************************************************************
/**
     A memory-based reader that reads from the initialized *data* passed during
     object creation.  This reader consumes the data as it is read.  It is not
     possible to reset the stream.  Create a new stream if this is required.
     Use the inherited *MStream* object's methods to perform reads.
     
     - Author:   Michael Griebling
     - Date:    11 May 2017
 
 ******************************************************************************** */
public class MemoryReaderStream : MStream {
    
    public init?(data: Data) {
        super.init(input: InputStream(data: data))
    }
    
}

//********************************************************************************
/**
     A file-based reader that reads from the *url* passed during
     object creation.  This reader consumes the data as it is read.  It is not
     possible to reset the stream.  Create a new stream if this is required.
     Use the inherited *MStream* object's methods to perform reads.
 
     - Author:   Michael Griebling
     - Date:    11 May 2017
 
 ******************************************************************************** */
public class FileReaderStream : MStream {

    public init?(url : URL) {
        if let input = InputStream(url: url) {
            super.init(input: input)
        } else {
            return nil
        }
    }
}

//********************************************************************************
/**
     A file-based writer that writes all its output to the *url* passed during
     object creation.  The writer stream must be closed using *close()* before 
     the file is created.  Enabling the *append* option during object creation
     allows appending to the end of an existing file.   
     Use the inherited *MStream* object's methods to perform writes.
     Default overwrites any existing file.
     
     - Author:   Michael Griebling
     - Date:    11 May 2017
 
 ******************************************************************************** */
public class FileWriterStream : MStream {
    
    public init?(url : URL, append: Bool = false) {
        if let output = OutputStream(url: url, append: append)  {
            super.init(output: output)
        } else {
            return nil
        }
    }
}

