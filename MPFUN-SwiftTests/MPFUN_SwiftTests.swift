//
//  MPFUN_SwiftTests.swift
//  MPFUN-SwiftTests
//
//  Created by Mike Griebling on 28 Apr 2019.
//  Copyright © 2019 Computer Inspirations. All rights reserved.
//

import XCTest
@testable import MPFUN_Swift

class MPFUN_SwiftTests: XCTestCase {

    override func setUp() {
        // Put setup code here. This method is called before the invocation of each test method in the class.
    }

    override func tearDown() {
        // Put teardown code here. This method is called after the invocation of each test method in the class.
    }
    
    private func checkdp (_ t1: MPFUN.MPRNumber, _ d1: Double, _ ndp: Int) {
//        let dtol = 1e-14
        
//        if (abs ((Double(t1) - d1) / d1) > dtol) {
//            write (6, '(a)') 'error:'
//            mpwrite (6, ndp + 20, ndp, t1)
//            write (6, '(1p,d30.16)') d1
//            stop
//        }
    }
    
    private func checkdc (_ z1: MPFUN.MPRComplex,  _ dc1: Complex64, _ ndp: Int) {
        
//        let dtol = 1e-14
        
//        if (abs ((Complex64(z1) - dc1) / dc1) > dtol) {
//            write (6, '("error")')
//            mpwrite (6, ndp + 20, ndp, z1)
//            write (6, '(1p,2d30.16)') dc1
//            stop
//        }
    }
    
    private func checkl (_ l1: Bool, _ l2: Bool) {
        XCTAssertTrue(l1 == l2, "error: \(l1) \(l2)")
    }

    func testExample() {
        // This is an example of a functional test case.
        // Use XCTAssert and related functions to verify your tests produce the correct results.
    }

    func testPerformanceExample() {
        // This is an example of a performance test case.
        self.measure {
            // Put the code you want to measure the time of here.
        }
    }

}
