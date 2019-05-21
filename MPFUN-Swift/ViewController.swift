//
//  ViewController.swift
//  MPFUN-Swift
//
//  Created by Mike Griebling on 28 Apr 2019.
//  Copyright © 2019 Computer Inspirations. All rights reserved.
//

import UIKit

class ViewController: UIViewController {

    override func viewDidLoad() {
        super.viewDidLoad()
        // Do any additional setup after loading the view.
        
        // check out MPFUN functions
        MPReal.digitPrecision = 500
        let t1 = MPReal.pi
        let t2 = -MPReal.log2
        let e1 = 3141.0/8192.0
        let e2 = 6931.0/8192.0
        let i1 = 5
//        let one = MPReal(1, 0)
        let half = MPReal(0.5, 0)
        let two = MPReal(2, 0)
        let c = half*t1
        let d = MPReal.exp(half)
        print(e1, e2, c, "\nt1+e2 = \(t1+e2)\n", "e1+t2 = \(e1+t2)\n", "t1-e2 = \(t1-e2)\n", "e1-t2 = \(e1-t2)\n")
        print("\nt1*t2 = \(t1*t2)\n", "t1*e2 = \(t1*e2)\n", "e1*t2 = \(e1*t2)\n")
        print("\nt1/t2 = \(t1/t2)\n", "t1/e2 = \(t1/e2)\n", "e1/t2 = \(e1/t2)\n")
        print("\nt1**i1 = \(t1**i1)\n", "t1**t2 = \(t1**t2)\n")
        print("\nabs(t1) = \(abs(t2))\n sqrt(2) = \(sqrt(two))\n exp(0.5) = \(d)")

    }
    

}

