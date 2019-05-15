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
        let a = MPReal.pi
        let b = -MPReal.log2
        let half = MPReal(0.5, 0)
        let c = half*a
        let d = MPReal.exp(half)
        print(a, b, c, d)

    }

}

