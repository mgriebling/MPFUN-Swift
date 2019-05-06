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
        var a = MPReal()
        var b = ""
        
        MPFUN.mpctomp("1.0e100", &a.number, 10)
        MPFUN.mpfformat(a.number, 140, 10, &b, 10)

    }


}

