# -*- coding: utf-8 -*-
import sys
sys.path.append(r'/Users/aitianyi/MPC')
from server import io_process

if __name__ == "__main__":

    print("-------The test of io-process is beginning!-------\n")
    a = [3, 5, 7, 8, 4, 9, 4]  # 即为密文a[]，由dealer生成并分发
    iotest = io_process.ioHandle()
    s1 = "0(2235523s+#1)(434135135s+#2)+(72342345s+23534645768)(434551341s+35222574745)(3524635239s+34592348626)"
    s2 = "1235246029526345$4095204936243$350238460284058"
    slist = iotest.choose(s1, a)
    print("slist = \n", slist)
    print("-------The test of io-process is finishing!-------\n")
