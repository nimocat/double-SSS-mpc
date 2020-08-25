# -*- coding: utf-8 -*-
import sys
sys.path.append(r'/Users/aitianyi/MPC')
from dealer import setup
from server import io_process
from server import response
from requester import recover
from verify import verify_cc
from verify import verify_res
import gmpy2


if __name__ == "__main__":

    print("-------The test of setup is beginning!-------\n")
    s = setup.setup()
    userid = dict()
    a = list()
    f = list()
    split = list()
    fk = dict()
    ck = dict()
    hk = dict()
    hksplit = dict()
    commitment_s = dict()
    commitment_fx_joint = str
    commitment_hk_joint = dict()
    commitment_ck_joint = dict()
    commitment_fr_k = dict()
    commitment_split = dict()
    commitment_hk_split = dict()
    userid, a, f, fk, ck, hk, split, hksplit, commitment_s, commitment_fx_joint, commitment_hk_joint, commitment_ck_joint, commitment_fr_k, commitment_split, commitment_hk_split = s.setup(
        3, 7, 10, 6)

    print("test userid = \n", userid)
    print("test a = \n", a)
    print("test f = \n", f)
    print("test fk = \n", fk)
    print("test ck = \n", ck)
    print("test hk = \n", hk)
    print("test split = \n", split)
    print("test hksplit = \n", hksplit)
    print("test commitment_s = \n", commitment_s)
    print("test commitment_fx_joint = \n", commitment_fx_joint)
    print("test commitment_hk_joint = \n", commitment_hk_joint)
    print("test commitment_ck_joint = \n", commitment_ck_joint)
    print("test commitment_fr_k = \n", commitment_fr_k)
    print("test commitment_split = \n", commitment_split)
    print("test commitment_hk_split = \n", commitment_hk_split)
    print("-------The test of setup is finishing!-------\n")

    print("-------The test of io-process is beginning!-------\n")
    iotest = io_process.ioHandle()
    s1 = "0(2235523s+#1)(242354s+543465436)(363635s+432523634)(34687356456s+4645334563462)(6457347s+4656346435)(467386234s+3785972365)(6347345743s+634745758)(3453463457s+89082894265)(83947689247s+523645724574)(346234s+64573444444)"
    sdic = iotest.choose(s1, a, 10)
    print("sdic = \n", sdic)
    print("-------The test of io-process is finishing!-------\n")

    print("-------The test of response is beginning!-------\n")
    res = response.res_generate()
    response_poly_mul = gmpy2.mpz()

    response_sum_comm, response_poly_mul = res.Gen_response(
        sdic, split[1], hksplit[1])
    print("response_sum_comm = ", response_sum_comm)
    print("response_poly_mul = ", response_poly_mul)

    print("-------The test of response is finishing!-------\n")

    print("-------The test of recover is beginning!-------\n")
    response_split = dict()
    for i in split:
        temp1, temp = res.Gen_response(sdic, split[i], hksplit[i])
        response_split[i] = temp
    print("response_split = ", response_split)
    rec = recover.recover()
    lag, secret_of_poly = rec.recover_poly(response_split, 3, userid)

    print("the recovered poly is :\n", lag)
    print("the secret of poly is :\n", secret_of_poly)

    print("-------The test of recover is finishing!-------\n")

    print("-----The test of verify commitment is beginning!-----\n")
    ver = verify_cc.verify_cc()
    ver_result = ver.verify(10, 3, commitment_s, commitment_fx_joint, commitment_fr_k,
                            commitment_ck_joint, commitment_hk_joint, commitment_split, commitment_hk_split, 1)
    if ver_result == True:
        print("Verification of all commitments is successful!\n")
    ver1 = verify_res.verify_response()

    response_coef_comm = dict()
    response_coef_comm = ver1.verify_response_comm(
        sdic, 3, 1, commitment_fx_joint, commitment_s, commitment_ck_joint, response_sum_comm)
    print("-----The test of verify commitment is finishing!-----\n")
