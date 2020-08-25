# -*- coding: utf-8 -*-
import gmpy2
import numpy as np
import global_variable
from pypbc import Element,G1
q = 730750818665451621361119245571504901405976559617

class res_generate():    
    def Gen_response(self,slist,split_ID,hksplit_ID):#t代表多项式的项数，n代表分片的数量
        pairing = global_variable.get_pairing()
        g1 = global_variable.get_g1()
        self.__response = gmpy2.mpz(0)
        coef_hk = dict()
        for i in range(len(slist)-2):
            temp = slist[i]
            temp1 = hksplit_ID[len(slist)-1-i]
            temp2 = gmpy2.mul(temp , temp1)
            temp2 = gmpy2.f_mod(temp2,q)
            coef_hk[i+2] = temp2
        #存m个f(i)+a[i]连乘的结果
        fk_ID = (np.polyval(slist,split_ID))
        fk_ID = fk_ID % q
        for l in coef_hk:
            coefhk = coef_hk[l]
            fk_ID = gmpy2.sub(fk_ID , coefhk)
            fk_ID = fk_ID % q
        self.__response = fk_ID

        
        self.__comm_response_ID = Element(pairing, G1)
        self.__comm_response_ID = g1 ** int(self.__response)
        return self.__comm_response_ID,self.__response