# -*- coding: utf-8 -*-
"""
本函数用于实现从收到的response分片中恢复出秘密
"""
import gmpy2
import numpy as np

q = 730750818665451621361119245571504901405976559617

class recover():
    def recover_poly(self,split,t,userid):
        poly_coef_up = dict()
        poly_coef_down = list()
        poly_coef = dict()
        qm = gmpy2.mpz(730750818665451621361119245571504901405976559617)
        for i in range(t):
            coef_of_variable = list()
            coef_of_constant = list()
            akcoef1 = list()
            temp1 = 1
            for j in range(t):  
                if j != i:
                    coef_of_variable.append(1)
                    coef_of_constant.append(-(userid[j+1]))
            coef = self.__F_x(len(coef_of_variable),coef_of_variable,coef_of_constant)
            poly_coef_up[i+1] = coef
            for m in range(t):
                akcoef1.append(userid[m+1])
            for m in range(t):
                if m != i:
                    temp1 = gmpy2.mul(temp1,(gmpy2.sub(akcoef1[i],akcoef1[m])))
                    temp1 = gmpy2.f_mod(temp1,q)
            poly_coef_down.append(temp1)
        for i in range(t):
            poly_coef_up_i = poly_coef_up[i+1]
            poly_coef_down_i = poly_coef_down[i]
            temp4 = list()
            for j in range(len(poly_coef_up_i)):
                poly_coef_up_i[0] = gmpy2.mpz(1)
                temp2 = poly_coef_up_i[j]
                temp3 = gmpy2.divm(temp2,poly_coef_down_i,qm)
                temp4.append(temp3)
            poly_coef[i+1] = temp4
        split_recover = list()#存用于恢复多项式的分片
        __recover_poly_coef = list()#恢复后的多项式系数
        split_mul_poly = dict()
        m = 0
        for i in range(t): 
            m+=1
            if m<=3:
                split_recover.append(split[akcoef1[i]])
        for i in range(t):
            temp7 = list()
            temp5 = poly_coef[i+1]
            for j in range(len(temp5)):
                temp6 = gmpy2.mul(split_recover[i],temp5[j])
                temp6 = gmpy2.f_mod(temp6,q)
                temp7.append(temp6)
            split_mul_poly[i] = temp7
        __recover_poly_coef = split_mul_poly[0]
        for i in range(t):
            if i !=0:
                temp = split_mul_poly[i]
                for j in range(len(__recover_poly_coef)):
                    __recover_poly_coef[j] = gmpy2.add(__recover_poly_coef[j],temp[j])
                    __recover_poly_coef[j] = gmpy2.f_mod(__recover_poly_coef[j],q)
        __secret_of_poly = np.polyval(__recover_poly_coef,0)#由恢复出来的多项式得到的秘密值
        __secret_of_poly = gmpy2.f_mod(__secret_of_poly,qm)
        return __recover_poly_coef,__secret_of_poly
        
    
    def __F_x(self,k,lsc,lc):#lsc是s的系数列表，lc是除s外的系数列表
        Fk = list()
        self.__Fkg = list()
        qm = gmpy2.mpz(730750818665451621361119245571504901405976559617)
        for i in range(k):
            temp1 = gmpy2.divm(1,lsc[i],qm)
            lc[i] = gmpy2.mul(lc[i],temp1)
            lsc[i] = gmpy2.mpz(1)
        for i in range(k):
            p = np.poly1d([1,lc[i]])
            Fk.append(p)
        temp = Fk[0]
        for j in range(1,k):
           temp = temp * Fk[j]
        self.__Fkg = list(temp.coef)
        for n in range(k):
            self.__Fkg[i] = gmpy2.mul(lsc[i],self.__Fkg[i])
            self.__Fkg[i] = np.mod(self.__Fkg[i],q)
        return self.__Fkg
