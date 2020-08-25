# -*- coding: utf-8 -*-
"""
本类用于实现对输入字符串的提取处理，使其转化为所需的系数字典

"""
import re
import gmpy2
import numpy as np

q = 730750818665451621361119245571504901405976559617

class ioHandle():
    
    #运用正则式方法实现对字符串的截断和提取
    def regular_prof(self,s,sym):
        sym = "("+sym+")"
        res = re.split(sym,s)
        return [x for x in res if x]
       
    def regular_proy(self,s,sym):
        sym = "["+sym+"]"
        res = re.split(sym,s)
        return [x for x in res if x]

    #循环调用上述两个函数，实现对多项式两个系数的提取，并对多项式乘式的次方的判断，确定是否有能力计算
    def handle0(self,s,a,k):
        #第一次切分
        sym1 = '\)\+\('
        res1 = self.regular_prof(s,sym1)
        sdic = dict()
        num = 0
        for i in range(len(res1)):
            if res1[i] != ')+(':
                #第二次切分
                si = res1[i]
                sym2 = ')(+s'
                res2 = self.regular_proy(si,sym2)
                lstemp = list()
                l1 = list()
                l2 = list()
                for j in range(len(res2)):
                    if j % 2 == 0:
                        l1.append(res2[j])
                    else:
                        l2.append(res2[j])
                lstemp = l1 + l2
                len1 = len(lstemp)
                len2 = int(len1/2)
                for j in range(len2,len1):
                    temp1 = lstemp[j]
                    # “#+i”表示密文a中的下标i对应的数据
                    if temp1[0] == '#':
                        sym3 = '#'
                        res3 = self.regular_proy(temp1,sym3)
                        temp2 = a[int(res3[0])]
                        lstemp[j] = gmpy2.mpz(temp2)
                for j in range(len1):
                    if type(lstemp[j]) == type('123'):
                        lstemp[j] = gmpy2.mpz(lstemp[j],0)
                sdic[num]=lstemp
                num = num + 1
        nv = 0
        for n in sdic:
            l = len(sdic[n])//2
            if l <= k:#这个循环就是对次方的判断，这里暂时设定为2
                nv = nv+1
        if nv == len(sdic):
            print("all of them is valid!")
            return sdic
        else:
            print("not valid!")

    # 若传输数据是多项式系数拼接形式，则调用该方法
    def handle1(self,s,k): 
        slist = list()
        sym = '$'
        res = self.regular_proy(s,sym)
        for i in range(len(res)):
            if res[i] != '$':
                temp = gmpy2.mpz(res[i],0)
                slist.append(temp)
        l = len(slist) - 1
        if l <= k:
            return slist
        else:
            print("list is not valid!")

    #对接收到的Fk(x)进行操作，使其趋于一般化
    def F_x(self,k,lsc,lc):#lsc是s的系数列表，lc是常数项列表
        Fk = list()
        Fc = gmpy2.mpz(1)
        self.Fkg = list()
        qm = gmpy2.mpz(730750818665451621361119245571504901405976559617)
        for i in range(k):
            Fc = Fc * lsc[i]
            temp1 = gmpy2.divm(gmpy2.mpz(1),lsc[i],qm)
            lc[i] = gmpy2.mul(lc[i],temp1)
            #lsc[i] = gmpy2.mpz(1)
        for i in range(k):
            p = np.poly1d([gmpy2.mpz(1),lc[i]])
            #print(p)
            Fk.append(p)
        temp = Fk[0]
        for j in range(1,k):
           temp = temp * Fk[j]
        self.Fkg = list(temp.coef)
        for n in range(k):
            self.Fkg[n] = self.Fkg[n] * Fc
            self.Fkg[i] = gmpy2.f_mod(self.Fkg[i],q)
        #return self.Fkg
        #print("the coefficient of Fk:",self.Fkg)
        
    #该类的接口函数
    def choose(self,s,a,k):      
        if s[0] == '0':
            s1 = s[1:]
            #print("s1 = \n",s1)
            llist = list()
            tdic = dict()
            sdic = self.handle0(s1,a,k)
            for key in sdic:
                lsc = list()
                lc = list()
                lf = sdic[key]
                len_lf = len(lf)
                k = int(len_lf/2)
                #这两步操作用于分别读取参数列表
                for i in range(k):
                    lsc.append(lf[i])
                for j in range(k):
                    lc.append(lf[j+k])
                self.F_x(k,lsc,lc)
                tdic[key] = self.Fkg
                llist.append(k+1)
            lm = max(llist)
            ltemp = [0]*lm
            for key in tdic:
                temp = tdic[key]
                l = len(temp)
                for i in range(l):
                    ltemp[i] = gmpy2.add(ltemp[i],temp[l-1-i])
            sl = list()
            ll = len(ltemp)
            for i in range(ll):
                sl.append(ltemp[ll-1-i])
            return sl
        
        if s[0] == '1':
            s2 = s[1:]
            slist = self.handle1(s2,k)
            if len(slist) != 0:
                return slist
            