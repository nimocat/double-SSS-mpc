# -*- coding: utf-8 -*-
"""
该类主要用以实现setup步骤所生成的所有公开参数的合法性验证，包括f(x)、hk(x)的分片和所有验证密钥。

"""

import gmpy2
from pypbc import Element,G1
from hashlib import sha256
import re
import sys
sys.path.append(r'/Users/aitianyi/MPC')
import global_variable

q = 730750818665451621361119245571504901405976559617

class verify_cc():
    
    def verify(self,k,t,commitment_s,commitment_fx_joint,commitment_fr_k,commitment_ck_joint,commitment_hk_joint,commitment_split,commitment_hk_split,ID):
        ver_result1 = self.__verify_com(k,t,commitment_s,commitment_fx_joint,commitment_fr_k,commitment_ck_joint,commitment_hk_joint)
        ver_result2 = self.__verify_core_share(k,t,commitment_split,commitment_hk_split,commitment_fx_joint,commitment_hk_joint,ID)
        if ver_result1 == True:
            if ver_result2 == True:
                return True
            
    
    def __verify_com(self,k,t,commitment_s,commitment_fx_joint,commitment_fr_k,commitment_ck_joint,commitment_hk_joint):
        #获取全局变量
        pairing = global_variable.get_pairing()
        g1 = global_variable.get_g1()
        
        #用来记录验证成功的次数
        success_num = 0
        
        #验证秘密 secret 的承诺
        count = 0
        temp1 = Element(pairing, G1)
        temp2 = Element(pairing, G1)
        for i in range(k):
            if i >= 1:
                temp1 = pairing.apply(commitment_s[1], commitment_s[i])
                temp2 = pairing.apply(commitment_s[i+1], g1)
                if (temp1 == temp2):
                    count = count + 1
                    print("the verify for commitment of secret is:^-^ ^-^ ^-^ ^-^ ^-^ ^-^success!^-^ ^-^ ^-^ ^-^ ^-^ ^-^")
                else:
                    print("the verify for commitment of secret is:>~< >~< >~< >~< >~< >~<fail!>~< >~< >~< >~< >~< >~<")
        if count == k-1:
            success_num = success_num + 1
        
        #生成 r'
        r = 0
        rt = sha256()
        rt.update(commitment_fx_joint.encode("utf8"))
        dr = rt.hexdigest()
        r = int(dr,16)
        r = gmpy2.mpz(r)
        r = gmpy2.f_mod(r,q)
        print("r' = \n",r)
                   
        #处理以字符串形式拼接起来的承诺 CM_fx,并生成CM_fr'
        temp = commitment_fx_joint
        str_fx = re.findall(r'.{130}', temp)
        #print("str_fx = \n",str_fx)
        gfr_list = list()
        for i in range(len(str_fx)):
            temp = Element(pairing,G1,value = str_fx[i])
            ri = gmpy2.powmod(r,i,q)
            temp = temp ** int(ri)
            gfr_list.append(temp)
        #print("gfr_list = \n",gfr_list)
        gfr = Element(pairing,G1)
        for i in range(len(gfr_list)):
            gfr = gfr * gfr_list[i]

        #验证秘密 fk(r) 的承诺
        count = 0
        temp3 = Element(pairing, G1)
        temp4 = Element(pairing, G1)
        for j in range(k):
            if j >= 1:
                temp3 = pairing.apply(gfr, commitment_fr_k[j])
                temp4 = pairing.apply(commitment_fr_k[j+1], g1)
                if (temp3 == temp4):
                    count = count + 1
                    print("the verify for commitment of fk(r) is:^-^ ^-^ ^-^ ^-^ ^-^ ^-^success!^-^ ^-^ ^-^ ^-^ ^-^ ^-^")
                else:
                    print("the verify for commitment of fk(r) is:>~< >~< >~< >~< >~< >~<fail!>~< >~< >~< >~< >~< >~<")
        if count == k-1:
            success_num = success_num + 1
        
        #处理以字符串形式拼接起来的承诺 CM_ck
        com_ck = dict() 
        for key in commitment_ck_joint:
            com_cki = list()
            temp = commitment_ck_joint[key]
            str_ck = re.findall(r'.{130}', temp)
            #print("str_ck = \n",str_ck)
            for i in range(len(str_ck)):
                temp1 = Element(pairing,G1,value = str_ck[i])
                com_cki.append(temp1)
            com_ck[key] = com_cki
        #print("com_ck = \n",com_ck)
        
        #生成 ck(r') 的承诺
        commitment_ck_r = dict()
        temp1 = list()
        for key in com_ck:
            comm_ck_r = list()
            temp1 = com_ck[key]
            temp2 = Element(pairing,G1)
            for i in range(t):
                temp2 = temp1[i] ** int(r ** i)
                comm_ck_r.append(temp2)
            temp3 = Element(pairing,G1)
            temp3 = comm_ck_r[0]
            for j in range(t-1):
                temp3 = temp3 * comm_ck_r[j+1]
            commitment_ck_r[key] = temp3
        #print("commitment_ck_r = \n",commitment_ck_r)    
    
        #生成 hk(r') 的承诺
        com_hk = dict()
        for key in commitment_hk_joint:
            com_hki = list()
            temp = commitment_hk_joint[key]
            str_hk = re.findall(r'.{130}', temp)
            for i in range(len(str_hk)):
                temp1 = Element(pairing,G1,value = str_hk[i])
                com_hki.append(temp1)
            com_hk[key] = com_hki
        #print("com_hk = \n",com_hk)
        
        temp1 = Element(pairing, G1)
        commitment_hk_r = dict()
        for key in com_hk:
            temp3 = list()
            temp2 = com_hk[key]
            for j in range(key*t-key+1):
                temp1 = temp2[j]**int(r**j)
                temp3.append(temp1)
            temp4 = Element(pairing, G1)
            temp4 = temp3[0]
            for k in range(key*t-key):
                temp4 = temp4 * temp3[k+1]
            commitment_hk_r[key] = temp4
        #print("commitment_hk_r = \n",commitment_hk_r)
    
        commitment_hk_r1 = dict()
        temp = Element(pairing, G1)
        for i in commitment_ck_r:
            temp = commitment_fr_k[i]-commitment_ck_r[i]-commitment_s[i]
            commitment_hk_r1[i] = temp
        #print ("commitment_hk_r1 = \n",commitment_hk_r1)
        
        #验证秘密 hk(r) 的承诺
        count = 0
        for key in commitment_hk_r1:
            if commitment_hk_r1[key] == commitment_hk_r[key]:
                count = count + 1
                print("the verify for commitment of hk(r) is:^-^ ^-^ ^-^ ^-^ ^-^ ^-^success!^-^ ^-^ ^-^ ^-^ ^-^ ^-^")
            else:
                print("the verify for commitment of hk(r) is:>~< >~< >~< >~< >~< >~<fail!>~< >~< >~< >~< >~< >~<")
        if count == len(commitment_hk_r1):
            success_num = success_num + 1
            
        #汇总验证结果
        if success_num == 3:
            return True



    def __verify_core_share(self,k,t,commitment_split,commitment_hk_split,commitment_fx_joint,commitment_hk_joint,ID):      
        #获取全局变量
        pairing = global_variable.get_pairing()

        #用来记录验证成功的次数
        success_num = 0
        
        #验证分片 fi 的承诺
        count = 0
        com_fx_coef = list()
        temp = commitment_fx_joint
        str_fx = re.findall(r'.{130}', temp)
        for i in range(len(str_fx)):
            temp = Element(pairing,G1,value = str_fx[i])
            com_fx_coef.append(temp)
        #print("com_fx_coef = \n",com_fx_coef)
        
        commitment_fID = Element(pairing, G1)
        comm_fi = Element(pairing, G1)  
        temp2 = Element(pairing, G1)
        commitment_fxcoef_ID = list()
        for j in range(t): 
            temp2 = com_fx_coef[j] ** int(ID**j)
            commitment_fxcoef_ID.append(temp2)
        #print("commitment_fxcoef_ID = ",commitment_fxcoef_ID)
        comm_fi = commitment_fxcoef_ID[0]
        for m in range(t-1):
            comm_fi = comm_fi * commitment_fxcoef_ID[m+1]
        commitment_fID = comm_fi
        #print("commitment_fID = ",commitment_fID)
        if commitment_fID == commitment_split[ID]:     
            count = count + 1
            print("the verify for commitment of fi is:^-^ ^-^ ^-^ ^-^ ^-^ ^-^success!^-^ ^-^ ^-^ ^-^ ^-^ ^-^")
        else:
            print("the verify for commitment of fi is:>~< >~< >~< >~< >~< >~<fail!>~< >~< >~< >~< >~< >~<")
        if count == len(commitment_fID):
            success_num = success_num + 1
            
        #验证分片 hki 的承诺
        count = 0
        temp1 = Element(pairing, G1)
        temp4 = Element(pairing, G1)
        commitment_hk_i = dict()
        com_hk = dict()
        for key in commitment_hk_joint:
            com_hki = list()
            temp = commitment_hk_joint[key]
            str_hk = re.findall(r'.{130}', temp)
            for i in range(len(str_hk)):
                temp1 = Element(pairing,G1,value = str_hk[i])
                com_hki.append(temp1)
            com_hk[key] = com_hki
        #print("com_hk = ",com_hk)
        
        for key in com_hk:
            com_hk_key = com_hk[key]
            comm_hk_ii = list()
            for j in range(key*t-key+1):
                temp1 = com_hk_key[j]**int(ID**j)
                comm_hk_ii.append(temp1)
            temp4 = comm_hk_ii[0]
            for k in range(key*t-key):
                temp4 = temp4 * comm_hk_ii[k+1]
            commitment_hk_i[key] = temp4
        print("commitment_hk_i = ",commitment_hk_i)

        if commitment_hk_split[ID] == commitment_hk_i:
            count = count + 1
            print("the verify for commitment of hksplit is:^-^ ^-^ ^-^ ^-^ ^-^ ^-^success!^-^ ^-^ ^-^ ^-^ ^-^ ^-^")
        else:
            print("the verify for commitment of hksplit is:>~< >~< >~< >~< >~< >~<fail!>~< >~< >~< >~< >~< >~<")
        if count == len(commitment_hk_split):
            success_num = success_num + 1
        
        if success_num == 2:
            return True
            
            
            
    def fun(self,coef_poly,k,t,x):
        ff = gmpy2.mpz(0)
        for j in range(k*t-k+1)[::-1]:   
           temp = gmpy2.mul(ff,x)
           ff = gmpy2.add(temp,coef_poly[j])
           if j == 1:
               break
        temp = gmpy2.mul(ff,x) #temp = at-1*xt-1+...+a1*x
        ff = gmpy2.add(temp,coef_poly[0]) #f = temp+a0
        ff = gmpy2.f_mod(ff,q)
        return ff