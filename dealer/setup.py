# -*- coding: utf-8 -*-
"""
该类主要实现了安全多方计算所需的所有公开参数的生成

"""
import gmpy2
import random
from pypbc import Element,G1
from hashlib import sha256
import sys
sys.path.append(r'/Users/aitianyi/MPC')
import global_variable

q = 730750818665451621361119245571504901405976559617

class setup():
    
    #该函数与该类同名，为该类的唯一对外的接口函数，其返回值为生成的所有公开参数和验证密钥。
    def setup(self,t,n,k,m):
        self.__Get_userid(n)  
        self.__Gen_polynomial(t,k)
        self.__Gen_core_share(t,n)
        self.__d_encryption(m)
        self.__Gen_vk(t,n,k)
        return self.__user_id,self.__a,self.__f,self.__fk,self.__ck,self.__hk,self.__split,self.__hksplit,self.__commitment_s,self.__commitment_fx_joint,self.__commitment_hk_joint,self.__commitment_ck_joint,self.__commitment_fr_k,self.__commitment_split,self.__commitment_hk_split
        

 
    def __Get_userid(self,n):
        self.__user_id = dict()
        print("Please input your id:\n")
        for i in range(n):
            Id = int(input())
            self.__user_id[i+1] = Id
        print("user_id = ",self.__user_id)


    def __Gen_polynomial(self,t,k):
        #生成 f(x)
        
        self.__f = list()#用于存储f(x)的系数
        secret = random.randint(0,2**255)
        secret = gmpy2.mpz(secret)
        #print("secret = \n",secret)
        self.__f.append(secret)
        for i in range(t-1):#f(x)多项式的其他系数
            fki = random.randint(0,2**255)
            fki = gmpy2.mpz(fki)
            self.__f.append(fki)
        for i in range(t):
            self.__f[i] = gmpy2.f_mod(self.__f[i],q)
        #print("f = \n",self.__f) 
               
        #生成 ck(x), k 取值为2，，，k，将某个 ck 的系数存在列表中，将生成的所有 ck(x) 存在字典中
        
        self.__ck = dict()
        for i in range(k+1):
            if i > 1:    
                cki = list()
                cki.append(gmpy2.mpz(0))
                for j in range(t-1):
                    ckii = random.randint(0,2**255)
                    ckii = gmpy2.mpz(ckii)
                    cki.append(ckii)
                for l in range(t):
                    cki[l] = gmpy2.f_mod(cki[l],q)
                self.__ck[i] = cki
        #print("ck = \n",self.__ck)
        
                      
        #生成 fk(x), k 取值为 2，，，k，将某个 fk 的系数存在列表中，将生成的所有 fk(x) 存在字典中
        self.__fk = dict()
        for m in range(k+1):#生成f^2(x)至f^k(x),m表示生成f(x)的m次方
            if m > 1:
                fki = [0] * (m*t-m+1)
                for i in range(t):
                    fki[i] = self.__f[i]
                for j in range(m-1):
                    fktemp = [0] * (m*t-m+1)
                    for p in range((m-1)*t-m+2):
                        for l in range(t):
                            temp = gmpy2.mul(fki[p] , self.__f[l]) 
                            fktemp[p+l] =  gmpy2.add(fktemp[p+l] , temp)
                    for w in range(m*t-m+1):
                        fki[w] = gmpy2.f_mod(fktemp[w],q)
                self.__fk[m] = fki
        #print("fk = \n",self.__fk)
    
        #计算 h^k(x)，并将系数存在字典中
        self.__hk = dict()
        fk1 = list()
        ck1 = list()
        for i in self.__fk:#取出字典中的列表进行h^k(x)的计算
            hki =  list()
            fk1 = self.__fk[i]
            ck1 = self.__ck[i]
            for j in range(t):
                temp = gmpy2.sub(fk1[j] , ck1[j])
                hki.append(temp)
            for l in range((i-1)*t-i+1):
                hki.append(fk1[t+l])
            hki[0] = gmpy2.sub(hki[0] , self.__f[0]**i)
            for m in range(i*t-i+1):
                hki[m] = gmpy2.f_mod(hki[m],q)
            self.__hk[i] = hki
        #print ("hk = \n",self.__hk)
                
        #return self.f,self.fk,self.ck,self.hk



    def __Gen_core_share(self,t,n):
        #生成秘密分片 f(i)
        self.__split = dict()
        for i in self.__user_id:
            temp = self.__fun(self.__f,1,t,self.__user_id[i])
            self.__split[self.__user_id[i]] = temp
        #print ("split = \n",self.__split)
    
    
        #生成 hk(x) 的分片
        
        self.__hksplit = dict()
        for i in range(n):
            hksplit_i = dict()
            for j in self.__hk:
                temp = self.__fun(self.__hk[j],j,t,self.__user_id[i+1])
                hksplit_i[j] = temp
            self.__hksplit[self.__user_id[i+1]] = hksplit_i
        #print("hksplit = \n",self.__hksplit)
        

     
    def __Gen_vk(self,t,n,k):
        #获取全局变量
        pairing = global_variable.get_pairing()
        g1 = global_variable.get_g1()

        #生成秘密 secret 的承诺
        secret_k = dict()
        secret = self.__f[0]
        temp = secret
        secret_k[1] = secret
        for i in range(k):
            if i >= 1:
                temp = gmpy2.mul(temp,secret)
                temp = gmpy2.f_mod(temp , q)
                secret_k[i+1] = temp
        commitment_s_i = Element(pairing, G1)
        self.__commitment_s = dict()
        for i in secret_k:
            commitment_s_i = g1 ** int(secret_k[i])
            self.__commitment_s[i] = commitment_s_i
        #print ("commitment_s = \n",self.__commitment_s)
        
        #生成 f(x) 的承诺
        commitment_fx_coef = self.__comm(self.__f,1,t,pairing,g1)
        print("the commitment of fx = \n",commitment_fx_coef)
        self.__commitment_fx_joint = str(commitment_fx_coef[0])
        for i in range(t-1):
            self.__commitment_fx_joint = self.__commitment_fx_joint+str(commitment_fx_coef[i+1])
        #print("the joint of fx = \n",self.__commitment_fx_joint)
                
        #生成 r
        self.__r = 0
        r = sha256()
        r.update(self.__commitment_fx_joint.encode("utf8"))
        dr = r.hexdigest()
        self.__r = int(dr,16)
        self.__r = gmpy2.mpz(self.__r)
        self.__r = gmpy2.f_mod(self.__r,q)
        #print("r = \n",self.__r)
        
        #生成 f(r) 和 f(r)^k 的承诺
        self.__commitment_fr_k = dict()
        fr = self.__fun(self.__f,1,t,self.__r)
        comm_fr = Element(pairing, G1)
        comm_fr = g1 ** int(fr)
        self.__commitment_fr_k[1] = comm_fr
        fr_i = fr
        for i in range(k-1):
            fr_i = fr_i * fr
            comm_fr_i = Element(pairing, G1)
            comm_fr_i = g1 ** int(fr_i)
            self.__commitment_fr_k[i+2] = comm_fr_i
        print("commitment_fr_k = \n",self.__commitment_fr_k)
                
        #生成 ck(x) 的承诺
        commitment_ck_coef = dict()
        self.__commitment_ck_joint = dict()
        for i in self.__ck:
            commitment_ck_i = list()  
            ck_i = self.__ck[i]
            commitment_ck_i = self.__comm(ck_i,1,t,pairing,g1)
            commitment_ck_coef[i] = commitment_ck_i
        #print("commitment_ck_coef = \n", commitment_ck_coef)
        for i in commitment_ck_coef:
            comm_ci_coef = commitment_ck_coef[i]
            str_comm_ci_coef = str(comm_ci_coef[0])
            str_comm_ci_coef = str_comm_ci_coef[:130]
            for j in range(t-1):
                str_comm_ci_coef = str_comm_ci_coef+str(comm_ci_coef[j+1])
            self.__commitment_ck_joint[i] = str_comm_ci_coef
        #print("commitment_ck_joint = \n",self.__commitment_ck_joint)

        
        #生成 hk(x) 的承诺
        commitment_hk_coef = dict()#用于存储基点与hk各个系数的点乘结果
        self.__commitment_hk_joint = dict()#用于存储拼接后的hk的承诺
        for i in self.__hk:
            commitment_hk_i = list()  
            hk_i = self.__hk[i]
            commitment_hk_i = self.__comm(hk_i,i,t,pairing,g1)
            commitment_hk_coef[i] = commitment_hk_i
        #print("commitment_hk_coef = \n",commitment_hk_coef)
        for i in commitment_hk_coef:
            comm_hi_coef = commitment_hk_coef[i]
            str_comm_hi_coef = str(comm_hi_coef[0])
            str_comm_hi_coef = str_comm_hi_coef[:130]
            for j in range(i*t-i):
                str_comm_hi_coef = str_comm_hi_coef + str(comm_hi_coef[j+1])
            self.__commitment_hk_joint[i] = str_comm_hi_coef
        #print("commitment_hk_joint = \n",self.__commitment_hk_joint)
      
        #生成 f(x) 分片 f(i) 的承诺,用于安全传输
        self.__commitment_split = dict()
        for i in self.__split:
            comm_split = Element(pairing, G1)
            comm_split = g1 ** int(self.__split[i])
            self.__commitment_split[i] = comm_split
        #print("commitment_split = \n",self.__commitment_split)

        #生成 hk(x) 分片 hk(i) 的承诺,用于安全传输
        self.__commitment_hk_split = dict()
        for i in self.__hksplit:
            hksplit_key = self.__hksplit[i]
            comm_hksplit_key = dict()
            for key in hksplit_key:
                comm_hksplit = Element(pairing, G1)
                comm_hksplit = g1 ** int(hksplit_key[key])
                comm_hksplit_key[key] = comm_hksplit
            self.__commitment_hk_split[i] = comm_hksplit_key
        #print("commitment_hk_split = ",self.__commitment_hk_split)
           
            
         


    def __d_encryption(self,m):
        d = list()
        self.__a = list()      
        for i in range(m):
            di = random.randint(0,2**255)
            di = gmpy2.mpz(di)
            d.append(di)
            d[i] = gmpy2.f_mod(d[i],q)          
            self.__a.append(gmpy2.sub(d[i], self.__f[0]))
            self.__a[i] = gmpy2.f_mod(self.__a[i],q)    
        #print("a = \n",self.__a)
        
        
    def __comm(self,f,k,t,pair,g):#生成承诺,g1与多项式系数进行点乘
        fg = list()
        for i in range(k*t-k+1):
            bg = Element(pair, G1)
            bg = g ** int(f[i]) 
            fg.append(bg)
        return fg
        
    
    def __fun(self,coef_poly,k,t,x):
        ff = gmpy2.mpz(0)
        for j in range(k*t-k+1)[::-1]: 
           x = gmpy2.mpz(x)
           temp = gmpy2.mul(ff,x)
           ff = gmpy2.add(temp,coef_poly[j])
           if j == 1:
               break
        temp = gmpy2.mul(ff,x) #temp = at-1*xt-1+...+a1*x
        ff = gmpy2.add(temp,coef_poly[0]) #f = temp+a0
        ff = gmpy2.f_mod(ff,q)
        return ff