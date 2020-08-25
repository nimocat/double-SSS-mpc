# -*- coding: utf-8 -*-


from __future__ import division
from pypbc import Element, G1
import sys
sys.path.append(r'/Users/aitianyi/MPC')
import global_variable
import re


class verify_response():

    def verify_response_comm(self, slist, t, ID, commitment_fx_joint, commitment_s, commitment_cx_joint, comm_response_ID):
        pairing = global_variable.get_pairing()
        g1 = global_variable.get_g1()

        self.__response_coef_comm = dict()
        commitment_cx = dict()  # 存储ck(ID)的承诺
        commitment_fx = dict()  # 存储f(ID)^k的承诺
        com_ck = dict()  # 存ck(x)系数的承诺
        for key in commitment_cx_joint:
            if key <= (len(slist)-1):
                com_cki = list()  # 中间变量
                temp = commitment_cx_joint[key]
                str_ck = re.findall(r'.{130}', temp)
                for i in range(len(str_ck)):
                    temp1 = Element(pairing, G1, value=str_ck[i])
                    com_cki.append(temp1)
                com_ck[key] = com_cki
                com_ck_key = com_ck[key]
                comm_ck_coef_i = list()
                for j in range(t):
                    temp2 = com_ck_key[j] ** (int(ID)**j)
                    comm_ck_coef_i.append(temp2)
                temp3 = comm_ck_coef_i[0]
                for j in range(t-1):
                    temp3 = temp3 * comm_ck_coef_i[j+1]
                commitment_cx[key] = temp3
        temp = commitment_fx_joint
        str_fx = re.findall(r'.{130}', temp)
        gfr_list = list()
        for i in range(len(str_fx)):
            temp1 = Element(pairing, G1, value=str_fx[i])
            gfr_list.append(temp1)
        comm_fk_coef_i = list()
        for j in range(t):
            temp2 = gfr_list[j] ** (int(ID)**j)
            comm_fk_coef_i.append(temp2)
        temp3 = comm_fk_coef_i[0]
        for j in range(t-1):
            temp3 = temp3 * comm_fk_coef_i[j+1]
        commitment_fx = temp3
        commck_mul_commsk = dict()  # 存储ck(ID)的承诺*s^k的承诺
        for key in commitment_cx:
            temp1 = commitment_cx[key]*commitment_s[key]
            commck_mul_commsk[key] = temp1
        if (len(slist)-1) > 1:
            response_commcksk = dict()
            for i in range(len(slist)-1):
                if i+2 <= len(slist)-1:
                    temp1 = commck_mul_commsk[i+2]
                    temp = temp1 ** int(slist[len(slist)-i-3])
                    response_commcksk[i+2] = temp
            temp2 = response_commcksk[2]
            for i in response_commcksk:
                if i != 2:
                    temp3 = response_commcksk[i]
                    temp2 = temp2 * temp3
            self.__response_coef_comm = temp2 * \
                (commitment_fx**int(slist[len(slist)-2])
                 )*(g1**int(slist[len(slist)-1]))
        else:
            self.__response_coef_comm = (
                commitment_fx**int(slist[len(slist)-2]))*(g1**int(slist[len(slist)-1]))
            #print("response_coef_comm = ",self.__response_coef_comm)

        if comm_response_ID == self.__response_coef_comm:
            print(
                "the verify for response is:^-^ ^-^ ^-^ ^-^ ^-^ ^-^success!^-^ ^-^ ^-^ ^-^ ^-^ ^-^")
        else:
            print(
                "the verify for response is:>~< >~< >~< >~< >~< >~<fali!>~< >~< >~< >~< >~< >~<")
