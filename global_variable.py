# -*- coding: utf-8 -*-
from __future__ import division
from pypbc import Parameters,Pairing,Element,G1

stored_params = """type a
q 8780710799663312522437781984754049815806883199414208211028653399266475630880222957078625179422662221423155858769582317459277713367317481324925129998224791
h 12016012264891146079388821366740534204802954401251311822919615131047207289359704531102844802183906537786776
r 730750818665451621361119245571504901405976559617
exp2 159
exp1 107
sign1 1
sign0 1
"""

class global_var(): 
    
    params = Parameters(param_string=stored_params)
    pairing = Pairing(params)
    g1 = Element.random(pairing, G1)

''''' 对于每个全局变量,都需要定义get_value和set_value接口'''''

''' 环境变量,只提供get函数'''
def get_pairing():
    return global_var.pairing
    
def get_g1():
    return global_var.g1

