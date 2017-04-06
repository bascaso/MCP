# Begona Ascaso

import numpy as N

#Multicompress (from BPZ)

def multicompress(condition,variables):
    lista=list(variables)
    n=len(lista)
    for i in range(n):	lista[i]=N.compress(condition,lista[i])
    return tuple(lista)

