from sympy import isprime
import math

def generate_factor_base(N):
    """Generate the factor base up to the bound B using SymPy."""
    B = math.exp(0.5 * math.sqrt(math.log(N) * math.log(math.log(N))))
    B = int(B)
    
    factor_base = []
    for p in range(2, B + 1):
        if isprime(p) and pow(N, (p - 1) // 2, p) == 1:
            factor_base.append(p)
    
    return factor_base

def lerEntradaArquivo(caminho_arquivo):
    with open(caminho_arquivo, 'r') as file:
        numero = int(file.readline().strip())
    return numero

caminho_arquivo = 'entrada.txt'
N = lerEntradaArquivo(caminho_arquivo)
factor_base = generate_factor_base(N)
print(factor_base)

l = len(factor_base)
print(l)

matriz = [[] for x in range(l+2)]
print(len(matriz))
xZero = math.floor(math.sqrt(N))
x = [] 
counterX = 0

while True:
    xAux = xZero
    if counterX % 2 == 0:
        X = xAux + 1
    else:
        X = xAux - 1

    xAux = xAux + 1
    funcX = pow((X + xZero), 2) - N

    # TEM QUE FATORAR funcX