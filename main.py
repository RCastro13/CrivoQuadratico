from sympy import isprime, factorint
import math

def factorize(n, factor_base):
    abs_n = abs(n)
    exponents = [0] * len(factor_base)
    
    for i, p in enumerate(factor_base):
        if abs_n == 1:
            break
        while abs_n % p == 0:
            abs_n //= p
            exponents[i] += 1
    
    # If n is not B-smooth, return None
    if abs_n != 1:
        return None
    
    return exponents

def generate_factor_base(N):
    """Generate the factor base up to the bound B using SymPy."""
    B = math.exp(0.5 * math.sqrt(math.log(N) * math.log(math.log(N))))
    B = int(B)
    
    factor_base = []
    for p in range(2, B + 1):
        if isprime(p) and pow(N, (p - 1) // 2, p) == 1:
            factor_base.append(p)

    factor_base.append(19)
    
    return factor_base, B

def lerEntradaArquivo(caminho_arquivo):
    with open(caminho_arquivo, 'r') as file:
        numero = int(file.readline().strip())
    return numero

caminho_arquivo = 'entrada.txt'
N = lerEntradaArquivo(caminho_arquivo)
primeList, B = generate_factor_base(N)
print(primeList)

l = len(primeList)
print("LIMITE SUPERIOR PARA OS PRIMOS DO CRIVO: ", B)
print("EXISTEM ", l, " PRIMOS PARA SEREM TESTADOS")
print("TAMANHO DOS VETORES: ", l+2)

matriz = [[] for x in range(l+2)]
#print(len(matriz))
xZero = math.ceil(math.sqrt(N))
print(xZero)
x = [] 
counterX = 0
vectorLine = 0
dist = 0
t = 0

functX = pow((xZero), 2) - N
if functX < 0:
    signal = 1
else:
    signal = 0
#print("functX: ", functX)
factor = functX % N
#print("fator: ", factor)
factors = factorize(factor, primeList)
if factors:
    vetor = [signal] + factors
    matriz[vectorLine] = vetor
    print(f"Fatores de {factor}: {factors}")
    print(f"Vetor adicionado: {vetor}")
    vectorLine = vectorLine + 1
# else:
#     print(f"{functX} não é B-smooth")

while True:
    if counterX % 2 == 0:
        t = dist + 1
    else:
        distNeg = -1*dist
        t = distNeg - 1

    #print("VALOR DA VEZ ", t + xZero)
    functX = pow((t + xZero), 2) - N
    if functX < 0:
        signal = 1
    else:
        signal = 0
    #print("functX: ", functX)
    factor = (functX) % N
    #print("fator: ", factor)
    factors = factorize(factor, primeList)
    if factors:
        vetor = [signal] + factors
        matriz[vectorLine] = vetor
        print(f"Fatores de {factor}: {factors}")
        print(f"Vetor adicionado: {vetor}")
        vectorLine = vectorLine + 1
    # else:
    #     print(f"{factor} não é B-smooth")

    # TEM QUE FATORAR funcX
    dist = dist + 1
    if(dist == 100):
        break

print(matriz)