from sympy import mod_inverse, isprime, Matrix
import math
from math import gcd
import numpy as np
from functools import reduce

def gauss_jordan_elimination(matrix):
    #faz a eliminacao de gauss-jordan em uma matriz mod 2
    m = np.array(matrix, dtype=int)
    rows, cols = m.shape
    
    for i in range(rows):
        if m[i, i] == 0:
            for j in range(i + 1, rows):
                if m[j, i] == 1:
                    m[[i, j]] = m[[j, i]]
                    break
        for j in range(rows):
            if j != i and m[j, i] == 1:
                m[j] = (m[j] + m[i]) % 2
    return m

def find_solution(matrix):
    #encontra solucao nao trivial para o sistema homogeneo mod 2
    m = np.array(matrix, dtype=int)
    rows, cols = m.shape
    pivot_columns = []
    
    #encontra o pivô
    for i in range(rows):
        for j in range(cols):
            if m[i, j] == 1:
                pivot_columns.append(j)
                break
    
    free_vars = [j for j in range(cols) if j not in pivot_columns]
    
    if not free_vars:
        return None  #nenhuma variavel livre implica em nenhuma solucao nao trivial
    
    #contrói o vetor de solucao
    solution = [0] * cols
    for free_var in free_vars:
        solution[free_var] = 1
        for i in range(rows):
            if m[i, free_var] == 1:
                for pivot_col in pivot_columns:
                    solution[pivot_col] = (solution[pivot_col] + m[i, pivot_col]) % 2
                break
    
    return solution

def factorize(n, factor_base):
    abs_n = abs(n)
    exponents = [0] * len(factor_base)
    
    for i, p in enumerate(factor_base):
        if abs_n == 1:
            break
        while abs_n % p == 0:
            abs_n //= p
            exponents[i] += 1
    
    #se n nao é B-smooth, retorna None
    if abs_n != 1:
        return None
    
    return exponents

def generate_factor_base(N):
    #gera a base de fatores até o limite B usando o SymPy
    B = math.exp(0.5 * math.sqrt(math.log(N) * math.log(math.log(N))))
    B = int(B)
    
    factor_base = []
    for p in range(2, B + 1):
        if isprime(p) and pow(N, (p - 1) // 2, p) == 1:
            factor_base.append(p)
    
    return factor_base, B

def lerEntradaArquivo(caminho_arquivo):
    with open(caminho_arquivo, 'r') as file:
        numero = int(file.readline().strip())
    return numero

def to_mod2(matrix):
    #converte todos os elementos da matriz para mod 2
    return [[element % 2 for element in row] for row in matrix]

def construct_xy(beta, x_vals, N):
    #controi as solucoes X e Y para o vetor beta e X_values
    x_selected = [x_vals[i] for i in range(len(beta)) if beta[i] == 1]
    X = 1
    Y = 1
    
    for x in x_selected:
        X *= x
        Y *= pow(x, 2) - N
    
    Y = int(math.sqrt(Y))
    
    return X, Y

def rref_mod2(matrix):
    rows, cols = matrix.shape
    A = matrix.copy()

    #transformando a matriz em uma forma escalonada reduzida sobre F2
    lead = 0
    for r in range(rows):
        if lead >= cols:
            return A
        i = r
        while A[i, lead] == 0:
            i += 1
            if i == rows:
                i = r
                lead += 1
                if lead == cols:
                    return A
        A[i], A[r] = A[r].copy(), A[i].copy()
        lv = A[r, lead]
        A[r] = (A[r] / lv) % 2
        for i in range(rows):
            if i != r:
                lv = A[i, lead]
                A[i] = (A[i] - lv * A[r]) % 2
        lead += 1
    return A

def is_integer_sqrt(x):
    if x < 0:
        return False
    root = math.isqrt(x)
    return root * root == x

caminho_arquivo = 'entrada.txt'
N = lerEntradaArquivo(caminho_arquivo)
primeList, B = generate_factor_base(N)
print(primeList)

l = len(primeList)
print("LIMITE SUPERIOR PARA OS PRIMOS DO CRIVO: ", B)
print("EXISTEM ", l, " PRIMOS PARA SEREM TESTADOS")
print("TAMANHO DOS VETORES: ", l+1)

matriz = np.zeros((l, l), dtype=int)
#matriz = np.zeros((l + 2, l + 1))
print(len(matriz))
xZero = math.floor(math.sqrt(N))
print(xZero)
x_vals = []
functX_vals = [] 
counterX = 0
vectorLine = 0
dist = 0
t = 0
factor1 = 0
factor2 = 0

functX = pow((xZero), 2) - N
if functX < 0:
    signal = 1
else:
    signal = 0
print("functX: ", functX)
originFunctX = functX
functX = abs(functX)
#factor = functX % N

if is_integer_sqrt(functX) and originFunctX >= 0:
    square = math.sqrt(functX)
    factor1 = (xZero) - square
    factor2 = (xZero) + square
else:
    factors = factorize(functX, primeList)
    if factors:
        #vetor = [signal] + factors
        vetor = factors
        matriz[vectorLine] = vetor
        x_vals.append(xZero)
        print(f"Fatores de {functX}: {factors}")
        print(f"Vetor adicionado: {vetor}")
        vectorLine = vectorLine + 1

tNeg = 0
while vectorLine < l:
    #print("VECTORLINE: ", vectorLine)
    t = dist + 1
    distNeg = -1*t
    tNeg = distNeg

    #print("VALOR DA VEZ ", t,  " + ", xZero)
    functX = pow((t + xZero), 2) - N
    if functX < 0:
        signal = 1
    else:
        signal = 0
    originFunctX = functX
    functX = abs(functX)
    #print("functX: ", functX)
    #factor = (functX) % N
    #print("fator: ", factor)
    if is_integer_sqrt(functX) and originFunctX >= 0:
        square = math.sqrt(functX)
        factor1 = (t + xZero) - square
        factor2 = (t + xZero) + square
        break
    else:
        factors = factorize(functX, primeList)
        if factors:
            print("X da vez: ", t + xZero)
            #vetor = [signal] + factors
            vetor = factors
            x_vals.append(t + xZero)
            matriz[vectorLine] = vetor
            print(f"Fatores de {functX}: {factors}")
            print(f"Vetor adicionado: {vetor}")
            vectorLine = vectorLine + 1
            print("VECTORLINE: ",vectorLine)

    # if vectorLine == l+1:
    #     break

    # #print("VALOR DA VEZ ", tNeg,  " + ", xZero)
    # functX = pow((tNeg + xZero), 2) - N
    # if functX < 0:
    #     signal = 1
    # else:
    #     signal = 0
    # originFunctX = functX
    # functX = abs(functX)
    # #print("functX: ", functX)
    # #factor = (functX) % N
    # #print("fator: ", factor)
    # if is_integer_sqrt(functX) and originFunctX >= 0:
    #     print("FUNCTX: ", functX)
    #     print("ORIGIN: ", originFunctX)
    #     square = math.sqrt(functX)
    #     print("SQUARE: ", square)
    #     factor1 = (tNeg + xZero) - square
    #     factor2 = (tNeg + xZero) + square
    #     break
    # else:
    #     factors = factorize(functX, primeList)
    #     if factors:
    #         print("X da vez: ", t + xZero)
    #         vetor = [signal] + factors
    #         #vetor = factors
    #         #vetor = apply_mod2(vetor)
    #         x_vals.append(t + xZero)
    #         #matriz[:, vectorLine] = vetor
    #         matriz[vectorLine] = vetor
    #         print(f"Fatores de {functX}: {factors}")
    #         print(f"Vetor adicionado: {vetor}")
    #         vectorLine = vectorLine + 1

    dist = dist + 1
    counterX = counterX + 1
# factor = quadratic_sieve(matriz, N)
# print(f"Found factor: {factor}")

if factor1 and factor2:
    print("FATORES TRIVIAIS: ")
    print(factor1)
    print(factor2)
else:
    matrizMod2 = to_mod2(matriz)
    matrizMod2 = np.array(matrizMod2)
    matrizMod2 = matrizMod2.T
    for linha in matrizMod2:
        print(linha)

    rref_matrix = gauss_jordan_elimination(matrizMod2)
    print("Matriz na forma escalonada reduzida:")
    print(rref_matrix)
    beta = find_solution(rref_matrix)
    if beta:
        print("Solução não trivial encontrada:")
        print(beta)
        X, Y = construct_xy(beta, x_vals, N)
        print("X:", X)
        print("Y:", Y)
        print("MDC(N, X-Y):", gcd(N, X - Y))
        print("MDC(N, X+Y):", gcd(N, X + Y))
    else:
        print("Nenhuma solução não trivial encontrada.")