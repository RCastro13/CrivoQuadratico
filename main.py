from sympy import mod_inverse, isprime
import math
from math import gcd
import numpy as np

def gauss_jordan_elimination(matrix):
    """Perform Gauss-Jordan elimination on the matrix mod 2."""
    m = np.array(matrix, dtype=int)
    rows, cols = m.shape
    
    for i in range(rows):
        # Make sure the pivot is 1
        if m[i, i] == 0:
            for j in range(i + 1, rows):
                if m[j, i] == 1:
                    m[[i, j]] = m[[j, i]]
                    break
        # Eliminate column entries above and below the pivot
        for j in range(rows):
            if j != i and m[j, i] == 1:
                m[j] = (m[j] + m[i]) % 2
    return m

def find_solution(matrix):
    """Find a non-trivial solution to the homogeneous system mod 2."""
    m = np.array(matrix, dtype=int)
    rows, cols = m.shape
    pivot_columns = []
    
    # Find the pivot columns
    for i in range(rows):
        for j in range(cols):
            if m[i, j] == 1:
                pivot_columns.append(j)
                break
    
    free_vars = [j for j in range(cols) if j not in pivot_columns]
    
    if not free_vars:
        return None  # No free variables, no non-trivial solution
    
    # Construct solution vector
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
    
    return factor_base, B

def lerEntradaArquivo(caminho_arquivo):
    with open(caminho_arquivo, 'r') as file:
        numero = int(file.readline().strip())
    return numero

def to_mod2(matrix):
    """Convert all elements of the matrix to mod 2."""
    return [[element % 2 for element in row] for row in matrix]

def solve_mod2(matrix):
    """
    Given a (l+2) x (l+1) matrix in mod 2, returns the vector beta such that matrix * beta = 0.
    """
    matrix = np.array(matrix, dtype=int)
    rows, cols = matrix.shape
    
    # Perform Gaussian elimination in mod 2
    for i in range(min(rows, cols)):
        # Find the pivot row
        if matrix[i, i] == 0:
            for j in range(i+1, rows):
                if matrix[j, i] == 1:
                    matrix[[i, j]] = matrix[[j, i]]
                    break
        
        # Eliminate below
        for j in range(i+1, rows):
            if matrix[j, i] == 1:
                matrix[j] = (matrix[j] + matrix[i]) % 2
    
    # Back-substitution to find solution
    beta = np.zeros(cols, dtype=int)
    
    for i in range(min(rows, cols)-1, -1, -1):
        if matrix[i, i] == 1:
            beta[i] = matrix[i, -1]
            for j in range(i+1, cols-1):
                beta[i] = (beta[i] + matrix[i, j] * beta[j]) % 2
    
    return beta.tolist()

def construct_xy(beta, x_vals, N):
    """Construct the solutions X and Y from the beta vector and x values."""
    x_selected = [x_vals[i] for i in range(len(beta)) if beta[i] == 1]
    X = 1
    Y = 1
    
    for x in x_selected:
        X *= x
        Y *= pow(x, 2) - N
    
    Y = int(math.sqrt(Y))
    
    return X, Y

caminho_arquivo = 'entrada.txt'
N = lerEntradaArquivo(caminho_arquivo)
primeList, B = generate_factor_base(N)
print(primeList)

l = len(primeList)
print("LIMITE SUPERIOR PARA OS PRIMOS DO CRIVO: ", B)
print("EXISTEM ", l, " PRIMOS PARA SEREM TESTADOS")
print("TAMANHO DOS VETORES: ", l+1)

matriz = [[] for x in range(4)]
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

functX = pow((xZero), 2) - N
if functX < 0:
    signal = 1
else:
    signal = 0
print("functX: ", functX)
functX = abs(functX)
#factor = functX % N
#print("fator: ", factor)
factors = factorize(functX, primeList)
if factors:
    vetor = [signal] + factors
    #vetor = factors
    #vetor = apply_mod2(vetor)
    matriz[vectorLine] = vetor
    #matriz[:, vectorLine] = vetor
    x_vals.append(xZero)
    print(f"Fatores de {functX}: {factors}")
    print(f"Vetor adicionado: {vetor}")
    vectorLine = vectorLine + 1

while vectorLine < 4:
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
    functX = abs(functX)
    #print("functX: ", functX)
    #factor = (functX) % N
    #print("fator: ", factor)
    factors = factorize(functX, primeList)
    if factors:
        print("X da vez: ", t + xZero)
        vetor = [signal] + factors
        #vetor = factors
        #vetor = apply_mod2(vetor)
        x_vals.append(t + xZero)
        matriz[vectorLine] = vetor
        #matriz[:, vectorLine] = vetor
        print(f"Fatores de {functX}: {factors}")
        print(f"Vetor adicionado: {vetor}")
        vectorLine = vectorLine + 1

    if vectorLine == l+2:
        break

    #print("VALOR DA VEZ ", tNeg,  " + ", xZero)
    functX = pow((tNeg + xZero), 2) - N
    if functX < 0:
        signal = 1
    else:
        signal = 0
    functX = abs(functX)
    #print("functX: ", functX)
    #factor = (functX) % N
    #print("fator: ", factor)
    factors = factorize(functX, primeList)
    if factors:
        print("X da vez: ", t + xZero)
        vetor = [signal] + factors
        #vetor = apply_mod2(vetor)
        x_vals.append(t + xZero)
        #matriz[:, vectorLine] = vetor
        matriz[vectorLine] = vetor
        print(f"Fatores de {functX}: {factors}")
        print(f"Vetor adicionado: {vetor}")
        vectorLine = vectorLine + 1

    dist = dist + 1
    counterX = counterX + 1

for linha in matriz:
    print(linha)

matrizMod2 = to_mod2(matriz)
for linha in matrizMod2:
    print(linha)

beta = solve_mod2(matrizMod2)
print("Beta vector:", beta)

X, Y = construct_xy(beta, x_vals, N)
print("X:", X)
print("Y:", Y)
print("GCD(N, X-Y):", gcd(N, X - Y))
print("GCD(N, X+Y):", gcd(N, X + Y))


# Resolver o sistema de congruências
#mod2_matriz_np = np.array(matrizMod2, dtype=int)
# rref_matrix = gauss_jordan_elimination(matrizMod2)
# print("Matriz na forma escalonada reduzida:")
# print(rref_matrix)
# beta = find_solution(rref_matrix)
# X=1
# if beta:
#     print("Solução não trivial encontrada:")
#     print(beta)
#     i = 0
#     # for value in x_vals:
#     #     X = X*pow(value, beta[i])
#     #     i += 1
#     # print("X: ", X)

# else:
#     print("Nenhuma solução não trivial encontrada.")