from sympy import isprime
from math import isqrt, sqrt, log, exp
from itertools import chain

#função para cálculo de mdc(a,b)
def gcd(a, b):
    if b == 0:
        return a
    elif a >= b:
        return gcd(b,a % b)
    else:
        return gcd(b,a)

#função para verificar se 'a' é um resíduo quadrado de n
def quadResidue(a, n):
    l = 1
    q = (n-1)//2
    x = q**l
    if x == 0:
        return 1
        
    a = a % n
    z = 1
    while x != 0:
        if x % 2 == 0:
            a = (a**2) % n
            x //= 2
        else:
            x -= 1
            z = (z*a) % n

    return z

#função que retorna k tal que b^k = 1 mod p
def order(p, q):
	if gcd(p, q) != 1:
		return -1
	k = 3
	while True:
		if pow(q, k, p) == 1:
			return k
		k += 1

#função que p-1 (= x de parâmetro) tal que x*2^e onde x é par 
def convertx2e(x):
	e = 0
	while x % 2 == 0:
		x /= 2
		e += 1
	return int(x), e

#função principal que aplica o algoritmo de tonelli-shanks para resolver x^2 = N mod p
def tonelliShanks(n, p):
    assert quadResidue(n, p) == 1, "não é um quadrado (mod p)"
    q = p - 1
    s = 0
    
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        r = pow(n, (p + 1) // 4, p)
        return r,p-r
    for z in range(2, p):
        if p - 1 == quadResidue(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i

    return (r, p - r)

#função que dado um vetor solução, calcula X e Y e retorna ambos e um dos fatores de N
def solve(solutionVec, smoothNums, xlist, N):
    solutionNums = [smoothNums[i] for i in solutionVec]
    xNums = [xlist[i] for i in solutionVec]
    aSquare = 1
    for num in solutionNums:
        aSquare *= num
        
    x = 1
    for num in xNums:
        x *= num
    y = isqrt(aSquare)
    factor = gcd(x - y, N)
    return factor, x, y

#função que executa eliminação gaussiana (encontra o espaço nulo a esquerda)
def gaussElim(matrix):
    marks = [False] * len(matrix[0])  
    for i in range(len(matrix)):
        row = matrix[i]
        #procura pelo pivô
        for num in row:
            if num == 1:
                j = row.index(num)
                marks[j] = True
                #procura por outros 1s na mesma coluna
                for k in chain(range(0, i), range(i+1, len(matrix))): 
                    if matrix[k][j] == 1:
                        for i in range(len(matrix[k])):
                            matrix[k][i] = (matrix[k][i] + row[i]) % 2
                break
    
    matrix = transpose(matrix)
    solRows = []
    #encontra colunas livre (são linhas agora)
    for i in range(len(marks)): 
        if marks[i] == False:
            freeRow = [matrix[i],i]
            solRows.append(freeRow)
    
    if not solRows:
        return(" ")
    return solRows, marks, matrix

#função que recebe as possíveis soluções e um index e retorna a solução
def solveRow(solRows, matrix, marks, K=0):
    solutionVec, indices = [], []
    freeRow = solRows[K][0]
    for i in range(len(freeRow)):
        if freeRow[i] == 1: 
            indices.append(i)
    #linhas com 1 na mesma coluna são dependentes
    for r in range(len(matrix)):
        for i in indices:
            if matrix[r][i] == 1 and marks[r]:
                solutionVec.append(r)
                break
            
    solutionVec.append(solRows[K][1])       
    return(solutionVec)

#função que transpõe uma matriz
def transpose(matrix):
    TMatrix = []
    for i in range(len(matrix[0])):
        newRow = []
        for row in matrix:
            newRow.append(row[i])
        TMatrix.append(newRow)
    return(TMatrix)

#função que constrói a matriz de expoentes dos números B-smooth
def buildMatrix(smoothNums, factorBase):
    def factor(n, factorBase):
        factors = []
        if n < 0:
            factors.append(-1)
        for p in factorBase:
            if p == -1:
                pass
            else:
                while n % p == 0:
                    factors.append(p)
                    n //= p
        return factors

    M = []
    factorBase.insert(0, -1)
    for n in smoothNums:
        expVector = [0] * (len(factorBase))
        nFactors = factor(n, factorBase)
        for i in range(len(factorBase)):
            if factorBase[i] in nFactors:
                expVector[i] = (expVector[i] + nFactors.count(factorBase[i])) % 2

        #procura por quadrados
        if 1 not in expVector:
            return True, n
        else:
            pass
        
        M.append(expVector)  

    return(False, transpose(M))

#função que encontra os números B-smooth usando o método do crivo
def findSmooth(factorBase, N, I, T, xZero):
    #gera o valor da função Q(x) = x^2 - N, começando em x = xZero
    def sievePrep(N, sieveInt, xZero):
        sieveSeq = [x**2 - N for x in range(xZero, xZero+sieveInt)]
        return sieveSeq

    sieveSeq = sievePrep(N, I, xZero)
    sieveList = sieveSeq.copy()
    if factorBase[0] == 2:
        i = 0
        while sieveList[i] % 2 != 0:
            i += 1
        #encontrei o 1º termo par, agora todos os outros também serão
        for j in range(i, len(sieveList), 2):
            while sieveList[j] % 2 == 0:
                sieveList[j] //= 2

    for p in factorBase[1:]:
        #encontro x tal que x^2 = n mod p com o algoritmo de Tonelli-Shanks
        residues = tonelliShanks(N, p) 
        for r in residues:
            for i in range((r - xZero) % p, len(sieveList), p):
                while sieveList[i] % p == 0:
                    sieveList[i] //= p
    xlistOriginal = []
    smoothNums = []
    indexes = []
    
    # dada a lista de candidatos, determina os números B-smooth e os valores para X associados
    for i in range(len(sieveList)):
        if len(smoothNums) >= len(factorBase) + T:
            break
        if sieveList[i] == 1: # found B-smooth number
            smoothNums.append(sieveSeq[i])
            xlistOriginal.append(i + xZero)
            indexes.append(i)

    return(smoothNums, xlistOriginal, indexes)

#função para gerar o B heurístico e os fatores primos com base nesse
def generateFactorBase(N):
    B = exp(0.5 * sqrt(log(N) * log(log(N))))
    B = int(B)
    factorBase = []
    for p in range(2, B + 1):
        if isprime(p) and pow(N, (p - 1) // 2, p) == 1:
            factorBase.append(p)
    
    return factorBase, B

#função que recebe um valor de B e gera os fatores primos com base nele
def generateFactorBaseWithB(N, B):
    factorBase = []
    for p in range(2, B + 1):
        if isprime(p) and pow(N, (p - 1) // 2, p) == 1:
            factorBase.append(p)
    
    return factorBase, B

#função principal para cálculo do crivo quadrático
def quadraticSieve(sieveIntervalMultiplicator, primeList, T, xZero):
    while(sieveIntervalMultiplicator <= 10000000):
        #encontrando os números B-smooth usando o método de sieve
        smoothNums, xlist, indexes = findSmooth(primeList, N, sieveIntervalMultiplicator, T, xZero)
        sieveIntervalMultiplicator = sieveIntervalMultiplicator * 10

        #Verifica se foram encontrados números smooth suficientes
        if len(smoothNums) < len(primeList):
            continue

        #montagem da matriz de expoentes
        isSquare, expMatrix = buildMatrix(smoothNums, primeList)

        #caso já tenha um quadrado imediato
        if isSquare == True:
            x = smoothNums.index(expMatrix)
            factor = gcd(xlist[x]+sqrt(expMatrix),N)
            print("X =", xlist[x])
            print("Y =", sqrt(expMatrix))
            print("Fator 1:", factor)
            print("Fator 2:", int(N / factor))
            return True

        #caso não, resolvo sol*matrix = 0 e testo todas as possíveis respostas
        else:
            solRows, marks, M = gaussElim(expMatrix)
            solutionVec = solveRow(solRows, M, marks, 0)
            factor, x, y = solve(solutionVec, smoothNums, xlist, N)

            #testo para todos os vetores solução até achar um que me dá a fatoração não trivial
            for K in range(1, len(solRows)):
                if (factor == 1 or factor == N):
                    solutionVec = solveRow(solRows, M, marks, K)
                    factor, x, y = solve(solutionVec, smoothNums, xlist, N)
                else:
                    print("X =", x)
                    print("Y =", y)
                    print("Fator 1:", factor)
                    print("Fator 2:", int(N / factor))
                    return True
                
    return False

#função para leitura do arquivo de entrada
def readFile(filePath):
    with open(filePath, 'r') as file:
        num = int(file.readline().strip())
    return num

#leitura da entrada
filePath = 'entrada.txt'
n = readFile(filePath)

#verificando se a entrada é primo
isPrime = isprime(n)
if isPrime:
    print("O número de entrada é primo!")
    print("Fator1: 1")
    print("Fator2:", n)
    exit(0)

#gerando B e a lista de primos heuristicamente
primeList, B = generateFactorBase(n)
print("Limite Superior para os primos usados no crivo: ", B)
print("A quantidade de primos que podem aparecer na fatoração após aplicação da heurística é: ", len(primeList))
print("O tamanho dos vetores incluídos na matriz é: ", len(primeList) + 1)

global N
N = n
T = 1
xZero = int(sqrt(n))
sieveIntervalMultiplicator = 10

resp = quadraticSieve(sieveIntervalMultiplicator, primeList, T, xZero)
if resp == 0:
    print("Não foi possível encontrar fatores com o B heurístico, digite um novo valor de B para ser aplicado: ")
    B = int(input())

while resp == 0:
    primeList, B = generateFactorBaseWithB(n, B)
    print("Limite Superior para os primos usados no crivo: ", B)
    print("A quantidade de primos que podem aparecer na fatoração após aplicação da heurística é: ", len(primeList))
    print("O tamanho dos vetores incluídos na matriz é: ", len(primeList) + 1)
    sieveIntervalMultiplicator = 10
    resp = quadraticSieve(sieveIntervalMultiplicator, primeList, T, xZero)
    if resp: break
    print("Não foi possível encontrar fatores com o B digitado, digite um novo valor de B para ser aplicado: ")
    B = int(input())