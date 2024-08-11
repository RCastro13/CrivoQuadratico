from sympy import isprime
from math import isqrt, sqrt, log, exp
from itertools import chain

#função para cálculo de mdc(a,b)
def gcd(a,b):
    if b == 0:
        return a
    elif a >= b:
        return gcd(b,a % b)
    else:
        return gcd(b,a)

#função para verificar se 'a' é um resíduo quadrado de n
def quad_residue(a,n):
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
    assert quad_residue(n, p) == 1, "não é um quadrado (mod p)"
    q = p - 1
    s = 0
    
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        r = pow(n, (p + 1) // 4, p)
        return r,p-r
    for z in range(2, p):
        if p - 1 == quad_residue(z, p):
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

    return (r,p-r)

#função que dado um vetor solução, calcula X e Y e retorna ambos e um dos fatores de N
def solve(solution_vec,smooth_nums,xlist,N):
    solution_nums = [smooth_nums[i] for i in solution_vec]
    x_nums = [xlist[i] for i in solution_vec]
    
    Asquare = 1
    for n in solution_nums:
        Asquare *= n
        
    x = 1
    for n in x_nums:
        x *= n

    y = isqrt(Asquare)
    
    factor = gcd(x-y,N)
    return factor, x, y

def gauss_elim(M):
#reduced form of gaussian elimination, finds rref and reads off the nullspace
    
    marks = [False]*len(M[0])
    for i in range(len(M)): #do for all rows
        row = M[i]
        
        for num in row: #search for pivot
            if num == 1:
                j = row.index(num) # column index
                marks[j] = True
                
                for k in chain(range(0,i),range(i+1,len(M))): #search for other 1s in the same column
                    if M[k][j] == 1:
                        for i in range(len(M[k])):
                            M[k][i] = (M[k][i] + row[i])%2
                break
    
    M = transpose(M)
        
    sol_rows = []
    for i in range(len(marks)): #find free columns (which have now become rows)
        if marks[i]== False:
            free_row = [M[i],i]
            sol_rows.append(free_row)
    
    if not sol_rows:
        #return("Nenhuma solução encontrada. É preciso de mais número smooth.")
        return
    #print("Foram encontradas {} soluções potenciais".format(len(sol_rows)))
    return sol_rows,marks,M

def solve_row(sol_rows,M,marks,K=0):
    solution_vec, indices = [],[]
    free_row = sol_rows[K][0] # may be multiple K
    for i in range(len(free_row)):
        if free_row[i] == 1: 
            indices.append(i)
    for r in range(len(M)): #rows with 1 in the same column will be dependent
        for i in indices:
            if M[r][i] == 1 and marks[r]:
                solution_vec.append(r)
                break
            
    solution_vec.append(sol_rows[K][1])       
    return(solution_vec)

def transpose(matrix):
#transpose matrix so columns become rows, makes list comp easier to work with
    new_matrix = []
    for i in range(len(matrix[0])):
        new_row = []
        for row in matrix:
            new_row.append(row[i])
        new_matrix.append(new_row)
    return(new_matrix)

def build_matrix(smooth_nums,factor_base):
# generates exponent vectors mod 2 from previously obtained smooth numbers, then builds matrix

    def factor(n,factor_base):#trial division from factor base
        factors = []
        if n < 0:
            factors.append(-1)
        for p in factor_base:
            if p == -1:
                pass
            else:
                while n % p == 0:
                    factors.append(p)
                    n //= p
        return factors

    M = []
    factor_base.insert(0,-1)
    for n in smooth_nums:
        exp_vector = [0]*(len(factor_base))
        n_factors = factor(n,factor_base)

        for i in range(len(factor_base)):
            if factor_base[i] in n_factors:
                exp_vector[i] = (exp_vector[i] + n_factors.count(factor_base[i])) % 2

        if 1 not in exp_vector: #search for squares
            return True, n
        else:
            pass
        
        M.append(exp_vector)  

    return(False, transpose(M))

def find_smooth(factor_base, N, I, T, xZero):
# tries to find B-smooth numbers in sieve_seq, using sieving

    def sieve_prep(N,sieve_int, xZero):
    # generates a sequence from Y(x) = x^2 - N, starting at x = root 
        sieve_seq = [x**2 - N for x in range(xZero,xZero+sieve_int)]
        #sieve_seq_neg = [x**2 - N for x in range(xZero,xZero-sieve_int,-1)]
        return sieve_seq

    sieve_seq = sieve_prep(N, I, xZero)
    sieve_list = sieve_seq.copy() # keep a copy of sieve_seq for later
    if factor_base[0] == 2:
        i = 0
        while sieve_list[i] % 2 != 0:
            i += 1
        for j in range(i,len(sieve_list),2): # found the 1st even term, now every other term will also be even
            while sieve_list[j] % 2 == 0: #account for powers of 2
                sieve_list[j] //= 2

    for p in factor_base[1:]: #not including 2
        residues = tonelliShanks(N,p) #finds x such that x^2 = n (mod p). There are two start solutions
        
        for r in residues:
            for i in range((r-xZero) % p, len(sieve_list), p): # Now every pth term will also be divisible
                while sieve_list[i] % p == 0: #account for prime powers
                    sieve_list[i] //= p
    xlist = [] #original x terms
    smooth_nums = []
    indices = [] # index of discovery
    
    for i in range(len(sieve_list)):
        if len(smooth_nums) >= len(factor_base) + T: #probability of no solutions is 2^-T
            break
        if sieve_list[i] == 1: # found B-smooth number
            smooth_nums.append(sieve_seq[i])
            xlist.append(i+xZero)
            indices.append(i)

    return(smooth_nums,xlist,indices)

#função para gerar o B heurístico e os fatores primos com base nesse
def generateFactorBase(N):
    B = exp(0.5 * sqrt(log(N) * log(log(N))))
    B = int(B)
    
    factor_base = []
    for p in range(2, B + 1):
        if isprime(p) and pow(N, (p - 1) // 2, p) == 1:
            factor_base.append(p)
    
    return factor_base, B

#função que recebe um valor de B e gera os fatores primos com base nele
def generateFactorBaseWithB(N, B):
    factor_base = []
    for p in range(2, B + 1):
        if isprime(p) and pow(N, (p - 1) // 2, p) == 1:
            factor_base.append(p)
    
    return factor_base, B

#função para leitura do arquivo de entrada
def lerEntradaArquivo(caminho_arquivo):
    with open(caminho_arquivo, 'r') as file:
        numero = int(file.readline().strip())
    return numero

#função principal para cálculo do crivo quadrático
def quadraticSieve(crivoIntervalMultiplicator, primeList, T, xZero):
    while(crivoIntervalMultiplicator <= 1000000):
        #encontrando os números B-smooth usando o método de sieve
        smooth_nums, xlist, indices = find_smooth(primeList, N, 10*crivoIntervalMultiplicator, T, xZero)
        crivoIntervalMultiplicator = crivoIntervalMultiplicator * 10

        #Verifica se foram encontrados números smooth suficientes
        if len(smooth_nums) < len(primeList):
            continue

        #montagem da matriz de expoentes
        is_square, t_matrix = build_matrix(smooth_nums,primeList)

        #caso já tenha um quadrado imediato
        if is_square == True:
            x = smooth_nums.index(t_matrix)
            factor = gcd(xlist[x]+sqrt(t_matrix),N)
            print("X =", xlist[x])
            print("Y =", sqrt(t_matrix))
            print("Fator 1:", factor)
            print("Fator 2:", int(N/factor))
            return True

        #caso não, resolvo sol*matrix = 0 e testo todas as possíveis respostas
        else:
            sol_rows, marks, M = gauss_elim(t_matrix)
            solution_vec = solve_row(sol_rows,M,marks,0)
            factor, x, y = solve(solution_vec,smooth_nums,xlist,N)

            #testo para todos os vetores solução até achar um que me dá a fatoração não trivial
            for K in range(1, len(sol_rows)):
                if (factor == 1 or factor == N):
                    solution_vec = solve_row(sol_rows,M,marks,K)
                    factor, x, y = solve(solution_vec,smooth_nums,xlist,N)
                else:
                    print("X =", x)
                    print("Y =", y)
                    print("Fator 1:", factor)
                    print("Fator 2:", int(N/factor))
                    return True

        print("Não foi possível encontrar fatores não triviais!")
    
    return False

#leitura da entrada
caminho_arquivo = 'entrada.txt'
n = lerEntradaArquivo(caminho_arquivo)

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
print("O tamanho dos vetores incluídos na matriz é: ", len(primeList)+1)

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
    print("O tamanho dos vetores incluídos na matriz é: ", len(primeList)+1)
    sieveIntervalMultiplicator = 10
    resp = quadraticSieve(sieveIntervalMultiplicator, primeList, T, xZero)
    if resp: break
    print("Não foi possível encontrar fatores com o B digitado, digite um novo valor de B para ser aplicado: ")
    B = int(input())