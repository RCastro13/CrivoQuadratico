from sympy import isprime
from math import isqrt, sqrt, log, exp
#from shanks import STonelli
from itertools import chain

def gcd(a,b):
    if b == 0:
        return a
    elif a >= b:
        return gcd(b,a % b)
    else:
        return gcd(b,a)
    

def quad_residue(a,n):
    #checks if a is quad residue of n
    l=1
    q=(n-1)//2
    x = q**l
    if x==0:
        return 1
        
    a =a%n
    z=1
    while x!= 0:
        if x%2==0:
            a=(a **2) % n
            x//= 2
        else:
            x-=1
            z=(z*a) % n

    return z

# Returns k such that b^k = 1 (mod p)
def order(p,q):
	if gcd(p,q) !=1:
		#print(p,' and ',q,'are not co-prime')
		return -1
	k=3
	while True:
		if pow(q,k,p)==1:
			return k
		k+=1

# function return p - 1 (= x argument) as x * 2^e,
# where x will be odd sending e as reference because
# updation is needed in actual e
def convertx2e(x):
	e=0
	while x%2==0:
		x/=2
		e+=1
	return int(x),e

# Main function for finding the modular square root
def STonelli(n, p): #tonelli-shanks to solve modular square root, x^2 = N (mod p)
    assert quad_residue(n, p) == 1, "not a square (mod p)"
    q = p - 1
    s = 0
    
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        r = pow(n, (p + 1) // 4, p)
        return r,p-r
    for z in range(2, p):
        #print(quad_residue(z, p))
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

def solve(solution_vec,smooth_nums,xlist,N):
    
    solution_nums = [smooth_nums[i] for i in solution_vec]
    x_nums = [xlist[i] for i in solution_vec]
    
    Asquare = 1
    for n in solution_nums:
        Asquare *= n
        
    b = 1
    for n in x_nums:
        b *= n

    a = isqrt(Asquare)
    
    factor = gcd(b-a,N)
    return factor, b, a

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

def find_smooth(factor_base,N,I):
# tries to find B-smooth numbers in sieve_seq, using sieving

    def sieve_prep(N,sieve_int):
    # generates a sequence from Y(x) = x^2 - N, starting at x = root 
        sieve_seq = [x**2 - N for x in range(root,root+sieve_int)]
        #sieve_seq_neg = [x**2 - N for x in range(root,root-sieve_int,-1)]
        return sieve_seq

    sieve_seq = sieve_prep(N,I)
    sieve_list = sieve_seq.copy() # keep a copy of sieve_seq for later
    if factor_base[0] == 2:
        i = 0
        while sieve_list[i] % 2 != 0:
            i += 1
        for j in range(i,len(sieve_list),2): # found the 1st even term, now every other term will also be even
            while sieve_list[j] % 2 == 0: #account for powers of 2
                sieve_list[j] //= 2

    for p in factor_base[1:]: #not including 2
        residues = STonelli(N,p) #finds x such that x^2 = n (mod p). There are two start solutions
        
        for r in residues:
            for i in range((r-root) % p, len(sieve_list), p): # Now every pth term will also be divisible
                while sieve_list[i] % p == 0: #account for prime powers
                    sieve_list[i] //= p
    xlist = [] #original x terms
    smooth_nums = []
    indices = [] # index of discovery
    
    for i in range(len(sieve_list)):
        if len(smooth_nums) >= len(factor_base)+T: #probability of no solutions is 2^-T
            break
        if sieve_list[i] == 1: # found B-smooth number
            smooth_nums.append(sieve_seq[i])
            xlist.append(i+root)
            indices.append(i)

    return(smooth_nums,xlist,indices)

def generate_factor_base(N):
    #gera a base de fatores até o limite B usando o SymPy
    B = exp(0.5 * sqrt(log(N) * log(log(N))))
    B = int(B)
    
    factor_base = []
    for p in range(2, B + 1):
        if isprime(p) and pow(N, (p - 1) // 2, p) == 1:
            factor_base.append(p)
    
    return factor_base, B

def generate_factor_base_with_B(N, B):
    #gera a base de fatores até o limite B usando o SymPy
    
    factor_base = []
    for p in range(2, B + 1):
        if isprime(p) and pow(N, (p - 1) // 2, p) == 1:
            factor_base.append(p)
    
    return factor_base, B

def lerEntradaArquivo(caminho_arquivo):
    with open(caminho_arquivo, 'r') as file:
        numero = int(file.readline().strip())
    return numero

def quadraticSieve(crivoIntervalMultiplicator, primeList):
    while(crivoIntervalMultiplicator <= 1000000):
        #print("Testando com o intervalo de crivo =", 10*crivoIntervalMultiplicator)
        
        #encontrando os números B-smooth usando o método de sieve
        smooth_nums, xlist, indices = find_smooth(primeList, N,10*crivoIntervalMultiplicator)
        crivoIntervalMultiplicator = crivoIntervalMultiplicator * 10
        #print("Foram encontrados {} números B-smooth.".format(len(smooth_nums)))
        #print(smooth_nums)

        if len(smooth_nums) < len(primeList):
            #print("Não foram encontrados números smooth suficientes..")
            continue

        #montagem da matriz de expoentes
        is_square, t_matrix = build_matrix(smooth_nums,primeList)
        #print(t_matrix)

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
            sol_rows,marks,M = gauss_elim(t_matrix)
            solution_vec = solve_row(sol_rows,M,marks,0)
            factor, b, a = solve(solution_vec,smooth_nums,xlist,N)

            #testo para todos os vetores solução até achar um que me dá a fatoração não trivial
            for K in range(1,len(sol_rows)):
                if (factor == 1 or factor == N):
                    solution_vec = solve_row(sol_rows,M,marks,K)
                    factor, b, a = solve(solution_vec,smooth_nums,xlist,N)
                else:
                    print("X =", b)
                    print("Y =", a)
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
primeList, B = generate_factor_base(n)
#print("Testando o B heurístico igual a ", B)
#print("Lista de Primos gerada: ", primeList)
print("Limite Superior para os primos usados no crivo: ", B)
print("A quantidade de primos que podem aparecer na fatoração após aplicação da heurística é: ", len(primeList))
print("O tamanho dos vetores incluídos na matriz é: ", len(primeList)+1)

global N
global root
global T #tolerance factor
N,root,K,T = n,int(sqrt(n)),0,1
crivoIntervalMultiplicator = 10

resp = quadraticSieve(crivoIntervalMultiplicator, primeList)
if resp == 0:
    print("Como o B heurístico falhou, digite um novo valor de B para ser aplicado: ")
    B = int(input())

while resp == 0:
    primeList, B = generate_factor_base_with_B(n, B)
    #print("Testando o B igual a ", B)
    #print("Lista de Primos gerada: ", primeList)
    print("Limite Superior para os primos usados no crivo: ", B)
    print("A quantidade de primos que podem aparecer na fatoração após aplicação da heurística é: ", len(primeList))
    print("O tamanho dos vetores incluídos na matriz é: ", len(primeList)+1)
    crivoIntervalMultiplicator = 10
    resp = quadraticSieve(crivoIntervalMultiplicator, primeList)
    if resp: break
    print("Como o B digitado falhou, digite um novo valor de B para ser aplicado: ")
    B = int(input())