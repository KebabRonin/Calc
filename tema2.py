import numpy, math, sys, os

# 1. Descompunere LU
# 2. det a = det L * det U
# 3. aprox. x = x_LU (a . x = b)
# 4. eroare euclidean_norm(a . x_LU - b)
#	 Restrictie memorie (Doar A_init si a = L si U)
# 5. a^-1_lib, solutia exacta Ax_lib = b,
#	 normele (x_LU - x_lib), (x_LU - a^-1_lib * b^init)
#	 Si pe sisteme > 100
# Bonus: vector de n(n+1)/2 pentru L si U, calc Ax = b


# Init: precizia si fisierul date ca argumente la linia de comanda
PRECISION = 10
IN_FILE = "t2.in"
norma = lambda x: numpy.linalg.norm(x, ord=2)

for arg in sys.argv[1:]:
	try:
		PRECISION = int(arg)
	except ValueError:
		if os.path.isfile(arg):
			IN_FILE = arg

EPS = 10 ** -PRECISION


def is_zero(x: float) -> bool:
	return math.fabs(x) < EPS


def calc_dets(A_init, a=None):
	if a is None:
		a = LU_decomp(A_init)
	return {"A": numpy.linalg.det(A_init),
			"L": math.prod([a[i, i] for i in range(a.shape[0])]),
			"U": 1}


def read_data(file_name: str):
	with open(file_name, "rt") as f:
		data = f.read()

	data = list(map(lambda a: float(a), data.split()))

	N = int(data[0])
	data = data[1:]

	if len(data) != N ** 2 + N:
		raise Exception("Invalid data format. Format: first number is N(matrix dimensions), next N * N are matrix values, next N are b(solutions) values")

	b = numpy.asarray(data[N ** 2:])
	A_init = numpy.matrix([[x for x in data[i*N:i*N+N]] for i in range(N)])
	return A_init, b


def LU_decomp(A_init):
	N = A_init.shape[0]
	a = numpy.zeros_like(A_init)

	for p in range(N):
		# L
		for i in range(p, N):
			a[i, p] = A_init[i, p] - sum([a[i, k] * a[k, p] for k in range(p)])

		# U
		if is_zero(a[p,p]):
			raise Exception(f"Matricea nu are derminant 0 pentru ca a[{p}, {p}] e 0")
		for i in range(p+1, N):
			a[p, i] = (A_init[p, i] - sum([a[p, k] * a[k, i] for k in range(p)])) / a[p, p]
	return a


def solve(A_init, b):
	a = LU_decomp(A_init)
	#Ax = b -> 1. Ly = b; 2. Ux = y;
	#Ly = b
	y = numpy.zeros_like(b)
	for i in range(b.shape[0]):
		y[i] = (b[i] - sum([a[i, j] * y[j] for j in range(i)])) / a[i, i]

	#Ux = y
	x = numpy.zeros_like(b)
	for i in reversed(range(b.shape[0])):
		x[i] = y[i] - sum([a[i, j] * x[j] for j in range(i+1, b.shape[0])])

	return x


print(read_data(IN_FILE))
A_init, b = read_data(IN_FILE)
A = LU_decomp(A_init)
print(f"{A_init=}\n{b=}")
print(f"\n----> Ex1: descompunerea LU intr-o singura matrice\n\n{A=}")
print(f"\n----> Ex2: determinanti:")
det = calc_dets(A_init, A)
print(f"{det['A']=} (numpy)")
print(f"{det['L']=}")
print(f"{det['U']=}")
print(f"det(A) numpy == det(L) * det(U) ? {is_zero(det['A'] - det['L'] * det['U'])}")

print(f"\n----> Ex3: Solutii:")
x = solve(A_init, b)
x_lib = numpy.linalg.solve(A_init, b)
A_inv_lib = numpy.linalg.inv(A_init)
err1 = norma(A_init * x.reshape((3, 1)) - b)
err2 = norma(x - x_lib)
err3 = norma(x - A_inv_lib * b.reshape((3, 1)))
print(f"Solutie calc : {x}")
print(f"Solutie numpy: {x_lib}")
print(f"{'norma(A_init * x - b)':<35}: {err1:<25}; {(err1 < 1e-8)=}")
print(f"{'norma(x - x_lib)':<35}: {err2:<25};")
print(f"{'norma(x - A_inv_lib * b)':<35}: {err3:<25};")
print(f"Inversa A^-1 numpy:\n{A_inv_lib}")

def d1(p, n):
	return sum([n - q for q in range(1, p+1)])

print([d1(i, 3) for i in range(3)])

def LU_decomp_bonus(A_init):
	N = A_init.shape[0]
	# L este retinut pe coloane, U pe linii
	L = numpy.ones((N*(N+1)//2,))
	U = numpy.ones((N*(N+1)//2,))

	for p in range(N):
		id_l = d1(p, N)
		id_u = d1(p, N)
		print(1, L, U)
		# L
		for i in range(p, N):
			L[id_l + i] = A_init[i, p] - sum([L[d1(i,N) + k] * U[k + d1(p,N)] for k in range(p)])

		# U
		if is_zero(L[id_l]):
			raise Exception(f"Matricea nu are derminant 0 pentru ca a[{p}, {p}] e 0")
		for i in range(p+1, N):
			U[id_u + i] = (A_init[p, i] - sum([L[d1(p,N) + k] * U[k +d1(i,N)] for k in range(p)])) / L[id_l]

	return (L, U)


def solve_bonus(A_init, b):
	L, U = LU_decomp_bonus(A_init)
	#Ax = b -> 1. Ly = b; 2. Ux = y;
	#Ly = b
	y = numpy.zeros_like(b)
	for i in range(b.shape[0]):
		y[i] = (b[i] - sum([a[i, j] * y[j] for j in range(i)])) / a[i, i]

	#Ux = y
	x = numpy.zeros_like(b)
	for i in reversed(range(b.shape[0])):
		x[i] = y[i] - sum([a[i, j] * x[j] for j in range(i+1, b.shape[0])])

	return x


# print(LU_decomp(A_init))
# print(LU_decomp_bonus(A_init))