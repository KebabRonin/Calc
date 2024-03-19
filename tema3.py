import numpy, math, sys, os, copy

# 1. Calc b
# 2. Descompunerea QR
# 3. Rezolvarea sistemului
# 4. Erori
# 5. Inversa A
# 6. Random init
# Bonus: limita A^(k)
# Restrictii: nr mic matrici, inmultirea cu partea superioara R, fara if-uri


# Init: precizia si fisierul date ca argumente la linia de comanda
PRECISION = 6
IN_FILE = "t3.in"
norma = lambda x: numpy.linalg.norm(x, ord=2)

for arg in sys.argv[1:]:
	try:
		PRECISION = int(arg)
	except ValueError:
		if os.path.isfile(arg):
			IN_FILE = arg

EPS = 10 ** -PRECISION

def read_data(file_name: str):
	with open(file_name, "rt") as f:
		data = f.read()

	data = list(map(lambda a: float(a), data.split()))

	n = int(data[0])
	data = data[1:]

	if len(data) != n ** 2 + n:
		raise Exception("Invalid data format. Format: first number is n(matrix dimensions), next n * n are matrix values, next n are s(vector) values")

	s = numpy.asarray(data[n ** 2:])
	A_init = numpy.matrix([[x for x in data[i*n:i*n+n]] for i in range(n)])
	return A_init, s

def rand_data(n):
	A_init = numpy.random.rand(n, n)
	s = numpy.random.rand(n)
	return A_init, s


def is_zero(x: float) -> bool:
	return math.fabs(x) < EPS

def calc_b(A_init, s):
	b = numpy.zeros_like(s)
	for i in range(s.shape[0]):
		b[i] = sum([s[j] * A_init[i, j] for j in range(s.shape[0])])
	return b

def QR_decomp(A_init, b_init):
	b = copy.deepcopy(b_init)
	Qt = numpy.identity(A_init.shape[0])
	R = copy.deepcopy(A_init)
	n = A_init.shape[0]
	for r in range(n-1):
		sigma = sum([R[i, r] ** 2 for i in range(r, n)])
		if is_zero(sigma):
			print("Matricea A e singulara")
			break
		k = math.sqrt(sigma) * (-1 if R[r, r] > 0 else 1)
		beta = sigma - k * R[r, r]

		u = numpy.zeros((n,))
		u[r] = R[r, r] - k
		for i in range(r+1, n):
			u[i] = R[i, r]

		# R = P * R
		for j in range(r+1, n):
			gamma = sum([u[i] * R[i, j] for i in range(r, n)]) / beta
			for i in range(r, n):
				R[i, j] -= gamma * u[i]
		R[r, r] = k
		for i in range(r+1, n):
			R[i, r] = 0

		# b = P * b
		gamma = sum([u[i] * b[i] for i in range(r, n)]) / beta
		for i in range(r, n):
			b[i] -= gamma * u[i]

		# Q = P * Q
		for j in range(n):
			gamma = sum([u[i] * Qt[i, j] for i in range(r, n)]) / beta
			for i in range(r, n):
				Qt[i, j] -= gamma * u[i]

	Q = numpy.transpose(Qt)
	return Q, R


def solve(A_init, b, Q, R):
	# Ceva de genul asta, nu e gata
	# a = LU_decomp(A_init)
	n = b.shape[0]
	#Ax = b -> 1. Ly = b; 2. Ux = y;
	#Ly = b
	y = numpy.zeros_like(b)
	for i in range(n):
		y[i] = (b[i] - sum([a[i, j] * y[j] for j in range(i)])) / a[i, i]

	#Ux = y
	x = numpy.zeros_like(b)
	for i in reversed(range(n)):
		x[i] = y[i] - sum([a[i, j] * x[j] for j in range(i+1, n)])

	return x

def sols(A_init, s):
	b_init = calc_b(A_init, s)
	n = s.shape[0]
	print(f"\n----> Ex3: Solutii:")
	xQR = numpy.linalg.qr(A_init)
	xHouseholder = QR_decomp(A_init, b_init)
	print("Lib:")
	print(xQR[0])
	print(xQR[1])
	print("Tema:")
	print(xHouseholder[0])
	print(xHouseholder[1])
	# err = norma(xQR - xHouseholder)
	# print(f"norma(xQR - xHouseholder): {err}")

def afisare(A_init, s):
	sols(A_init, s)

A_init, s = read_data(IN_FILE)
Q, R = QR_decomp(A_init, calc_b(A_init, s))
print(calc_b(A_init, s))
print(Q)
print(R)
afisare(A_init, s)