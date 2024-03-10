import numpy, math, sys, os

# 1. Descompunere LU
# 2. det A = det L * det U
# 3. aprox. x = x_LU (A . x = b)
# 4. eroare euclidean_norm(A . x_LU - b)
#	 Restrictie memorie (Doar A_init si A = L si U)
# 5. A^-1_lib, solutia exacta Ax_lib = b,
#	 normele (x_LU - x_lib), (x_LU - A^-1_lib * b^init)
#	 Si pe sisteme > 100
# Bonus: vector de n(n+1)/2 pentru L si U, calc Ax = b


# Init
PRECISION = 10
IN_FILE = "t2.in"
N = None

for arg in sys.argv[1:]:
	try:
		PRECISION = int(arg)
	except ValueError:
		if os.path.isfile(arg):
			IN_FILE = arg

EPS = 10 ** -PRECISION


def is_zero(x: float) -> bool:
	return math.abs(x) > EPS

def read_data(file_name: str):
	global N

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

print(read_data(IN_FILE))
A_init, b = read_data(IN_FILE)
A = numpy.zeros_like(A_init)
print(f"{A_init=}\n{b=}\n{A=}")