EPS = 1e-4
lr = 0.01
kmax = 30_000
h = 1e-6
import random, numpy

norma = lambda x: numpy.linalg.norm(x, ord=2)

exemple = (
	{
		'f':lambda x,y: x**2 + y**2 - 2*x - 4*y - 1,
		'grad':lambda _,x,y:(2*x-2, 2*y-4)
	},
	{
		'f':lambda x,y: 3*(x**2) -12*x + 2*(y**2) +16*y -10,
		'grad':lambda _,x,y:(6*x-12, 4*y+16)
	},
	{
		'f':lambda x,y: x**2 - 4*x*y + 5*(y**2) - 4*y + 3,
		'grad':lambda _,x,y:(2*x-4*y, -4*x+10*y-4)
	},
	{
		'f':lambda x,y: (x**2)*y - 2*x*(y**2) + 3*x*y + 4,
		'grad':lambda _,x,y:(2*x*y-2*(y**2)+3*y, x**2-4*x*y+3*x)
	}
)

ex = exemple[1]

f,g = ex['f'], ex['grad']

def lr_beta(beta, f, grad, x, y):
	lr=1
	p=1
	while f(x-grad(f,x,y)[0],y-grad(f,x,y)[1])>f(x,y)-lr/2*norma(grad(f,x,y))**2 and p < 8:
		lr *= beta
		p += 1
	return lr


def aprox(f, grad, x, y, iters, lr):
	# x, y = x_init, y_init
	print(f"Init: {x=}, {y=}")
	k=0
	# hist = [(x, y)]
	g = grad(f, x, y)
	while EPS <= lr * norma(g) <= 1e11 and k <= kmax:
		g = grad(f, x, y)
		lr = lr #lr_beta(0.8, f, grad, x, y)
		x -= lr * g[0] # dx
		y -= lr * g[1] # dy
		# hist += [(x, y)]
		k += 1
	if lr * norma(g) <= EPS:
		print("Solutie in limita epsilon:", x, y)
	else:
		print("Divergenta")
	print("Iters:", k)


def grad(f, x: float, y: float):
	def g1(f, x: float, y: float):
		global h
		return (3*f(x, y) - 4*f(x-h, y) + f(x-2*h, y)) / (2 * h)
	def g2(f, x: float, y: float):
		global h
		return (3*f(x, y) - 4*f(x, y-h) + f(x, y-2*h)) / (2*h)
	return (g1(f, x, y), g2(f, x, y))


x, y = random.random(), random.random()

print("Aprox:")
aprox(f=f, grad=grad, x=x, y=y, iters=kmax, lr=lr)
print("Analitic:")
aprox(f=f, grad=g   , x=x, y=y, iters=kmax, lr=lr)