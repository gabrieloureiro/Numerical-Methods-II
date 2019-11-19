#EQUIPE 10
#GABRIEL FRANÇA LOUREIRO  //#MATRÍCULA: 388835
#ABNER DE LIMA ARAÚJO    //#MATRÍCULA: 398067
#MN2 2019
from __future__ import print_function
from math import sqrt, ceil

def ft(t):
    if ( (t < 0) or (t > 1)):
        return 0;
    elif (t <= 0.5):
        return 4*t
    elif (t <= 1):
        return 4*(1-t)
    else:
        return 0

um = 1
zero = 0
while zero < um:
    matricula = (input("\n>>> Insira a matrícula a ser utilizada: "))
    if len(matricula)>6 or len(matricula)<6: 
        print("\nMatrícula inválida!")
    else:
        zero = 1

A, B, C, D, E, F = map(int, list(matricula))

m = 1 + (A + B + C + D + E + F) % 4
k = 4
w = sqrt(k/m)
z = 0.05

tn = 1.2

v0 = 0
x0 = (E + F) % 3

h = 0.6 # Delta T


e = 0.0001


class S:
    def __init__(self, v, x):
        self.v = v
        self.x = x
        
    def __add__(self, other):
        return S(self.v+other.v, self.x+other.x)
    def __sub__(self, other):
        return S(self.v-other.v, self.x-other.x)
    def __mul__(self, other):
        if type(other) is not int and type(other) is not float:
            raise Exception('MultiplicationError')
        return S(self.v*other, self.x*other)
    def __truediv__(self, other):
        if type(other) is not int and type(other) is not float:
            raise Exception('DivisionError')
        return S(self.v/other, self.x/other)
    
    def __radd__(self, value):
        return S(self.v+value, self.x+value)
    def __rsub__(self, value):
        return S(self.v-value, self.x-value)
    def __rmul__(self, value):
        if type(value) is not int and type(value) is not float:
            raise Exception('MultiplicationError')
        return S(self.v*value, self.x*value)
    def __rtruediv__(self, value):
        if type(value) is not int and type(value) is not float:
            raise Exception('DivisionError')
        return S(self.v/value, self.x/value)
    
    def __str__(self):
        return "(v={}, x={})".format(self.v, self.x)




def Function(Si, t):
    return S( (ft(t)/m - 2*z*w*Si.v - pow(w,2)*Si.x), Si.v)


# Programa principal
stopCondition = 0
erroRelativo = 1
s = S(v0, x0)
while stopCondition < 2:
    #valores iniciais
    ############# Range Kutta ############
    x_old = s.x
    s = S(v0, x0)
    t0 = 0
    tf = 1.2
    while t0 < tf:

        k1 = Function(s, t0)

        k2 = Function( (s + k1 * (h/2)) , (t0 + h/2) )
   
        k3 = Function( (s + k2 * (h/2)) , (t0 + h/2) )

        k4 = Function(s + k3 * h, t0 + h)

        s = s + h* (k1 + 2*k2 + 2*k3 + k4)/6
        t0 += h
    #print(s)
    ################# FIM #################

    # Verifica o erro
    erroRelativo = abs(s.x - x_old)/abs(s.x)
    #print("Erro relativo: {}".format(erroRelativo))

    if (erroRelativo < e):
        stopCondition += 1
    else:
        stopCondition = 0
    h = h/2


print("\n⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥ >>> AP3 QUESTÃO ➀ <<< ⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤\n")
print("\nDensenvolvendo a questão, temos que:\n>>> dx/dt = x' = v\n>>> (d^2)x(t)/dt^2 = v'\n>>> x(t) = x\n\nLogo, temos: v' = f(t)/m - 2*ζ*ω*v - ω^2*x\n")
print("\nCriando a estrutura S = (v  x), temos que dS/dt = F(S,t) = (v' x')\n")
print("\nPortanto, se S_0 = (v_0 x_0), então S_0 = (0  2)\n")
print("\nCom t = 1.2s, aplicaremos Range-Kutta de 4ª Ordem com ∆t = 0.6\n")    
print("v(1.2) = {}".format(s.v))
print("x(1.2) = {}".format(s.x))
print("\n")