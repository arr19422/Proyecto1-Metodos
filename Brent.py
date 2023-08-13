import math
import numpy as np

def brent(func, lower_bound=0.0, upper_bound=1.0, tolerance=1e-12, max_iterations=100):
    """
    Encuentra la raíz de la función f(x) en el intervalo [a, b] utilizando el método de Brent.
    Parámetros
    ----------
    func : callable
        Función a la que se desea encontrar su raíz.
    lower_bound : float, opcional (valor por defecto = 0.0)
        Extremo inferior del intervalo.
    upper_bound : float, opcional (valor por defecto = 1.0)
        Extremo superior del intervalo.
    tolerance : float, opcional (valor por defecto = 1e-12)
        Tolerancia para el error en la aproximación de la raíz.
    max_iterations : int, opcional (valor por defecto = 100)
        Número máximo de iteraciones permitidas.
    Devuelve
    --------
    root : float
        Aproximación de la raíz de f(x).
    """
    fa = func(lower_bound)
    fb = func(upper_bound)
    
    if fa * fb >= 0:
        raise ValueError("La función debe tener signos opuestos en los extremos del intervalo.")
    
    if abs(fa) < abs(fb):
        lower_bound, upper_bound = upper_bound, lower_bound
        fa, fb = fb, fa
    
    c = lower_bound
    fc = fa
    d = None
    interpolation_flag = True
    iteration = 0
    
    while iteration < max_iterations and abs(upper_bound - lower_bound) > tolerance:
        if fa != fc and fb != fc:
            
            # Interpolación cuadrática inversa
            s = lower_bound * fb * fc / ((fa - fb) * (fa - fc)) + \
                upper_bound * fa * fc / ((fb - fa) * (fb - fc)) + \
                c * fa * fb / ((fc - fa) * (fc - fb))
        elif fa != fb:
            
            # Método de la secante
            s = upper_bound - fb * (upper_bound - lower_bound) / (fb - fa)
        
        if (s < (3 * lower_bound + upper_bound) / 4 or s > upper_bound) or \
           (interpolation_flag and abs(s - upper_bound) >= abs(upper_bound - c) / 2) or \
           (not interpolation_flag and abs(s - upper_bound) >= abs(c - d) / 2) or \
           (interpolation_flag and abs(upper_bound - c) < tolerance) or (not interpolation_flag and abs(c - d) < tolerance):
               
            # Método de la bisección
            s = (lower_bound + upper_bound) / 2
            interpolation_flag = True
        else:
            interpolation_flag = False
        
        fs = func(s)
        d = c
        c = upper_bound
        fc = fb
        
        if fa * fs < 0:
            upper_bound = s
            fb = fs
        else:
            lower_bound = s
            fa = fs
        
        if abs(fa) < abs(fb):
            lower_bound, upper_bound = upper_bound, lower_bound
            fa, fb = fb, fa
        
        iteration += 1
    
    if iteration == max_iterations:
        raise RuntimeError("El método de Brent no converge")
    else:
        return upper_bound
    

f1 = lambda x: (3*math.cos(2*x) - (2*math.sin(3*x)))/(x**2 + math.e**-x)
f2 = lambda x: 7*(math.e**(-x/7)) * math.sin(2*x) - 1
f3 = lambda x: 4 - (23*x) + (21*(x**2)) - (3*(x**3))

a = 0
b = 2
print ("\n\nRaiz para la funcion: (3*math.cos(2*x) - (2*math.sin(3*x)))/(x**2 + math.e**-x), para el intervalo [",a,"]", b, ", es : ", brent(f1, a, b))

a = 1
b = 2
print ("\n\nRaiz para la funcion: 7*(math.e**(-x/7)) * math.sin(2*x) - 1, para el intervalo [",a,"]", b, ", es : ", brent(f2, a, b))

a = 1
b = 2
print ("\n\nRaiz para la funcion: 4 - (23*x) + (21*(x**2)) - (3*(x**3)), para el intervalo [",a,"]", b, ", es : ", brent(f3, a, b))