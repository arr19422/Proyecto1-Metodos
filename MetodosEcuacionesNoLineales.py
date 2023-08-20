import math
import numpy as np
import matplotlib.pyplot as plt
import tabulate as tb
from scipy.optimize import fsolve
from sympy.abc import x, y
import sympy as sp

##########################################################################################
##########################################################################################


def MetodoBiseccion(func, a, b, iters=8, tabla='N', tol=1E-16):
    # Evalúe los extremos del intervalo
    fa = func(a)
    fb = func(b)

    if tabla == 'S' or tabla == 'Y':  # Imprima la tabla si selecciona la opción
        print(' ')
        print('Resultados del Método de Bisección')  # Encabezado Tabla
        print(' iter        a          c           b        f(a)       f(c)       f(b)     Error aproximado')

    for i in range(iters):  # Empiece el ciclo for de la Bisección
        # Verifique si hay un cambio de signo entre a y b
        if fa * fb > 0:
            break  # Salga del ciclo for si no hay un cambio de signo
        xc = 0.5 * (a + b)  # Punto Medio
        ErrorAprox = abs(b-xc)
        fc = func(xc)  # Valor funcional en el punto medio
        if tabla == 'S' or tabla == 'Y':  # Imprima los resultados de cada iteración
            print('{:3d}    {:9.7f}   {:9.7f}  {:9.7f}   {:8.4f}   {:8.4f}   {:8.4f}     {:5.3e}'. format(
                i+1, a, xc, b, func(a), func(xc), func(b), ErrorAprox))
        if fc == 0:
            break  # Salga del ciclo for si se encuentra la raíz exacta
        # Actualice el nuevo intervalo
        if fa * fc < 0:  # El nuevo intervalo es [a xc]
            b = xc
        else:           # El nuevo intervalo es [xc xb]
            a = xc
            fa = fc  # Actualice el valor funcional del extremo izquierdo
        # fin if-else
        # Termine la iteración si el Error Aproximado es menor a una tolerancia
        if ErrorAprox < tol:
            print(' ')
            print('Se satisface la tolerancia especificada después de',
                  i + 1, 'iteraciones.')
            break

    # fin for

    # La raiz aproximada es el punto medio del último intervalo
    xr = 0.5 * (a + b)

    # Imprima los resultados del método de Bisección cuando No se satisface y se satisface la tolerancia
    if fa * fb > 0:
        print(
            'El método de Bisección no encontró una raíz aproximada después de {:3d} iteraciones'. format(i+1))
        print('El cambio de signo f(a)f(b)<0 no se satisface')
        print('--------------------------------------------------------------------')
        print(' ')
    else:
        print(' ')
        print('--------------------------------------------------------------------')
        print('Método de Bisección: Después de {:3d} iteraciones'. format(i+1))
        print(
            'Una raíz de f(x) se encuentra encerrada entre [{:7.5f} ,  {:7.5f}] . \n'. format(a, b))
        print('La raíz aproximada de f(x) es   {:2.7f} . '. format(xr))
        print('--------------------------------------------------------------------')
        print(' ')
    # Fin del if con el resultado final del Método

    return xr  # El único output es la raíz aproximada

# Fin de la función MetodoBisección

##########################################################################################
##########################################################################################


def MetodoSecante(funcion, x1, x2, iters=4, tabla='N', tol=1E-16):
    xi0 = x1  # El 1er estimado inicial es la iteracion x0
    xi1 = x2  # El 1er estimado inicial es la iteracion x0

    if tabla == 'S' or tabla == 'Y':  # Imprima la tabla si selecciona la opción
        print(' ')
        print('Iteraciones del Método de la Secante')  # Encabezado Tabla
        print(' iter         xi             f(xi)       Error Aproximado ')
        #print(' iter     xi        f(xi)      df(xi)    xi+1')
    for i in range(iters):
        if funcion(xi1)-funcion(xi0) == 0:  # Termine la iteración si la secante es horizontal
            break
        # Evalúe la función y la derivada en la estimación inicial
        fxi0 = funcion(xi0)  # valor funcional en x0
        fxi1 = funcion(xi1)  # valor funcional en x1
        msec = (fxi1-fxi0)/(xi1-xi0)  # Pendiente de la recta secante
        xi0 = xi1  # Actualice x0 a x1
        xi1 = xi1 - fxi1/msec  # Nueva estimación de la raíz
        ErrorAprox = abs(xi1-xi0)  # Error aproximado
        # Imprima la tabla de Resultados
        if tabla == 'S' or tabla == 'Y':  # Imprima la tabla si selecciona la opción
            print('{:3d}    {:13.10f}      {:8.3e}      {:5.3e} '. format(
                i+1, xi1, funcion(xi1), ErrorAprox))
            #print('{:3d}    {:7.4f}    {:7.4f}    {:7.4f}    {:7.4f}'.format(i+1, xi0, fxi1, msec, xi1))
        # Termine la iteración si el Error Aproximado es menor a una tolerancia
        if ErrorAprox < tol:
            print(' ')
            print('Se satisface la tolerancia especificada después de',
                  i+1, 'iteraciones.')
            break
        # fin if
    # fin ciclo for

    # Imprima los resultados del método de la secante
    if funcion(xi1)-funcion(xi0) == 0:  # Termine la iteración si la secante es horizontal
        print(
            'El método de la Secante no encontró una raíz aproximada después de {:3d} iteraciones'.format(i+1))
        print('La secante tiene pendiente horizontal entre xi+1 y xi')
        print('--------------------------------------------------------------------')
        print(' ')
    else:
        print(' ')
        print('--------------------------------------------------------------------')
        print('Método de la Secante: Después de {:3d} iteraciones'.format(i+1))
        print('La raíz aproximada de f(x) es   {:2.12f} . '.format(xi1))
        print('--------------------------------------------------------------------')
        print(' ')

        # Fin if-else

    return xi1  # El único output es la raíz aproximada

# Fin de la función MetodoSecante
