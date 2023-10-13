import numpy as np
import matplotlib.pyplot as plt

def gen_array_espines():
    """Función que genera aleatoriamente un array bidimensional de espines 

    Returns:
        [np.array]: [Array bidimensional de espines de tamaño (long_lado, longlado)]
    """
    s = np.random.uniform(0, 1, (long_lado, long_lado))

    ii = (s < 0.5)
    s[ii] = -1
    s[~ii] = 1

    return s

def pos_vecinos(i, j):
    """Función que determina los indices i-1, i+1, j-1, j+1 para una posición (i, j) del array 
    bidimensional de espines a partir de las condiciones de frontera impuestas

    Args:
        i (int): [Posicion i en el array bidimensional de espines]
        j (int): [Posicion j en el array bidimensional de espines]

    Returns:
        [np.array]: [Array con las posciones i-1, i+1, j-1, j+1]
    """
    pos = np.array([i-1, i+1, j-1, j+1])
    ii_1 = (pos == -1)
    ii_2 = (pos == long_lado)

    pos[ii_1] = long_lado - 1
    pos[ii_2] = 0

    return pos

def var_termo(s):
    """Función que calcula las variables termodinámicas del sistema (Energía y Magnetización)
    para una configuración dada del sistema

    Args:
        s (np.array): [Array bidimensional de espines de tamaño (long_lado, longlado)]

    Returns:
        [tuple(float, float)]: [Valores de la Energía E y Magnetización M]
    """
    E = 0
    M = 0

    for i in range(0, long_lado):
        for j in range(0, long_lado):
            pos = pos_vecinos(i, j)

            M += s[i, j]
            E += -J*s[i, j]*(s[pos[0], j] + s[pos[1], j] + s[i, pos[2]] + s[i, pos[3]])

    E /= 2
    M /= N
    return E, M

def cambio_espin(s, beta):
    """Función que genera un cambio aleatorio de un espin y evalua si lo acepta o lo rechaza

    Args:
        s (np.array): [Array bidimensional de espines de tamaño (long_lado, longlado)]
        beta (float): [Temperatura inversa beta del sistema]

    Returns:
        [tuple(float, float)]: [Valores de la Energía E_f y Magnetización M_f finales después de haber aceptado o rechazado un cambio de espin]
    """
    pos_x = np.random.randint(0, long_lado)
    pos_y = np.random.randint(0, long_lado)

    E_0, M_0 = var_termo(s)
    s[pos_x, pos_y] = -s[pos_x, pos_y]
    E_f, M_f = var_termo(s)

    Delta_E = E_f - E_0

    if Delta_E > 0:
        p = np.random.uniform(0, 1)
        if p >= np.exp(-beta*Delta_E):
            s[pos_x, pos_y] = -s[pos_x, pos_y]

    return E_f, M_0

def evolucion(tiempo_sim, beta, s):
    """Función que evoluciona el sistema a partir de cambios aleatorios de espin en función del
    tiempo de evolución

    Args:
        tiempo_sim (np.array): [Array con los tiempos de simulación]
        beta (float): [Temperatura inversa beta del sistema]
        s (np.array): [Array bidimensional de espines de tamaño (long_lado, longlado)]

    Returns:
        [tuple(np.array, np.array)]: [Arrays de energía y magnetización para todos los tiempos de simulacion]
    """
    E_0, M_0 = var_termo(s)

    array_E = np.ones(cant_iter)
    array_M = np.ones(cant_iter)

    array_E[0] = E_0
    array_M[0] = M_0

    for i in tiempo_sim[1:]:
        E_f, M_f = cambio_espin(s, beta)
        array_E[i] = E_f
        array_M[i] = M_f

    return array_E, array_M

def promedios(array_E, array_M, t_estable):
    """Función que calcula los promedios de la Energía y la Magnetización a partir de un tiempo en el
    que se estabillizan los valores

    Args:
        array_E (np.array): [Array de energía para los tiempos de simulacion]
        array_M (np.array): [Array de magnetizacion para los tiempos de simulacion]
        t_estable (int): [Tiempo de simulación en que se estabilizan los valores de energía y magnetización]

    Returns:
        [tuple(float, float)]: [Tupla de energía y magnetización promedio]
    """
    E_prom = np.mean(array_E[t_estable:])
    M_prom = np.mean(array_M[t_estable:])

    return E_prom, M_prom

def graficas_EyM(tiempo_sim, array_E, array_M, T, j):
    """Función que realiza las graficas de la Energía y Magnetización en función del tiempo de simulación

    Args:
        tiempo_sim (np.array): [Array con los tiempos de simulación]
        array_E (np.array): [Array de energía para los tiempos de simulacion]
        array_M (np.array): [Array de magnetizacion para los tiempos de simulacion]
        T (float): [Temperatura en el que se encuentra el sistema]
        j (int): [Numero j que se asocia con los nombres de los archivos de las graficas]

    """
    razon_T = T/T_critico

    fig = plt.figure()
    plt.plot(tiempo_sim, array_E)
    plt.xlabel("$t$")
    plt.ylabel("$E$")
    plt.title("Energía vs Tiempo simulación para $T/T_c = {:.2f}$".format(razon_T))
    plt.savefig("Grafica_Evst_{}.png".format(j))
    plt.close(fig)

    fig = plt.figure()
    plt.plot(tiempo_sim, array_M)
    plt.xlabel("$t$")
    plt.ylabel("$M$")
    plt.title("Magnetización vs Tiempo simulación para $T/T_c = {:.2f}$".format(T/T_critico))
    plt.savefig("Grafica_Mvst_{}.png".format(j))
    plt.close(fig)

    return None

def grafica_MvsT(array_M_prom, array_T):
    """Función que realiza la grafica de Magnetización en función de la Temperatura

    Args:
        array_M_prom (np.array): [Array de los valores promedio de la Magnetización]
        array_T (np.array): [Array de los valores de Temperatura]
    
    """
    fig = plt.figure()
    plt.scatter(array_T/T_critico, array_M_prom, s=20)
    plt.xlabel("$T/T_c$")
    plt.ylabel("$M$")
    plt.title("Magnetización vs Temperatura")
    plt.savefig("Grafica_MvsT.png")
    plt.close(fig)
    return None

def main():
    """Funcion principal del código

    Returns:
        [int]: [Retorna 0 si el codigo se ejecuta correctamente]
    """
    array_T = np.linspace(0.01, 3, cant_Temp)*T_critico
    array_betas = 1/array_T

    array_E_prom = np.zeros(cant_Temp)
    array_M_prom = np.zeros(cant_Temp)

    j = 0

    for i in range(cant_Temp):
        beta = array_betas[i]

        s = gen_array_espines()
        array_E, array_M = evolucion(tiempo_sim, beta, s)
        
        if i == 0 or i == int(cant_Temp/4) or i == cant_Temp-1:
            j += 1
            graficas_EyM(tiempo_sim, array_E, array_M, array_T[i], j)

        E_prom, M_prom = promedios(array_E, array_M, t_estable)

        array_E_prom[i] = E_prom
        array_M_prom[i] = M_prom

    grafica_MvsT(array_M_prom, array_T)
    return 0


cant_iter = 8000
t_estable = 6000
tiempo_sim = np.arange(cant_iter)

N = 9**2
long_lado = int(np.sqrt(N))
J = 1

beta_critico = np.arcsinh(1)/(2*J)
T_critico = 1/beta_critico
cant_Temp = 100

main()
