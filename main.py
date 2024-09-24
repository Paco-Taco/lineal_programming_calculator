def mostrar_menu():
    print("Seleccione el método que desea utilizar:")
    print("1. Método Simplex")
    print("2. Método de la M Grande")
    print("3. Método de Dos Fases")
    opcion = input("Ingrese el número de la opción: ")
    return opcion

def ingresar_problema():
    tipo = input("¿Desea maximizar o minimizar la función Z? (max/min): ").strip().lower()
    num_vars = int(input("Ingrese el número de variables de decisión: "))
    coeficientes_z = []
    print("Ingrese los coeficientes de la función objetivo Z:")
    for i in range(num_vars):
        coef = float(input(f"Coeficiente de x{i+1}: "))
        coeficientes_z.append(coef if tipo == 'max' else -coef)
    num_restricciones = int(input("Ingrese el número de restricciones: "))
    restricciones = []
    print("Ingrese las restricciones:")
    for i in range(num_restricciones):
        print(f"Restricción {i+1}:")
        coeficientes = []
        for j in range(num_vars):
            coef = float(input(f"Coeficiente de x{j+1}: "))
            coeficientes.append(coef)
        signo = input("Ingrese el signo de la restricción (<=, >=, =): ").strip()
        derecho = float(input("Ingrese el término independiente: "))
        restricciones.append({'coeficientes': coeficientes, 'signo': signo, 'derecho': derecho})
    return tipo, num_vars, coeficientes_z, restricciones

def inicializar_tabla(num_vars, coeficientes_z, restricciones):
    num_filas = len(restricciones) + 1
    num_columnas = num_vars + len(restricciones) + 1
    tabla = [[0 for _ in range(num_columnas)] for _ in range(num_filas)]
    for i, restriccion in enumerate(restricciones):
        for j in range(num_vars):
            tabla[i][j] = restriccion['coeficientes'][j]
        tabla[i][-1] = restriccion['derecho']
    for j in range(num_vars):
        tabla[-1][j] = -coeficientes_z[j]
    return tabla

def agregar_variables_holgura(tabla, restricciones, num_vars):
    indice_columna = num_vars
    for i, restriccion in enumerate(restricciones):
        if restriccion['signo'] == '<=':
            tabla[i][indice_columna] = 1
        elif restriccion['signo'] == '>=':
            tabla[i][indice_columna] = -1
        indice_columna += 1
    return tabla

def mostrar_tabla(tabla, variables_basicas):
    print("\nTabla Simplex:")
    encabezados = ["VB"] + [f"x{j+1}" for j in range(len(tabla[0])-1)] + ["LD"]
    print(" | ".join(f"{h:>6}" for h in encabezados))
    for i, fila in enumerate(tabla):
        vb = variables_basicas[i] if i < len(variables_basicas) else ""
        print(f"{vb:>6} | " + " | ".join(f"{v:6.2f}" for v in fila))
    print()

def simplex(tabla, num_vars):
    variables_basicas = [f"x{num_vars+i+1}" for i in range(len(tabla)-1)]
    iteracion = 0
    while True:
        mostrar_tabla(tabla, variables_basicas)
        # Verificar si la solución es óptima
        if all(c >= 0 for c in tabla[-1][:-1]):
            print("Se ha alcanzado la solución óptima.")
            break
        # Identificar la columna pivote
        columna_pivote = tabla[-1][:-1].index(min(tabla[-1][:-1]))
        # Identificar la fila pivote
        razones = []
        for i in range(len(tabla)-1):
            elemento = tabla[i][columna_pivote]
            if elemento > 0:
                razones.append(tabla[i][-1] / elemento)
            else:
                razones.append(float('inf'))
        fila_pivote = razones.index(min(razones))
        if razones[fila_pivote] == float('inf'):
            print("El problema no tiene solución óptima finita.")
            return
        # Mostrar variables que entran y salen
        print(f"Iteración {iteracion+1}:")
        print(f"Variable que entra: x{columna_pivote+1}")
        print(f"Variable que sale: {variables_basicas[fila_pivote]}")
        variables_basicas[fila_pivote] = f"x{columna_pivote+1}"
        # Realizar pivoteo
        pivote = tabla[fila_pivote][columna_pivote]
        tabla[fila_pivote] = [v / pivote for v in tabla[fila_pivote]]
        for i in range(len(tabla)):
            if i != fila_pivote:
                factor = tabla[i][columna_pivote]
                tabla[i] = [tabla[i][j] - factor * tabla[fila_pivote][j] for j in range(len(tabla[0]))]
        iteracion += 1
    # Mostrar resultados finales
    mostrar_tabla(tabla, variables_basicas)
    solucion = [0] * num_vars
    for i, vb in enumerate(variables_basicas):
        if vb.startswith('x') and int(vb[1:]) <= num_vars:
            indice = int(vb[1:]) - 1
            solucion[indice] = tabla[i][-1]
    print("Valor óptimo de Z:", tabla[-1][-1])
    for i, valor in enumerate(solucion):
        print(f"x{i+1} = {valor}")

def metodo_simplex():
    tipo, num_vars, coeficientes_z, restricciones = ingresar_problema()
    tabla = inicializar_tabla(num_vars, coeficientes_z, restricciones)
    tabla = agregar_variables_holgura(tabla, restricciones, num_vars)
    simplex(tabla, num_vars)

def metodo_m_grande():
    print("Método de la M Grande aún no implementado.")

def metodo_dos_fases():
    print("Método de Dos Fases aún no implementado.")

def main():
    opcion = mostrar_menu()
    if opcion == '1':
        metodo_simplex()
    elif opcion == '2':
        metodo_m_grande()
    elif opcion == '3':
        metodo_dos_fases()
    else:
        print("Opción inválida.")

if __name__ == "__main__":
    main()
