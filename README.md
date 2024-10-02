# Calculadora de Programación Lineal

Esta aplicación es una calculadora interactiva para resolver problemas de programación lineal utilizando los métodos Simplex, de la M Grande y de las Dos Fases. Está diseñada con una interfaz gráfica (GUI) creada en `Tkinter`.

## Características

- **Métodos de solución**: 
  - Simplex
  - Método de la M Grande
  - Método de las Dos Fases
- **Optimización**: 
  - Maximización
  - Minimización
- **Entrada de usuario**:
  - Coeficientes de la función objetivo
  - Restricciones y sus signos (<=, >=, =)
  - Variables no negativas por defecto
- **Salida**: 
  - Tabla de iteraciones
  - Solución óptima, en caso de existir
  - Detección de casos ilimitados o no factibles

## Requisitos

- Python 3.x
- tkinter

## Instalación

1. Clona o descarga este repositorio.
2. Asegúrate de tener `tkinter` instalado:

3. Ejecuta el archivo principal:
    ```bash
    python main.py
    ```

## Uso

1. Al abrir la aplicación, selecciona el **método de solución** en el primer menú desplegable.
2. Define el **tipo de optimización** (Maximizar o Minimizar).
3. Introduce los coeficientes de la **función objetivo** y las **restricciones** en los campos correspondientes.
4. Haz clic en el botón `Resolver` para obtener la solución del problema.

### Métodos de Solución

- **Método Simplex**: Utilizado para resolver problemas de maximización con restricciones tipo `<=`. Si hay restricciones `>=` o `=`, el método mostrará un mensaje de error.
- **Método de la M Grande**: Soluciona problemas con restricciones `>=` o `=` agregando variables artificiales.
- **Método de las Dos Fases**: Resuelve problemas con restricciones `>=` o `=`, trabajando en dos fases para encontrar una solución factible antes de optimizar.

## Ejemplo

1. Selecciona el **Método Simplex**.
2. Introduce la siguiente función objetivo para maximizar:
   - `Z = 3x1 + 5x2`
3. Añade las siguientes restricciones:
   - `x1 + 2x2 <= 10`
   - `x1 + x2 <= 6`
   - `x2 <= 4`
4. Presiona `Resolver`. La solución se mostrará en la columna de la derecha.
