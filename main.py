import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

class LinearProgrammingSolver:
    def __init__(self, root):
        self.root = root
        self.root.title("Solver de Programación Lineal")

        # Hacer que la ventana sea redimensionable
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        # Variables
        self.num_variables = 2
        self.num_constraints = 3
        self.method = tk.StringVar(value="Método Simplex")
        self.optimization_type = tk.StringVar(value="Maximizar")

        self.create_main_frame()

    def create_main_frame(self):
        # Main frame
        main_frame = ttk.Frame(self.root)
        main_frame.grid(row=0, column=0, sticky="nsew")
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        # Configurar expansión del main_frame
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(6, weight=1)

        # Método de Solución (con Combobox)
        method_frame = ttk.LabelFrame(main_frame, text="Método de Solución")
        method_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        method_frame.columnconfigure(0, weight=1)

        method_combobox = ttk.Combobox(
            method_frame,
            textvariable=self.method,
            values=["Método Simplex", "Método de la M grande", "Método de las dos fases"],
            state="readonly",
        )
        method_combobox.grid(row=0, column=0, padx=10, pady=5, sticky="ew")
        method_combobox.current(0)  # Seleccionar por defecto el primer método

        # Optimización
        opt_frame = ttk.LabelFrame(main_frame, text="Tipo de Optimización")
        opt_frame.grid(row=1, column=0, padx=10, pady=10, sticky="ew")

        opt_frame.columnconfigure(0, weight=1)
        opt_frame.columnconfigure(1, weight=1)

        ttk.Radiobutton(
            opt_frame, text="Maximizar", variable=self.optimization_type, value="Maximizar"
        ).grid(row=0, column=0, padx=5, pady=5, sticky="w")
        ttk.Radiobutton(
            opt_frame, text="Minimizar", variable=self.optimization_type, value="Minimizar"
        ).grid(row=0, column=1, padx=5, pady=5, sticky="w")

        # Función objetivo
        obj_frame = ttk.LabelFrame(main_frame, text="Función Objetivo (Z)")
        obj_frame.grid(row=2, column=0, padx=10, pady=10, sticky="ew")
        obj_frame.columnconfigure('all', weight=1)

        self.obj_coeffs = []
        for i in range(self.num_variables):
            entry = ttk.Entry(obj_frame, width=5)
            entry.grid(row=0, column=i * 2, padx=5, pady=5, sticky="ew")
            lbl = ttk.Label(obj_frame, text=f"x{i + 1}")
            lbl.grid(row=0, column=i * 2 + 1, padx=5, pady=5, sticky="w")
            self.obj_coeffs.append(entry)

        # Restricciones
        constr_frame = ttk.LabelFrame(main_frame, text="Restricciones")
        constr_frame.grid(row=3, column=0, padx=10, pady=10, sticky="ew")
        constr_frame.columnconfigure('all', weight=1)

        self.constraints_coeffs = []
        self.constraints_signs = []
        self.constraints_rhs = []

        for i in range(self.num_constraints):
            constr_frame.rowconfigure(i, weight=1)
            row_entries = []
            for j in range(self.num_variables):
                entry = ttk.Entry(constr_frame, width=5)
                entry.grid(row=i, column=j * 2, padx=5, pady=5, sticky="ew")
                lbl = ttk.Label(constr_frame, text=f"x{j + 1}")
                lbl.grid(row=i, column=j * 2 + 1, padx=5, pady=5, sticky="w")
                row_entries.append(entry)
            self.constraints_coeffs.append(row_entries)

            sign_cb = ttk.Combobox(constr_frame, values=["<=", ">=", "="], width=5, state="readonly")
            sign_cb.grid(row=i, column=self.num_variables * 2, padx=5, pady=5, sticky="ew")
            sign_cb.current(0)
            self.constraints_signs.append(sign_cb)

            rhs_entry = ttk.Entry(constr_frame, width=5)
            rhs_entry.grid(row=i, column=self.num_variables * 2 + 1, padx=5, pady=5, sticky="ew")
            self.constraints_rhs.append(rhs_entry)

        # Restricciones de variables
        var_constraints_frame = ttk.LabelFrame(
            main_frame, text="Restricciones de Variables (activadas por defecto)"
        )
        var_constraints_frame.grid(row=4, column=0, padx=10, pady=10, sticky="ew")
        var_constraints_frame.columnconfigure('all', weight=1)

        for i in range(self.num_variables):
            lbl = ttk.Label(var_constraints_frame, text=f"x{i + 1} ≥ 0")
            lbl.grid(row=0, column=i, padx=5, pady=5, sticky="w")

        # Botón para resolver
        solve_button = ttk.Button(main_frame, text="Resolver", command=self.solve)
        solve_button.grid(row=5, column=0, pady=20, sticky="ew")

        # Salida
        output_frame = ttk.LabelFrame(main_frame, text="Solución")
        output_frame.grid(row=6, column=0, padx=10, pady=10, sticky="nsew")
        output_frame.columnconfigure(0, weight=1)
        output_frame.rowconfigure(0, weight=1)

        self.output_text = tk.Text(output_frame, wrap="word", height=10)
        self.output_text.grid(row=0, column=0, sticky="nsew")

        # Añadir una barra de desplazamiento para el texto de salida
        scrollbar = ttk.Scrollbar(output_frame, orient="vertical", command=self.output_text.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.output_text.configure(yscrollcommand=scrollbar.set)

    def solve(self):
        # Recolectar datos
        try:
            c = [float(entry.get()) for entry in self.obj_coeffs]
            if self.optimization_type.get() == "Minimizar":
                c = [-coeff for coeff in c]
            A = []
            b = []
            signs = []
            for i in range(self.num_constraints):
                row = [float(entry.get()) for entry in self.constraints_coeffs[i]]
                A.append(row)
                b.append(float(self.constraints_rhs[i].get()))
                signs.append(self.constraints_signs[i].get())
            # Restricciones de variables
            # var_constraints = [var.get() for var in self.var_constraints]
        except ValueError:
            messagebox.showerror("Error de Entrada", "Por favor, ingrese valores numéricos válidos.")
            return
        
        # Implementar el método seleccionado
        if self.method.get() == "Método Simplex":
            self.simplex_method(c, A, b, signs)
        elif self.method.get() == "Método de la M grande":
            self.big_m_method(c, A, b, signs)
        elif self.method.get() == "Método de las dos fases":
            self.two_phase_method(c, A, b, signs)
        else:
            messagebox.showerror("Error de Método", "Por favor, seleccione un método válido.")
            return
        
    def simplex_method(self, c, A, b, signs):
        self.output_text.delete(1.0, tk.END)
        
        # Verificar si hay restricciones '>=' o '='
        if '>=' in signs or '=' in signs:
            self.output_text.insert(tk.END, "El Método Simplex estándar no puede resolver restricciones '≥' o '='.\n")
            return
        
        # Convertir a forma estándar
        num_variables = len(c)
        num_constraints = len(A)
        
        # Extender coeficientes de la función objetivo
        c_extended = c[:]
        
        # Agregar variables de holgura
        tableau = []
        basis = []
        slack_var_index = num_variables
        for i in range(num_constraints):
            row = A[i][:]
            slack = [0] * num_constraints
            slack[i] = 1
            c_extended.append(0)
            basis.append(slack_var_index)
            slack_var_index +=1
            row.extend(slack)
            row.append(b[i])
            tableau.append(row)
        
        # Mostrar el tableau inicial
        self.output_text.insert(tk.END, "Tabla Inicial:\n")
        self.display_tableau(tableau, c_extended, basis)
        
        # Iniciar iteraciones
        iteration = 0
        max_iterations = 100  # Para evitar ciclos infinitos
        while iteration < max_iterations:
            iteration +=1
            self.output_text.insert(tk.END, f"\nIteración {iteration}:\n")
            # Calcular las ganancias relativas (cj - zj)
            zj = [0]*len(c_extended)
            for i in range(len(basis)):
                for j in range(len(c_extended)):
                    zj[j] += c_extended[basis[i]] * tableau[i][j]
            cj_zj = [c_extended[j] - zj[j] for j in range(len(c_extended))]
            
            # Mostrar el tableau con las ganancias relativas
            self.display_tableau(tableau, c_extended, basis, cj_zj=cj_zj)
            
            # Verificar optimalidad
            if all(value <= 1e-5 for value in cj_zj):
                self.output_text.insert(tk.END, "Se encontró la solución óptima.\n")
                break
            
            # Variable entrante (mayor cj - zj)
            entering_candidates = [j for j in range(len(cj_zj)) if cj_zj[j] > 1e-5]
            if not entering_candidates:
                self.output_text.insert(tk.END, "No hay variables entrantes. El problema puede tener soluciones múltiples.\n")
                break
            entering = max(entering_candidates, key=lambda j: cj_zj[j])
            self.output_text.insert(tk.END, f"Variable entrante: x{entering+1}\n")
            
            # Calcular razones (b_i / a_i)
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][entering] > 1e-5:
                    ratios.append(tableau[i][-1] / tableau[i][entering])
                else:
                    ratios.append(float('inf'))
            
            # Mostrar las razones en el tableau
            self.display_tableau(tableau, c_extended, basis, ratios=ratios, entering=entering)
            
            # Verificar si es ilimitado
            if all(r == float('inf') for r in ratios):
                self.output_text.insert(tk.END, "La solución es ilimitada.\n")
                return
            
            # Variable saliente (mínima razón positiva)
            leaving = ratios.index(min(ratios))
            self.output_text.insert(tk.END, f"Variable saliente: x{basis[leaving]+1}\n")
            
            # Pivoteo
            pivot_element = tableau[leaving][entering]
            tableau[leaving] = [x / pivot_element for x in tableau[leaving]]
            
            for i in range(len(tableau)):
                if i != leaving:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][j] - factor * tableau[leaving][j] for j in range(len(tableau[i]))]
            
            # Actualizar base
            basis[leaving] = entering
        else:
            self.output_text.insert(tk.END, "Se alcanzó el número máximo de iteraciones. El problema puede no tener solución óptima.\n")
            return
        
        # Extraer solución
        solution = [0]*len(c_extended)
        for i in range(len(basis)):
            if basis[i] < len(solution):
                solution[basis[i]] = tableau[i][-1]
        
        z = sum(c_extended[i]*solution[i] for i in range(len(c_extended)))
        if self.optimization_type.get() == "Minimizar":
            z = -z
        self.output_text.insert(tk.END, f"\nSolución Óptima:\n")
        for i in range(num_variables):
            self.output_text.insert(tk.END, f"x{i+1} = {solution[i]:.2f}\n")
        self.output_text.insert(tk.END, f"Z = {z:.2f}\n")
        
    def display_tableau(self, tableau, c_extended, basis, cj_zj=None, ratios=None, entering=None):
        num_vars = len(c_extended)
        num_decision_vars = self.num_variables
        num_slack_surplus_vars = num_vars - num_decision_vars
        # Construir encabezados
        header = "+------+"
        header += "-----------" * num_decision_vars + "+"
        header += "-----------" * num_slack_surplus_vars + "+"
        header += "----------+\n"
        
        title_row = "| Base |"
        for i in range(num_decision_vars):
            title_row += f"   X{i+1}   |"
        for i in range(num_slack_surplus_vars):
            title_row += f"   S{i+1}   |"
        title_row += " Solución |\n"
        
        header += title_row
        header += "+------+"
        header += "-----------" * num_decision_vars + "+"
        header += "-----------" * num_slack_surplus_vars + "+"
        header += "----------+\n"
        
        self.output_text.insert(tk.END, header)
        
        # Mostrar filas del tableau
        for i in range(len(tableau)):
            row = tableau[i]
            base_var = basis[i]
            row_str = f"|  x{base_var+1:<3}|"
            for val in row[:-1]:
                row_str += f" {val:>8.2f} |"
            row_str += f" {row[-1]:>8.2f} |"
            # Añadir razón si corresponde
            if ratios:
                if ratios[i] != float('inf'):
                    ratio_str = f"{row[-1]:.2f}/{row[entering]:.2f}={ratios[i]:.2f}"
                else:
                    ratio_str = "Inf"
                row_str += f" {ratio_str} |"
            row_str += "\n"
            self.output_text.insert(tk.END, row_str)
            self.output_text.insert(tk.END, "+------+")
            self.output_text.insert(tk.END, "-----------" * num_decision_vars + "+")
            self.output_text.insert(tk.END, "-----------" * num_slack_surplus_vars + "+")
            self.output_text.insert(tk.END, "----------+\n")
        
        # Mostrar fila de Z
        if cj_zj:
            z_value = sum(c_extended[basis[i]] * tableau[i][-1] for i in range(len(basis)))
            row_str = f"|  Z   |"
            for val in cj_zj:
                row_str += f" {val:>8.2f} |"
            row_str += f" {z_value:>8.2f} |\n"
            self.output_text.insert(tk.END, row_str)
            self.output_text.insert(tk.END, "+------+")
            self.output_text.insert(tk.END, "-----------" * num_decision_vars + "+")
            self.output_text.insert(tk.END, "-----------" * num_slack_surplus_vars + "+")
            self.output_text.insert(tk.END, "----------+\n")
    
    def big_m_method(self, c, A, b, signs):
        self.output_text.delete(1.0, tk.END)
        
        M = 1e6  # Valor grande para M
        num_variables = len(c)
        num_constraints = len(A)
        
        c_extended = c[:]
        tableau = []
        basis = []
        artificial_vars = []
        var_index = num_variables
        
        for i in range(num_constraints):
            row = A[i][:]
            slack_surplus = [0]*num_constraints
            artificial = [0]*num_constraints
            if signs[i] == "<=":
                # Agregar variable de holgura
                slack_surplus[i] = 1
                c_extended.append(0)
                basis.append(var_index)
                var_index += 1
            elif signs[i] == ">=":
                # Agregar variable de superávit y variable artificial
                slack_surplus[i] = -1
                c_extended.append(0)
                artificial[i] = 1
                c_extended.append(M)
                artificial_vars.append(var_index + 1)
                basis.append(var_index + 1)
                var_index += 2
            elif signs[i] == "=":
                # Agregar variable artificial
                artificial[i] = 1
                c_extended.append(M)
                artificial_vars.append(var_index)
                basis.append(var_index)
                var_index += 1
            else:
                self.output_text.insert(tk.END, "Signo de restricción no reconocido.\n")
                return
            row.extend(slack_surplus)
            row.extend(artificial)
            row.append(b[i])
            tableau.append(row)
        
        # Ajustar longitud de c_extended
        while len(c_extended) < len(tableau[0])-1:
            c_extended.append(0)
        
        # Ajustar los coeficientes de la función objetivo para las variables artificiales
        for idx in artificial_vars:
            if self.optimization_type.get() == "Minimizar":
                c_extended[idx] = -M
            else:
                c_extended[idx] = -M  # Cambiar a -M para maximizar
        
        # Mostrar el tableau inicial
        self.output_text.insert(tk.END, "Tabla Inicial (Método de la M Grande):\n")
        self.display_tableau(tableau, c_extended, basis)
        
        # Iniciar iteraciones
        iteration = 0
        max_iterations = 100
        while iteration < max_iterations:
            iteration += 1
            self.output_text.insert(tk.END, f"\nIteración {iteration}:\n")
            
            # Calcular zj y cj - zj
            zj = [0]*(len(tableau[0])-1)
            for i in range(len(basis)):
                for j in range(len(zj)):
                    zj[j] += c_extended[basis[i]] * tableau[i][j]
            cj_zj = [c_extended[j] - zj[j] for j in range(len(zj))]
            
            # Mostrar el tableau con las ganancias relativas
            self.display_tableau(tableau, c_extended, basis, cj_zj=cj_zj)
            
            # Verificar optimalidad
            if all(value <= 1e-5 for value in cj_zj):
                self.output_text.insert(tk.END, "Se encontró la solución óptima.\n")
                break
            
            # Variable entrante
            entering_candidates = [j for j in range(len(cj_zj)) if cj_zj[j] > 1e-5]
            if not entering_candidates:
                self.output_text.insert(tk.END, "No hay variables entrantes. El problema puede tener soluciones múltiples.\n")
                break
            entering = entering_candidates[0]  # Puede usarse alguna estrategia
            self.output_text.insert(tk.END, f"Variable entrante: x{entering+1}\n")
            
            # Calcular razones
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][entering] > 1e-5:
                    ratios.append(tableau[i][-1] / tableau[i][entering])
                else:
                    ratios.append(float('inf'))
            
            # Mostrar las razones en el tableau
            self.display_tableau(tableau, c_extended, basis, ratios=ratios, entering=entering)
            
            # Verificar ilimitado
            if all(r == float('inf') for r in ratios):
                self.output_text.insert(tk.END, "La solución es ilimitada.\n")
                return
            
            # Variable saliente
            leaving = ratios.index(min(ratios))
            self.output_text.insert(tk.END, f"Variable saliente: x{basis[leaving]+1}\n")
            
            # Pivoteo
            pivot_element = tableau[leaving][entering]
            tableau[leaving] = [x / pivot_element for x in tableau[leaving]]
            
            for i in range(len(tableau)):
                if i != leaving:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][j] - factor * tableau[leaving][j] for j in range(len(tableau[i]))]
            
            # Actualizar base
            basis[leaving] = entering
        else:
            self.output_text.insert(tk.END, "Se alcanzó el número máximo de iteraciones.\n")
            return
        
        # Verificar si hay variables artificiales en la base
        if any(var in basis for var in artificial_vars):
            self.output_text.insert(tk.END, "El problema no tiene solución factible.\n")
            return
        
        # Extraer solución
        solution = [0]*(len(tableau[0])-1)
        for i in range(len(basis)):
            solution[basis[i]] = tableau[i][-1]
        
        z = sum(c_extended[i]*solution[i] for i in range(len(c_extended)))
        if self.optimization_type.get() == "Minimizar":
            z = -z
        self.output_text.insert(tk.END, f"\nSolución Óptima:\n")
        for i in range(num_variables):
            self.output_text.insert(tk.END, f"x{i+1} = {solution[i]:.2f}\n")
        self.output_text.insert(tk.END, f"Z = {z:.2f}\n")
        
    def two_phase_method(self, c, A, b, signs):
        self.output_text.delete(1.0, tk.END)
        
        num_variables = len(c)
        num_constraints = len(A)
        
        c_extended = c[:]
        tableau = []
        basis = []
        artificial_vars = []
        var_index = num_variables
        
        # Fase 1: Construir problema auxiliar
        c_phase1 = [0]*(num_variables)
        for i in range(num_constraints):
            row = A[i][:]
            slack_surplus = [0]*num_constraints
            artificial = [0]*num_constraints
            if signs[i] == "<=":
                slack_surplus[i] = 1
                c_phase1.append(0)
                basis.append(var_index)
                var_index += 1
            elif signs[i] == ">=":
                slack_surplus[i] = -1
                c_phase1.append(0)
                artificial[i] = 1
                c_phase1.append(1)
                artificial_vars.append(var_index + 1)
                basis.append(var_index + 1)
                var_index += 2
            elif signs[i] == "=":
                artificial[i] = 1
                c_phase1.append(1)
                artificial_vars.append(var_index)
                basis.append(var_index)
                var_index += 1
            else:
                self.output_text.insert(tk.END, "Signo de restricción no reconocido.\n")
                return
            row.extend(slack_surplus)
            row.extend(artificial)
            row.append(b[i])
            tableau.append(row)
        
        # Ajustar longitud de c_phase1
        while len(c_phase1) < len(tableau[0])-1:
            c_phase1.append(0)
        
        # Mostrar el tableau inicial de la Fase 1
        self.output_text.insert(tk.END, "Fase 1: Tabla Inicial\n")
        self.display_tableau(tableau, c_phase1, basis)
        
        # Iniciar iteraciones de la Fase 1
        iteration = 0
        max_iterations = 100
        while iteration < max_iterations:
            iteration += 1
            self.output_text.insert(tk.END, f"\nFase 1 - Iteración {iteration}:\n")
            
            # Calcular zj y cj - zj
            zj = [0]*(len(tableau[0])-1)
            for i in range(len(basis)):
                for j in range(len(zj)):
                    zj[j] += c_phase1[basis[i]] * tableau[i][j]
            cj_zj = [c_phase1[j] - zj[j] for j in range(len(zj))]
            
            # Mostrar el tableau con las ganancias relativas
            self.display_tableau(tableau, c_phase1, basis, cj_zj=cj_zj)
            
            # Verificar optimalidad
            if all(value >= -1e-5 for value in cj_zj):
                self.output_text.insert(tk.END, "Fase 1 completada.\n")
                break
            
            # Variable entrante
            entering_candidates = [j for j in range(len(cj_zj)) if cj_zj[j] < -1e-5]
            if not entering_candidates:
                self.output_text.insert(tk.END, "No hay variables entrantes en Fase 1.\n")
                break
            entering = entering_candidates[0]
            self.output_text.insert(tk.END, f"Variable entrante: x{entering+1}\n")
            
            # Calcular razones
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][entering] > 1e-5:
                    ratios.append(tableau[i][-1] / tableau[i][entering])
                else:
                    ratios.append(float('inf'))
            
            # Mostrar las razones en el tableau
            self.display_tableau(tableau, c_phase1, basis, ratios=ratios, entering=entering)
            
            # Verificar ilimitado
            if all(r == float('inf') for r in ratios):
                self.output_text.insert(tk.END, "El problema no tiene solución factible.\n")
                return
            
            # Variable saliente
            leaving = ratios.index(min(ratios))
            self.output_text.insert(tk.END, f"Variable saliente: x{basis[leaving]+1}\n")
            
            # Pivoteo
            pivot_element = tableau[leaving][entering]
            tableau[leaving] = [x / pivot_element for x in tableau[leaving]]
            
            for i in range(len(tableau)):
                if i != leaving:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][j] - factor * tableau[leaving][j] for j in range(len(tableau[i]))]
            
            # Actualizar base
            basis[leaving] = entering
        else:
            self.output_text.insert(tk.END, "Se alcanzó el número máximo de iteraciones en Fase 1.\n")
            return
        
        # Verificar si valor óptimo es cero
        z0 = sum(c_phase1[basis[i]] * tableau[i][-1] for i in range(len(basis)))
        if abs(z0) > 1e-5:
            self.output_text.insert(tk.END, "El problema no tiene solución factible.\n")
            return
        
        # Eliminar columnas de variables artificiales
        artificial_vars = sorted(artificial_vars, reverse=True)
        for var in artificial_vars:
            for row in tableau:
                del row[var]
            del c_extended[var]
        
        # Actualizar base
        basis = [b for b in basis if b not in artificial_vars]
        
        # Ajustar c_extended
        while len(c_extended) < len(tableau[0])-1:
            c_extended.append(0)
        
        # Fase 2: Resolver el problema original
        self.output_text.insert(tk.END, "\nFase 2: Resolver el problema original\n")
        self.display_tableau(tableau, c_extended, basis)
        
        # Iniciar iteraciones de la Fase 2
        iteration = 0
        while iteration < max_iterations:
            iteration += 1
            self.output_text.insert(tk.END, f"\nFase 2 - Iteración {iteration}:\n")
            
            # Calcular zj y cj - zj
            zj = [0]*(len(tableau[0])-1)
            for i in range(len(basis)):
                for j in range(len(zj)):
                    zj[j] += c_extended[basis[i]] * tableau[i][j]
            cj_zj = [c_extended[j] - zj[j] for j in range(len(zj))]
            
            # Mostrar el tableau con las ganancias relativas
            self.display_tableau(tableau, c_extended, basis, cj_zj=cj_zj)
            
            # Verificar optimalidad
            if all(value <= 1e-5 for value in cj_zj):
                self.output_text.insert(tk.END, "Se encontró la solución óptima.\n")
                break
            
            # Variable entrante
            entering_candidates = [j for j in range(len(cj_zj)) if cj_zj[j] > 1e-5]
            if not entering_candidates:
                self.output_text.insert(tk.END, "No hay variables entrantes. El problema puede tener soluciones múltiples.\n")
                break
            entering = entering_candidates[0]
            self.output_text.insert(tk.END, f"Variable entrante: x{entering+1}\n")
            
            # Calcular razones
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][entering] > 1e-5:
                    ratios.append(tableau[i][-1] / tableau[i][entering])
                else:
                    ratios.append(float('inf'))
            
            # Mostrar las razones en el tableau
            self.display_tableau(tableau, c_extended, basis, ratios=ratios, entering=entering)
            
            # Verificar ilimitado
            if all(r == float('inf') for r in ratios):
                self.output_text.insert(tk.END, "La solución es ilimitada.\n")
                return
            
            # Variable saliente
            leaving = ratios.index(min(ratios))
            self.output_text.insert(tk.END, f"Variable saliente: x{basis[leaving]+1}\n")
            
            # Pivoteo
            pivot_element = tableau[leaving][entering]
            tableau[leaving] = [x / pivot_element for x in tableau[leaving]]
            
            for i in range(len(tableau)):
                if i != leaving:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][j] - factor * tableau[leaving][j] for j in range(len(tableau[i]))]
            
            # Actualizar base
            basis[leaving] = entering
        else:
            self.output_text.insert(tk.END, "Se alcanzó el número máximo de iteraciones en Fase 2.\n")
            return
        
        # Extraer solución
        solution = [0]*(len(tableau[0])-1)
        for i in range(len(basis)):
            solution[basis[i]] = tableau[i][-1]
        
        z = sum(c_extended[i]*solution[i] for i in range(len(c_extended)))
        if self.optimization_type.get() == "Minimizar":
            z = -z
        self.output_text.insert(tk.END, f"\nSolución Óptima:\n")
        for i in range(num_variables):
            self.output_text.insert(tk.END, f"x{i+1} = {solution[i]:.2f}\n")
        self.output_text.insert(tk.END, f"Z = {z:.2f}\n")
 


if __name__ == "__main__":
    root = tk.Tk()
    app = LinearProgrammingSolver(root)
    root.mainloop()
