import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

class calc_simplex_mBig_dosFases:
    def __init__(self, root):
        self.root = root
        self.root.title("Calculadora de Programación Lineal")

        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        self.num_variables = 2
        self.numero_restricciones = 3
        self.method = tk.StringVar(value="Método Simplex")
        self.optimization_type = tk.StringVar(value="Maximizar")

        self.interfaz()

    def interfaz(self):

        interfaz_princpial = ttk.Frame(self.root)
        interfaz_princpial.grid(row=0, column=0, sticky="nsew")
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        interfaz_princpial.columnconfigure(0, weight=1)
        interfaz_princpial.columnconfigure(1, weight=2)  
        interfaz_princpial.rowconfigure(6, weight=1)

        metodo_interfaz = ttk.LabelFrame(interfaz_princpial, text="Método: ")
        metodo_interfaz.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        metodo_interfaz.columnconfigure(0, weight=1)

        method_combobox = ttk.Combobox(
            metodo_interfaz,
            textvariable=self.method,
            values=["Método Simplex", "Método de la M grande", "Método de las dos fases"],
            state="readonly",
        )
        method_combobox.grid(row=0, column=0, padx=10, pady=5, sticky="ew")
        method_combobox.current(0)  

        opt_frame = ttk.LabelFrame(interfaz_princpial, text="Tipo de Optimización")
        opt_frame.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        opt_frame.columnconfigure(0, weight=1)
        opt_frame.columnconfigure(1, weight=1)

        ttk.Radiobutton(
            opt_frame, text="Maximizar", variable=self.optimization_type, value="Maximizar"
        ).grid(row=0, column=0, padx=5, pady=5, sticky="w")
        ttk.Radiobutton(
            opt_frame, text="Minimizar", variable=self.optimization_type, value="Minimizar"
        ).grid(row=0, column=1, padx=5, pady=5, sticky="w")

        obj_frame = ttk.LabelFrame(interfaz_princpial, text="Función Objetivo (Z)")
        obj_frame.grid(row=2, column=0, padx=10, pady=10, sticky="ew")
        obj_frame.columnconfigure('all', weight=1)

        self.obj_coeffs = []
        for i in range(self.num_variables):
            entry = ttk.Entry(obj_frame, width=5)
            entry.grid(row=0, column=i * 2, padx=5, pady=5, sticky="ew")
            lbl = ttk.Label(obj_frame, text=f"x{i + 1}")
            lbl.grid(row=0, column=i * 2 + 1, padx=5, pady=5, sticky="w")
            self.obj_coeffs.append(entry)

        constr_frame = ttk.LabelFrame(interfaz_princpial, text="Restricciones")
        constr_frame.grid(row=3, column=0, padx=10, pady=10, sticky="ew")
        constr_frame.columnconfigure('all', weight=1)

        self.restricciones_coeffs = []
        self.restricciones_signs = []
        self.restricciones_rhs = []

        for i in range(self.numero_restricciones):
            constr_frame.rowconfigure(i, weight=1)
            row_entries = []
            for j in range(self.num_variables):
                entry = ttk.Entry(constr_frame, width=5)
                entry.grid(row=i, column=j * 2, padx=5, pady=5, sticky="ew")
                lbl = ttk.Label(constr_frame, text=f"x{j + 1}")
                lbl.grid(row=i, column=j * 2 + 1, padx=5, pady=5, sticky="w")
                row_entries.append(entry)
            self.restricciones_coeffs.append(row_entries)

            sign_cb = ttk.Combobox(constr_frame, values=["<=", ">=", "="], width=5, state="readonly")
            sign_cb.grid(row=i, column=self.num_variables * 2, padx=5, pady=5, sticky="ew")
            sign_cb.current(0)
            self.restricciones_signs.append(sign_cb)

            rhs_entry = ttk.Entry(constr_frame, width=5)
            rhs_entry.grid(row=i, column=self.num_variables * 2 + 1, padx=5, pady=5, sticky="ew")
            self.restricciones_rhs.append(rhs_entry)

        var_restricciones_frame = ttk.LabelFrame(
            interfaz_princpial, text="Restricciones de Variables (activadas por defecto)"
        )
        var_restricciones_frame.grid(row=4, column=0, padx=10, pady=10, sticky="ew")
        var_restricciones_frame.columnconfigure('all', weight=1)

        for i in range(self.num_variables):
            lbl = ttk.Label(var_restricciones_frame, text=f"x{i + 1} ≥ 0")
            lbl.grid(row=0, column=i, padx=5, pady=5, sticky="w")

        solucion_button = ttk.Button(interfaz_princpial, text="Solucionar", command=self.solucion)
        solucion_button.grid(row=5, column=0, pady=20, sticky="ew")

        output_frame = ttk.LabelFrame(interfaz_princpial, text="Solución")
        output_frame.grid(row=0, column=1, rowspan=6, padx=10, pady=10, sticky="nsew")
        output_frame.columnconfigure(0, weight=1)
        output_frame.rowconfigure(0, weight=1)

        self.output_text = tk.Text(output_frame, wrap="word", height=10)
        self.output_text.grid(row=0, column=0, sticky="nsew")

        scrollbar = ttk.Scrollbar(output_frame, orient="vertical", command=self.output_text.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.output_text.configure(yscrollcommand=scrollbar.set)


    def solucion(self):
        try:
            c = [float(entry.get()) for entry in self.obj_coeffs]
            if self.optimization_type.get() == "Minimizar":
                c = [-coeff for coeff in c]
            A = []
            b = []
            signs = []
            for i in range(self.numero_restricciones):
                row = [float(entry.get()) for entry in self.restricciones_coeffs[i]]
                A.append(row)
                b.append(float(self.restricciones_rhs[i].get()))
                signs.append(self.restricciones_signs[i].get())
        except ValueError:
            messagebox.showerror("Error de Entrada", "Por favor, ingrese valores numéricos válidos.")
            return

        if self.method.get() == "Método Simplex":
            self.simplex(c, A, b, signs)
        elif self.method.get() == "Método de la M grande":
            self.m_grande(c, A, b, signs)
        elif self.method.get() == "Método de las dos fases":
            self.dos_fases(c, A, b, signs)
        else:
            messagebox.showerror("Error de Método", "Por favor, seleccione un método válido.")
            return
        
    def simplex(self, c, A, b, signs):
        self.output_text.delete(1.0, tk.END)
        
        if '>=' in signs or '=' in signs:
            self.output_text.insert(tk.END, "El Método Simplex estándar no puede resolucionr restricciones '≥' o '='.\n")
            return

        num_variables = len(c)
        numero_restricciones = len(A)
        
        c_extended = c[:]
        tableau = []
        basis = []
        slack_var_index = num_variables
        
        variable_names = [f"X{i+1}" for i in range(num_variables)]
        
        for i in range(numero_restricciones):
            row = A[i][:]
            slack = [0] * numero_restricciones
            slack[i] = 1
            c_extended.append(0)
            variable_names.append(f"S{slack_var_index - num_variables +1}")
            basis.append(slack_var_index)
            slack_var_index +=1
            row.extend(slack)
            row.append(b[i])
            tableau.append(row)

        self.output_text.insert(tk.END, "Tabla Inicial:\n")
        self.mostrar_tabla(tableau, c_extended, basis, variable_names)

        iteracion = 0
        max_iteracions = 100 
        while iteracion < max_iteracions:
            iteracion +=1
            self.output_text.insert(tk.END, f"\nIteración {iteracion}:\n")

            zj = [0]*len(c_extended)
            for i in range(len(basis)):
                for j in range(len(c_extended)):
                    zj[j] += c_extended[basis[i]] * tableau[i][j]
            cj_zj = [c_extended[j] - zj[j] for j in range(len(c_extended))]
            
            self.mostrar_tabla(tableau, c_extended, basis, variable_names, cj_zj=cj_zj)
            
            if all(value <= 1e-5 for value in cj_zj):
                self.output_text.insert(tk.END, "Se encontró la solución óptima.\n")
                break
            
            entering_candidates = [j for j in range(len(cj_zj)) if cj_zj[j] > 1e-5]
            if not entering_candidates:
                self.output_text.insert(tk.END, "No hay variables entrantes. El problema puede tener soluciones múltiples.\n")
                break
            entering = max(entering_candidates, key=lambda j: cj_zj[j])
            self.output_text.insert(tk.END, f"Variable entrante: {variable_names[entering]}\n")
            
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][entering] > 1e-5:
                    ratios.append(tableau[i][-1] / tableau[i][entering])
                else:
                    ratios.append(float('inf'))

            self.mostrar_tabla(tableau, c_extended, basis, variable_names, ratios=ratios, entering=entering)
            
            if all(r == float('inf') for r in ratios):
                self.output_text.insert(tk.END, "La solución es ilimitada.\n")
                return
            
            leaving = ratios.index(min(ratios))
            self.output_text.insert(tk.END, f"Variable saliente: {variable_names[basis[leaving]]}\n")

            pivot_element = tableau[leaving][entering]
            tableau[leaving] = [x / pivot_element for x in tableau[leaving]]
            
            for i in range(len(tableau)):
                if i != leaving:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][j] - factor * tableau[leaving][j] for j in range(len(tableau[i]))]
            
            basis[leaving] = entering
        else:
            self.output_text.insert(tk.END, "Se alcanzó el número máximo de iteraciones. El problema puede no tener solución óptima.\n")
            return
        
        solution = [0]*len(c_extended)
        for i in range(len(basis)):
            if basis[i] < len(solution):
                solution[basis[i]] = tableau[i][-1]
        
        z = sum(c_extended[i]*solution[i] for i in range(len(c_extended)))
        if self.optimization_type.get() == "Minimizar":
            z = -z
        self.output_text.insert(tk.END, f"\nSolución Óptima:\n")
        for i in range(num_variables):
            self.output_text.insert(tk.END, f"{variable_names[i]} = {solution[i]:.2f}\n")
        self.output_text.insert(tk.END, f"Z = {z:.2f}\n")
        
    def mostrar_tabla(self, tableau, c_extended, basis, variable_names, cj_zj=None, ratios=None, entering=None):
        num_vars = len(c_extended)
        header = "+------+"
        header += "-----------" * num_vars + "+"
        header += "----------+\n"
        
        title_row = "| Base |"
        for var_name in variable_names:
            title_row += f"   {var_name:<5}|"
        title_row += " Solución |\n"
        
        header += title_row
        header += "+------+"
        header += "-----------" * num_vars + "+"
        header += "----------+\n"
        
        self.output_text.insert(tk.END, header)

        for i in range(len(tableau)):
            row = tableau[i]
            base_var = basis[i]
            base_var_name = variable_names[base_var]
            row_str = f"| {base_var_name:<4}|"
            for val in row[:-1]:
                row_str += f" {val:>8.2f} |"

            if ratios:
                if ratios[i] != float('inf'):
                    ratio_str = f"{row[-1]:.2f}/{row[entering]:.2f}={ratios[i]:.2f}"
                else:
                    ratio_str = "Inf"
                row_str += f" {ratio_str} |"
            row_str += "\n"
            self.output_text.insert(tk.END, row_str)
            self.output_text.insert(tk.END, "+------+")
            self.output_text.insert(tk.END, "-----------" * num_vars + "+")
            self.output_text.insert(tk.END, "----------+\n")
        
        if cj_zj:
            z_value = sum(c_extended[basis[i]] * tableau[i][-1] for i in range(len(basis)))
            row_str = f"|   Z  |"
            for val in cj_zj:
                row_str += f" {val:>8.2f} |"
            row_str += f" {z_value:>8.2f} |\n"
            self.output_text.insert(tk.END, row_str)
            self.output_text.insert(tk.END, "+------+")
            self.output_text.insert(tk.END, "-----------" * num_vars + "+")
            self.output_text.insert(tk.END, "----------+\n")
    
    def m_grande(self, c, A, b, signs):
        self.output_text.delete(1.0, tk.END)
        
        M = 1e6  
        num_variables = len(c)
        numero_restricciones = len(A)
        
        c_extended = c[:]
        tableau = []
        basis = []
        varaibles_artificiales = []
        var_index = num_variables
        
        variable_names = [f"X{i+1}" for i in range(num_variables)]
        
        for i in range(numero_restricciones):
            row = A[i][:]
            slack_surplus = [0]*numero_restricciones
            artificial = [0]*numero_restricciones
            if signs[i] == "<=":
        
                slack_surplus[i] = 1
                c_extended.append(0)
                variable_names.append(f"S{var_index - num_variables +1}")
                basis.append(var_index)
                var_index += 1
            elif signs[i] == ">=":

                slack_surplus[i] = -1
                c_extended.append(0)
                variable_names.append(f"S{var_index - num_variables +1}")
                var_index += 1
                artificial[i] = 1
                c_extended.append(-M)
                variable_names.append(f"A{var_index - num_variables +1}")
                varaibles_artificiales.append(var_index)
                basis.append(var_index)
                var_index += 1
            elif signs[i] == "=":

                artificial[i] = 1
                c_extended.append(-M)
                variable_names.append(f"A{var_index - num_variables +1}")
                varaibles_artificiales.append(var_index)
                basis.append(var_index)
                var_index += 1
            else:
                self.output_text.insert(tk.END, "Signo de restricción no reconocido.\n")
                return
            row.extend(slack_surplus)
            row.extend(artificial)
            row.append(b[i])
            tableau.append(row)
        
        while len(c_extended) < len(tableau[0])-1:
            c_extended.append(0)
        
        self.output_text.insert(tk.END, "Tabla Inicial (Método de la M Grande):\n")
        self.mostrar_tabla(tableau, c_extended, basis, variable_names)
        
        iteracion = 0
        max_iteracions = 100
        while iteracion < max_iteracions:
            iteracion += 1
            self.output_text.insert(tk.END, f"\nIteración {iteracion}:\n")
            
            zj = [0]*(len(tableau[0])-1)
            for i in range(len(basis)):
                for j in range(len(zj)):
                    zj[j] += c_extended[basis[i]] * tableau[i][j]
            cj_zj = [c_extended[j] - zj[j] for j in range(len(zj))]
            
            self.mostrar_tabla(tableau, c_extended, basis, variable_names, cj_zj=cj_zj)
            
            if all(value <= 1e-5 for value in cj_zj):
                self.output_text.insert(tk.END, "Se encontró la solución óptima.\n")
                break
            
            entering_candidates = [j for j in range(len(cj_zj)) if cj_zj[j] > 1e-5]
            if not entering_candidates:
                self.output_text.insert(tk.END, "No hay variables entrantes. El problema puede tener soluciones múltiples.\n")
                break
            entering = entering_candidates[0] 
            self.output_text.insert(tk.END, f"Variable entrante: {variable_names[entering]}\n")
            
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][entering] > 1e-5:
                    ratios.append(tableau[i][-1] / tableau[i][entering])
                else:
                    ratios.append(float('inf'))
            
            self.mostrar_tabla(tableau, c_extended, basis, variable_names, ratios=ratios, entering=entering)
            
            if all(r == float('inf') for r in ratios):
                self.output_text.insert(tk.END, "La solución es ilimitada.\n")
                return
            
            leaving = ratios.index(min(ratios))
            self.output_text.insert(tk.END, f"Variable saliente: {variable_names[basis[leaving]]}\n")
            
            pivot_element = tableau[leaving][entering]
            tableau[leaving] = [x / pivot_element for x in tableau[leaving]]
            
            for i in range(len(tableau)):
                if i != leaving:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][j] - factor * tableau[leaving][j] for j in range(len(tableau[i]))]
            
            basis[leaving] = entering
        else:
            self.output_text.insert(tk.END, "Se alcanzó el número máximo de iteraciones.\n")
            return

        if any(var in basis for var in varaibles_artificiales):
            self.output_text.insert(tk.END, "El problema no tiene solución factible.\n")
            return
        
        solution = [0]*(len(tableau[0])-1)
        for i in range(len(basis)):
            solution[basis[i]] = tableau[i][-1]
        
        z = sum(c_extended[i]*solution[i] for i in range(len(c_extended)))
        if self.optimization_type.get() == "Minimizar":
            z = -z
        self.output_text.insert(tk.END, f"\nSolución Óptima:\n")
        for i in range(num_variables):
            self.output_text.insert(tk.END, f"{variable_names[i]} = {solution[i]:.2f}\n")
        self.output_text.insert(tk.END, f"Z = {z:.2f}\n")
    
    def dos_fases(self, c, A, b, signs):
        self.output_text.delete(1.0, tk.END)
        
        num_variables = len(c)
        numero_restricciones = len(A)
        
        variable_names = [f"X{i+1}" for i in range(num_variables)]
        c_phase1 = [0]*(num_variables)
        tableau = []
        basis = []
        varaibles_artificiales = []
        var_index = num_variables
        
        for i in range(numero_restricciones):
            row = A[i][:]
            slack_surplus = [0]*numero_restricciones
            artificial = [0]*numero_restricciones
            if signs[i] == "<=":
                slack_surplus[i] = 1
                c_phase1.append(0)
                variable_names.append(f"S{var_index - num_variables + 1}")
                basis.append(var_index)
                var_index += 1
            elif signs[i] == ">=":
                slack_surplus[i] = -1
                c_phase1.append(0)
                variable_names.append(f"S{var_index - num_variables + 1}")
                var_index += 1
                artificial[i] = 1
                c_phase1.append(1)
                variable_names.append(f"A{var_index - num_variables + 1}")
                varaibles_artificiales.append(var_index)
                basis.append(var_index)
                var_index += 1
            elif signs[i] == "=":
                artificial[i] = 1
                c_phase1.append(1)
                variable_names.append(f"A{var_index - num_variables + 1}")
                varaibles_artificiales.append(var_index)
                basis.append(var_index)
                var_index += 1
            else:
                self.output_text.insert(tk.END, "Signo de restricción no reconocido.\n")
                return
            row.extend(slack_surplus)
            row.extend(artificial)
            row.append(b[i])
            tableau.append(row)
        
        while len(c_phase1) < len(tableau[0])-1:
            c_phase1.append(0)
        
        self.output_text.insert(tk.END, "Fase 1: Tabla Inicial\n")
        self.mostrar_tabla(tableau, c_phase1, basis, variable_names)
        
        iteracion = 0
        max_iteracions = 100
        while iteracion < max_iteracions:
            iteracion += 1
            self.output_text.insert(tk.END, f"\nFase 1 - Iteración {iteracion}:\n")
            
            zj = [0]*(len(tableau[0])-1)
            for i in range(len(basis)):
                for j in range(len(zj)):
                    zj[j] += c_phase1[basis[i]] * tableau[i][j]
            cj_zj = [c_phase1[j] - zj[j] for j in range(len(zj))]
            
            self.mostrar_tabla(tableau, c_phase1, basis, variable_names, cj_zj=cj_zj)
            
            if all(value >= -1e-5 for value in cj_zj):
                self.output_text.insert(tk.END, "Fase 1 completada.\n")
                break
            
            entering_candidates = [j for j in range(len(cj_zj)) if cj_zj[j] < -1e-5]
            if not entering_candidates:
                self.output_text.insert(tk.END, "No hay variables entrantes en Fase 1.\n")
                break
            entering = entering_candidates[0]
            self.output_text.insert(tk.END, f"Variable entrante: {variable_names[entering]}\n")
            
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][entering] > 1e-5:
                    ratios.append(tableau[i][-1] / tableau[i][entering])
                else:
                    ratios.append(float('inf'))

            self.mostrar_tabla(tableau, c_phase1, basis, variable_names, ratios=ratios, entering=entering)
            
            if all(r == float('inf') for r in ratios):
                self.output_text.insert(tk.END, "El problema no tiene solución factible.\n")
                return

            leaving = ratios.index(min(ratios))
            self.output_text.insert(tk.END, f"Variable saliente: {variable_names[basis[leaving]]}\n")
            
            pivot_element = tableau[leaving][entering]
            tableau[leaving] = [x / pivot_element for x in tableau[leaving]]
            
            for i in range(len(tableau)):
                if i != leaving:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][j] - factor * tableau[leaving][j] for j in range(len(tableau[i]))]

            basis[leaving] = entering
        else:
            self.output_text.insert(tk.END, "Se alcanzó el número máximo de iteraciones en Fase 1.\n")
            return

        z0 = sum(c_phase1[basis[i]] * tableau[i][-1] for i in range(len(basis)))
        if abs(z0) > 1e-5:
            self.output_text.insert(tk.END, "El problema no tiene solución factible.\n")
            return

        varaibles_artificiales = sorted(varaibles_artificiales)
        variables_to_keep = [i for i in range(len(variable_names)) if i not in varaibles_artificiales]
        variable_names = [variable_names[i] for i in variables_to_keep]

        basis = [b for b in basis if b not in varaibles_artificiales]
        old_to_new_indices = {old_idx: new_idx for new_idx, old_idx in enumerate(variables_to_keep)}
        basis = [old_to_new_indices[b] for b in basis]

        for i in range(len(tableau)):
            row = tableau[i]
            new_row = [row[old_idx] for old_idx in variables_to_keep] + [row[-1]]  
            tableau[i] = new_row

        c_extended_phase2 = [0]*len(variable_names)
        for idx, var_name in enumerate(variable_names):
            if var_name.startswith("X"):
                original_idx = int(var_name[1:]) - 1
                c_extended_phase2[idx] = c[original_idx]
   
        self.output_text.insert(tk.END, "\nFase 2: Solucionar el problema original\n")
        self.mostrar_tabla(tableau, c_extended_phase2, basis, variable_names)
        
        iteracion = 0
        while iteracion < max_iteracions:
            iteracion += 1
            self.output_text.insert(tk.END, f"\nFase 2 - Iteración {iteracion}:\n")
            
            zj = [0]*(len(tableau[0])-1)
            for i in range(len(basis)):
                for j in range(len(zj)):
                    zj[j] += c_extended_phase2[basis[i]] * tableau[i][j]
            cj_zj = [c_extended_phase2[j] - zj[j] for j in range(len(zj))]
            
            self.mostrar_tabla(tableau, c_extended_phase2, basis, variable_names, cj_zj=cj_zj)
            
            if all(value <= 1e-5 for value in cj_zj):
                self.output_text.insert(tk.END, "Se encontró la solución óptima.\n")
                break
            
            entering_candidates = [j for j in range(len(cj_zj)) if cj_zj[j] > 1e-5]
            if not entering_candidates:
                self.output_text.insert(tk.END, "No hay variables entrantes. El problema puede tener soluciones múltiples.\n")
                break
            entering = entering_candidates[0]
            self.output_text.insert(tk.END, f"Variable entrante: {variable_names[entering]}\n")
            
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][entering] > 1e-5:
                    ratios.append(tableau[i][-1] / tableau[i][entering])
                else:
                    ratios.append(float('inf'))
            
            self.mostrar_tabla(tableau, c_extended_phase2, basis, variable_names, ratios=ratios, entering=entering)
            
            if all(r == float('inf') for r in ratios):
                self.output_text.insert(tk.END, "La solución es ilimitada.\n")
                return
            
            leaving = ratios.index(min(ratios))
            self.output_text.insert(tk.END, f"Variable saliente: {variable_names[basis[leaving]]}\n")
            
            pivot_element = tableau[leaving][entering]
            tableau[leaving] = [x / pivot_element for x in tableau[leaving]]
            
            for i in range(len(tableau)):
                if i != leaving:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][j] - factor * tableau[leaving][j] for j in range(len(tableau[i]))]
            
            basis[leaving] = entering
        else:
            self.output_text.insert(tk.END, "Se alcanzó el número máximo de iteraciones en Fase 2.\n")
            return
        
        solution = [0]*(len(variable_names))
        for i in range(len(basis)):
            solution[basis[i]] = tableau[i][-1]
        
        z = sum(c_extended_phase2[i]*solution[i] for i in range(len(c_extended_phase2)))
        if self.optimization_type.get() == "Minimizar":
            z = -z
        self.output_text.insert(tk.END, f"\nSolución Óptima:\n")
        for idx, var_name in enumerate(variable_names):
            if var_name.startswith("X"):
                self.output_text.insert(tk.END, f"{var_name} = {solution[idx]:.2f}\n")
        self.output_text.insert(tk.END, f"Z = {z:.2f}\n")
        


if __name__ == "__main__":
    root = tk.Tk()
    app = calc_simplex_mBig_dosFases(root)
    root.mainloop()
