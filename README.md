# Linear Programming Solver

This application is an interactive calculator for solving linear programming problems using the Simplex, Big M, and Two-Phase methods. It is built with a graphical user interface (GUI) using `Tkinter`.

## Features

- **Solution Methods**:
  - Simplex
  - Big M Method
  - Two-Phase Method
- **Optimization**:
  - Maximization
  - Minimization
- **User Input**:
  - Coefficients for the objective function
  - Constraints and their signs (<=, >=, =)
  - Non-negative variable constraints by default
- **Output**:
  - Iteration table
  - Optimal solution, if available
  - Detection of unbounded or infeasible cases

## Requirements

- Python 3.x
- tkinter

## Installation

1. Clone or download this repository.
2. Ensure `tkinter` is installed:
3. Run the main file:
    ```bash
    python main.py
    ```

## Usage

1. Open the application and select the **solution method** from the first dropdown menu.
2. Define the **type of optimization** (Maximize or Minimize).
3. Enter the coefficients for the **objective function** and the **constraints** in the corresponding fields.
4. Click the `Solve` button to obtain the solution to the problem.

### Solution Methods

- **Simplex Method**: Used to solve maximization problems with `<=` constraints. If there are `>=` or `=` constraints, the method will display an error message.
- **Big M Method**: Solves problems with `>=` or `=` constraints by adding artificial variables.
- **Two-Phase Method**: Handles problems with `>=` or `=` constraints, working in two phases to find a feasible solution before optimization.

## Example

1. Select the **Simplex Method**.
2. Enter the following objective function to maximize:
   - `Z = 3x1 + 5x2`
3. Add the following constraints:
   - `x1 + 2x2 <= 10`
   - `x1 + x2 <= 6`
   - `x2 <= 4`
4. Press `Solve`. The solution will appear in the output column on the right.
