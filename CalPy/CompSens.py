import os
import shutil
import csv
import argparse
def CompSens(InputFileName, penalty, radius, delta=0.0001):
    """
    Perform sensitivity analysis and optimization using a specified input file.
    
    This function automates the process of running a simulation with a given input file, 
    modifying the density values in small increments (controlled by `delta`), 
    re-running the simulation, and calculating the sensitivity of the system 
    by comparing the resulting compliance values.

    Parameters:
    ----------
    InputFileName : str
        The name (without the extension) of the input file used in the simulation. 
        For example, if the input file is "T.inp", pass "T" as the argument.

    penalty : float
        A penalty value used in the simulation. Default is 2. (Not currently used in this function.)

    radius : float
        The radius value used in the simulation command. Default is 0.001.
    
    delata : float, optional
        The stepsize relative to the current density in a particular cell to compute compliance gradient         using FD scheme. The default is 0.0001
    Returns:
    -------
    None
        This function does not return any value but writes results to a `result.csv` file.
        The file contains the sensitivity analysis results, with columns for the adjoint values,
        finite difference values (FD), and the difference between them.

        argument can be parsed using argparse: For instance you could do python3 CompSens.py T --p 2 --r         0.01 --delta 0.00001  where T.inp is the input file of interest.

        """
    # Full name for the input file
    fullname = InputFileName + '.inp'
    
    # 1. Read the density.dat and store it in an array called density
    os.system(f"calTop.exe {InputFileName} -p {penalty} -r {radius}")  # Run the simulation
    with open('density.dat', 'r') as file:
        density = [float(line.strip()) for line in file.readlines()]

    # 2. Make a new directory called Sensitivity and copy over the input file and density.dat
    if not os.path.exists('Sensitivity'):
        os.makedirs('Sensitivity')

    # Copy the dynamically passed input file to the new directory
    shutil.copy(f'{InputFileName}.inp', 'Sensitivity/')  # Copy the input file
    shutil.copy('density.dat', 'Sensitivity/density.dat')  # Copy density.dat to the new directory
    
    # Change to the Sensitivity directory
    os.chdir('Sensitivity')
    
    # 3. Run "CalTop.exe T -p 2" in the Sensitivity directory
    os.system(f"calTop.exe {InputFileName} -p {penalty} -r {radius}")
    
    # 4. Read the first column of compliance_sens.csv and store it in an array called adjoint
    adjoint = []
    with open('compliance_sens.csv', mode='r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip the header row
        for row in reader:
            adjoint.append(float(row[0]))  # Store the first column in adjoint

    # 5. Read the objective.csv and store the first value in the second row as comp_old
    with open('objectives.csv', mode='r') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
        comp_old = float(rows[1][0])  # Store the first value of the second row as comp_old

    # 6. Iterate through the density array, change values to 0.99, re-run the simulation, and compute FD
    FD = []
    for i in range(len(density)):
        # Change the i-th density value to 0.99
        modified_density = density.copy()
        modified_density[i] = modified_density[i] - delta * density[i]

        # Write the modified density array back to the density.dat file in the Sensitivity directory
        with open('density.dat', 'w') as file:
            for value in modified_density:
                file.write(f"{float(value)}\n")

        # Re-run "CalTop.exe T -p 2" in the Sensitivity directory after modifying the density
        os.system(f"calTop.exe {InputFileName} -p {penalty} -r {radius}")

        # Read objective.csv again to get the new compliance value (comp_new)
        with open('objectives.csv', mode='r') as csvfile:
            reader = csv.reader(csvfile)
            rows = list(reader)
            comp_new = float(rows[1][0])  # Store the first value of the second row as comp_new

        # Calculate FD and append to the FD array
        FD.append(-(comp_new - comp_old) / (delta * density[i]))
        print(f"Iteration {i + 1} completed")

    # 7. Write adjoint, FD, and difference to a result.csv file
    with open('result.csv', mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write header
        writer.writerow(['Adjoint', 'FD', 'Difference'])

        # Write each row: adjoint, FD, and the difference (FD - adjoint)
        for a, f in zip(adjoint, FD):
            difference = f - a  # Calculate the difference
            writer.writerow([a, f, difference])  # Write the values to the CSV file

    print("Results written to result.csv")
def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run sensitivity analysis and optimization.")
    
    # Define command-line arguments
    parser.add_argument('InputFileName', type=str, help="Name of the input file (without extension).")
    parser.add_argument('--delta', type=float, default=0.0001, help="Small change in density for FD calculation.")
    parser.add_argument('--p', type=float, default=2, help="Penalty value used in the simulation.")
    parser.add_argument('--r', type=float, default=0.001, help="Radius value used in the simulation.")
    
    # Parse the arguments
    args = parser.parse_args()

    # Call the function with the parsed arguments
    CompSens(args.InputFileName, penalty=args.p, radius=args.r, delta=args.delta)

if __name__ == "__main__":
    main()
