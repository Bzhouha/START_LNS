# START_LNS

## Installation
To install the solver, a MPI implementation and a fortran compiler are needed, PETSc has follwing configuration:
```
./configure --with-debugging=no --with-scalar-type=complex --with-fortran-bindings=1 --with-precision=double --download-slepc --download-hdf5
```
Once you have these prerequisites installed, you can follow these steps to install the software:
1. Clone the repository to your local machine.
2. Change directories to the project directory:'cd START_LNS'
3. Run the make command to build the software:'make'
4. Run the make install command to install the software:'make install'
5. This will install the software to the ~/bin directory. You can then run the software by typing the following command in a terminal:'START+LNS -f config.ini'

## Usage
## License

This project is licensed under the GNU Lesser General Public License (LGPL v3). You are free to use, modify, and distribute the code in this project, but you must comply with the terms of the LGPL license.

You can get the text of the LGPL license from the following link: https://www.gnu.org/licenses/lgpl-3.0.html .

