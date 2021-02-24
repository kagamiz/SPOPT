SPOPT
====

SPOPT is a **S**parse **P**olynomial **OPT**imizing package written in C++.

---

SPOPT numerically gains a lower bound for the polynomial optimization programs

```
minimize        f(x)
subject to      x in K
```

where K is a semi-algebraic set, by using Lasserre Hierarchy [1].

---

### Requirements and Dependencies

SPOPT uses [CMake](http://www.cmake.org) to build the entire package. Please download and install [CMake](http://www.cmake.org).

Also, following C++ libraries are used in the source code. Please install them, and place the header file (or, its symbolic link) to `ext/include/` and the static library files to `ext/lib/` (or any other place where your C++ compilier can recognize).

- [Eigen](http://eigen.tuxfamily.org/)
- [yaml-cpp (ver 0.6.3)](https://github.com/jbeder/yaml-cpp)
- [Boost (ver 1.74.0)](https://www.boost.org/)
    - SPOPT only uses `boost::program_options`.
- [AA](https://github.com/cvxgrp/aa/tree/645cd5dd3970020ad78cd2837260725dbd433e23)
- [MOSEK fusion API](https://docs.mosek.com/9.2/cxxfusion/api-reference.html) (Optional)
    - Change the parameter `MOSEK_EXTENSION` in `CMakeLists.txt` to `true` to use MOSEK to solve SDPs. You need to obtain an appropriate license to use MOSEK.
- [SCS](https://github.com/cvxgrp/scs/tree/681b049f8f942c6f1a81014deaef92e5c869f354) (Optional)
    - Change the parameter `SCS_EXTENSION` in `CMakeLists.txt` to `true` to use SCS to solve SDPs.
    Add the directory `include` as `scs` under the directory `ext/include`, and also add the directory `linsys` under the directory `ext/include`.
    Also, add the files made by `make` command to `ext/lib/scs`, i.e., add an symbolic link from
    `(scs_install_dir_on_your_file_system)/out` to `ext/lib/scs`.
- BLAS
- LAPACK

---

### Compilation

By typing following commands, `cmake` will automatically generate `Makefile` and `make` command makes it able to compile the code and create SPOPT executable program in the `build` folder.

```sh
cd build
cmake ..
make
```

---

### Basic Usage

See the file `problem_data.yml` and `solver.yml` in `example` directory.
If you have correctly compiled SPOPT, then by running the following command,

```sh
(SPOPT_install_dir)/build/SPOPT -p (SPOPT_install_dir)/example/problem_data.yml -s (SPOPT_install_dir)/example/solver.yml
```

SPOPT will solve the following POP problem :

```
minimize        x_1 + x_2 + x_3 + x_4
subject to      x_1^2 + x_2^2 + x_3^2 + x_4^2 = 1.
```

It should show the optimal value `-0.5` and its corresponding solution vector.

Type

```sh
(SPOPT_install_dir)/build/SPOPT --help
```

to obtain further information.


---

### References

- [1] J. B. Lasserre. Global optimization with polynomials and the problem of moments. *SIAM Journal on Optimization*, 11 (2001), pp. 796-817.