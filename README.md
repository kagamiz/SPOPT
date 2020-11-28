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

### Requirements and dependencies

SPOPT uses [CMake](http://www.cmake.org) to build the entire package. Please download and install [CMake](http://www.cmake.org).

Also, following C++ libraries are used in the source code. Please install them, and place the header file (or, its symbolic link) to `ext/include/` and the static library files to `ext/lib/`.

- [Eigen](http://eigen.tuxfamily.org/)
- [yaml-cpp (ver 0.6.3)](https://github.com/jbeder/yaml-cpp)
- [Boost (ver 1.74.0)](https://www.boost.org/)
    - SPOPT only uses `boost::program_options`.

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

To be written...

---

### References

- [1] J. B. Lasserre. Global optimization with polynomials and the problem of moments. *SIAM Journal on Optimization*, 11 (2001), pp. 796-817.