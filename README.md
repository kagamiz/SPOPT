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

Also, following C++ libraries are used in the source code. Please install them, and place the header file and the library file inside the search path used by your system.

- [Eigen](http://eigen.tuxfamily.org/)
- [yaml-cpp (ver 0.6.3)](https://github.com/jbeder/yaml-cpp)

---

### Compilation

Typing `cmake` at the command line will compile the code and create SPOPT executable program in the `bin` folder.

```sh
cmake
```

---

### Basic Usage

To be written...

---

### References

- [1] J. B. Lasserre. Global optimization with polynomials and the problem of moments. *SIAM Journal on Optimization*, 11 (2001), pp. 796-817.