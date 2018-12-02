# matvev

Matrice vector product for tensor train tensor

## Instructions:

### Compilation
Type make.

### Creating a ttmat:
Here is an example for creating a 3-dimensional TT-matrix of sizes (10x15), (12x12), (14x8) and ranks (1, 5, 8, 1).

```
create-ttmat -f ttmat.bin -d 3 -m 10,12,14 -n 15,12,8 -r 5,8

```

### Creating a ttvec:
Here is an example for creating a 3-dimensional TT-vector of sizes (15, 12, 8) and ranks (1, 4, 2, 1).

```
create-ttvec -f ttvec.bin -d 3 -m 15,12,8 -r 4,2

```

### Running the skeleton code
After having craeted a TT-matrix and a TT-vector in a file (say "mata.bin" and "vecx.bin"), you can run the skeleton code that computes the matrix-vector multiplication and saves it in the file "vecy.bin" as follows

```
./ttmatvec -a mata.bin -x vecx.bin -y vecy.bin
```
