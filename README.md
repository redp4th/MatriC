# MatriC

MatriC is a simple matrix library written in C++ originally for course experiment.

# Usage
## Declaring a Vector of dimension $n$ with all entries initialized to 0.
```
Vector v(n);
```

## Declaring an elementary vector with $1$ in $j$th entry.
```
Vector ej =  elementary(n, j);
```

## Declaring a $m\cross n$ matrix.
```
Matrix matrix(m, n);
```

## Declaring an identity matrix.
```
Matrix I = eye(n);
```
