# capdmathieu
Computer-assisted proofs in Hill's equations

In the folder "scripts" there are Newton's interval method, bisection interval method and Krawczyk method for different cases of the Mathieu equation. The C++ code uses the library CAPD::DynSys for rigorous numerics.

In the folder "enclosures" there are the data enclosures for the different results of the Hill's equation. The data is organized in four columns, where each column represents one of the following values $a_{left}$, $a_{right}$, $b_{left}$ and $b_{right}$, of the enclosures. If there is a fifth column, corresponds to the type of boundary: type 1 if $0\in\boldsymbol P_{21}$; type 2 if $0\in\boldsymbol P_{12}$; and type 3 if both things happen.
