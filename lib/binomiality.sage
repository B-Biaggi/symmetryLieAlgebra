def test_binomiality_random(flist, B, epsilon = 0):
    """
    Tests if there exists a g in GL_n(CC) such that the vanishing ideal of g  cdot overline{im(flist)} is binomial.
    
    Input: 
        flist (list): A list of rational functions over QQ[x_1,...,x_m]
        B (list of Matrices): Input the basis of the symmetriy Lie algebra given by symmetry_lie_algebra_random
        epsilon (rational, optional): Enter epsilon in (0,1) to have probablity as in the paper
        
    Output:
        boolean
    """
    d = max([max_degree(fi) for fi in flist])
    n = len(flist)
    m = len(flist[0].parent().gens())

    if 0 < epsilon < 1:
        small_entries = False
        eps = epsilon
    elif epsilon == 0:
        small_entries = True
        eps = 1
    else:
        raise ValueError("epsilon must be 0 or a number strictly between 0 and 1")

    M = n^2 / eps

    dim_lie_alg = len(B)
    A = sum( [ZZ.random_element(1,M+1) * B_mat for B_mat in B])
    As = semisimple_part_hensel_lift(A)
    AsDList = [As * b - b * As for b in B]
    AsDMat = matrix(QQ, [list(vector(AsD)) for AsD in AsDList]).transpose()
    sols = AsDMat.right_kernel()
    C = commutator_kernel_basis(As, B)
    Cs = []

    for D in C:
        # print("Computing the semsisimple part of elements in D: " + str(C.index(D)) +" out of  " + str(len(C)) + " computed")
        Ds = semisimple_part_hensel_lift(D)
        Cs.append(Ds)

    N = n*d*(1+2*m+n) /eps
    p = [ZZ.random_element(1,N+1)  for i in range(m)]
    x = get_generators(flist)

    # Coerce all fi into the fraction field of R
    F = x[0].parent().fraction_field()
    flist = [F(fi) for fi in flist]
    J = jacobian(flist, x)(*p)

    qlist = []
    for fi in flist:
        try:
            qlist.append(fi(*p))
        except (ZeroDivisionError, ValueError, ArithmeticError):
            return "fail"

    q = vector(qlist)
    Csq = [Ds * q for Ds in Cs]

    DqAsRows = matrix(QQ, [list(Dq) for Dq in Csq])
 
    try:
        DqAsRows.solve_left(J.transpose())
        return (True)
    except (ZeroDivisionError, ValueError, ArithmeticError):
        return (False)


def semisimple_part_using_JNF_ober_QQbar(A):
    """
    Computes the semisimple part of a matrix A using the JNF of A over QQbar.
    This version is the old slow verion and is now only used to test the code for the new version semisimple_part_hensel_lift().
    
    Input: 
        A (matrix): Matrix A over QQ
    Output: 
        matrix: Semisimple part of A over QQ
    """
    n = A.nrows()
    MQQbar = MatrixSpace(QQbar,n)
    if MQQbar(A).is_diagonalizable():
        As = A
    else:
        Jf = MQQbar(A).jordan_form(transformation=True)
        D = copy(Jf[0])
        T = copy(Jf[1])
        for i1 in range(n):
            for i2 in range(i1+1,n):
                D[i1, i2] = QQbar(0)
        As = T * D * T.inverse()
    MQQ = MatrixSpace(QQ,n)
    return MQQ(As)    

def semisimple_part_hensel_lift(A):
    """
    Computes the semisimple part of a matrix A using Hensel lift of the minimalpolynomial and its squarefree factor
    This version is faster and should be used.
    
    Input: 
        A (matrix): Matrix A over QQ
    Output: 
        matrix: Semisimple part of A over QQ
    """
    p_min = A.minpoly()

    #square-free part of p_min
    f = prod( [ factor for (factor,power) in p_min.factor() ] )
    z = find_zi(p_min, f)
    
    return z(A)

def find_zi(p_min, f):
    """
    Computes the zi such that f(zi) = 0 mod f^i , zi = t mod f and i is large enough such that f^i is divisible by p_min

    Input: 
        p_min (polynomial in ring QQ[x] with one generator): minimalpolynomial of a matrix A
        f (polynomial in ring QQ[x]): squarefreepart of p_min

    Ouput:
        polynomial in ring QQ[x]
    """
    z = f.parent().gen()
    if p_min.divides(f):
        return z

    i = 1
    df = f.derivative()
    a = inverse_mod(df,f)
    q = f
    
    while not p_min.divides(f^i):
        z = z -a(z) * f(z)
        i += 1
        q *= f
    return z
        
        

def commutator_kernel_basis(As, B):
    """
    Compute a basis C for { D = sum_i c_i * B[i] : As*D - D*As = 0 }. 
    Input: 
        As (matrix): matrix over QQ
        B (list): list of matrices over QQ
    
    Output:
        list pof matrices:basis for the kernel subspace 
    """

    AsDList = [As * b - b * As for b in B]

    AsDMat = matrix(QQ, [list(vector(M)) for M in AsDList]).transpose()

    sols = AsDMat.right_kernel().basis_matrix()
    C = []
    for row in sols.rows():
        D = sum(coeff * b for coeff, b in zip(row, B))

        C.append(D)

    return C


