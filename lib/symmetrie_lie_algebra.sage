import random


def symmetry_lie_algebra_random(flist , epsilon = 0):
    """
    Computes the Symmetry lie algebra
    
    Input
        flist (list): A list of rational functions over QQ[x_1,...,x_m]
        epsilon (rational, optional): Enter epsilon in (0,1) to have probablity as in the paper, if left empty reduce random choices to [0,1,-1,2,-2,1/2,-1/2, 1/3, 2/3, -1/3, -2/3, 3,-3]
    
    Output: 
        list of Matrices:  which is a basis of the space of the symmetry Lie alebra
    """

    if 0 < epsilon < 1:
        small_entries = False
        eps = epsilon
    elif epsilon == 0:
        small_entries = True
        eps = 1
    else:
        raise ValueError("epsilon must be 0 or a number strictly between 0 and 1")
        
    d = max([max_degree(fi) for fi in flist])
    n = len(flist)
    m = len(flist[0].parent().gens())
    N = 2 * n^3 * d * (2 * m + 1) / eps
    L = list(MatrixSpace(QQ, n).basis())
    x = get_generators(flist)
    # Coerce all fi into the fraction field of R
    F = x[0].parent().fraction_field()

    flist = [F(fi) for fi in flist]
    J = jacobian(flist, x)
    i = 0
    spaces_are_not_equal = True

    rank = rank_guess(J, m, 30, N)

    # print("Rank J: " + str(rank))


    #In case that L = Lnew to early, we repeat for the same L 5 times.
    stuck = 0
    while i < n^2 and spaces_are_not_equal and stuck < 5:
        p = [random.choice([0,1,-1,2,-2,1/2,-1/2, 1/3, 2/3, -1/3, -2/3, 3,-3]) if small_entries else ZZ.random_element(1,N+1) for i in range(m)]
        qlist = []
        for fi in flist:
            try:
                qlist.append(fi(*p))
            except (ZeroDivisionError, ValueError, ArithmeticError):
                print("Try: " + str(i) + ": ZeroDivisionError")
                i += 1
                continue
        q = vector(qlist)
        Jp = J(*p)
        if rank > Jp.rank():
            # print("Try: " + str(i) + ": rankexception")
            i += 1
            continue
        
        Lnew = find_linear_subspace(L,q,Jp)
        if len(Lnew) == len(L) and stuck <15:
            stuck += 1
        elif len(Lnew) == len(L) and stuck == 15:
            spaces_are_not_equal = False
        elif len(Lnew) < len(L):
            stuck = 0
        # print("Max bit-size of numerators in L:", max_entry_bits(L))
        L = Lnew
        i += 1
        # print("Try: " + str(i) )
        # print("Dimension of L: " +str(len(Lnew)))
    return L


def rank_guess(J, m, tries, max_int):
    """
    Guesses the rank of J: Evaluate J at random point p and returns the largest value.

    Input:
        J (matrix): Jacobian with entries in QQ[x_1,...,x_m]
        m (int): number of variables
        tries (int): number of tries to evaluate the Jacobian
        max_int (int): random entries in p ly in the intervall (1,max_int)

    Output:
        int: the rank of J
    """
    rank_guess = 0
    for i in range(tries):
        p = [ZZ.random_element(1,max_int+1) for i in range(m)]
        try:
            rank_new = J(*p).rank()
        except (ZeroDivisionError, ValueError, ArithmeticError):
            continue
        rank_new = J(*p).rank()
        if rank_new > rank_guess:
            rank_guess = rank_new
    return rank_guess
    

def max_degree(f):
    """
    Returns the maximum degree of numerator and denominator of f
    
    Input: 
        f: a rational function over QQ[x_1,...,x_m].

    Output:
        int: maximal degree
    """
    # Try to get the base polynomial ring
    if hasattr(f.parent(), "ring"):    # case: fraction field
        R = f.parent().ring()
    elif f in f.parent().fraction_field():  # case: polynomial
        R = f.parent()
    else:
        # fallback: build polynomial ring from numerator
        R = f.numerator().parent()

    # Coerce numerator and denominator into R
    num_deg = R(f.numerator()).degree()
    den_deg = R(f.denominator()).degree()
    return max(num_deg, den_deg)

def get_generators(flist):
    """
    Input: 
        flist: list of rational functions/polynomials/integers over QQ
        
    Output: 
        list of generators (of the form (x1, ..., xm))
    """
    for f in flist:
        P = f.parent()
        if hasattr(P, "gens"):  # polynomial ring
            return P.gens()
        if hasattr(P, "ring"):  # fraction field -> go down to ring
            return P.ring().gens()
    raise ValueError("No element in list has polynomial generators.")

def find_linear_subspace(L,q,J):
    """
    Computes the linear subspace L' of L sucht that all Elements of A in L' satisfy Aq is contained in the image of J
    
    Input:  
        L (list of matrices): list of  Basis of linear subspace 
        q (vector)
        J (matrix)
    Ouptut: 
        list of matrices: Basis of the linear subspace L' of L satisfying Aq in im(J) for all A in L'
    """
    Aq_list = [list(B * q) for B in L]
    Aq_matrix = transpose(matrix(QQ, Aq_list))
    J_null = J.transpose().right_kernel().basis_matrix()
    sols = (J_null * Aq_matrix).right_kernel()
    Lnew_basis = []
    for row in sols.basis():
        lnew = sum(coeff * A for coeff, A in zip(row, L))
        Lnew_basis.append(lnew)
    return Lnew_basis

def max_entry_bits(L):
    """
    Outputs the maximal entry size of elements of L, can be used to analyize the code
    """
    sizes = []
    for A in L:
        for e in A.list():  # all entries of the matrix
            if e != 0:
                num = QQ(e).numerator()
                sizes.append(num.nbits())
    return max(sizes) if sizes else 0

def normalize_matrix(A):
    """
    Returns a new matrix with smaller integer entries.
    """
    # clear denominators
    A = A.apply_map(QQ)   # ensure rational entries
    den = lcm([a.denominator() for a in A.list() if a != 0]) or 1
    A = (den * A).change_ring(ZZ)  # integer matrix

    # divide out gcd of all entries
    g = gcd([a for a in A.list() if a != 0]) or 1
    return (A / g)
