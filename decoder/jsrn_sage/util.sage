#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Various utility functions used in codinglib"""
import sys, re
import sage.sets
import sage.sets.cartesian_product

######
###### Plotting and Presenting
######


# Short list of various colors
some_colors = ["blue","red","green","orange","pink","black",
               "yellow","darkgreen","turquoise"]

def plotc(funs, varRange, *args, **kwargs):
    """Plot the list of funs with different colors.
       If each element of funs is a tuple, the second argument will be
       interpreted as the legend label."""
    if isinstance(funs[0], tuple):
        return sum([ plot(funs[i][0], varRange, *args,
                        color=some_colors[i], legend_label=funs[i][1],
                        **kwargs)
                for i in range(0, len(funs)) ])
    else:
        return sum([ plot(funs[i], varRange, *args,
                        color=some_colors[i], **kwargs)
                for i in range(0, len(funs)) ])

def int_plot(fun, varRange, log_yscale=False, **kwargs):
    """Plot the function in the range of the variable. Range is specified
    as in plot: (var, minimal, max)"""
    (var, minx, maxx) = varRange
    if not log_yscale:
        funPoints = [ (x,fun(x)) for x in range(minx, maxx+1) ]
    else:
        funPoints = [ (x, log(fun(x))) for x in range(minx, maxx+1) ]
    #Use our version of point plot to support diff markers
    import codinglib.plot_util
    return codinglib.plot_util.point(funPoints, **kwargs)

def multi_plot(funs, varRange, shift_equals=True, log_yscale=False, **kwargs):
    """Plot all functions in the range of the variable in different colors.
    These are plotted as list_plot for the integer values. Range is specified
    as in plot: (var, minimal, max)"""
    nfuns = len(funs)
    (var, minx, maxx) = varRange
    numx = maxx-minx+1
    isLabeled = False
    if isinstance(funs[0], tuple):
        isLabeled = True
        (funs, funLabels) = zip(*funs)
        
    if not log_yscale:
        funVals = [[ fun(x) for x in range(minx, maxx+1) ] for fun in funs ]
    else:
        funVals = [[ log(fun(x)) for x in range(minx, maxx+1) ]
                        for fun in funs ]

    if shift_equals:
        # Get a scale for how much "a little" is for the raising
        shiftLeft = 1.5/400*numx
        allVals = flatten(funVals)
        minVal = min(allVals)
        maxVal = max(allVals)
        shiftUp = 1.5/400*(maxVal-minVal)
    else:
        shiftLeft = 0
        shiftUp = 0

    # For each var value, construct a map N -> N that tells how many functions
    # so far have a certain function value.
    atVal = [ dict() for i in range(minx, maxx+1) ] 
            
    # Make the points to place the functions
    # Go through all function values for each function (in reverse order) and
    # slightly raise them to the up and left depending on the no. of lower
    # functions which have this function value
    funi = range(0,nfuns)
    funi.reverse()
    funPoints = [ [None for i in range(0,numx)] for j in range(0, nfuns) ]
    for i in funi:
        for x in range(0, numx):
            lower = atVal[x].get(funVals[i][x], 0)
            funPoints[i][x] = (x + minx - lower*shiftLeft,
                                funVals[i][x] + lower*shiftUp)
            atVal[x][funVals[i][x]] = lower + 1
    g = Graphics()
    #Use our version of point plot to support diff markers
    import codinglib.plot_util
    for i in range(0, nfuns):
        if isLabeled:
            g += codinglib.plot_util.point(funPoints[i], color=some_colors[i], 
                    legend_label=funLabels[i], **kwargs)
        else:
            g += codinglib.plot_util.point(funPoints[i], color=some_colors[i], **kwargs)
    return g

def ltext(s, surroundDollar=True, returnEndState=False):
    r"""Wrap the string s in LaTeX markup so that it will be rendered as LaTeX
    text in e.g. the legend of a plot. The text is treated as inline LaTeX
    beginning in normal mode; thus, usual text is rendered until a $ is met, at
    which point it changes to math until the end $ is met.
    Some uses (e.g. when used in the legend of a plot) needs the string to
    begin with $ for it to be recognised as LaTeX; if this is not the case set
    surroundDollar to False.
    If ''returnEndState'' is True a tuple is return with the first element
    being the typeset string and the second being True if the input string
    terminated in math mode, otherwise False. This might be useful if treating
    a longer string, each line being handled by ''ltext''. Note that no matter
    whether the input string terminated in math mode or not, ''ltext'' always
    returns a balanced string."""
    i = 0
    last = 0
    if surroundDollar:
        ss="$"
    else:
        ss=""
    ss += r"{\rm "
    inText = True
    for i in range(0,len(s)):
        if s[i] == '$':
            ss += s[last:i]
            if inText:
                ss += r"}"
            else:
                ss += r"{\rm "
            last = i+1
            inText = not inText
        if s[i] == ' ' and inText:
            ss += s[last:i]
            ss += r"\ "
            last = i+1
    ss += s[last:]
    if inText:
        ss += "}"
    if surroundDollar:
        ss += "$"
    if returnEndState:
        return (ss, not inText)
    else:
        return ss

def write(s, *args):
    r"""For the notebook: similar to ''view'' but wraps the string in ''ltext'', thus rendering
    text nicer and treating text in $$ as LaTeX mathematics. The string will be
    split at newlines, each one being printed separately. ''
    
    EXAMPLE:
        sage: res = 1/2
        sage: write("Result are $n + m = %s$\nTook %s secs", res, timing)
        <html><span class="math">{\rm Result\ are\ }n + m = \frac{1}{2}{\rm }</span></html>
        <html><span class="math">{\rm Took\ 2\ secs}</span></html>'
    """
    def to_latex(obj):
        str = latex(obj)
        return str.replace(r"\Bold", r"\mathbf") #Stupid Sage crap
    tag='span'
    lines = s.split('\n')
    htmlLines = []
    for line in lines:
        line = "$"+line+"$"
        typed =ltext(line, surroundDollar=False)
        htmlLines.append( '<html><%s class="math">%s</%s></html>' % (tag, typed, tag) )
    print( ("\n".join(htmlLines)) % tuple( to_latex(a) for a in args ) )

def printlist(l):
    r"""Print a list with one element per line"""
    elems = [ latex(e) for e in l ]
    return LatexExpr(r"\begin{array}{l}[ \\" + r", \\".join(elems) + r"\\ ]\end{array}")

def result_table(titles, results, types=None, col_width=16, float_precision=2, exp_precision=6):
    """Writes a table of results from e.g. testing.

    - titles: a list of strings, zeroth element is title of results' keys
    - results: a dictionary of lists, each list of length len(titles)
    Sorts the dictionary's keys canonically or using sort (a cmp function), and
    writes results on separate lines.

    Types can be omitted in which case they are guessed.
    Types are given as a list of strings, the zeroth element being the type of
    the keys of `results`. Allowed types are:
    - "str"
    - "int"
    - "float"
    - "exp"
    - "%"
    
    TODO: Specify sort function
    """
    type_to_format = {
        "str"   : "{:>%d}" % col_width ,
        "int"   : "{:>%dd}" % col_width ,
        "float" : "{:>%d.%d}" % (col_width, float_precision) ,
        "exp"   : "{:>%d.%de}" % (col_width, exp_precision) ,
        "%"     : "{:>%d.2%%}" % col_width, 
    }
    if types is None:
        raise NotImplementedError("Type guessing for result table")
    def get_type(col, x):
        typ = types[col]
        if not typ in type_to_format:
            raise ValueError("Unrecognised type: %s" % typ)
        tx = None
        if typ=="int":
            tx = int(x)
        elif typ=="float" or typ=="%" or typ=="exp":
            tx = float(x)
        else:
            tx = x
        return (type_to_format[typ], tx)
    title_str_s = "\t".join([type_to_format["str"]]*len(titles))
    title_str = title_str_s.format(*titles)
    print title_str
    print "-"*len(title_str.expandtabs())
    res_keys = results.keys()
    res_keys.sort()
    for key in res_keys:
        res = [ get_type(i, x) for (i, x) in zip(range(len(results[key])+1), [key]+results[key]) ]
        print "\t".join([ f for (f,x) in res ]).format(*( x for (f,x) in res ))

def search_dir(obj, pattern):
    """Find the methods and fields of `obj` which matches the search pattern
    (with ".*" added at either end).
    This is similar to %psearch obj.*<pattern>*
    """
    reg = re.compile(".*"+ pattern + ".*")
    return [ m for m in dir(obj) if reg.match(m) ] 
    
    
######
###### Various testing
######


def find_minimal_satisfiable(f, startn=1, contiguous=True):
    """Find the minimal integral n>0 such that f(n) == true. If the interval
    for which f is true is contiguous and open towards infinity, a logarithmic
    algorithm is used, otherwise linear. startn can be given as a hint to a
    value that might be true."""
    if not contiguous:
        n=startn
        if f(n):
            while f(n) and n>0:
                n=n-1
            return n+1
        else:
            while not f(n):
                n=n+1
            return n
    else:
        maxn = startn
        minn = 1
        # Keep doubling n to find one that works and then binary
        while not f(maxn):
            minn = maxn+1
            maxn *= 2
        while minn < maxn:
            tryn = minn + floor((maxn-minn)/2)
            if f(tryn):
                maxn=tryn
            else:
                minn=tryn+1
        return maxn

def find_best_monte_carlo(new, order=None, score=None):
    """Generate random solutions for something using function new()
    indefinitely, printing the best solution whenever one is found. "Best" can
    be determined in two ways:
    1) by a total ordering of solutions by setting function order, where
    order(a,b) = True if a is a better solution b.
    2) by a scoring using function score(a) which should return an element of a
    type orderable by ">".
    """
    if order is None and score is None:
        raise Exception("Supply either order or score function")
    if not order is None:
        best = new()
        while True:
            cand = new()
            if order(cand, best):
                print "New best: %s" % (cand,)
                best = cand
    if not score is None:
        best = new()
        bestScore = score(best)
        while True:
            cand = new()
            candScore = score(cand)
            if candScore > bestScore:
                print "New best %s having score %s" % (cand,candScore)
                best = cand
                bestScore = candScore



print_progress_counter=1000
def print_progress(cur_iter, total_iters, updates=50):
    """Prints a progress bar for long-running computations.  Call each iteration of the computation, giving current iteration and total number of iterations."""
    global print_progress_counter
    if print_progress_counter > cur_iter:
        print("|"+"-"*updates+"|")
        sys.stdout.write("|")
        sys.stdout.flush()
    print_progress_counter=cur_iter
    #TODO: When total_iters<50 more dots should be printed at a time
    if total_iters < 50 or cur_iter % floor(total_iters/50) == 0:
        sys.stdout.write(".")
        sys.stdout.flush()
    if cur_iter >= total_iters-1:
        print("|")


    
######
###### Comparison testing framework
######
def comparison_test(workers, parameter_iter, problem_builder, iters=1,
                    worker_allowed=None, continue_results=None, check_results=False, debug=None):
    r"""Framework for running timing comparisons between different algorithms
    for solving the same problem. Problems are generated from give parameters
    and solved by the supplied worker algorithms, and results are harvested.

    ARGUMENTS::
    'workers' - a dictionary from name into a function. Each function has signature Problems -> Answers

    'parameter_iter' - An iterator for producing parameters. Should return
                       Parameter * Key tuples. Keys identify a coherent set of
                       parameters for a given test. For instance, Parameters
                       might be 5-tuples of intricately dependent values, but in
                       the given test, they are calculated deterministcally,
                       based on only one value. This value is then the key.
                       Testing stops when this returns None or raises
                       StopIteration.
                       To use comparison_test_plot, Key must be a non-zero
                       number and it should denote the "size" of the problem.
                       
    'problem_builder' - A function from Parameters -> Problems.

    'iters=1' - How many problems are created for each set of parameters.

    'worker_allowed' - A function (Worker-name * Key) -> Boolean for
                       fine-grained control over which workers should be applied
                       to which sets of parameters (e.g. if some workers are
                       known to be very slow on large input)

    'continue_results=None' - Can be set to the results of a previous run in
                              order to continue on the aggregated results there.
                              Already tested sets of parameters will then be
                              skipped.

    'check_results=False' - If True, then check that the outputs of the workers
                            are all the same.
    

    OUTPUT::
    A dictionary results from Worker-Name into dictionary of timings. The inner
    dictionaries are from Key to list of timings, each timing being a float
    counting the number of seconds the worker took on a given problem. The
    problems themselves are not saved.
    """
    import gc, timeit

    if continue_results is None:
        results = dict({ name: dict() for name in workers })
    else:
        results = continue_results

    try:
        params_key = parameter_iter()
        while params_key:
            params, key = params_key
            if debug > 0:
                print "Parameters given for key %s" % key
            for name in workers:
                if not key in results[name]:
                    results[name][key] = []
            if any( len(results[name][key]) < iters for name in workers if worker_allowed is None or worker_allowed(name, key)):
                iters_left = iters - min(len(results[name][key]) for name in workers if worker_allowed is None or worker_allowed(name, key))
                for i in range(iters_left):
                    problem  = problem_builder(params)
                    if debug:
                        if debug > 1:
                            print "Comparison problem: %s" % problem
                        else:
                            print "\tCreated problem"
                    old_answer = None
                    for name in workers:
                        if len(results[name][key]) < iters and (worker_allowed is None or worker_allowed(name, key)):
                            alg = workers[name]
                            before = timeit.default_timer()
                            answer = alg(problem)
                            clock = timeit.default_timer() - before
                            if debug > 1:
                                print "\t\t%s on key %s took %s secs" % (name, key, clock)
                            results[name][key].append(clock)
                            if check_results and old_answer:
                                assert answer == old_answer , "Answer from %s did not correspond to earlier answer.\n\tProblem: %s" % (name, problem)
                            gc.collect()
            for name in workers:
                results[name][key].sort()
            params_key = parameter_iter()
        if debug:
            print "Simulation done"
    except StopIteration:
        return results
    return results

def comparison_test_save(test, results):
    """Save the results of a comparison test"""
    import pickle
    afile = open(test + ".data", 'wb')
    pickle.dump(results, afile)
    afile.close()

def comparison_test_load(test):
    """Load the results of a comparison test"""
    import pickle
    afile = open(test + ".data", 'r')
    results = pickle.load(afile)
    afile.close()
    return results

def comparison_test_plot(workers, results, **options):
    """Plot the results of a comparison test using log-scales. On the x-axis
    will be the log of parameter keys, and y-axis will be the log of time spent.
    Only the median value over all iterations is plotted.
    """
    #Choose colors
    colors = dict()
    for name, color in zip(workers.keys(), rainbow(len(workers.keys()))):
        colors[name] = color

    def options_set(name, val):
        if not name in options:
            options[name] = val
        
    options_set("marker", "o")
    options_set("linestyle", "-")
    options_set("plotjoined", True)
    if not "scale" in options:
        options["scale"] = "loglog"
        def get_log_interval(vals):
            return floor(log(min(vals))), ceil(log(max(vals)))
        def get_log_ticks(log_interval):
            return list([ 10^i for i in range(log_interval[0], log_interval[1]+1) ])
        def get_tick_formatter(log_interval):
            return list([ "$10^{%s}$" % i for i in range(log_interval[0], log_interval[1]+1) ])
        time_interval = get_log_interval(list(flatten_once(
            [ (results[name][key][len(results[name][key])//2])
              for key in results[name] if results[name][key] ] for name in workers)))
        key_interval = get_log_interval(list(flatten_once(results[name].keys() for name in workers)))
        options_set("ticks", [ get_log_ticks(key_interval), get_log_ticks(time_interval) ])
        options_set("tick_formatter", [ get_tick_formatter(key_interval), get_tick_formatter(time_interval) ])

    options_set("axes_labels", [ "", ltext("secs.") ])
    options_set("fontsize", 15)
    #options_set("legend_handlelength", 2)
        
    g = Graphics()
    for name in workers:
        pts = [ (key, results[name][key][len(results[name][key])//2]) for key in results[name] if results[name][key] ]
        pts.sort(key=lambda (key,t): key)
        options['color']=colors[name]
        options['legend_label'] = ltext("%s" % (name,))
        g += list_plot(pts, **options)
    return g



        

######
###### Number theory
######

def ligt(x):
    """Least integer greater than x"""
    return int(x+1)

def gilt(x):
    """Greatest integer less than x"""
    if x == int(x):
        return x-1
    else:
        return int(x)

def find_integral_max(xmax, f):
    """Given non-integral xmax, find integral xmax locally maximising f(x), i.e.
    either f(floor(x)) or f(ceil(x)). Returns (x,f(x))"""
    if xmax==int(xmax):
        return (xmax,f(xmax))
    else:
        xmaxf = floor(xmax)
        xmaxc = xmaxf+1
        fxmaxf, fxmaxc = f(xmaxf), f(xmaxc)
        return (xmaxf,fxmaxf) if fxmaxf >= fxmaxc else (xmaxc,fxmaxc)
    
def find_modulo_power(n, delta, q):
    """Find (by search) minimal positive power of q such that q^m+delta % n ==
    0. Return None if no such exists."""
    m = 1
    qpow = q
    rem = (q+delta) % n
    seen = set()
    while not rem in seen:
        seen.add(rem)
        m += 1
        qpow = q*qpow
        rem = (qpow+delta) % n
        if rem == 0:
            return m
    return None

def solve2deg_int(a,b,c):
    """Utility function. Returns the greates integer range [i1, i2] such that
    i1 > x1 and i2 < x2 where x1,x2 are the two zeroes of the equation in x:
    ax^2+bx+c=0. If no solution, returns an empty, range w. negative coeff."""
    D = b^2 - 4*a*c
    if D < 0:
        return (-2,-1)
    sD = float(sqrt(D))
    minx , maxx = (-b-sD)/2/a , (-b+sD)/2/a
    mini , maxi = ( ligt(minx), gilt(maxx)  )
    if mini > maxi:
        return (-2,-1)
    else:
        return (mini,maxi)

def power_discr(a, p):
    """The power function a^p with the convention 0^0 = 1

    TODO: I think this is unnecessary in modern Sage and should be killed.
    """
    if type(p)==int:
        return a^p if not p==0 else 1
    else:
        return a^p if not p.is_zero() else 1

######
###### is_* functions
######
# Most of these were found by inspecting the list
#   <obj_of_interested_type>.__class__.__mro__
# and seeing what the most general super class in that list having the
# interesting type was
#####
def is_polynomial(p):
    return isinstance(p, sage.rings.polynomial.polynomial_element.Polynomial)

def is_matrix(M):
    return isinstance(M, sage.structure.element.Matrix)

def is_vector(v):
    return isinstance(v, sage.structure.element.Vector)

######
###### Algebra
######
def GFnice(q, names=None):
    """Return, if possible, a GF whose elements output in correct and nice log
    form. The pure GF function with my patch for repr=log has the bug that the
    generator used for expressing elements is not always the one used
    internally. This function fixes that changing which generator is
    returned."""
    if names==None:
        # No need to do anything, repr is fine as integers if q is prime
        # If q is not a prime, this will throw an exception
        return GF(q) 
    else:
        try:
            F = GF(q,names,repr='log')
            # If q=2^m then we don't need to change anything
            # Otherwise, overwrite the gens function in F to return the
            # internal multiplicative generator
            if not q.is_power_of(2):
                gnom = F.multiplicative_generator()
                ginternal = gnom^(1/gnom._element() % (q-1))
                def gen(self=None):
                    return ginternal
                def gens():
                    return [gen()]
                F.gen  = gen
                F.gens = gens
        except TypeError: #If repr is unsupported
            F = GF(q,names)
        return F

def list2poly(l, var):
    """Helper-function converting a list of coefficients to a polynomial in
    variable [var]"""
    p = 0
    for i in range(0,len(l)):
        p = p*var + l[len(l)-i-1]
    return p

def poly2list(p, len):
    """Convert the polynomial p into a list of coefficients of length len"""
    return list(p) + [0]*max(0, len-p.degree()-1)

def dft(p, w):
    """Calculate the Discrete Fourier Transform of polynomial p over primitive
    element w. Returns the list of the spectrum. Straightforward, slow
    algorithm"""
    return [ p(w^i) for i in range(0,w.multiplicative_order())]

def poly_collect(exp, collector, field=QQ):
    """If both exp and collector is a symbolic expression polynomial in all
    variables, this returns an expression where exp is described as a
    polynomial in the collector expression. The return value is a polynomial
    expression in a variable T which represent the collector expression."""
    orig_vars = list( set(exp.variables()).union(collector.variables()) )
    Nv = len(orig_vars)
    ovar_str = ",".join([ str(v) for v in orig_vars ])
    # Create the polynomial ring over the variables as well as a new one T
    P = PolynomialRing(field, ovar_str+",T",
            order = TermOrder('degrevlex',Nv) + TermOrder('degrevlex', 1))
    T = P.gens()[Nv] # The last variable must be our new one
    # Express the exp and collector in this ring
    pExp, pCollector = P(exp), P(collector)
    # Calculate the reduction of the expression by using an eliminator ideal
    J = Ideal(pCollector - T)
    gred = J.reduce(pExp)
    # Present this as a polynomial in T by coercing to a ring with same
    # variables but different ordering
    Pend = PolynomialRing(PolynomialRing(field, ovar_str), "T")
    return Pend(gred)

def matrix_subfield(M):
    """Expand each row in the matrix M to its row-wise representation in the
    prime subfield of M.base_ring()"""
    Md = [ [ M[i,j]._vector_() for j in range(0, M.ncols()) ]
                               for i in range(0, M.nrows()) ]
    Fbig = M.base_ring()
    F = Fbig.prime_subfield()
    deg = Fbig.degree()
    return matrix(F, M.nrows()*deg, M.ncols(),
            lambda i,j: Md[i//deg][j][i%deg])

def cyclotomic_coset(n, r, q):
    """The q-cyclotomic coset of r modulo n."""
    r = r%n
    cyc = set([r])
    rem = (r*q) % n
    while not rem in cyc:
        cyc.add(rem)
        rem = (rem*q) % n
    return list(cyc)

def nonzero_word(vec):
    """Returns True iff there is a nonzero element in vec. Faster than calling
    weight(vec)!=0."""
    for e in vec:
        if not e.is_zero():
            return True
    return False

def weight(vec):
    """Return the Hamming weight of the vector (iterable)"""
    i = 0
    for e in vec:
        if not e.is_zero():
            i = i+1
    return i

def support(vec):
    """Return a set consisting of the i such that vec[i] != 0"""
    s = set()
    for i in range(len(vec)):
        if not vec[i].is_zero():
            s.add(i)
    return s

def subset_word(vec, tau):
    """Return a random word which equals vec but where $weight(vec)-tau$ of the
    non-zero positions have been set to zero."""
    F = vec[0].parent()
    nonZero = [ i for i in range(0,len(vec)) if not vec[i].is_zero() ]
    shuffle(nonZero)
    subvec = [F.zero()]*len(vec)
    for i in nonZero[:tau]:
        subvec[i] = vec[i]
    return vector(subvec)

def poly_degs(T):
    """If T is a list/list of list/vector/matrix of univariate polynomials,
    return the degrees of these elements"""
    if type(T)==list or is_vector(T):
        if type(T[0])==list or is_vector(T[0]) \
           or (hasattr(T[0], "parent") and isinstance(T[0].parent(), sage.sets.cartesian_product.CartesianProduct)):
            action = 1
        else:
            action = 0
    else:
        # test whether T is a matrix, otherwise assume it is a somehow differently iterable
        if is_matrix(T):
            action = 2
        else:
            action = 0
    if action==0:
        return [ T[i].degree() for i in range(0, len(T))]
    elif action==1:
        return Matrix(ZZ, len(T), len(T[0]), lambda i,j: T[i][j].degree() if j < len(T[i]) else -2)
    elif action==2:
        return Matrix(ZZ, T.nrows(), T.ncols(), lambda i,j: T[i,j].degree())
    return None

def maxdeg(T, weights=None):
    """Return the maximal degree of elements in T. T can be any structure allowed by poly_degs"""
    degs = poly_degs(T)
    maxd = -1
    if weights is None:
        for s in degs:
            try:
                for ss in s:
                    if ss > maxd:
                        maxd = ss
            except:
                if s > maxd:
                    maxd = s
    else:
        if len(weights) != len(degs):
            raise ValueError("Given weights do not match dimensions of input")
        for i in range(len(weights)):
            s = degs[i]
            ws = weights[i]
            try:
                for j in range(len(s)):
                    if s[j] + ws[j] > maxd:
                        maxd = s[j] + ws[j]
            except:
                if s + ws > maxd:
                    maxd = s + ws
    return maxd

def LP(v, weights=None):
    """If v is a vector of polynomials, return the leading position of v using
    <_w where w is the weights vector (0 is assumed as all weights if none
    given). In case of tie, the highest position is given"""
    if not weights:
        weights=[0]*len(v)
    best = None
    bestp = None
    for p in range(0,len(v)):
        if not v[p].is_zero():
            vpdeg = v[p].degree() + weights[p]
            if vpdeg >= best:
                best = vpdeg
                bestp = p
    if best == None:
        return -1
    else:
        return bestp

def LT(v, weights=None):
    """If v is a vector of polynomials, return the leading term of v using <_w
    where w is the weights vector (0 is assumed as all weights if none
    given)."""
    return v[LP(v, weights=weights)]

def reverse_poly(p, deg):
    """Reverse the polynomial p according to the upper bound on the degree,
    deg. Note that you must promise deg >= p.degree()"""
    #TODO: Rename to poly_reverse
    ls = p.list()
    ls.reverse()
    return p.parent()([0]*(deg-len(ls)+1) + ls)

def vandermonde(xs, ncols, power=power_discr):
    r"""Return the matrix of dimension $|xs| \times ncols$ and whose $i,j$'th
    entry is $xs[i]^j$, defaulting 0^0=1."""
    return matrix(len(xs), ncols, lambda i,j: power(xs[i],j))

def hankel_matrix(S, cols):
    """Return the list of ring elements S as a hankel matrix with cols
    columns."""
    rows = len(S) - cols+1
    F = S[0].parent()
    return matrix(F, rows, cols, lambda i,j: S[i+j])

def irreducible_polynomial(PF, gdeg):
    """Find a random irreducible polynomial in polynomial ring PF of degree
    exactly gdeg"""
    # There seems to be a bug in polynomial ring random element with degree
    # bound, so it's not always the right degree
    g = None
    while (not g) or (not g.is_irreducible()):
        g = PF.random_element(degree=gdeg)
    return g

def poly_evaluate_build_M(points, PF):
    """Build M-table for use for multipoint evaluation. Gathen & Gerhard 10.3, p. 297"""
    (x,) = PF.gens()
    M = []
    M.append( [ x-a for a in points ] )
    lasth = len(points)
    h = lasth
    i = 0
    while lasth > 2:
        h = lasth - lasth//2 # half rounded up
        Mm = []
        for j in range(h): 
            m = M[i][2*j]
            if 2*j + 1 < lasth:
                m *= M[i][2*j + 1]
            Mm.append(m)
        M.append(Mm)
        lasth = h
        i += 1
    return M
    
def poly_evaluate(p, points, algorithm=None, M=None):
    """Evaluate (by pure Python) the polynomial p in all points points.
    High-degree polynomials are handled using Gathen & Gerhard 10.5,
    p. 298. However, this is only faster than trivial evaluation if fast
    remainder has been implemented (which it isn't for generic polynomials). For
    this, one can optionally supply the M table (using poly_evaluate_build_M) if
    repeatedly evaluating polynomials in the same points.
    """
    if algorithm=='naive' or (algorithm is None and p.degree() < 50):
        return [ p(a) for a in points ]
    else:
        F = p.base_ring()
        PF = p.parent()
        (x,) = PF.gens()
        if M is None:
            M = poly_evaluate_build_M(points, PF)
        evals = [ 0 for i in range(len(points)) ]
        def eval(p, i, j):
            if i<0:
                evals[j] = F(p)
            elif 2*j+1 < len(M[i]):
                r0, r1 = p % M[i][2*j] , p % M[i][2*j+1]
                eval(r0, i-1, 2*j)
                eval(r1, i-1, 2*j+1)
            else:
                eval(p, i-1, 2*j)
        eval(p, len(M)-1, 0)
        return evals


######
###### Functional list/set/dict
######
def flatten_once(lstlst):
    """Given a list of lists, flatten that list once, i.e. return one list with
    the elements being the total of the elements of the sublists, in that
    order.
    Returns a generator."""
    for lst in lstlst:
        for e in lst:
            yield e

def pad(lst, n, padwith):
    """Pads the list lst in-place with padwith until it has length n. If its
    length is <= n, then it is unchanged"""
    if n > len(lst):
        lst.extend([padwith]*(n-len(lst)))


