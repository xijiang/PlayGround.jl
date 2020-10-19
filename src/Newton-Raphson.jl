"""
    Newton_Raphson(f, fp, x; tol = 1e-14, MITR = 100)
---
Given function `f` and its derivative `fp`, this function finds a root for `f = 0` with
a starter `x`
"""
function Newton_Raphson(f, fp, x; tol = 1e-14, MITR = 100)
    xnew, xold = x, Inf
    fn, fo = f(xnew), Inf

    ctr = 1

    while (ctr < 100) && (abs(xnew - xold) > tol) && ( abs(fn - fo) > tol )
	x = xnew - f(xnew)/fp(xnew) # update step
	xnew, xold = x, xnew
        fn, fo = f(xnew), fn
	ctr = ctr + 1
    end

    if ctr == MITR
	error("Did not converge in $MITR steps")
    else
	xnew, ctr
    end
end

"""
    Newton_Raphson(f, x; tol = 1e-14, MITR = 100, inc = 1e-8)
---
Given function `f`, this function numerically finds a root for `f = 0` with a starter `x` using `secant` method.
"""
function Newton_Raphson(f, x; tol = 1e-14, MITR = 100, inc = 1e-8)
    fp(x1, x2) = (f(x1) - f(x2))/(x1 - x2)

    fn, fo = f(x1), Inf

    ctr = 1

    x2, x1 = x + inc, x
    while (ctr < 100) && (abs(x2 - x1) > tol) && ( abs(fn - fo) > tol )
        x2, x1 = x2 - f(x2)/fp(x1, x2), x2
        fn, fo = f(x2), fn
	ctr = ctr + 1
    end

    if ctr == 100
	error("Did not converge in 100 steps")
    else
	x2, ctr
    end
end
