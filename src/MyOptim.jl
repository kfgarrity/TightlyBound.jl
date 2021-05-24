module MyOptim

using LinearAlgebra

function conjgrad(fn, grad, x0; maxstep=2.0, niters=50, conv_thr = 1e-2)

#    println("CG START")
    x = deepcopy(x0)
    storage = zeros(size(x))

    xold = deepcopy(x)

#    println("x0 ", x)

#    println("CG RUN GRAD1")
    #iter 1
    f, g = grad(storage, x)

#    println("f ", f, " x   ",  x)

    dx = deepcopy(-g)

    dx1 = deepcopy(dx)

    s = deepcopy(dx)
    
    step_size = maxstep/10.0

#    println("CG RUN LS1")

    x, step_size = linesearch(x, dx, fn, f, step_size)
    step_size = min(step_size, maxstep)



    for i = 1:niters
        

        if sum(abs.(xold - x)) > 1e-7  #only update grad if x changes
            f, g = grad(storage, x)
            xold = deepcopy(x)
        end

        if sqrt(sum(g.^2)) < conv_thr && f < 0.1
            println("yes conv")
            return x, f, g
        end

#        println("MY CONJGRAD $i ", f, " sg ", sqrt(sum(g.^2)), " step_size $step_size --------------------")


        dx1[:] = dx[:]
        dx[:] = -g[:]

#        println("dx ", dx)
#        println("dx1 ", dx1)

        beta = ( dx' * (dx - dx1)) / (dx1' * dx1)  #PR
#        println("beta $beta")
        beta = max(0, beta)

        s = dx + beta * s

        x, step_size, good = linesearch(x, s, fn, f, step_size)
        if good
            step_size =min(step_size, maxstep)
        else
            beta = 0.0
        end
        println("MY CG $i ", f, " sg ", sqrt(sum(g.^2)), " step_size $step_size good $good ")


    end

    println("conv not achieved")
    return x, f, g

end

function linesearch(x, dx, fn, f0, step_size)
#    println("CG FN0 $f0")
#    println("CG dx $dx")
#    println("ls $x $dx $f0 $step_size")
#    println("CG FN1 ss $step_size")
    f0r, flag0r = fn( x + dx * step_size * 0.0)
#    println("CG real $f0r $flag0r")

    f1, flag1 = fn( x + dx * step_size * 0.5)
#    println("CG FN1 $f1 $flag1")
    if f1 > f0 || flag1 == false
#        println("MY LS uphill $f1 > $f0 flag $flag1, reduce step")
        return x, step_size/3.1, false
    end
#    println("CG FN2")    
    f2, flag2 = fn( x + dx * step_size)
#    println("CG FN2 $f2 $flag2")

    if flag2 == false
#        println("MY LS flag2 $flag2")
        return x + dx * step_size * 0.5, step_size/1.75, false
    end

    c = f0
    a = 2 * (f2-2*f1+f0)
    b = f2 - a - c

#    println("a b c $a $b $c ", -b/(2*a))

    if a < 0
#        println("MY LS a negative $a, take max step increase stepsize")
        return x + dx * step_size * 1.0, step_size * 2.1, true
    else
        step = min(-b/(2*a), 1.0)
#        println("MY LS NORMAL $step")
        return x + dx * step * step_size, step * 2.0, true
    end
    

end


end #end module

#=
function myf(x)

    return -0.5*x'*A*x + 0.1*sum(x.^4)

end

using ForwardDiff



function grad(storage, x)
    sleep(0.1)

    return myf(x), ForwardDiff.gradient(myf, x)

end

function grad_only(storage, x)
    storage[:] = ForwardDiff.gradient(myf, x)
    sleep(0.1)
    return storage
end
=#
