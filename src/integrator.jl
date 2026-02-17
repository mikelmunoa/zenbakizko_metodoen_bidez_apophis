# src/integratzaileak.jl

"""
    RK4(u0, t0, T, h, p, f!, m = 1)
Runge-Kutta 4. ordenako integratzailea.
u0: hasierako egoera, t0: hasiera denbora, T: amaiera denbora,
h: integrazio pausoa, p: parametroak, f!: deribatu funtzioa, m: barne-pausoak.
"""
function RK4(u0, t0, T, h, p, f!, m = 1)
    dt_out = h * m
    N_steps = round(Int, (T - t0) / dt_out)
    
    uu = [copy(u0) for _ in 1:(N_steps + 1)]
    tt = Vector{Float64}(undef, N_steps + 1)
    
    tt[1] = t0
    
    uk = copy(u0)
    utmp = copy(u0)
    k1 = similar(u0)
    k2 = similar(u0)
    k3 = similar(u0)
    k4 = similar(u0)
    
    tk = t0
    
    for i in 1:N_steps
        for _ in 1:m
            f!(k1, uk, p, tk)
            
            @. utmp = uk + 0.5 * h * k1
            f!(k2, utmp, p, tk + h/2)
            
            @. utmp = uk + 0.5 * h * k2
            f!(k3, utmp, p, tk + h/2)
            
            @. utmp = uk + h * k3
            f!(k4, utmp, p, tk + h)
            
            @. uk = uk + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
            tk += h
        end
        tt[i+1] = tk
        uu[i+1] = copy(uk)
    end
    return tt, uu
end
