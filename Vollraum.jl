using SpecialFunctions
using LinearAlgebra

function theta(t, sigma)
    sqrt.(pi * 4e-7 * sigma / 4 ./ t)
end

function Et(r, m, t, sigma)
    mu = 4e-7 * pi
    x = r[1]
    y = r[2]
    z = r[3]
    R = norm(r)
    
    theta = sqrt(mu * sigma / 4.0 / t)
    tr = theta * R
    t2r2 = tr^2
    alpha = -mu / (4 * pi)
    
    E = alpha / R^3 * cross(m, r) * 
        2.0 * tr / (sqrt(pi) * t) * t2r2 * exp(-t2r2)

    return E
end

function Ht(r, m, t, sigma)
    mu = 4e-7 * pi
    x = r[1]
    y = r[2]
    z = r[3]
    R = norm(r)
    
    theta = sqrt(mu * sigma / 4.0 / t)
    tr = theta * R
    t2r2 = tr^2
    t3r3 = t2r2 * tr
    sqpi = sqrt(pi)
    alpha = 1.0 / (4 * pi * R^3)

    H = alpha * 
        (dot(m, r) / R^2  * r
        *
        (3 * SpecialFunctions.erfc(tr) +
            (6 * tr / sqpi + 4 * t3r3 / sqpi) * exp(-t2r2)) 
            .-
        (SpecialFunctions.erfc(tr) +
            (4 * t3r3 / sqpi + 2 * tr / sqpi) * exp(-t2r2)) 
        * m)

    return H
end

function dHdt(r, m, t, sigma)
    mu = 4e-7 * pi
    x = r[1]
    y = r[2]
    z = r[3]
    R = norm(r)
    
    theta = sqrt(mu * sigma / 4.0 / t)
    tr = theta * R
    t2r2 = tr^2
    sqpi = sqrt(pi)

    dHdt = theta^3 / sqpi^3 / t * exp(-t2r2) *
        (
            t2r2 * dot(m, r) * r / R^2 .+
            (1 - t2r2) * m
        )

    return dHdt
end