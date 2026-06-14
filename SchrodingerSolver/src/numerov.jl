"""
    solve_Numerov(V, x, m, planck, E_guess; tol=1e-10, max_iter=200)

Rozwiązuje jednowymiarowe, niezależne od czasu równanie Schrödingera metodą Numerova
(shooting method) z dokładnością O(h^4).

# Argumenty:
- `V`: Wektor wartości potencjału w węzłach siatki przestrzennej.
- `x`: Wektor punktów siatki przestrzennej (równomierna).
- `m`: Masa cząstki (`Float64`).
- `planck`: Zredukowana stała Plancka ℏ (`Float64`).
- `E_guess`: Przybliżenie szukanej energii własnej (`Float64`).
- `tol`: Tolerancja bisekcji (domyślnie `1e-10`).
- `max_iter`: Maksymalna liczba iteracji bisekcji (domyślnie `200`).

# Zwraca:
- Krotkę `(E, psi)`:
  - `E`: Obliczona energia własna (`Float64`).
  - `psi`: Znormalizowana funkcja falowa (`Vector{Float64}`).
"""
function solve_Numerov(V::Vector{Float64}, x::Vector{Float64},
                       m::Float64, planck::Float64, E_guess::Float64;
                       tol::Float64=1e-10, max_iter::Int=200)

    N  = length(x)
    dx = x[2] - x[1]

    # f(x, E) = (2m/ℏ²)(V(x) - E)  —  prawa strona równania ψ'' = f·ψ
    f_vec(E) = (2.0 * m / planck^2) .* (V .- E)

    # Propagacja wzorem Numerova: ψ(x_1)=0, ψ(x_2)=małe ziarno
    function numerov_shoot(E::Float64)
        f = f_vec(E)
        psi = zeros(N)
        psi[1] = 0.0
        psi[2] = dx * 1e-6

        for n in 2:(N-1)
            num = 2.0 * (1.0 + (5.0/12.0) * dx^2 * f[n])   * psi[n] -
                         (1.0 - (1.0/12.0) * dx^2 * f[n-1]) * psi[n-1]
            den = 1.0 - (1.0/12.0) * dx^2 * f[n+1]
            psi[n+1] = num / den
        end
        return psi
    end

    # Wartość brzegowa ψ(x_N) — ma być zerem dla poprawnego E
    boundary(E) = numerov_shoot(E)[N]

    # Wyznaczenie przedziału bisekcji wokół E_guess
    dE   = abs(E_guess) * 0.5 + 0.5
    E_lo = E_guess - dE
    E_hi = E_guess + dE
    bv_lo, bv_hi = boundary(E_lo), boundary(E_hi)

    for _ in 1:50
        bv_lo * bv_hi < 0 && break
        E_lo -= dE;  E_hi += dE
        bv_lo, bv_hi = boundary(E_lo), boundary(E_hi)
    end

    bv_lo * bv_hi >= 0 && error(
        "Nie znaleziono zmiany znaku wokół E_guess=$E_guess. Zmień E_guess lub zakres x.")

    # Bisekcja
    for _ in 1:max_iter
        E_mid  = (E_lo + E_hi) / 2.0
        abs(E_hi - E_lo) < tol && (E_lo = E_mid; break)
        bv_mid = boundary(E_mid)
        if bv_lo * bv_mid < 0
            E_hi = E_mid;  bv_hi = bv_mid
        else
            E_lo = E_mid;  bv_lo = bv_mid
        end
    end

    E_final  = (E_lo + E_hi) / 2.0
    psi_raw  = numerov_shoot(E_final)

    # Normalizacja trapezami
    norm2    = sum((psi_raw[1:end-1].^2 .+ psi_raw[2:end].^2) .* (dx / 2.0))
    psi_norm = psi_raw ./ sqrt(norm2)

    return E_final, psi_norm
end


"""
    E_analytic(n, m, planck, k)

Analityczna energia własna n-tego stanu oscylatora harmonicznego V(x) = ½kx².

    E_n = ℏ·ω·(n - ½),  ω = √(k/m)

# Argumenty:
- `n`: Numer stanu (n = 1, 2, 3, …).
- `m`: Masa cząstki.
- `planck`: ℏ.
- `k`: Stała sprężystości potencjału.
"""
function E_analytic(n::Int, m::Float64, planck::Float64, k::Float64)
    return planck * sqrt(k / m) * (n - 0.5)
end
