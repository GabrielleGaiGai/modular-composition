
currentdir(worksheetdir);
read "MultiModular_Prime_Field.mpl";

Primes_Less_than_N := proc(N) local i; with(NumberTheory); return [seq(ithprime(i), i = 1 .. PrimeCounting(N))]; end proc;
compute_d := proc(f, X) local d, L, i; L := []; for i to numelems(X) do L := [op(L), degree(f, X[i])]; end do; return max(L) + 1; end proc;

MultiModular_Z_rZ := proc(f, alpha, r, t) local X, m, d, N, fBar, `&alpha;Bar`, l, primes, k, L_f, L_alpha, alpha_i, i, j, h, f_alpha, P; X := convert(indets(f), list); m := numelems(X); d := compute_d(f, X); N := numelems(alpha); if f = 0 then return [seq(0, i = 1 .. numelems(alpha))]; end if; fBar := f mod r; `&alpha;Bar` := alpha mod r; l := ceil(16*log(d^m*(r - 1)^(d*m))); primes := Primes_Less_than_N(l); k := numelems(primes); L_f := []; for h to k do L_f := [op(L_f), fBar mod primes[h]]; end do; L_alpha := []; for i to k do L_alpha := [op(L_alpha), `&alpha;Bar` mod primes[i]]; end do; f_alpha := []; if t = 1 then for h to k do f_alpha := [op(f_alpha), MultiModular_Prime_Field(L_f[h], primes[h], L_alpha[h])]; end do; else for h to k do f_alpha := [op(f_alpha), MultiModular_Z_rZ(L_f[h], L_alpha[h], primes[h], t - 1)]; end do; end if; P := chrem(f_alpha, primes); return P mod r; end proc;

