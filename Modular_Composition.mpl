
currentdir(worksheetdir);
read "MultiModular_Extension_Ring.mpl";
Select_N_Distinct_Elements := proc(N, r, E) local total, returnL, i, j, newL, newElementL, element, deg; total := r^degree(E, Z); if total < N then return false; else returnL := [seq(i, i = 0 .. r - 1)]; while numelems(returnL) < N + 1 do deg := degree(returnL[-1], Z) + 1; newL := []; for i to numelems(returnL) do element := returnL[i]; newElementL := []; for j from 0 to r - 1 do newElementL := [op(newElementL), element + j*Z^deg]; end do; newL := [op(newL), op(newElementL)]; end do; returnL := newL; end do; return returnL[2 .. N + 1]; end if; end proc;
Kronecker_SubStitution_Inverse := proc(f, h, l) local g, deg, a, coefficient, i, j, Y_term, Y_list; g := 0; deg := degree(f, X); Y_list := [seq(Y[i], i = 0 .. l - 1)]; for i to deg do a := convert(i, base, h); Y_term := 1; for j to numelems(a) do Y_term := Y_term*Y_list[j]^a[j]; end do; g := g + coeff(f, X, i)*Y_term; end do; return g + coeff(f, X, 0); end proc;
Modular_Composition := proc(f, g, h, E, r) local fPrime, d, d0, N, l, mPrime, NPrime, gPrime, generate_d0, Y_list, g_list, betas, alpha_i_k, alpha_i, f_alpha_list, i, k, interpolation, data; d := degree(f, X) + 1; generate_d0 := rand(2 .. d - 1); d0 := generate_d0(); N := max(degree(g, X), degree(h, X)) + 1; l := ceil(log[d0](d)); mPrime := l; Y_list := [seq(Y[i], i = 0 .. l - 1)]; fPrime := Kronecker_SubStitution_Inverse(f, d0, l); fPrime := rem(fPrime, E, Z) mod r; g_list := []; for i from 0 to l - 1 do if Y_list[i + 1] in indets(fPrime) then g_list := [op(g_list), rem(rem(g^(d0^i), h, X), E, Z) mod r]; end if; end do; if Z in indets(fPrime) then NPrime := N*(numelems(indets(fPrime)) - 1)*d0; else NPrime := N*numelems(indets(fPrime))*d0; end if; betas := Select_N_Distinct_Elements(NPrime, r, E); alpha_i_k := []; for k to NPrime do alpha_i := []; for i to numelems(g_list) do alpha_i := [op(alpha_i), eval(g_list[i], X = betas[k])]; end do; alpha_i_k := [op(alpha_i_k), alpha_i]; end do; f_alpha_list := MultiModular_Extension_Ring(fPrime, alpha_i_k, 1, E, r); interpolation := []; for i to NPrime do data := [betas[i], f_alpha_list[i]]; interpolation := [op(interpolation), data]; end do; return rem(rem(CurveFitting[PolynomialInterpolation](interpolation, X), h, X), E, Z) mod r; end proc;
