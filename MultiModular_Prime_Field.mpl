

Coeffs_from_zero_to_p_minus_one := proc(f, p, x) local L, i; L := []; for i from 0 to p - 1 do L := [op(L), coeff(f, 'x', i)]; end do; return L; end proc;
Kth_element_in_each_list := proc(L, k) local len, outputL, i; len := numelems(L); outputL := []; for i to len do outputL := [op(outputL), L[i][k]]; end do; return outputL; end proc;
MultidimentionFFT := proc(f, m, p, X) local i, j, x, fcoeffs, fValues, f_i, g, outputL; with(PolynomialTools); if X = [] then return [seq(f), i = 1 .. p^m]; end if; outputL := []; x := X[1]; if m = 1 then for i from 0 to p - 1 do outputL := [op(outputL), eval(f, x = i) mod p]; end do; return outputL; else fcoeffs := Coeffs_from_zero_to_p_minus_one(f, p, x); fValues := []; for i to numelems(fcoeffs) do f_i := fcoeffs[i]; fValues := [op(fValues), MultidimentionFFT(f_i, m - 1, p, X[2 .. -1])]; end do; for i to p^(numelems(X) - 1) do g := FromCoefficientList(Kth_element_in_each_list(fValues, i), x); for j from 0 to p - 1 do outputL := [op(outputL), eval(g, x = j) mod p]; end do; end do; return outputL; end if; end proc;
MultiModular_Prime_Field := proc(f, p, alpha) local X, m, x, M, fBar, L, outputL, entry, i, j, alpha_i, alpha_i_j; with(PolynomialTools); if f mod p = 0 then return [seq(0, i = 1 .. numelems(alpha))]; end if; X := convert(indets(f), list); m := numelems(X); fBar := f; for i to m do x := X[i]; fBar := rem(fBar, x^p - x, x) mod p; end do; L := MultidimentionFFT(fBar, m, p, X); outputL := []; for i to numelems(alpha) do alpha_i := alpha[i]; entry := 1; for j to m do alpha_i_j := alpha_i[j] mod p; entry := entry + alpha_i_j*p^(j - 1); end do; outputL := [op(outputL), L[entry]]; end do; return outputL; end proc;
