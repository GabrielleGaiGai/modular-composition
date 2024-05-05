
currentdir(worksheetdir);
read "MultiModular_Z_rZ.mpl";

alpha_mod_poly := proc(alpha, quo, Z, m) local outputL, N, i, j, `&alpha;_i`; N := numelems(alpha); outputL := []; for i to N do `&alpha;_i` := []; for j to m do `&alpha;_i` := [op(`&alpha;_i`), rem(alpha[i][j], quo, Z)]; end do; outputL := [op(outputL), `&alpha;_i`]; end do; return outputL; end proc;
MultiModular_Extension_Ring := proc(f, alpha, t, E, r) local e, d, X, m, N, M, rPrime, fBar, `&alpha;Bar`, fTilde, `&alpha;Tilde`, betas, Q, result, i; X := indets(f); if 'Z' in X then X := convert(X minus {Z}, list); else X := convert(X, list); end if; e := degree(E); d := compute_d(f, X); m := numelems(X); N := numelems(alpha); M := d^m*(e*(r - 1))^((d - 1)*m + 1) + 1; rPrime := M^((e - 1)*d*m + 1); fTilde := f mod r; `&alpha;Tilde` := alpha mod r; fTilde := rem(fTilde, E, Z); `&alpha;Tilde` := alpha_mod_poly(`&alpha;Tilde`, E, Z, m); fBar := fTilde mod rPrime; `&alpha;Bar` := `&alpha;Tilde` mod rPrime; fBar := rem(fBar, Z - M, Z); `&alpha;Bar` := alpha_mod_poly(`&alpha;Bar`, Z - M, Z, m); betas := MultiModular_Z_rZ(fBar, `&alpha;Bar`, rPrime, t); result := []; for i to N do Q := FromCoefficientList(convert(betas[i], base, M), Z); result := [op(result), rem(Q, E, Z) mod r]; end do; return result; end proc;

