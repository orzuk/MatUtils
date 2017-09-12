function DrawRLMMOneSNP(RLMM, x, mult_ratio, legend_flag, color_vec)


if(x == -1)
    M=(vec_into_mat(RLMM.MuMuVec, 2))'
    S{1} = [RLMM.SigmaMuVec(1) RLMM.SigmaMuVec(3); RLMM.SigmaMuVec(3) RLMM.SigmaMuVec(2)];
    S{2} = [RLMM.SigmaMuVec(4) RLMM.SigmaMuVec(6); RLMM.SigmaMuVec(6) RLMM.SigmaMuVec(5)];
    S{3} = [RLMM.SigmaMuVec(7) RLMM.SigmaMuVec(9); RLMM.SigmaMuVec(9) RLMM.SigmaMuVec(8)];
    color_vec = '---';
else
    M=[RLMM.MuMats.AA(x,:)' RLMM.MuMats.AB(x,:)' RLMM.MuMats.BB(x,:)']';
    S{1} = mult_ratio*mult_ratio*[RLMM.SigmaMats.AA(x,1), RLMM.SigmaMats.AA(x,3); RLMM.SigmaMats.AA(x,3), RLMM.SigmaMats.AA(x,2)];
    S{2} = mult_ratio*mult_ratio*[RLMM.SigmaMats.AB(x,1), RLMM.SigmaMats.AB(x,3); RLMM.SigmaMats.AB(x,3), RLMM.SigmaMats.AB(x,2)];
    S{3} = mult_ratio*mult_ratio*[RLMM.SigmaMats.BB(x,1), RLMM.SigmaMats.BB(x,3); RLMM.SigmaMats.BB(x,3), RLMM.SigmaMats.BB(x,2)];
end

if(nargin < 4)
    legend_flag = 1;
end
if(legend_flag)
    legends_vec = {'AA', 'AB', 'BB'};
else
    legends_vec = [];
end

axes_vec = [0 4 0 4];

if(nargin < 5)
    MixtureOfGaussiansDraw2dGaussians(M, S, {'A Int.', 'B Int.'}, legends_vec, [], axes_vec);
else
    MixtureOfGaussiansDraw2dGaussians(M, S, {'A Int.', 'B Int.'}, legends_vec, color_vec, axes_vec);
end
