function [pL,pR,kappa,A,info] = get_p_opt( theta,beta,lambda,varargin )
%% Optimal CHAMP Robin weights p_L, p_R (order-2 / BDF2 scheme) from the
%% precomputed table: p = r/(1-r), with r from get_r_opt (same arguments and
%% optional outputs). p = Inf is the Dirichlet limit (r = 1), p = 0 Neumann;
%% NaN where no family is stable or the query is outside the table.
%%
%%   [pL,pR]               = get_p_opt( theta,beta,lambda )
%%   [pL,pR,kappa,A,info]  = get_p_opt( theta,beta,lambda,method,file )
%%
%% See get_r_opt for the inputs, the per-family interpolation and the outputs
%% kappa, A (tracked peaks) and info.

  [rL,rR,kappa,A,info] = get_r_opt( theta,beta,lambda,varargin{:} );
  pL = rL./(1-rL); pR = rR./(1-rR);
end
