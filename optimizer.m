function[f_hat,mincut] = optimizer( mode, partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j, lambda)

%function: to sample two labels from the two sets randomly
%Author:Shaona Ghosh
%Date:10.01.2014


%Check for number of arguments
%error(nargchk(1, 10, nargin)');
error(nargchk(1, 11, nargin)');

switch lower(mode)
    case 'harmonicinv'
      %Harminic energy minimization by the inverse method
      f = harmonicinverse(partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j);
    case 'harmonicquad'
      %Energy minimization using quadratic programming
      f = quadraticprog(partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j);
    case 'maxflowlinear'
      [f,mincut] = linearprog(partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j);
    case 'harplusmaxflowcons'
       [f, mincut] = linearprog(partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j);
       f =  quadraticprogflowcons(partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j, mincut);
      %f =  quadraticprogflowcons(partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j);
    case 'harplusmaxflowobj'
      %Energy minimization with max flow as the linear part of the
      %objective
      f = quadraticprogflowobj(partlab1, partlab2, nodelabi, nodelabj, A, D, L, num, f_j, lambda);
    case 'maxflowFord'
      %TODO
end
% Prediction
f_hat = sign(f);
idzeros = find(f_hat == 0);
if ~isempty(idzeros)
    f_hat(idzeros) = 1; %+1 label for evens
end
end





