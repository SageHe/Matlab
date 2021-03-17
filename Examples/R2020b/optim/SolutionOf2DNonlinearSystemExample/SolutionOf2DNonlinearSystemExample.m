%% Solution of 2-D Nonlinear System
% This example shows how to solve two nonlinear equations in two variables.
% The equations are
%
% $$ \begin{array}{c}
% {e^{ - {e^{ - ({x_1} + {x_2})}}}} = {x_2}\left( {1 + x_1^2} \right)\\
% {x_1}\cos \left( {{x_2}} \right) + {x_2}\sin \left( {{x_1}} \right) = \frac{1}{2}.
% \end{array} $$
%%
% Convert the equations to the form $F(x) = \bf{0}$.
%
% $$\begin{array}{c}
% {e^{ - {e^{ - ({x_1} + {x_2})}}}} - {x_2}\left( {1 + x_1^2} \right) = 0\\
% {x_1}\cos \left( {{x_2}} \right) + {x_2}\sin \left( {{x_1}} \right)
% - \frac{1}{2} = 0. \end{array} $$
%%
% Write a function that computes the left-hand side of these two equations.
%
% <include>root2d.m</include>
%
% Save this code as a file named |root2d.m| on your MATLAB(R) path.
%%
% Solve the system of equations starting at the point |[0,0]|.
fun = @root2d;
x0 = [0,0];
x = fsolve(fun,x0)

%% 
% Copyright 2012 The MathWorks, Inc.