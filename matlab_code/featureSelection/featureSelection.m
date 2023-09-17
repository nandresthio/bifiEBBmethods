% -------------------------------------------------------------------------
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------

% Copyright 2023, Nicolau Andrés-Thió

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% The original creator is Mario Andres Munoz Acosta. To see repository
% from which this code was obtained, visit https://github.com/andremun/InstanceSpace

% -------------------------------------------------------------------------

rootdir = './';

opts.sifted.flag = true;            % Automatic feature selectio on. Set to false if you don't want any feature selection.
opts.sifted.rho = 0.3;              % Minimum correlation value acceptable between performance and a feature. Between 0 and 1
opts.sifted.K = 9;                 % Number of final features. Ideally less than 10.
opts.sifted.NTREES = 50;            % Number of trees for the Random Forest (to determine highest separability in the 2-d projection)
opts.sifted.MaxIter = 1000;
opts.sifted.Replicates = 100;

% Saving all the information as a JSON file
fid = fopen([rootdir 'options.json'],'w+');
fprintf(fid,'%s',jsonencode(opts));
fclose(fid);

try
    model = buildIS(rootdir);
catch ME
    disp('EOF:ERROR');
    rethrow(ME)
end