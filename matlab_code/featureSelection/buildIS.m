function model = buildIS(rootdir)
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


startProcess = tic;
scriptdisc('buildIS.m');
% -------------------------------------------------------------------------
% Collect all the data from the files
disp(['Root Directory: ' rootdir]);
datafile = [rootdir 'metadata.csv'];
optsfile = [rootdir 'options.json'];
if ~isfile(datafile) || ~isfile(optsfile)
    error(['Please place the datafiles in the directory ''' rootdir '''']);
end
opts = jsondecode(fileread(optsfile));
disp('-------------------------------------------------------------------------');
disp('-> Listing options to be used:');
optfields = fieldnames(opts);
for i = 1:length(optfields)
    disp(optfields{i});
    disp(opts.(optfields{i}));
end

disp('-------------------------------------------------------------------------');
disp('-> Loading the data.');
Xbar = readtable(datafile);

varlabels = Xbar.Properties.VariableNames;
isname = strcmpi(varlabels,'instances');
isfeat = strncmpi(varlabels,'feature_',8);
isalgo = strncmpi(varlabels,'algo_',5);
issource = strcmpi(varlabels,'source');

model.data.instlabels = Xbar{:,isname};
if isnumeric(model.data.instlabels)
    model.data.instlabels = num2cell(model.data.instlabels);
    model.data.instlabels = cellfun(@(x) num2str(x),model.data.instlabels,'UniformOutput',false);
end
if any(issource)
    model.data.S = categorical(Xbar{:,issource});
end
model.data.X = Xbar{:,isfeat};
model.data.Y = Xbar{:,isalgo};

model.data.featlabels = varlabels(isfeat);
model.data.algolabels = varlabels(isalgo);



% -------------------------------------------------------------------------
% Removing the template data such that it can be used in the labels of
% graphs and figures.
model.data.featlabels = strrep(model.data.featlabels,'feature_','');
model.data.algolabels = strrep(model.data.algolabels,'algo_','');
% -------------------------------------------------------------------------
% Find best algorithms and binary performance
[model.data.Ybest,~] = max(model.data.Y,[],2);
model.data.Ybin = model.data.Y >= 0.5;

disp('=========================================================================');
disp('-> Calling SIFTED for auto-feature selection.');
disp('=========================================================================');
[model.data.X, model.sifted] = SIFTED(model.data.X, model.data.Y, model.data.Ybin, opts.sifted, model.data.featlabels);
model.data.featlabels = model.data.featlabels(model.sifted.selvars);
writecell(model.data.featlabels, "features")

% -------------------------------------------------------------------------
disp(['-> Completed! Elapsed time: ' num2str(toc(startProcess)) 's']);
disp('EOF:SUCCESS');
end
