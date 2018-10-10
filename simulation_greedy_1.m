function X = simulation_main

% Kajsa Mollersen (kajsa.mollersen@uit.no), October 9th 2018

% Requires: structure_matrix.m, cell_gene_effect.m, observation_matrix.m

% Creates a structure matrix, the cell and/or gene effects, and finally the
% observed binary gene expression matrix.

%% For the figures
% The size
n = 1000; 
d = 5000;

rng('default')

fig_nr = 2;
clf(figure(2))
IM = ones(n,d);
xlab = 'Cells';
ylab = 'Genes';

%% The structure matrix

C = cell(1,3);
C{1} = 1: n;
C{2} = 1: ceil(n/3);
C{3} = setdiff(C{1},C{2});

G = cell(1,3);
G{1} = 1: 1500;
G{2} = 1501: 2000;
G{3} = 2001: 2750;

S = structure_matrix(n,d,C,G);

figure(1), colormap(gray)
imagesc((IM - S)', [0 1])
title('S')
set(gca,'xaxisLocation','top')
xlabel(xlab)
ylabel(ylab)
drawnow

%% The Bernoulli parameter matrix. 
% Each entry gives defines the Bernoulli distribution from which x_{ij} is drawn. 

p_const = 0.8; % An overall probability of X_{ij} = s_{ij}

p = 0; % The cell and/or gene effect

block = 3;
distr = 'Normal';
param = [0 0.15];

p = cell_gene_effect(distr, param, block, p_const, d);

figure(fig_nr), subplot(1,3,3) , histogram(p,100), title(strcat('Truncated ',{' '}, distr))

%% The simulated observation matrix
X = observation_matrix(S, p_const, p, fig_nr, xlab, ylab);




