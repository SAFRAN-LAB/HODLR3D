clc;
clearvars;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code finds the singular values of various interaction matrices
% INPUT: 'choice' that decides what kernel is to be used
% OUTPUT: singular values of the face sharing interaction, edge sharing
% interaction, vertex sharing interaction, well-separated interaction.
% The OUTPUT is written to file "output_file_%d_%d.mat". The first argument
% is 'N' (the size of the matrix) and the second argument is 'choice'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmt = 'output_file_%d_%d.mat';

choice = input('enter 1 for 1/r, 2 for 1/r^4, 3 for cos(r)/r: ');

% n = [5,10,15,20,25,30,35,39]; % vector of n's. The matrix size for which the rank is to be computed is n^3. 
n = 5;

L = 1;
for nnum = 1:length(n)
N = n(nnum)*n(nnum)*n(nnum);
fprintf('\n')
fprintf('==================================================================\n')
fprintf('N = %d \n',N);

[Txs,Tys,Tzs] = gen_points(0,0,0,L,n(nnum),n(nnum),n(nnum)); % SOURCE Points
[Txv,Tyv,Tzv] = gen_points(2*L,2*L,2*L,L,n(nnum),n(nnum),n(nnum)); % DESTINATION points (Vertex Sharing)
[Txf,Tyf,Tzf] = gen_points(2*L,0,0,L,n(nnum),n(nnum),n(nnum)); % DESTINATION points (Face Sharing)
[Txe,Tye,Tze] = gen_points(2*L,2*L,0,L,n(nnum),n(nnum),n(nnum)); % DESTINATION points (Edge Sharing)
[Txw,Tyw,Tzw] = gen_points(4*L,0,0,L,n(nnum),n(nnum),n(nnum)); % DESTINATION points (Well seperated, One box away)



%%%%%    Vertex Shared domain Intraction %%%%%
Txd = Txv;
Tyd = Tyv;
Tzd = Tzv;
Kxy = zeros(N,N);
fprintf('\n')
disp('Vertex Sharing')
tic
% r^2
for i = 1:N
    Kxy(:,i)     =  (Txs - Txd(i)).^2 + (Tys - Tyd(i)).^2 + (Tzs - Tzd(i)).^2;
end
% r
Kxy = sqrt(Kxy);
fprintf('Interaction matrix formed in\n')
toc

switch choice
    case 1
        svd_v = svd(1./Kxy);
    case 2
        svd_v = svd(1./Kxy.^4);
    case 3
        svd_v = svd(cos(Kxy)./Kxy);
    otherwise
        return;
end

fprintf('singular values of Vertex sharing done in\n')
toc

%%%%%    Edge Shared domain Intraction %%%%%
Txd = Txe;
Tyd = Tye;
Tzd = Tze;
Kxy = zeros(N,N);
fprintf('\n')
disp('Edge Sharing')
tic
% r^2
for i = 1:N
    Kxy(:,i)     =  (Txs - Txd(i)).^2 + (Tys - Tyd(i)).^2 + (Tzs - Tzd(i)).^2;
end
% r
Kxy = sqrt(Kxy);
fprintf('Interaction matrix formed in\n')
toc
switch choice
    case 1
        svd_e = svd(1./Kxy);
    case 2
        svd_e = svd(1./Kxy.^4);
    case 3
        svd_e = svd(cos(Kxy)./Kxy);
    otherwise
        return;
end

fprintf('singular values of Edge sharing done in\n')
toc

%%%%%    Face Shared domain Intraction %%%%%
Txd = Txf;
Tyd = Tyf;
Tzd = Tzf;
Kxy = zeros(N,N);
fprintf('\n')
disp('Face Sharing')
tic
% r^2
for i = 1:N
    Kxy(:,i)     =  (Txs - Txd(i)).^2 + (Tys - Tyd(i)).^2 + (Tzs - Tzd(i)).^2;
end
% r
Kxy = sqrt(Kxy);
fprintf('Interaction matrix formed in\n')
toc

switch choice
    case 1
        svd_f = svd(1./Kxy);
    case 2
        svd_f = svd(1./Kxy.^4);
    case 3
        svd_f = svd(cos(Kxy)./Kxy);
    otherwise
        return;
end

fprintf('singular values of Face sharing done in\n')
toc

%%%%%    Well seperated domain Intraction %%%%%
Txd = Txw;
Tyd = Tyw;
Tzd = Tzw;
Kxy = zeros(N,N);
fprintf('\n')
disp('Vertex Sharing')
tic
% r^2
for i = 1:N
    Kxy(:,i)     =  (Txs - Txd(i)).^2 + (Tys - Tyd(i)).^2 + (Tzs - Tzd(i)).^2;
end
% r
Kxy = sqrt(Kxy);
fprintf('Interaction matrix formed in\n')
toc
switch choice
    case 1
        svd_w = svd(1./Kxy);
    case 2
        svd_w = svd(1./Kxy.^4);
    case 3
        svd_w = svd(cos(Kxy)./Kxy);
    otherwise
        return;
end
fprintf('singular values of Well seperated done in\n')
toc
fprintf('\n')

fname = sprintf(fmt,N,choice)
save(fname,'svd_v','svd_e','svd_f','svd_w','N')
end
