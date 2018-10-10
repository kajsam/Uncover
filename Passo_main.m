function Passo_main(X,K,tau, ass, fig_nr)

% Input:    X - binary matrix of gene expression
%           K - maximum rank
%           tau - threshold

addpath('C:\Users\kajsam\Documents\MATLAB\mdl4bmf')
addpath('C:\Users\kajsam\Documents\MATLAB\Simulation')

[n,d] = size(X);
IM = ones(n,d);
if ass
  tau_vec = 0.1:0.1:1;
  % [k, t] = mdl4bmf(double(X), tau_vec,'errorMeasure', 'all'); % the actual work
 
  % K = k.DtMtypedXor 
  % tau = t.DtMtypedXor
  tic
  [W, H] = asso(double(X'), K, tau); % notice the transpose of A!
  toc
 
  W = W';
  H = H';
  As = logical(W*H);
  
  imK = max(K,3);

  figure(fig_nr), subplot(3,imK,1),imagesc(IM'-X'), colormap(gray), title('X')
  subplot(3,imK,2), imagesc(IM'-As'), title(strcat('Asso cover K = ',' ',num2str(K)))
  drawnow
  
  figure(fig_nr)
  for k = 1: K
    Ask = logical(W(:,k)*H(k,:));
    subplot(3,imK,imK+k),imagesc(IM'-Ask'), colormap(gray), xlabel(k)
    if k == 1
      ylabel('Asso components')
    end
  end
end

%% First time is always special

% Association matrix

tic
[Z, sum0] = passociation_matrix(X,X,tau);
toc
X0 = X(sum0,:);
n0 = size(X0,1);
mask = false(n0,d);

W = false(n0,K); 
H = false(K,d);
tic
for k = 1: K
  [w, h, Z] = select_column(X0,Z,mask);
  W(:,k) = w;
  
%   if k < K
%     mask_tmp = w*h;
%     left_ma = mask_tmp(:,h) & ~X(:,h);
%   
%     [w_l,h_l] = select_column(~X(:,h),~Z, false(size(left_ma)));
%     h_left = false(1,d);
%     h_left(h) = h_l;
%     h = h & ~h_left;
%     figure(4), subplot(1,2,1), imagesc(w_l*h_l)
%   end
  
  H(k,:) = h;
  mask = mask | w*h;
  % figure(4), subplot(1,2,2), imagesc(mask)
  
end  
A = false(n,d);
A(sum0,:) = logical(W*H);

toc

if ass && sum(sum(As == A))==n*d
  'equality!'
end

figure(fig_nr),subplot(3,imK,3), imagesc(IM'-A'), title(strcat('Uncover K = ',' ',num2str(K)))

figure(fig_nr)
for k = 1: K
  Ask = logical(W(:,k)*H(k,:));
  subplot(3,imK,2*imK+k),imagesc(IM'-Ask'), colormap(gray), xlabel(k)
  if k == 1
    ylabel('Uncover components')
  end
end


