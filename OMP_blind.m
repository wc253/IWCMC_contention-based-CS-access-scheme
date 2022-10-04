function [B_est, H_est] = OMP_blind(Y, P,SNR)
max_itr_OMP = max(size(P));
B = [];
Y_r = Y;
t1 = 1;
[~,n]=size(P);
H = zeros(n,size(Y,2));
H_est= zeros(n,size(Y,2));
norm_save(t1) = norm(Y_r,'fro');
err=5*10^(-SNR/10)*norm_save(1);
while 1
    t1 = t1 + 1;
    B_last = B;    
    [~,k] = max(sum(abs(P'*Y_r).^2,2));
    B = union(k,B);
    H(B,:) = P(:,B) \ Y;
    Y_r = Y - P(:,B)*H(B,:);
    norm_save(t1) = norm(Y_r,'fro');
    
    if norm_save(t1) < err
        break;
    end    
       
    if norm_save(t1)/norm_save(t1-1) >= 1 
        B = B_last';
        H(B,:) = pinv(P(:,B))*Y;
        break;
    end
    
    if t1 >max_itr_OMP
        break;
    end
    
end
erre=10^(-SNR/10);
[H_sort,pos]=sort(sum(abs(H).^2,2),'descend');
B_est=pos(H_sort>erre);
H_est(B_est,:)= P(:,B_est) \ Y; 
