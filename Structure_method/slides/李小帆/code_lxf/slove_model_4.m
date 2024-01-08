%this function is to solve model-3
%multiple regios
%multiple sectors with  IO-linkage
%without migration
%%
function eqlm=slove_model_4(data);

v2struct(data);

%initial guess of w_i
w_n=ones(N,1);
P_nj=ones(N,J);
X_nj=ones(N,J);

Ga_j=gamma((theta_j+1-Sigma_j)./theta_j).^(1./(1-Sigma_j));

%%%%%%
%iteration on wage w_i
%%%%%%
diff3=1;
while diff3>1e-7
    
    %%%%%%
    %iteration on price X_i
    %%%%%%
    diff2=1;
    while diff2>1e-7
        
        %%%%%%
        %iteration on price X_i
        %%%%%%
        diff1=1;
        while diff1>1e-7
            %%%calculate price (given wage and price guess)
            P_nk=P_nj;
            c_nj=w_n.^gamma_nj.*permute(prod(P_nk.^gamma_nkj,2),[1,3,2]); %unit cost
            c_ij=c_nj;
            h_nij=(reshape(c_ij,[1,N,J]).*d_nij).^reshape(-theta_j,[1,1,J]);
            P_nj_new=Ga_j'.*permute(sum(reshape(T_ij,[1,N,J]).*h_nij,2),[1,3,2]).^(-1./theta_j');
            diff1=max(max(abs(P_nj-P_nj_new)));
            P_nj=P_nj_new;
        end
        %trade share
        pi_nij=(reshape(T_ij,[1,N,J]).*h_nij)./sum(reshape(T_ij,[1,N,J]).*h_nij,2);
        
        %implied expenditure X_i
        Alpha_ij=Alpha_nj;
        gamma_ikj=gamma_nkj;
        w_i=w_n;
        v_ikj=gamma_ikj.*permute(sum(pi_nij.*reshape(X_nj,[N,1,J])),[2,1,3]);
        v_ik=permute(sum(permute(v_ikj,[3,1,2])),[2,3,1]);
        v_ij=v_ik;
        X_ij=Alpha_ij.*w_i.*L_i+v_ij;
        X_nj_new=X_ij;
        diff2=max(max(abs(X_nj-X_nj_new)));
        X_nj=0.5*X_nj+0.5.*X_nj.*(X_nj_new./X_nj).^0.5;
    end
    expend_i=w_i.*L_i;
    gamma_ij=gamma_nj;
    income_i=sum(gamma_ij.*permute(sum(pi_nij.*reshape(X_nj,[N,1,J])),[2,3,1]),2);
    diff3=max(abs(income_i-expend_i));
    w_i=0.5*w_i+0.5*w_i.*(income_i./expend_i).^(0.5);
    w_i=w_i/w_i(1); %normalization
    w_n=w_i;
end









eqlm=v2struct;
end


