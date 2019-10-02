function [imgR,count]=SENSE_CG_recon(sMaps,U1,kU)
% this is a Conjugate Gradient SENSE based on sensitivity map
% input sMaps: Sensitivity Map
%       U1: undersampling matrix
%       kU: undersampled k-space zero filled
%       opt: other iteration parameters

% 2013.04.22 Yi Guo
% simplified 01/10/2014

[kx,ky,coils]=size(kU);
x0=zeros(kx*ky,1);
b=conj(sMaps).*ifft2c(kU);
b=sum(b,3);
b=b(:);

count=1; % set initial parameters
xk=x0;
rk=b-compQ(x0);
pk=rk;
rk1=rk;
itermax=5;

while count<itermax && norm(rk1)>1e-6 % CG iteration
    % disp(['Residue:' num2str(norm(rk1))]);
    alpha=(rk'*rk)/(pk'*compQ(pk)); % modified CG equations
    xk1=xk+alpha*pk;
    rk1=rk-alpha*compQ(pk);
    beta=(rk1'*rk1)/(rk'*rk);
    
    pk1=rk1+beta*pk;
    count=count+1;
    
    xk=xk1; % update variable
    pk=pk1;
    rk=rk1;
    
    if (norm(rk) < 1e-6)
        disp('Iteration stops due to small residue');
    elseif count == itermax
        disp('Iteration stops due to max iteration');
    end
end

imgR=reshape(xk,[kx,ky]);

    function Qx=compQ(x) % calculate AL'*AL*x, eg Qx
        x_k=reshape(x,kx,ky);
        % calculate S'*DFT'*U'*U*DFT*S*x_k
        QxC=conj(sMaps).*ifft2c(U1.*fft2c(sMaps.*repmat(x_k,[1 1 coils])));
        Qx=sum(QxC,3);
        Qx=Qx(:);
    end
end
