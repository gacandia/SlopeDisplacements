function[lambdaD]=BT2007_cdmM(Ts_, ky_, param,im,M,dPm,Haz50)

% Bray, J. D., & Travasarou, T. (2007). Simplified procedure for estimating
% earthquake-induced deviatoric slope displacements. Journal of Geotechnical
% and Geoenvironmental Engineering, 133(4), 381-392.

%%
d       = param.d;
realD   = param.realD;
Nd      = length(d);     % N° of displacement values of lambdaD
Nim     = length(im);    % N° of acceleration values of lambdaSa
dlambda = -[diff(Haz50,1);0];
pd      = makedist('Normal');
t       = truncate(pd,-2,2);
xrnd    = random(t,1,realD);
dPm(isnan(dPm))=0;

% ky and Ts realizations
Nky  = 100; realky = trlognpdf_psda(ky_,Nky);
NTs  = 100; realTs = trlognpdf_psda(Ts_,NTs);
[ky,Ts]  = meshgrid(realky,realTs);

ky       = ky(:)';
Ts       = Ts(:)';
Nmag     = length(M);
lnd      = zeros(Nim,Nmag,Nky*NTs);
Tlow     = Ts<0.05;
sigmaD   = 0.67;

for m=1:Nmag
    for j=1:Nim
        lnd(j,m,Tlow)  = -0.22-2.83*log(ky( Tlow))-0.333*(log(ky( Tlow))).^2+0.566*log(ky( Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts( Tlow)+ 0.278*(M(m)-7);
        lnd(j,m,~Tlow) = -1.10-2.83*log(ky(~Tlow))-0.333*(log(ky(~Tlow))).^2+0.566*log(ky(~Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts(~Tlow)+ 0.278*(M(m)-7);
    end
end

mean_d  = mean(lnd   , 3);
std_d   = std( lnd,[], 3);

logd    = log(d(:));

switch param.integration
    case 'MC'
        lambdaD = zeros(realD,Nd);
        for i=1:realD
            for j=1:Nd
                for m = 1:Nmag
                    PSA          = 1-normcdf((logd(j)-(mean_d(:,m) + std_d(:,m)*xrnd(i)))/sigmaD);
                    lambdaD(i,j) = lambdaD(i,j) + dPm(m,:)*(PSA.*dlambda);
                end
            end
        end
        
    case 'PC'
        PC0_array = zeros(Nim,Nd);
        PC1_array = zeros(Nim,Nd);
        PC2_array = zeros(Nim,Nd);
        PC3_array = zeros(Nim,Nd);
        PC4_array = zeros(Nim,Nd);
        
        PC0 = zeros(1, Nd);
        PC1 = zeros(1, Nd);
        PC2 = zeros(1, Nd);
        PC3 = zeros(1, Nd);
        PC4 = zeros(1, Nd);
        
        for m=1:length(M)
            for i = 1:Nim
                A = - std_d(i,m)^2/(2*sigmaD^2) - 1/2;
                B = (logd - mean_d(i,m))*std_d(i,m)/(sigmaD^2);
                C = -(logd - mean_d(i,m)).^2 * 1/(2*sigmaD^2);
                PC0_array(i,:) = 1/1 * (1 - normcdf((logd - mean_d(i,m))/(sqrt(sigmaD^2 + std_d(i)^2))))*dlambda(i);
                PC1_array(i,:) = 1/1 * (std_d(i,m)/(sigmaD *2*pi) .* (sqrt(pi)/sqrt(-A))* exp(C -B.^2/(4*A)))*dlambda(i);
                PC2_array(i,:) = 1/2 * (std_d(i,m)/(sigmaD *2*pi) .* (sqrt(pi)*B/(2*(-A)^(3/2))).* exp(C -B.^2/(4*A)))*dlambda(i);
                PC3_array(i,:) = 1/6 * (std_d(i,m)/(sigmaD *2*pi) .* (sqrt(pi)*(-2*A*(1 + 2*A) + B.^2)/(4*(-A)^(5/2))).* exp(C -B.^2/(4*A)))*dlambda(i);
                PC4_array(i,:) = 1/24* (std_d(i,m)/(sigmaD *2*pi) .* (sqrt(pi)*(-B).*(6*A*(1 + 2*A) - B.^2)/(8*(-A)^(7/2))).* exp(C -B.^2/(4*A)))*dlambda(i);
            end
            
            PC0 = PC0 + dPm(m,:)*PC0_array; % sums over the rows
            PC1 = PC1 + dPm(m,:)*PC1_array;
            PC2 = PC2 + dPm(m,:)*PC2_array;
            PC3 = PC3 + dPm(m,:)*PC3_array;
            PC4 = PC4 + dPm(m,:)*PC4_array;
        end
        
        HD      = [H(0,xrnd);H(1,xrnd); H(2,xrnd);H(3,xrnd);H(4,xrnd)];
        PC      = [PC0;PC1;PC2;PC3;PC4];
        lambdaD = HD'*PC;
        lambdaD(lambdaD<0)=0; %careful here. 
end



