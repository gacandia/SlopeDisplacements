function[lambdaD]=BMT2017_cdmM(Ts_, ky_, param,im,M,dPm,Haz50)
                                                     
% im, M,dPm,Cz,realSa,realD
% Bray, J. D., Macedo, J., & Travasarou, T. (2017). Simplified procedure 
% for estimating seismic slope displacements for subduction zone 
% earthquakes. Journal of Geotechnical and Geoenvironmental Engineering, 
% 144(3), 04017124.

%% 
d       = param.d;
realD   = param.realD;
Nd      = length(d);     % N° of displacement values of lambdaD
Nim     = length(im);    % N° of acceleration values of lambdaSa
dlambda = -[diff(Haz50,1);0];
pd      = makedist('Normal');
t       = truncate(pd,-2,2);
xrnd    = random(t,1,realD);

% ky and Ts realizations
realky = unique(trlognpdf_psda(ky_,100)); Nky=length(realky);
realTs = unique(trlognpdf_psda(Ts_,100)); NTs=length(realTs);
[ky,Ts]  = meshgrid(realky,realTs);

ky       = ky(:)';
Ts       = Ts(:)';
Nmag     = length(M);
lnd      = zeros(Nim,Nmag,Nky*NTs);
TL       = find(Ts < 0.1);
TH       = find(Ts>= 0.1);
sigmaD   = 0.73;

for m=1:Nmag
    for j=1:Nim
        lnd(j,m,TL) = -5.864-3.353*log(ky(TL))-0.390*(log(ky(TL))).^2+0.538*log(ky(TL))*log(im(j))+3.060*log(im(j))-0.225*(log(im(j))).^2-9.421*Ts(TL)+0.55*M(m);
        lnd(j,m,TH) = -6.896-3.353*log(ky(TH))-0.390*(log(ky(TH))).^2+0.538*log(ky(TH))*log(im(j))+3.060*log(im(j))-0.225*(log(im(j))).^2+3.081*Ts(TH)+0.55*M(m)-0.803*Ts(TH).^2;
    end
end

mean_d  = mean(lnd  , 3);
std_d   = std(lnd,[], 3);
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
end
