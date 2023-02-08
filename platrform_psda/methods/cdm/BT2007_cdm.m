function[lambdaD]=BT2007_cdm (Ts_, ky_, param, im,MRE,Cz,NMmin)

% Bray, J. D., & Travasarou, T. (2007). Simplified procedure for estimating
% earthquake-induced deviatoric slope displacements. Journal of Geotechnical
% and Geoenvironmental Engineering, 133(4), 381-392.

d      = param.d;
realSa = param.realSa;
realD  = param.realD;
Nd     = length(d);     % N° of displacement values of lambdaD
Nim    = length(im);    % N° of acceleration values of lambdaSa
pd     = makedist('Normal');
t      = truncate(pd,-2,2);  % changed here
xrnd   = random(t,1,realD);
zrnd   = random(t,1,realSa);

% ky and Ts realizations
Nky  = 100; realky = trlognpdf_psda(ky_,Nky);
NTs  = 100; realTs = trlognpdf_psda(Ts_,NTs);
[ky,Ts]  = meshgrid(realky,realTs);

ky       = ky(:)';
Ts       = Ts(:)';

lnd  = zeros(Nim,Nky*NTs);
Tlow = Ts<0.05;
for j=1:Nim
    lnd(j, Tlow) = -0.22-2.83*log(ky( Tlow))-0.333*(log(ky( Tlow))).^2+0.566*log(ky( Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts( Tlow);
    lnd(j,~Tlow) = -1.10-2.83*log(ky(~Tlow))-0.333*(log(ky(~Tlow))).^2+0.566*log(ky(~Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts(~Tlow);
end
mean_d  = mean(lnd, 2);
std_d   = std(lnd,[], 2);
sigmaD  = 0.67;

st      = sprintf('%s-%s',param.integration,param.hazard);
switch st
    case 'MC-full'
        MRE     = permute(MRE,[1 3 2]);
        lambdaD = zeros(realD*realSa, Nd);
        for i = 1:Nd
            for j = 1:realD
                xhat = (log(d(i)) - (mean_d + std_d* xrnd(j)))/sigmaD;
                Pd   = 0.5*(1-erf(xhat/sqrt(2)));
                IND  = (1:realSa)+realSa*(j-1);
                lambdaD(IND,i) = trapz(Pd,MRE)';   % simpson rule
            end
        end
        
    case 'PC-full'
        lambdaD   = zeros(realD*realSa, Nd);
        dlambdaPC = - NMmin*[diff(Cz,1);zeros(1,5)]';
        HD        = [H(0,xrnd);H(1,xrnd); H(2,xrnd);H(3,xrnd);H(4,xrnd)];
        HSa       = [H(0,zrnd);H(1,zrnd); H(2,zrnd);H(3,zrnd);H(4,zrnd)];
        PC0       = zeros(Nd,5,Nim);
        PC1       = zeros(Nd,5,Nim);
        PC2       = zeros(Nd,5,Nim);
        PC3       = zeros(Nd,5,Nim);
        PC4       = zeros(Nd,5,Nim);
        logd      = log(d(:));
        for k = 1:5
            for i = 1:Nim
                A_s = - std_d(i)^2/(2*sigmaD^2) - 1/2;
                B_s = (logd - mean_d(i))*std_d(i)/(sigmaD^2);
                C_s = -(logd - mean_d(i)).^2 * 1/(2*sigmaD^2);
                PC0(:,k,i) =  1/1 * (1 - normcdf((logd - mean_d(i))/(sqrt(sigmaD^2 + std_d(i)^2))))*dlambdaPC(k, i);
                PC1(:,k,i) =  1/1 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)/sqrt(-A_s))* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(k, i);
                PC2(:,k,i) =  1/2 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*B_s/(2*(-A_s)^(3/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(k, i);
                PC3(:,k,i) =  1/6 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*(-2*A_s*(1 + 2*A_s) + B_s.^2)/(4*(-A_s)^(5/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(k, i);
                PC4(:,k,i) = 1/24 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*(-B_s).*(6*A_s*(1 + 2*A_s) - B_s.^2)/(8*(-A_s)^(7/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(k, i);
            end
        end % End loop over PC terms Hazard!!
        
        PC0 = sum(PC0,3);
        PC1 = sum(PC1,3);
        PC2 = sum(PC2,3);
        PC3 = sum(PC3,3);
        PC4 = sum(PC4,3);
        for i = 1:Nd
            cont = 1;
            CPC = [PC0(i,:);PC1(i,:);PC2(i,:);PC3(i,:);PC4(i,:)];
            for j = 1:realD
                for k = 1:realSa
                    lambdaD(cont, i) = HD(:,j)'*CPC*HSa(:,k);
                    cont = cont + 1;
                end
            end
        end

    case 'MC-average'
        lambdaD = zeros(realD, Nd);
        Percent = 50;
        MRE     = permute(MRE,[1 3 2]);
        MRE     = prctile(MRE,Percent,2);
        
        for i = 1:Nd
            for j = 1:realD
                xhat = (log(d(i)) - (mean_d + std_d* xrnd(j)))/sigmaD;
                Pd   = 0.5*(1-erf(xhat/sqrt(2)));
                lambdaD(j,i) = trapz(Pd,MRE);   % simpson rule
            end
        end
        
    case 'PC-average'
        Haz50      = prctile(MRE,50,3);
        dlambdaPCP = -[diff(Haz50,1);0]';
        HD         = [H(0,xrnd);H(1,xrnd); H(2,xrnd);H(3,xrnd);H(4,xrnd)];
        
        PC0 = zeros(Nim,Nd);
        PC1 = zeros(Nim,Nd);
        PC2 = zeros(Nim,Nd);
        PC3 = zeros(Nim,Nd);
        PC4 = zeros(Nim,Nd);
        logd = log(d(:));
        
        for i = 1:Nim
            A_s = - std_d(i)^2/(2*sigmaD^2) - 1/2;
            B_s = (logd - mean_d(i))*std_d(i)/(sigmaD^2);
            C_s = -(logd - mean_d(i)).^2 * 1/(2*sigmaD^2);
            PC0(i,:) =  1/1 * (1 - normcdf((logd - mean_d(i))/(sqrt(sigmaD^2 + std_d(i)^2))))*dlambdaPCP(i);
            PC1(i,:) =  1/1 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)/sqrt(-A_s))* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPCP(i);
            PC2(i,:) =  1/2 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*B_s/(2*(-A_s)^(3/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPCP(i);
            PC3(i,:) =  1/6 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*(-2*A_s*(1 + 2*A_s) + B_s.^2)/(4*(-A_s)^(5/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPCP(i);
            PC4(i,:) = 1/24 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*(-B_s).*(6*A_s*(1 + 2*A_s) - B_s.^2)/(8*(-A_s)^(7/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPCP(i);
        end
        
        PC0 = sum(PC0,1);
        PC1 = sum(PC1,1);
        PC2 = sum(PC2,1);
        PC3 = sum(PC3,1);
        PC4 = sum(PC4,1);
        
        PC      = [PC0;PC1;PC2;PC3;PC4];
        lambdaD = HD'*PC;
end
