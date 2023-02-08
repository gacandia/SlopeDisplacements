function[lambdaD]=BT2007_cdm_benchmark(Ts_param, ky_param, psda_param, im,HAZ,Cz,varargin)

% Bray, J. D., & Travasarou, T. (2007). Simplified procedure for estimating
% earthquake-induced deviatoric slope displacements. Journal of Geotechnical
% and Geoenvironmental Engineering, 133(4), 381-392.

%%
d      = psda_param.d;
realSa = psda_param.realSa;
realD  = psda_param.realD;
Nd     = length(d);     % N° of displacement values of lambdaD
Nim    = length(im);    % N° of acceleration values of lambdaSa

if nargin==6
    pd    = makedist('Normal');
    t     = truncate(pd,-2,2);
    xrnd  = random(t,1,realD);
    zrnd  = random(t,1,realSa);
    irule = 'simpson';
else
    xrnd  = varargin{1}; realD  = length(xrnd);
    zrnd  = varargin{2}; realSa = length(zrnd);
    irule = varargin{3};
end

if length(Ts_param)==2, Ts_param = [Ts_param,100]; end; NTs = Ts_param(3);
if length(ky_param)==2, ky_param = [ky_param,100]; end; Nky = ky_param(3);

[ky,Ts]  = meshgrid(trlognpdf_psda(ky_param),trlognpdf_psda(Ts_param));
ky       = ky(:)';
Ts       = Ts(:)';

if strcmp(psda_param.imhazard,'average')
    Percent = 50;
    HAZ     = repmat(prctile(HAZ,Percent,2),1,realSa);
end

if strcmp(irule,'rectangular')
    dHAZ  = diff(HAZ,1); dHAZ(end+1,:)=dHAZ(end,:);
    dHAZ  = dHAZ';
end

switch psda_param.method
    case 'MC'
        lnd  = zeros(Nim,Nky*NTs);
        Tlow = Ts<0.05;
        for j=1:Nim
            lnd(j, Tlow) = -0.22-2.83*log(ky( Tlow))-0.333*(log(ky( Tlow))).^2+0.566*log(ky( Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts( Tlow);
            lnd(j,~Tlow) = -1.10-2.83*log(ky(~Tlow))-0.333*(log(ky(~Tlow))).^2+0.566*log(ky(~Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts(~Tlow);
        end
        sigmaD  = 0.67;
        mean_d  = mean(lnd, 2);
        std_d   = std(lnd,[], 2);
        lambdaD = zeros(realD * realSa, length(d));
        
        switch irule
            case 'simpson'
                for i = 1:Nd
                    for j = 1:realD
                        xhat = (log(d(i)) - (mean_d + std_d* xrnd(j)))/sigmaD;
                        Pd   = 0.5*(1-erf(xhat/sqrt(2)));
                        IND  = (1:realSa)+realSa*(j-1);
                        lambdaD(IND,i) = trapz(Pd,HAZ)';   % simpson rule
                    end
                end
                
            case 'rectangular'
                for i = 1:Nd
                    for j = 1:realD
                        xhat = (log(d(i)) - (mean_d + std_d* xrnd(j)))/sigmaD;
                        Pd   = 0.5*(1-erf(xhat/sqrt(2)));
                        IND  = (1:realSa)+realSa*(j-1);
                        lambdaD(IND,i) = -dHAZ*Pd;   % rectangular rule
                    end
                end
        end
        
    case 'PC'
        lnd  = zeros(Nim,Nky*NTs);
        Tlow = Ts<0.05;
        for j=1:Nim
            lnd(j, Tlow) = -0.22-2.83*log(ky( Tlow))-0.333*(log(ky( Tlow))).^2+0.566*log(ky( Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts( Tlow);
            lnd(j,~Tlow) = -1.10-2.83*log(ky(~Tlow))-0.333*(log(ky(~Tlow))).^2+0.566*log(ky(~Tlow))*log(im(j))+3.04*log(im(j))-0.244*(log(im(j))).^2+1.5*Ts(~Tlow);
        end
        mean_d  = mean(lnd, 2);
        std_d   = std(lnd,[], 2);
        sigmaD  = 0.67;
        
        PC_term_0_GM = Cz(1:Nim, 1)';
        PC_term_1_GM = Cz(1:Nim, 2)';
        PC_term_2_GM = Cz(1:Nim, 3)';
        PC_term_3_GM = Cz(1:Nim, 4)';
        PC_term_4_GM = Cz(1:Nim, 5)';
        
        Hazard_PC_samples = zeros(realSa, Nim);
        
        for i = 1:Nim
            Hazard_PC_samples(:, i) = PC_term_0_GM(i) + ...
                PC_term_1_GM(i) * zrnd + ...
                PC_term_2_GM(i) * (zrnd.^2-1) + ...
                PC_term_3_GM(i) * (zrnd.^3 - 3*zrnd) + ...
                PC_term_4_GM(i) * (zrnd.^4 - 6*zrnd.^2 + 3);
        end
        
        dlambdaPC = zeros(5, Nim);
        dlambdaPC(1, 1:end-1) = diff(-PC_term_0_GM(1, :),1,2);
        dlambdaPC(2, 1:end-1) = diff(-PC_term_1_GM(1, :),1,2);
        dlambdaPC(3, 1:end-1) = diff(-PC_term_2_GM(1, :),1,2);
        dlambdaPC(4, 1:end-1) = diff(-PC_term_3_GM(1, :),1,2);
        dlambdaPC(5, 1:end-1) = diff(-PC_term_4_GM(1, :),1,2);
        dlambdaPC(:,end) = dlambdaPC(:,end);
        
        PC_term_0_array = zeros(Nd,5,Nim);
        PC_term_1_array = zeros(Nd,5,Nim);
        PC_term_2_array = zeros(Nd,5,Nim);
        PC_term_3_array = zeros(Nd,5,Nim);
        PC_term_4_array = zeros(Nd,5,Nim);
        logd = log(d(:));
        for k = 1:5
            for i = 1:Nim
                A_s = - std_d(i)^2/(2*sigmaD^2) - 1/2;
                B_s = (logd - mean_d(i))*std_d(i)/(sigmaD^2);
                C_s = -(logd - mean_d(i)).^2 * 1/(2*sigmaD^2);
                PC_term_0_array(:,k,i) = 1/1 * (1 - normcdf((logd - mean_d(i))/(sqrt(sigmaD^2 + std_d(i)^2))))*dlambdaPC(k, i);
                PC_term_1_array(:,k,i) = 1/1 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)/sqrt(-A_s))* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(k, i);
                PC_term_2_array(:,k,i) = 1/2 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*B_s/(2*(-A_s)^(3/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(k, i);
                PC_term_3_array(:,k,i) = 1/6 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*(-2*A_s*(1 + 2*A_s) + B_s.^2)/(4*(-A_s)^(5/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(k, i);
                PC_term_4_array(:,k,i) = 1/24 * (std_d(i)/(sigmaD *2*pi) .* (sqrt(pi)*(-B_s).*(6*A_s*(1 + 2*A_s) - B_s.^2)/(8*(-A_s)^(7/2))).* exp(C_s -B_s.^2/(4*A_s)))*dlambdaPC(k, i);
            end
            
        end % End loop over PC terms Hazard!!
        PC0 = sum(PC_term_0_array,3);
        PC1 = sum(PC_term_1_array,3);
        PC2 = sum(PC_term_2_array,3);
        PC3 = sum(PC_term_3_array,3);
        PC4 = sum(PC_term_4_array,3);
        
        HD      = zeros(5,realD);
        HSa     = zeros(5,realSa);
        for i=1:5
            HD (i,:)= H(i-1,xrnd);
            HSa(i,:)= H(i-1,zrnd');
        end
        
        lambdaD = zeros(realD*realSa, Nd);
        for i = 1:Nd
            cont = 1;
            for j = 1:realD
                HD1=HD(1,j);
                HD2=HD(2,j);
                HD3=HD(3,j);
                HD4=HD(4,j);
                HD5=HD(5,j);
                for k = 1:realSa
                    for l = 1:5 % Loop over PC terms hazard
                        Hlk = HSa(l,k);
                        lambdaD(cont, i) =lambdaD(cont, i)+Hlk*(PC0(i, l)*HD1+PC1(i,l)*HD2+PC2(i,l)*HD3+PC3(i,l)*HD4+ PC4(i,l)*HD5);
                    end
                    cont = cont + 1;
                end
            end
        end
end
