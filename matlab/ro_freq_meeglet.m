function [out] = ro_freq_meeglet(dat,cfg)

% meeglet (ver 1.0)
% Denis Engemann & Jörg Hipp, Feb 2024
%
% Log-spaced frequency transform using Morlet wavelets for continous
% electrophysiological signals. Derive power, covariance and various
% conncetivity metrics
%
% % Input
% - dat ... [channels x samples], set invalid data sections to NaN;
%           Input signal needs to have in µV
% - cfg ... struct, see script itself for default parameters, cfg.fsample
%           has to be set by the user
% 
% % Output
% Metrics indicated in cfg.output (default: all)
% Values returned as fields in struct out:
% - foi .. frequencies for which below metrics are derived in Hz
% - n .. number of time-windows used to compute below metrics
% - cfg .. the cfg struct containing all parameters
% - unit .. unit of the power spectral density (pow)
% - pow .. Power spectrum. The total power in µV² for a channel (ich) over
%          the frequency range [cfg.f_start,cfg.f_end] can be derived as:
%          if cfg.density is 'linear'         -> total_power = trapz(foi,pow(ich,:))
%          if cfg.density is 'log2' (default) -> total_power = trapz(log2(foi),pow(ich,:))
% - csd .. Cross-spectral density matrix
% - cov .. Covariance matrix, same as real(csd)
% - coh .. Coherence, compley values use coh.*conj(coh) to derive magnitude squared coherence
% - icoh .. Imaginary coherence (Nolte et al., Clin Neurophys, 2004); 
%           function returns signed values, use abs(icoh) to be line with the original definition
%           is the same as the cross-sepectral density of the orthogonalized signal
% - plv .. Phase-locking value (Lachaux et al., HBM, 1999)
% - pli .. Phase lag index (Stam et al., HBM, 2007); function returns signed values,
%          use abs(pli) to be in line with the original definition;
%          is the same as the phase-locking value of the orthognalized signal.
% - dwpli .. Debiased weighted phase lag index (Vinck et al., NeuroImage 2011)
% - r_plain .. Power correlation
% - r_orth .. Orthogonalized power correlation (Hipp et al., Nature Neurosci, 2012)
% - gim .. Global interaction measure (Ewald et al., NeuroImage, 2012)
%
% %% Example 1:
% % Single channel with background noise and sine wave at 10 Hz
% fsample = 250;
% T = 60; % [s]
% n_chan=1;
% n_sample = T*fsample;
% t=(0:n_sample-1)/fsample; % [s]
% f = 10; % Hz
% rng(0)
% x=randn(n_chan,n_sample)*10+repmat(sin(2*pi*f*t)*10,[n_chan,1]); % [µV]
% %
% clear cfg
% cfg.fsample = fsample;
% cfg.delta_oct = 0.05;
% cfg.foi_start = 0.1; % given parameterization cannot start at 0
% cfg.foi_end = fsample/2;
% cfg.density='Hz';
% % µV²/Hz density
% freq = ro_freq_meeglet(x,cfg);
% figure('Color','w')
% plot(freq.foi,freq.pow,'r')
% xlabel('Frequency (Hz)'), ylabel(sprintf('%s',freq.unit))
% % µV²/oct density
% cfg.density='oct';
% freq2 = ro_freq_meeglet(x,cfg);
% %
% % Test Parseval's theoreme, i.e. power should be identical in time and
% % frequency domain (note: expected not to numerically match given the parameterization)
% Pt = mean(x.^2,2);    % time domain; same as sum(x.^2*dt)/T, where dt=1/fsample, T=n_sample*dt
% Pf_lin = sum(freq.pow(1:end-1).*diff(freq.foi)); % Frequency domain with linear frequency spaceing
% Pf_log = sum(freq2.pow(1:end-1).*diff(log2(freq2.foi))); % Frequency domain with log frequency spaceing
% fprintf('Power derived in time and frequency domains:\nPt = %.1f µV²\nPf_lin = %.1f µV²\nPf_log = %.1f µV²\n',Pt,Pf_lin,Pf_log)
%
% %% Example 2:
% % Artificial signal with background noise on 19 channels and two coherent sine-waves at two channels
% fsample = 250; % [Hz]
% T = 60; % [s]
% n_chan = 19;
% n_sample = T*fsample;
% t=(0:n_sample-1)/fsample; % [s]
% f = 10; % [Hz]
% rng(0)
% x=cumsum(randn(n_chan,n_sample),2); x=5*x./repmat(std(x,[],2),[1,size(x,2)]); % pink noise
% x(1,:) = x(1,:)+sin(2*pi*f*t);
% x(10,:) = x(10,:)+sin(2*pi*f*t+pi/4);
% %
% clear cfg
% cfg.output={'pow','coh','icoh'};
% cfg.fsample = fsample;
% freq = ro_freq_meeglet(x,cfg);
% %
% figure('Color','w')
% h=subplot(2,1,1); plot(log2(freq.foi),mean(freq.pow)), set(h,'XTick',log2(freq.foi(mod(log2(freq.foi),1)==0)),'XTickLabel',freq.foi(mod(log2(freq.foi),1)==0))
% title('Power spetral density'), xlabel('Frequency [Hz]'), ylabel(freq.unit), a=axis; a(3)=0; axis(a)
% subplot(2,2,3), imagesc(abs(freq.coh(:,:,20)),[0,1]), axis square, colorbar, title(sprintf('abs(COH), f=%.1fHz',freq.foi(20)))
% subplot(2,2,4), imagesc(freq.icoh(:,:,20),[-1,1]), axis square, colorbar, title(sprintf('iCOH, f=%.1fHz',freq.foi(20)))

% 2024, Joerg Hipp, F. Hoffmann La-Roche Ltd

% sanity checks
if ~isfield(cfg,'fsample'), error('specify cfg.fsample!'), end
if isfield(cfg,'bw_oct') && isfield(cfg,'qt');
    error('Please provide either bw_oct or qt in cfg but not both!'); 
end
n_sens=size(dat,1); n_sample=size(dat,2);
if n_sens>n_sample, fprintf('Warning, number of channles (%i) larger than number of samples (%i), check orientation of data matrix\n',n_sens,n_sample), end

% input parameter definition / default values
cfg.ver_meeglet = 1.0; % Version number
cfg.ver_matlab = ver; % Store information about the matlab environment
if ~isfield(cfg,'verbose'),            cfg.verbose=1;                        end % if true, plot progress in command line
if ~isfield(cfg,'bw_oct') && ~isfield(cfg,'qt')                              
                                       cfg.bw_oct = 0.5;                         % [oct] spectral smoothing
                                       cfg.qt = bw2qt(cfg.bw_oct);               % characteristic wavelet parameter
elseif isfield(cfg,'bw_oct')
                                       cfg.qt = bw2qt(cfg.bw_oct);               % characteristic wavelet parameter
else isfield(cfg,'qt');
                                       cfg.bw_oct = qt2bw(cfg.qt);               % [oct] spectral smoothing
end
if ~isfield(cfg,'delta_oct'),          cfg.delta_oct = cfg.bw_oct/4;         end % [oct] spectral spacing between neighboring wavelets
if ~isfield(cfg,'foi_start'),          cfg.foi_start = 2;                    end % [Hz] start frequency
if ~isfield(cfg,'foi_end'),            cfg.foi_end = 32;                     end % [Hz] end frequency
if ~isfield(cfg,'window_shift'),       cfg.window_shift=0.25;                end % fraction the window is shifted of total window size
if ~isfield(cfg,'kernel_width'),       cfg.kernel_width=5;                   end % length of the Kernel in units of the temporal standarad deviation
if ~isfield(cfg,'allow_fraction_nan'), cfg.allow_fraction_nan=0;             end % number between 0 and 1, if >0 windows a fractio of NaNs up to this number are allowed, solution derived using a minimum norm approach
if ~isfield(cfg,'density');            cfg.density='oct';                    end % 'Hz', 'oct'
if ~isfield(cfg,'freq_shift_factor');  cfg.freq_shift_factor=1;              end % the whole log-frequency axis is "shifted" by this amount, this can be used to e.g. align all spectra to individual alpha peak frequencies (see e.g. Hipp et al., Sci Reports, 2021), no effect if set to 1 (default)
if ~isfield(cfg,'n_dim');              cfg.n_dim=size(dat,1);                end % rank of the dataset (used to derive the global interaction measure)
if ~isfield(cfg,'output');             cfg.output={'pow'};                   end % type of connectivity measures to derive 'csd','cov','gim','dwpli','r_orth','pli','icoh','r_plain','coh','plv'; power metrics are always returned

%% spectral parameter
foi       = 2.^(log2(cfg.foi_start):cfg.delta_oct:log2(cfg.foi_end)); % [Hz] center frequencies
foi       = foi*cfg.freq_shift_factor;
foi_min = 2*foi/(2^cfg.bw_oct+1);                         % formular for arithmetic mean
foi_max = 2*foi/(2^-cfg.bw_oct+1);
sigma_freq = (foi_max-foi_min)/(2*sqrt(2*log(2)));        % standard deviation in the frequency domain
sigma_time = 1./(2*pi*sigma_freq);                        % standard deviation in the time domain

if cfg.verbose, tic, fprintf('Morlet wavelet transform ['), end
for ifoi=1:length(foi)
    if cfg.verbose, fprintf('.'), end
    
    % convolution kernel
    n_win = int64(ceil(cfg.kernel_width*sigma_time(ifoi)*cfg.fsample+1));
    n_shift = int64(ceil(double(n_win)*cfg.window_shift));
    t = ((1:double(n_win))-double(n_win)/2-0.5)/cfg.fsample;
    z = t./sigma_time(ifoi);
    TAPER = exp(-(1/2)*z.^2);
    TAPER = TAPER/sqrt(sum(abs(TAPER).^2));
    iEXP = exp(1i*2*pi*foi(ifoi)*t);
    KERNEL = (TAPER.*iEXP).';
    switch cfg.density
        case 'Hz'
            SCALING = sqrt(2/cfg.fsample);
            unit='µV²/Hz';
        case 'oct'
            SCALING = sqrt(2/cfg.fsample) * sqrt(log(2)*foi(ifoi));
            unit='µV²/oct';
    end
    
    % collect info on nan sections
    idx_up = find(diff([0,isnan(sum(dat))])==1);
    idx_down = find(diff([0,isnan(sum(dat))])==-1);
    nan_width = zeros(1,size(dat,2));
    for icnt=1:length(idx_up)-1
        nan_width(idx_up(icnt):idx_down(icnt)) = idx_down(icnt)-idx_up(icnt);
    end
    if length(idx_up)>length(idx_down)
        nan_width(idx_up(end):end) = length(nan_width)+1-idx_up(end);
    end
    
    % memory allocation
    DAT      = nan(n_sens,length(1:n_shift:n_sample-n_win+1));
    frac_nan = nan(1,size(DAT,2));
    
    % convolution
    cnt = 0;
    for isection = 1:n_shift:size(dat,2)-n_win+1
        section = double(dat(:,isection:isection+n_win-1));
        nan_width_section = nan_width(isection:isection+n_win-1);
        n_nan = sum(isnan(section(1,:)));
        cnt=cnt+1;
        frac_nan(cnt) = n_nan/size(section,2);
        if n_nan==0
            DAT(:,cnt) = section*flip(KERNEL,1)*SCALING; % mirror image kernel for convolution
        elseif n_nan<size(section,2)*cfg.allow_fraction_nan && ...
                max(nan_width_section)<size(section,2)*cfg.allow_fraction_nan
            idx_valid = find(~isnan(section(1,:)));
            KERNEL_tmp = flip(KERNEL,1); % mirror image kernel for convolution
            KERNEL_tmp = KERNEL_tmp(idx_valid)/sqrt(sum(abs(KERNEL_tmp(idx_valid)).^2));
            DAT(:,cnt) = section(:,idx_valid)*KERNEL_tmp*SCALING;
        else
            DAT(:,cnt) = nan(size(section,1),1);
        end % valid section
    end
    
    % derive metrics for frequency-transformed data
    idx_valid = find(~isnan(DAT(1,:)));
    n_valid = length(idx_valid);
    DAT = DAT(:,idx_valid);
    frac_nan = frac_nan(idx_valid);
    
    if n_valid>0
        %% Power (compute by default)
        if ifoi==1; pow = zeros(n_sens,length(foi)); pow_geo = pow; pow_median = pow; pow_var = pow; n = zeros(1,length(foi)); end
        pow(:,ifoi)        = mean(abs(DAT).^2,2);
        pow_geo(:,ifoi)    = exp(mean(log(abs(DAT).^2),2));
        pow_median(:,ifoi) = median(abs(DAT).^2,2);
        pow_var(:,ifoi)    = var(abs(DAT).^2,[],2);
        n(ifoi)            = n_valid;
        
        %% Cross-spectral density matrix
        if any(ismember({'csd','cov','icoh','coh'},cfg.output))
            if ifoi==1, csd = zeros(n_sens,n_sens,length(foi)); end
            csd(:,:,ifoi) = DAT*DAT'/n_valid;
        end
        
        %% Covariance matrix
        if ismember('cov',cfg.output)
            if ifoi==1, cov = zeros(n_sens,n_sens,length(foi)); end
            cov(:,:,ifoi)=real(csd(:,:,ifoi));
        end
        
        %% Coherence matrix
        if any(ismember({'icoh','coh'},cfg.output))
            if ifoi==1, coh = zeros(n_sens,n_sens,length(foi)); end
            coh(:,:,ifoi)=csd(:,:,ifoi)./sqrt(diag(csd(:,:,ifoi))*diag(csd(:,:,ifoi))');
        end
        
        %% Imaginary coherence matrix
        if ismember('icoh',cfg.output)
            if ifoi==1, icoh = zeros(n_sens,n_sens,length(foi)); end
            icoh(:,:,ifoi)=imag(coh(:,:,ifoi));
        end

        %% Global interaction measure, % Ewald et al., Neuroimage, 2012, eq. 15, note, the pinv instead of inv reduces the
        if ismember('gim',cfg.output)
            if ifoi==1, gim = zeros(1,length(foi)); end
            C = DAT*DAT'/n_valid;
            if cfg.n_dim < size(C, 1)
                C_inv = ro_pinv(real(C),cfg.n_dim);
            else
                C_inv = pinv(real(C));
            end
            gim(ifoi) = 1/2*trace(C_inv*imag(C)*C_inv*imag(C)');
        end
        %% Phase-locking value matrix
        if ismember('plv',cfg.output)
            if ifoi==1, plv = zeros(n_sens,n_sens,length(foi)); end
            DATN=DAT./abs(DAT);
            plv(:,:,ifoi) = DATN*DATN'/n_valid;
        end

        %% Phase-lag index matrix
        if ismember('pli',cfg.output)
            if ifoi==1, pli = zeros(n_sens,n_sens,length(foi)); end
            DATN=DAT./abs(DAT);
            for icnt=1:n_sens
                for jcnt = (icnt+1):n_sens % upper right triangle
                    pli(icnt,jcnt,ifoi) = mean(sign(imag(DATN(icnt,:).*conj(DATN(jcnt,:)))));
                end
            end
            pli(:,:,ifoi)=pli(:,:,ifoi)+pli(:,:,ifoi)';
        end

        %% dwPLI (debiased weighted phase-lag index)
        if ismember('dwpli',cfg.output)
            if ifoi==1, dwpli = zeros(n_sens,n_sens,length(foi)); end
            for icnt=1:n_sens
                for jcnt = (icnt+1):n_sens % upper right triangle
                    cdi = imag(DAT(icnt,:).*conj(DAT(jcnt,:)));
                    imagsum      = sum(cdi,2);
                    imagsumW     = sum(abs(cdi),2);
                    debiasfactor = sum(cdi.^2,2);
                    dwpli(icnt,jcnt,ifoi)  = (imagsum.^2 - debiasfactor)./(imagsumW.^2 - debiasfactor);
                end
            end
            dwpli(:,:,ifoi)=dwpli(:,:,ifoi)+dwpli(:,:,ifoi)';
        end

        %% Power correlations
        if ifoi==1
            if ismember('r_plain',cfg.output), r_plain = zeros(n_sens,n_sens,length(foi)); end
            if ismember('r_orth',cfg.output), r_orth = zeros(n_sens,n_sens,length(foi)); end
        end
        for isens=1:n_sens
            seed = DAT(isens,:);
            seed_logpow = log(seed.*conj(seed));
            src = DAT;
            src_logpow = log(src.*conj(src));
            % plain
            if ismember('r_plain',cfg.output)
                r_plain(isens,:,ifoi) = ro_corrcoef(repmat(seed_logpow,[n_sens,1]),src_logpow,2);
            end  
            % orthogonalized
            if ismember('r_orth',cfg.output)
                src_orth = imag(DAT.*conj(repmat(seed./abs(seed),[n_sens,1]))) .* repmat(sqrt(-1)*seed./abs(seed),[n_sens,1]);
                src_logpow_orth = log(src_orth.*conj(src_orth));
                r_orth(isens,:,ifoi) = ro_corrcoef(repmat(seed_logpow,[n_sens,1]),src_logpow_orth,2);
            end
        end
    else
        n(ifoi)                 = 0;
        pow(:,ifoi)             = nan(1,n_sens);
        pow_var(:,ifoi)         = nan(1,n_sens);
        pow_geo(:,ifoi)         = nan(1,n_sens);
        pow_median(:,ifoi)      = nan(1,n_sens);
        if ismember('gim',cfg.output), gim(ifoi)             = nan; end
        if ismember('csd',cfg.output), csd(:,:,ifoi)         = nan(n_sens); end
        if ismember('cov',cfg.output), cov(:,:,ifoi)         = nan(n_sens); end
        if ismember('dwpli',cfg.output), dwpli(:,:,ifoi)     = nan(n_sens); end
        if ismember('r_plain',cfg.output), r_plain(:,:,ifoi) = nan(n_sens); end
        if ismember('coh',cfg.output), coh(:,:,ifoi)         = nan(n_sens); end
        if ismember('plv',cfg.output), plv(:,:,ifoi)         = nan(n_sens); end
        if ismember('r_orth',cfg.output), r_orth(:,:,ifoi)   = nan(n_sens); end
        if ismember('icoh',cfg.output), icoh(:,:,ifoi)       = nan(n_sens); end
        if ismember('pli',cfg.output), pli(:,:,ifoi)         = nan(n_sens); end
    end
end % loop foi

if cfg.verbose, fprintf('] %.1f sec\n',toc), end

% Output struct
out.cfg = cfg;
out.foi = foi;
out.pow = pow;
out.pow_var = pow_var;
out.unit = unit;
out.pow_median = pow_median;
out.pow_geo = pow_geo; % geometic mean across power values from each window
out.n = n;
out.qt = cfg.qt;
out.bw_oct = cfg.bw_oct;
for icnt=1:length(cfg.output)
    eval(sprintf('out.%s = %s;',cfg.output{icnt},cfg.output{icnt}));
end

%%-------------------------------------
function qt = bw2qt(bw)

L = sqrt(2*log(2));
qt = (2^(bw)+2^(-bw)+2)/(2^(bw)-2^(-bw))*L;


%%-------------------------------------
function bw = qt2bw(qt)

L = sqrt(2*log(2));
bw = log2( sqrt((L^2/(L-qt)^2)-((L+qt)/(L-qt))) - (L/(L-qt)));


%%-------------------------------------
function X = ro_pinv(A,r)

% X = ro_pinv(A,r)
%
% A covariance matrix 
% r dimension to use
% X pseudoinverse
%
% Computes the pseudoinverse of a matrix using svd with only r dimensions

if nargin<2
    r = size(A,1);
end
[U,S,V] = svd(A,0);
s = diag(S);
s = diag(ones(r,1)./s(1:r));
X = V(:,1:r)*s*U(:,1:r)';

%%-------------------------------------
function [r,p,z,r_f,r_op] = ro_corrcoef(x,y,dim)

% instead of using corrcoeff
% dimension dim contains the data
%
% r_f .. see: Fisher, R.A. (1915). Frequency distribution of the values of the correlation coefficient in samples from an indefinitely large population. Biometrika, 10, 507-521.
% http://digital.library.adelaide.edu.au/dspace/bitstream/2440/15166/1/4.pdf
% r_op .. see Olkin, I., Pratt, J.W. (1958). Unbiased estimation of certain correlation coefficients. Annals of Mathematical Statistics, 29, 201-211.
% http://digital.library.adelaide.edu.au/dspace/bitstream/2440/15166/1/4.pdf
% see also: Zimmerman, D. W., Zumbo, B. D, Richard H. Williams, R. H. (2003). Bias in Estimation and Hypothesis Testing of Correlation. Psicolï¿œgica (2003), 24, 133-158.

if nargin<3
    dim = 1;
end

if ndims(x)==2 && min(size(x))==1
    x = x(:);
    y = y(:);
    dim = 1;
end

r = squeeze((mean(x.*y,dim)-mean(x,dim).*mean(y,dim)) ./ sqrt(mean(x.^2,dim)-mean(x,dim).^2) ./ sqrt(mean(y.^2,dim)-mean(y,dim).^2));

n = size(x,dim);
if nargout>1
    dof  = n-2;
    t    = r.*sqrt(dof./(1-r.^2));
    p = (1-cdf('t',abs(t),dof))*2;
end
if nargout >2
    z = sign(r).*norminv(1-p/2,0,1);
    %   z    = spm_t2z(t,dof);
    %p = normcdf(z,0,1)*2;
end
if nargout>3
    r_f = r*(1+(1-r.^2)/(2*n)); % note, assumes independent samples!
end
if nargout>4
    r_op = r*(1+(1-r.^2)/(2*(n-3))); % note, assumes independent samples!
end
