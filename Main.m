function [n_mean, ne_mean, xr_mean] = Main(varargin)
%MAIN Simulate trait and regulator evolution under sex-specific selection.
%   [n_mean, ne_mean, xr_mean] = Main('name', value, ...) returns the
%   generation-wise means of [n, ne, xr] over the last min(Tcd,1000)
%   generations. n counts distinct trait values with positive abundance;
%   ne is their inverse Simpson diversity; xr is the retained trait range.
%   Figures are displayed during/following a run; no files are saved.
%   See README.md for usage and output definitions.

% --- Simulation parameters --- %

% -- Demographic parameters -- %
Nmax      =    5e4; % Census population size (individuals)
T         =    2e5; % Generations with mutation input, before cooldown
Tevod     =      0; % Regulator mutation starts when t > Tevod
x0        =      1; % Initial displacement from the singular point
mux       =   5e-6; % Trait mutation rate per allele copy
mud       =   5e-5; % Regulator mutation rate per allele copy
delta_x   =   0.08; % Expected mutational step size
delta_d   = log(2); % Expected mutational step size in log-affinity (natural log)
rho       =      0; % Correlation of female and male regulator-mutation steps

% -- Evolutionary parameters -- %
sif     =      1; % Female selection strength (inverse sigma)
sim     =      1; % Male selection strength (inverse sigma)
dr      = [4 4]; % Female/male log10 affinity ranges; evolving d is bounded by +/-abs(dr)/2
XL      = false; % true for X-linked inheritance; false for autosomal


% -- Load input variables -- %
loadinput(nargin)

% -- Fixed runtime and plotting settings -- %
prpl     =     1; % Update the figure during the run
prog     =     1; % Display simulation progress
t0       =   1e2; % Logarithmic time offset and initial live-plot delay
Nlim_rel =  1e-2; % Allele-history pruning threshold relative to Nmax
Tcd      = min(2e5,round(T/2)); % Mutation-free cooldown generations
pl_res   = min(1e3,T+Tcd+1); % Target sampling times, capped for short runs
fs       =     7; % Figure font size

% -- Check parameters before initialization -- %
validateattributes(T,       {'numeric'}, {'real','finite','scalar','integer','positive'}, mfilename,'T');
validateattributes(Tevod,   {'numeric'}, {'real','finite','scalar','integer','nonnegative'}, mfilename,'Tevod');
validateattributes(Nmax,   {'numeric'}, {'real','finite','scalar','positive'}, mfilename,'Nmax');
validateattributes(dr,     {'numeric'}, {'real','finite','size',[1 2]}, mfilename,'dr');
validateattributes(sif,    {'numeric'}, {'real','finite','scalar','positive'}, mfilename,'sif');
validateattributes(sim,    {'numeric'}, {'real','finite','scalar','positive'}, mfilename,'sim');
validateattributes(mux,    {'numeric'}, {'real','finite','scalar','>=',0,'<=',1}, mfilename,'mux');
validateattributes(mud,    {'numeric'}, {'real','finite','scalar','>=',0,'<=',1}, mfilename,'mud');
validateattributes(rho,    {'numeric'}, {'real','finite','scalar','>=',-1,'<=',1}, mfilename,'rho');
validateattributes(delta_x,{'numeric'}, {'real','finite','scalar','nonnegative'}, mfilename,'delta_x');
validateattributes(delta_d,{'numeric'}, {'real','finite','scalar','nonnegative'}, mfilename,'delta_d');
validateattributes(x0,     {'numeric'}, {'real','finite','scalar'}, mfilename,'x0');
validateattributes(XL,     {'numeric','logical'}, {'scalar','binary'}, mfilename,'XL');
if round(Nmax/2)*2 < 200
    error('Main:PopulationSize','Nmax must round to at least 200 individuals for initialization.')
end
if mux+mud > 1
    error('Main:MutationRate','The sum of mux and mud must not exceed 1.')
end
if Tevod> T                      ; dr   =               [0 0]; end % If d never evolves, set range to 0
if ~any(dr)                      ; mud  =                   0; end % No mutations in d, if no dominance

% --- Model functions --- %

% -- Van Dooren dominance -- %
z_f = @(x1,x2,d1,d2) (x1.*exp(d1*log(10)) + x2.*exp(d2*log(10)))./(exp(d1*log(10))+exp(d2*log(10)));

% --- Simulation preparations --- %
T0 = T    ; % Time of mutation input
T  = T+Tcd; % Total time, including mutation cooldown period

% --- Variables shared by nested functions --- %
ll=[];
d_f=[]; d_m=[]; d=[];
n=[];ne=[];xr=[];xs=[];
N=[];Nf=[];Nm=[]; nxv=[]; nsi=[]; Tsamp=[]; nXstore=[]; Nth=[];
oi=[];   tplv=[];
sifsqo2=[]; simsqo2=[];
x=[]; xind=[]; xpih=[];
o_f=[];o_m=[];sf_mean=[];sm_mean=[];
t_last_report=0;time_left=inf; ism_str=[];
figh=[]; tocv_plot=[];
alive=[];xx=[];

tic_tot=tic; tic_last=tic; tic_plot=tic; % Timers for elapsed runtime, progress reports, and plotting
firstplot = true;

% --- Simulation initialization --- %
t=0;
sim_prework

% --- The simulation loop --- %
while t < T
    t = t+1; % Advance one generation

    % --- 1. Fitness and reproduction --- %
    Ni = Nf >0 | Nm >0;
    [Nf(Ni),Nm(Ni)]  =  fitness_and_reproduction(Nf(Ni), Nm(Ni)); %#ok<AGROW> % Reproduction

    if (length(Nf)>1 || length(Nm)>1) && (min(Nf(Ni)+Nm(Ni)) == 0)
        remove_extinct_alleles; % Release extinct simulation slots and prune small allele histories
    end

    % --- 2. Mutation --- %
    muti = mutation;
    Nth = storeAlleles(Nf+Nm,xind,muti,Nth);

    % --- 3. Plot and display progress --- %
    plotprog
end

% Mean manuscript statistics over the final sampling window.
nxv = nxv(1:oi,:);
mean_nxv = mean(nxv,1);
n_mean  = mean_nxv(1);
ne_mean = mean_nxv(2);
xr_mean = mean_nxv(3);
if prog
    disp('Simulation done!')
end

% ======================== %
% --- Nested functions --- %
% ======================== %

%   % ------------------------------------------ %
%   % --- Function: fitness_and_reproduction --- %
%   % ------------------------------------------ %
    function [Nf1_x, Nm1_x] = fitness_and_reproduction(Nf_x,Nm_x)
        n = length(Nf_x); % Number of active alleles

        Nf_sum = sum(Nf_x) ; % Total female allele copies
        Nm_sum = sum(Nm_x) ; % Total male allele copies

        pfm = Nf_x*Nm_x'/Nf_sum/Nm_sum;      % Frequencies of main trait genotypes among offspring
        pfm = triu(pfm) + tril(pfm,-1)'; % i,j = j,i
        ind = pfm>0; ind = ind(:);     % Indices of present genotypes among offspring

        xv = x(Ni,:);
        dv = d(Ni,:);

        d_f = dv(:,1);
        d_m = dv(:,2);

        if XL
            pfm_xm = Nf_x/Nf_sum;% Frequencies of main trait genotypes among male offspring
            ind_xm = pfm_xm>0;
        end

        xind_1 = repmat( 1:n  ,n,1); xind_1 = xind_1(ind);
        xind_2 = repmat((1:n)',n,1); xind_2 = xind_2(ind);

        zf_ind         = z_f(xv(xind_1),xv(xind_2), d_f(xind_1),d_f(xind_2)); % Realized female trait z_f among offspring
        if ~XL; zm_ind = z_f(xv(xind_1),xv(xind_2), d_m(xind_1),d_m(xind_2)); % Realized   male trait z_m among offspring
        else  ; zm_ind = xv(ind_xm);                                          % Realized   male trait z_m among offspring
        end

        wf_ind = exp(-(o_f-zf_ind).^2.*sifsqo2);   % Survival of females
        wm_ind = exp(-(o_m-zm_ind).^2.*simsqo2);   % Survival of   males
        pfm1     = pfm(ind);    % Frequency of all offspring genotypes, before survival
        p_sf     = wf_ind.*pfm1/2;  % Female offspring genotypes frequencies after survival
        if ~XL
            p_sm = wm_ind.*pfm1/2;  % Male   offspring genotypes frequencies after survival
        else
            pfm1m = pfm_xm(ind_xm); % Frequency of all male offspring genotypes, before survival
            p_sm = wm_ind.*pfm1m/2;  % Male   offspring genotypes frequencies after survival
        end

        sf_mean = sum(p_sf*2);  % Sum --> Mean female survival
        sm_mean = sum(p_sm*2);  % Sum --> Mean   male survival
        if n>1
            s_mean = (sf_mean + sm_mean)/2; % Population mean survival
            Uf   = p_sf/s_mean;             % Relative female offspring genotypes frequencies after survival = Genotype relative fitness
            Um   = p_sm/s_mean;             % Relative   male offspring genotypes frequencies after survival = Genotype relative fitness

            if abs(sum(Uf(:)) + sum(Um(:)) - 1) > 1e6*eps % Check that Uf and Um sums to 1
                error('Main:FitnessNormalization','%sOffspring frequencies sum to %.17g instead of 1.',errinf,sum(Uf(:))+sum(Um(:)))
            end

            % --- Main trait --- %
            % -- Expected allele population size to next generations size -- %
            % - Females - %
            Nfet1        = zeros(n  );  % Preallocate expected frequencies after survival, female allele pairs
            if ~XL;Nmet1 = zeros(n  );  % Preallocate expected frequencies after survival, male allele pairs
            else  ;Nmet1 = zeros(n,1);  % Preallocate expected frequencies after survival, X-linked male alleles
            end

            Nfet1(ind)         = Uf; % Expected allele-pair frequency after survival, females
            if ~XL; Nmet1(ind) = Um; % Expected allele-pair frequency after survival, males
            else  ; Nmet1(ind_xm)= Um; % Expected allele frequency after survival, X-linked males
            end

            % Adjusting sum from numeric error
            p_et = (sum(Nfet1(:)) + sum(Nmet1(:))); % expected total frequency
            if abs(1-p_et)>0.01
                error('Main:OffspringNormalization','%sOffspring frequency sum differs from 1 by %.17g.',errinf,1-p_et)
            end
            if abs(1-p_et)>numel(Nfet1)*eps
                Nfet1 = Nfet1/p_et;
                Nmet1 = Nmet1/p_et;
                p_et = (sum(Nfet1(:)) + sum(Nmet1(:))); % expected total frequency
                if abs(1-p_et)>numel(Nfet1)*2*eps
                    error('Main:OffspringNormalization','%sOffspring normalization left a residual of %.17g.',errinf,1-p_et)
                end
            end

            N1 = zeros(n,n,2); % Preallocate array of genotype numbers of the next generation
            N1(:,:,1)     = Nfet1; % Allocated expected female frequencies
            if ~XL
                N1(:,:,2) = Nmet1; % Allocated expected   male frequencies
            else
                N1(:,1,2) = Nmet1; % Allocated expected   male frequencies
            end

            N1(:) = mnrnd(Nmax,N1(:)); % Number of genotypes, after the sampling of Nmax surviving individuals

            Nf1_x     = sum(N1(:,:,1)+N1(:,:,1)',2); % one allele from the i and one from the j position
            if ~XL
                Nm1_x = sum(N1(:,:,2)+N1(:,:,2)',2); % one allele from the i and one from the j position
            else
                Nm1_x = N1(:,1,2); % only one allele from the i position (no second allele)
            end
        else
            Nf1_x = 2*binornd(Nmax,sf_mean/(sf_mean+sm_mean)); % Female allele copies after sampling Nmax surviving individuals
            Nm1_x = (Nmax*2 - Nf1_x)/(1+XL);                   % Male allele copies after sampling Nmax surviving individuals
        end
    end

%   % -------------------------- %
%   % --- Function: mutation --- %
%   % -------------------------- %
    function muti = mutation
        if t>T0 % No more mutations
            muti = {[],[]};
        else
            Nv{1} = Nf;
            Nv{2} = Nm;
            dv = d;
            xv = x;

            if t>Tevod
                mud2 = mud;
            else
                mud2 = 0;
            end

            muTotal = mux + mud2; % Total mutation rate
            muti = {[],[]};

            for k = 2:-1:1
                Nmut{k}          = zeros(size(Nv{k}));
                % Approximate binomial mutation counts using 50 randomized rounding terms.
                Nmut{k}(Nv{k}>1) = sum(floor(rand(sum(Nv{k}>1),50) + Nv{k}(Nv{k}>1)*(muTotal)/50),2);
                if sum(Nmut{k})>0
                    Nmut{k}(Nv{k}>1) = min(Nmut{k}(Nv{k}>1),Nv{k}(Nv{k}>1)-1); % mutants can't remove alleles.
                    Nv{k}(Nv{k}>1) = Nv{k}(Nv{k}>1) - Nmut{k}(Nv{k}>1);

                    while sum(Nmut{k})>sum(isnan(Nv{k}))
                        % -- Info string to display info about increased storage -- %
                        if prog
                            ism_str = [ism_str, strrep(sprintf(', matrix sizes (%3.0f): t = %6.3g',[length(Nv{k})*2,t]),'e+0','e')]; %#ok<AGROW>
                        end

                        xv   = [xv  ;  nan(size(xv  ))]; %#ok<AGROW>
                        dv   = [dv  ;  nan(size(dv  ))]; %#ok<AGROW>
                        Nv{1} = [Nv{1};  nan(size(Nv{1}))];
                        Nv{2} = [Nv{2};  nan(size(Nv{2}))];
                        xind = [xind;  nan(size(xind))]; %#ok<AGROW>
                    end

                    while sum(Nmut{k}) > sum(isnan(xpih(:,1)))
                        xpih = [xpih; nan( 100,size(xpih,2))]; %#ok<AGROW>
                        Nth  = [Nth ; cell(100,size(Nth ,2))]; %#ok<AGROW>

                        % -- Info string to display info about increased storage -- %
                        if prog
                            ism_str = [ism_str, strrep(sprintf(', storage sizes (%3.0f): t = %6.3g',[size(xpih,1),t]),'e+0','e')]; %#ok<AGROW>
                        end
                    end
                    % placing new mutants in empty slots
                    muti{k} = find(isnan(Nv{k}),sum(Nmut{k}));  % Indices of empty simulation slots

                    Nv{k  }(muti{k}) = 1; % place mutant in sex k
                    Nv{3-k}(muti{k}) = 0; % no mutant in other sex

                    [xv(muti{k}), dv(muti{k},:)]  = getmut(xv,dv,Nmut{k},mud2);

                    te = find(isnan(xpih(:,1)));
                    xind(muti{k},:)          = [te(1:sum(Nmut{k})),xind(reprow(find(Nmut{k}),Nmut{k}(Nmut{k}>0)),1)];
                    xpih(xind(muti{k},1),:)  = [xv(muti{k},:),dv(muti{k},:),xind(muti{k},2)];
                    for xi = xind(muti{k},1)'
                        if sum(~any(xpih(xi,:) - xpih(~isnan(xpih(:,1)),:),2))>1
                            if ~any(abs(xpih(xi,2:3)) == dr/2)
                                error('Main:DuplicateMutation','%sDuplicate allele entry at index %d.',errinf,xi)
                            end
                        end
                    end
                end
            end
            if sum(Nmut{1})>0 || sum(Nmut{2})>0

                Nf = Nv{1};
                Nm = Nv{2};
                x = xv;
                d = dv;
            end
        end
    end

%   % ------------------------ %
%   % --- Function: getmut --- %
%   % ------------------------ %
    function [x, logd] = getmut(x,logd,Nmut,mud2)
        xd = [x,logd];
        if any(Nmut>1) % if there is a type with more than one mutant
            xd = reprow(xd(Nmut>0,:),Nmut(Nmut>0));
        else % if there is only one type per mutant
            xd = xd(Nmut>0,:);
        end
        NmutT = sum(Nmut);
        % -- Mutate either x or the regulators -- %
        mutxi = rand(NmutT,1) < mux/(mux + mud2); % Mutations on x
        mutdi = ~mutxi;                         % Remaining mutations on regulators

        if mud % If d evolves independently of x
            dmut = zeros(NmutT,3);
            if any(mutdi); dmut(mutdi,2:3) = mvnrnd([0,0], Sigma_d, sum(mutdi)  ); end % Regulator mutations
            if any(mutxi); dmut(mutxi,1)   = normrnd(0   , sigx   , sum(mutxi),1); end % Trait mutations
            xd = xd + dmut;
        else % If d is a function of x
            xd(:,1) = xd(:,1) + normrnd(zeros(NmutT,1),delta_x*sqrt(pi)/sqrt(2));
            xd(:,2:3) = xd(:,1).*dr;
        end

        % -- Return mutated values -- %
        x = xd(:,  1);
        if mud
            logd = min(abs(dr)/2,max(-abs(dr)/2, xd(:,2:3)));
        else
            logd = xd(:,2:3);
        end
    end

%   % ---------------------------------------- %
%   % --- Function: remove_extinct_alleles --- %
%   % ---------------------------------------- %
    function remove_extinct_alleles
        xv=x(Ni,:);dv=d(Ni,:);Nf1=Nf(Ni);Nm1=Nm(Ni);
        x1ind=xind(Ni,:);
        for i = find(Nf1+Nm1==0)'
            ei = x1ind(i,1); %erase index
            if~any(xpih(:,end)==ei) && max(Nth{ei}(1,:))<Nlim_rel*Nmax
                % Erase the allele's plotting history
                [xpih,Nth,ll] = eraseAl(ei,xpih,Nth,ll,1);
            elseif  max(Nth{ei}(1,:))<Nlim_rel*Nmax
                te = find(xpih(:,end)==ei);
                ll = [ll;[repmat(ei,length(te),1) te] ]; %#ok<AGROW>
            end
            if t <= T-Tsamp
                xv   (i,:) = nan;
                dv   (i,:) = nan;
                x1ind(i,:) = nan;
            end
            Nf1  (i  ) = nan;
            Nm1  (i  ) = nan;
        end

        if any(Nf1+Nm1==0)
            error('Main:ExtinctionCleanup','%sZero-count allele remained after extinction cleanup.',errinf)
        end

        % --- Function output --- %
        x(Ni,:) = xv;
        d(Ni,:) = dv;
        Nf(Ni)  = Nf1;
        Nm(Ni)  = Nm1;
        xind(Ni,:) = x1ind;
    end

%   % ------------------------- %
%   % --- Function: eraseAl --- %
%   % ------------------------- %
    function [xpih,Nth,ll] = eraseAl(erind,xpih,Nth,ll,lev)
        xpih(erind,:) = nan;
        Nth{ erind  } = [] ;
        if ~isempty(ll)
            erind2 = (ll(:,2)==erind);
            if any(erind2)
                te = ll(erind2,1);
                ll2 = ll;
                ll = ll(~erind2,:);
                if sum(ll2(:,1) == ll2(erind2,1))==1
                    [xpih,Nth,ll] = eraseAl(te,xpih,Nth,ll,lev+1);
                end
            end
        end
    end

%   % ------------------------------ %
%   % --- Function: storeAlleles --- %
%   % ------------------------------ %
    function Nth = storeAlleles(N,xind,muti,Nth)
        muti = [muti{1};muti{2}]; % indices of all mutations in either x or d
        for j = xind(N>0,1)'
            Nt   = N(xind(:,1)==j); % Population size N at time t
            if size(Nth{j},2)>=2
                Ntm12 = Nth{j}(1,end-1:end); % Population size N at the last two time slots
            end
            % Compress allele histories when both recent log-count differences are <0.07
            % or count differences are <5, and both differences are <0.01*(2*Nmax).
            if ~any(j==xind(muti,1)) && size(Nth{j},2) > 1 && (max(abs(log(Nt./Ntm12))) < 0.07 || max(abs(Nt - Ntm12))<5) && max(abs(Nt - Ntm12)) < Nmax*2*0.01
                Nth{j} = [Nth{j}(:,1:end-1), [Nt;t]];

            else
                Nth{j} = [Nth{j}, [Nt;t]];
            end
        end

    end

%   % ---------------------- %
%   % --- Function: getH --- %
%   % ---------------------- %
    function H = getH(N)
        H = 1-sum((N/sum(N)).^2);
    end

%   % --------------------------- %
%   % --- Function: loadinput --- %
%   % --------------------------- %
    function loadinput(nargin)
        if mod(nargin,2)
            error('Main:NameValuePairs','Supply parameters as name/value pairs.')
        end
        for i = 1:2:nargin
            if ~((ischar(varargin{i}) && isrow(varargin{i})) || (isstring(varargin{i}) && isscalar(varargin{i})))
                error('Main:ParameterName','Each parameter name must be a character row or string scalar.')
            end
            switch varargin{i}
                case 'delta_d'  ; delta_d   = varargin{i+1};
                case 'delta_x'  ; delta_x   = varargin{i+1};
                case 'dr'       ; dr        = varargin{i+1};
                case 'mux'      ; mux       = varargin{i+1};
                case 'mud'      ; mud       = varargin{i+1};
                case 'Nmax'     ; Nmax      = varargin{i+1};
                case 'rho'      ; rho       = varargin{i+1};
                case 'sif'      ; sif       = varargin{i+1};
                case 'sim'      ; sim       = varargin{i+1};
                case 'T'        ; T         = varargin{i+1};
                case 'Tevod'    ; Tevod     = varargin{i+1};
                case 'x0'       ; x0        = varargin{i+1};
                case 'XL'       ; XL        = varargin{i+1};
                otherwise; error('Main:UnknownParameter','Unknown parameter: %s.',varargin{i})
            end
        end
    end

%   % ---------------------- %
%   % --- Function: setn --- %
%   % ---------------------- %
    function setn
        % Group allele copies by identical x, ignoring differences in d.
        if t <= T-Tsamp
            alive = ~isnan(N);
        end
        [xx, ~, idx]  = unique(x(alive));
        dat = N(alive);
        dat(isnan(dat)) = 0;
        Nx = accumarray(idx, dat, [], @nansum); %#ok<NANSUM>

        % -- Retained allele statistics -- %
        n  = sum(Nx >0);
        ne = 1/(1 - getH(Nx));
        xr = max(xx)-min(xx);
    end

%   % ----------------------------- %
%   % --- Function: sim_prework --- %
%   % ----------------------------- %
sigx=[];sigdf=[];sigdm=[];Sigma_d=[];
    function sim_prework

        % -- Mutational variances and covariances -- %
        sigx  = delta_x*sqrt(pi)/sqrt(2);                                 % sqrt of mutational variance in x
        sigdf = delta_d*sqrt(pi)/sqrt(2)/log(10)*abs(dr(1))/max(abs(dr)); % sqrt of mutational variance in df
        sigdm = delta_d*sqrt(pi)/sqrt(2)/log(10)*abs(dr(2))/max(abs(dr)); % sqrt of mutational variance in dm
        Sigma_d = [sigdf^2, rho*sigdf*sigdm;   rho*sigdf*sigdm, sigdm^2]; % mutational covariance matrix of (df,dm)

        Nmax = round(Nmax/2)*2; % Round census size to an even integer
        sifsqo2 = sif^2/2;
        simsqo2 = sim^2/2;

        Tsamp = min(Tcd,1e3); % Final generations used for the three returned means
        nxv = zeros(Tsamp+1,3); % Columns: [n, ne, xr]
        oi = 0;

        % --- Preallocating variables --- %
        ng0s    = max( 64,1); % starting number of types for storage
        x       = nan(ng0s, 1); % trait values of alleles (x)
        d       = nan(ng0s, 2); % trait values of alleles (df,dm)
        xind    = nan(ng0s,2); % Columns: allele index and parent allele index

        xpih=  nan(100,2+2) ; % Columns: [x, d_f, d_m, parent allele index]
        Nth = cell(100,1)  ; % Each allele history is [allele counts; generations]

        Nf  = nan(ng0s,    1); % Female allele-copy counts
        Nm	= nan(ng0s,    1); % Male allele-copy counts

        Nf(1) = 2*binornd(Nmax,0.5);     % Number of female alleles, sampling p=0.5 out of Nmax
        Nm(1) = (Nmax*2 - Nf(1))/(1+XL); % Number of   male alleles, remaining population

        xs = (1/sif^2 - 1/sim^2)/(2/sif^2 + 2/sim^2); % Singular point

        % --- Set optima and dominance --- %
        % - VD-dom - %
        o_f = -0.5; o_m =  0.5;
        % -- Set x0 to far end of singular point -- %
        if sif>sim; x0 =  abs(x0);
        else      ; x0 = -abs(x0);
        end

        % -- Adjust the starting trait to avoid an extreme expected sex ratio -- %
        rat_fm = 0; % ratio of females and males
        c = 0;  % ratio iteration counter

        while rat_fm*Nmax<100 % Adjust x deterministically until the expected count is >=100 per sex
            if c>0 && (prog ||c>95)
                disp(['rat_fm = ' num2str(rat_fm) ', c = ' num2str(c) ', x0 = ' num2str(x0)])
                drawnow
            end
            x(1) = xs + x0 - c*0.01*sign(x0);
            if mud
                d(1,:) = 0;
            else
                d(1,:) = x(1).*dr;
            end
            sf = exp(-(o_f-z_f(x(1),x(1),d(1,1),d(1,1)) ).^2.*sifsqo2);
            sm = exp(-(o_m-z_f(x(1),x(1),d(1,2),d(1,2)) ).^2.*simsqo2);
            rat_fm = min(sf,sm)/(sf+sm);
            c = c +1;
        end

        xind(1,:) = [1 0];
        xpih(1,:) = [x(1,:) d(1,:) 0];
        Nth{1}    = [Nf(1)+Nm(1);t];

        if prog
            fprintf('Start: %s\n', char(datetime('now','Format','yyyy-MM-dd HH.mm')))
            if XL
                fprintf('Inheritance: X-linked\n')
            else
                fprintf('Inheritance: Autosomal\n')
            end
            fprintf('Demographic parameters:\n')
            fprintf('  Nmax = %.8g individuals\n', Nmax)
            fprintf('  T = %.8g generations with mutation; cooldown = %.8g generations\n', T0, Tcd)
            fprintf('  Tevod = %.8g (regulator mutation starts when t > Tevod)\n', Tevod)
            fprintf('  |x0| = %.8g (initial displacement); starting trait = %.8g\n', abs(x0), x(1))
            fprintf('  mu_x = %.8g; mu_d = %.8g (per allele copy)\n', mux, mud)
            fprintf('  delta_x = %.8g (expected mutational step size)\n', delta_x)
            fprintf('  delta_d = %.8g (expected step size in natural-log affinity)\n', delta_d)
            fprintf('  rho = %.8g (female/male regulator-mutation correlation)\n', rho)
            fprintf('Evolutionary parameters:\n')
            fprintf('  sigma_f^-1 = %.8g; sigma_m^-1 = %.8g\n', sif, sim)
            fprintf('  dr = [%.8g %.8g] (female/male log10 affinity ranges)\n\n', dr)
        end

        % --- Plot Prework --- %
        c = t0;
        while length(tplv)<pl_res
            c = c + 1;
            tplv =unique(round(exp(linspace(log(0+c),log(T+c),pl_res+c)))-c);
        end
        pl_res = length(tplv);
        nXstore = zeros(pl_res,3); % Columns: [n, ne, xr]
        nsi = 1;

    end

%   % ------------------------ %
%   % --- Function: fixfig --- %
%   % ------------------------ %
    function fixfig
        fw = 90; % mm
        fh = 80 ; % mm
        figh = figure('visible','on','WindowStyle','normal', 'Units', 'centimeter', 'Position', [0.5 1.5 fw*0.1 fh*0.1]);
        theme light
        drawnow
    end

% --- Helper functions --- %

%   % ------------------------ %
%   % --- Function: reprow --- %
%   % ------------------------ %
    function A = reprow(A,v)
        p = accumarray(cumsum([1;v]),1);
        A = A(cumsum(p(1:end-1)),:);
    end

%   % ------------------------ %
%   % --- Function: errinf --- %
%   % ------------------------ %
    function s = errinf
        s = sprintf('time=%d, ',t);
    end

%   % -------------------------- %
%   % --- Function: plotprog --- %
%   % -------------------------- %
    function plotprog
        N = Nf+Nm; % Total allele count
        setn
        if     t/T0<1e-3; if time_left>60; prog_res = 5000; else; prog_res = 2000; end
        elseif t/T0<0.01; if time_left>60; prog_res = 2000; else; prog_res = 1000; end
        elseif t/T0< 0.1; if time_left>60; prog_res = 1000; else; prog_res =  200; end
        else            ; if time_left>60; prog_res =  200; else; prog_res =  100; end
        end

        toc_tot = toc(tic_tot);
        if prog && (t==T || ceil((t+1)/T0*prog_res)>ceil(t/T0*prog_res)) && (toc_tot>1 || t>=t0)
            toc_last = toc(tic_last);
            tic_last = tic;
            generations_elapsed = t - t_last_report;
            gen_p_sec = generations_elapsed/max(toc_last,eps);
            time_left2 = (T-t)/gen_p_sec/60; % Remaining minutes, including cooldown

            if t==T
                time_left = 0;
            elseif isinf(time_left) || ~time_left
                time_left = time_left2;
            else
                time_left = exp(0.75*log(time_left) + 0.25*log(time_left2)); % Gradually update estimated time left
            end
            t_last_report = t;

            % Round only the displayed rate, after computing remaining time.
            if gen_p_sec>=10
                gen_p_sec = round(gen_p_sec);
            else
                gen_p_sec = round(gen_p_sec,2,'significant');
            end

            per = t/T*100; % Include the mutation-free cooldown

            if toc_tot<60
                pltimeS = sprintf('treefig: %4.0f/%4.0f sec', [sum(tocv_plot) toc_tot]      );
            elseif toc_tot<2*60*60
                pltimeS = sprintf('treefig: %4.2g/%4.3g min', [sum(tocv_plot) toc_tot]/60   );
            else
                pltimeS = sprintf('treefig: %4.2g/%4.3g h  ', [sum(tocv_plot) toc_tot]/60/60);
            end
            dtime_digits = max(0,2-max(0,floor(log10(toc_last)))+(toc_last<=0.01)); % digits for presenting delta time
            pstr1 = ['%6.2f%%, t=%-6d, g/sec=%4g, (n,ne)=(%2d,%5.2f), ' char(916) 'T=%-5.' num2str(dtime_digits) 'fsec, Tl=%3.' num2str(time_left<9.95) 'f min, '];
            pvar1 = {per, t, gen_p_sec, n, ne, toc_last, time_left};
            pstr2 = '%s [%s]%s';
            pvar2 = {pltimeS, char(datetime('now','Format','HH.mm')), ism_str};
            ism_str=[]; % Reset info string for increased storage and matrix size

            printstr = [pstr1 pstr2];
            printvar = [pvar1 pvar2];

            prog_str = strrep(strrep(strrep(sprintf(printstr,printvar{:}),'e-0','e-'),'e+0','e'),'.0e','e');
            disp(prog_str);


        end
        if t > T-Tsamp
            oi = oi+1;
            nxv(oi,:) = [n, ne, xr];
        end
        if sum(t==tplv)
            nsi = find(t==tplv(:));
            nXstore(nsi,:) = [n, ne, xr];
        end

        % ---- Plot figures ---- %
        plot_toc = toc(tic_plot);
        plot_time_ok = isempty(tocv_plot) || plot_toc*0.025*(6-5*t/T) > min(tocv_plot(max(1,end-3):end)); % Spend no more than 2.5-15% on plotting

        if t==T || any(t==tplv) && prpl && t>t0 && plot_time_ok
            tic_plot = tic;

            % --- Plot allele tree --- %
            tlines =    [T-Tcd Tevod];
            if firstplot
                fixfig
                firstplot = false;
            end
            SimPlot(dr,figh,fs,nsi,Nth,nXstore,tlines,tplv,xpih,t0)
            drawnow

            tocv_plot = [tocv_plot, toc(tic_plot)];

        end
    end

% === Change log === %
% 260216: Fix dr range to consistently be (-dr/2,dr/2)

end
