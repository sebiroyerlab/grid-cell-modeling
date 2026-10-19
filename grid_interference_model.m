%function [thetaF, grid_score, grid_size] = grid_interference_model()
%
% thetaF: Cells' theta frequencies based on inter-peak-intervals within 50-300ms time range.
% grid_score: Cells' grid scores based on peak/through-ratio of ACGs (using ACGs of firing rate envelops)
% grid_size: Cells' grid scores based on inter-peak-intervals of firing rate envelops (beat size)
%
% To select frequency settings for PI, Cue or OF, uncomment the appropriate line in the 'parameters' section
% Notice: each model run will output a slightly different figure because of the varying random seeds 

function [thetaF, grid_score, grid_size] = grid_interference_model()




% --- Parameters ---
n_OC1 = 100; n_OC2 = 100; n_OC = n_OC1+n_OC2;       % Number of oscillatory cells (OCs)
n_GC = 100;                                         % Number of grid cells

%f1 = 7.1; f2 = 9.3;                                 % PI theta frequencies for the two populations of OCs
%f1 = 7.2; f2 = 8.7;                                 % Cue theta frequencies for the two populations of OCs
f1 = 8.4; f2 = 9.5;                                 % OF theta frequencies for the two populations of OCs

time_step = 0.001;                                  % time_step (s)
t = 0:time_step:100;                                % Time vector (s)
npt = length(t);

Initial_gamma = 0.005;                              % Gamma shape parameter for the random distribution of OC-to-GC synaptic weights 
alpha = 1;                                          % Global inhibition





% --- OC-to-GC synaptic weights ---
GCxOC_Wi = random('gam', Initial_gamma, ones(n_GC, n_OC));
GCxOC_W  = GCxOC_Wi ./ repmat(sum(GCxOC_Wi,2), 1, n_OC);        % Synaptic weight normalization 


% --- Activity of OCs versus t ---
sigma = 0.001;                                      % set a small frequency variability within a group of OC. As the value is very small, there is practically no variability within a group.
delays = repmat(0.01:0.01:0.1,1,100);               % time shifts (s) for individual oscillatory cells

OC1xt = zeros(n_OC1, npt);
for ii = 1:n_OC1
    f = f1 * (1 + sigma*randn(1));
    d = delays(ii);
    OC1xt(ii,:) = sin((t + d) .* f * 2*pi());
end
OC2xt = zeros(n_OC2, npt);
for ii = 1:n_OC2
    f = f2 * (1 + sigma*randn(1));
    d = delays(ii);
    OC2xt(ii,:) = sin((t + d) .* f * 2*pi());
end
OCxt = [OC1xt;OC2xt];


% --- Activity of GC versus t ---
EE  = GCxOC_W * OCxt;                                           % Excitation received by GCs (n_GC x npt)
thr = zeros(1, npt);
for tt = 1:npt                                                  % Numerical computation of inhibition threshold
    I_test_range = 0 : max(EE(:,tt))/1000 : max(EE(:,tt));
    E_I_mat = repmat(EE(:,tt),1,length(I_test_range)) - repmat(I_test_range,n_GC,1);
    tmp = I_test_range - alpha*sum(E_I_mat.*(E_I_mat>0));
    [~,imin] = min(abs(tmp));
    thr(tt) = I_test_range(imin);
end
inhibition_thr = repmat(thr, n_GC, 1);
GCxt = (EE - inhibition_thr) .* (EE > inhibition_thr);





% --- Computation of theta frequency, grid score, grid size ---
thetaF       = zeros(1, n_GC);
grid_score   = zeros(1, n_GC);
grid_size    = zeros(1, n_GC);
grid_envelop = zeros(npt, n_GC);
grid_acg     = zeros(2*npt-1,n_GC);

for ii = 1:n_GC
    % Theta frequency from inter-peak intervals in the 0.05–0.3s range
    tmp = GCxt(ii,:).*(GCxt(ii,:)>max(GCxt(ii,:))/10);
    [peaks,t_peaks] = findpeaks(tmp,t);
    Dt_peaks = diff(t_peaks);
    thetaF(ii) = 1/median(Dt_peaks(Dt_peaks>0.05 & Dt_peaks<0.3));

    % grid score, grid size
    if size(peaks,2)>=6
        
        %computation of cells' firing rate envelops
        ndx = find(Dt_peaks>0.3);
        new_t = t_peaks(ndx) + Dt_peaks(ndx)/2;
        t_peaks = [t_peaks new_t];
        peaks = [peaks zeros(1,length(new_t))];
        [t_peaks,isort] = sort(t_peaks);
        peaks = peaks(isort);
        envelop = interp1(t_peaks,peaks,t,'spline');  
        grid_envelop(:,ii) = envelop;
        
        %spatial autocorrelograms of envelops
        acg = xcorr(envelop,envelop,'normalized'); 
        grid_acg(:,ii) = acg; 

        %computation of grid size and grid score
        [grid_peaks,i_grid_peaks] = findpeaks(acg,1:length(acg));        
        if length(grid_peaks) >= 5
            Di_grid_peaks = diff(i_grid_peaks);
            grid_size_i = round(median(Di_grid_peaks));
            grid_size(ii) = grid_size_i*time_step;                      %grid size

            ndx1 = (length(acg)-1)/2;
            ndx2 = ndx1 - grid_size_i;
            ndx3 = ndx2 + round(grid_size_i/2);
            grid_score(ii) = acg(ndx2) - acg(ndx3)*(acg(ndx3)>0);       %grid score    
        end

    end

end


% --- make figures ---

grid_score_thr = 0.2;

figure;
[~,isort] = sort(grid_score,'descend');
for ii = 1:3
    subplot(6,1,ii);plot(t,GCxt(isort(ii),:),'b',t,grid_envelop(:,isort(ii)),'r');xlim([10 16])
    if ii == 1; title('For 3 cell examples, firing rate (blue) and envelop (red)');end
end
ylabel('firing rate (a.u.)');xlabel('time (s)')


subplot(2,4,5)
tmp = grid_envelop(10000:15000,grid_score>grid_score_thr);
[~,imax] = max(tmp(1:1000,:),[],1);
[~,isort] = sort(imax);
tmp = tmp(:,isort);
imagesc(matnorm(tmp',2)); 
set(gca,'XTick')
xlabel('time (ms)'); ylabel('grid cell #, sorted by peak position');title('grid envelop')

subplot(2,4,6)
gridC = grid_acg(:,grid_score>grid_score_thr);
gridC_size = grid_size(grid_score>grid_score_thr);
[~,isort] = sort(gridC_size);
imagesc(matnorm(gridC(:,isort)',2)); 
xlim([size(grid_acg,1)/2-2000 size(grid_acg,1)/2+2000]); 
xlabel('time (ms)'); ylabel('grid cell #, sorted by grid size');title('envelop ACG')

subplot(2,4,7)
thetaF = thetaF(grid_score>grid_score_thr);
bar(1,mean(thetaF),'FaceColor','none');hold on
errorbar(1,mean(thetaF),std(thetaF)/sqrt(numel(thetaF)),'k.');
bar(2,std(thetaF),'k')
title('theta frequency');xlabel('mean std');ylabel('Hz');

subplot(2,4,8)
grid_size = grid_size(grid_score>grid_score_thr);
bar(1,mean(grid_size),'FaceColor','none');hold on
errorbar(1,mean(grid_size),std(grid_size)/sqrt(numel(grid_size)),'k.');
ylim([0 1.2])
title('grid size (s)');ylabel('s');


end
