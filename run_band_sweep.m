%RUN_BAND_SWEEP  Sensitivity of the finalized design to the dual-mode hand-off band.
%
%   The finalized design uses acq_band = 2.5 deg (turret_params.m), chosen from the
%   geometry of the saturated-rate unwind (paper Sec. IV-D).  This script checks that
%   choice: it re-runs scenarios A-C and 10 randomized seeds with other bands, every
%   other parameter frozen, and reports success, interception time, post-acquisition
%   miss and the number of times the controller re-enters the slew mode after the first
%   hand-off (capture chatter).  It is a sensitivity check ONLY -- the band is not tuned
%   per case and the default design is unchanged.
%
%   Writes results/band_sweep.txt
clear; close all; clc;
root = fileparts(mfilename('fullpath'));
resdir = fullfile(root,'results'); if ~isfolder(resdir), mkdir(resdir); end
bands_deg = [1.0 1.5 2.0 2.5 3.5 5.0];
cases = [{'A','B','C'}, arrayfun(@(s) struct('seed',s), 101:110, 'uni',0)];
fid = fopen(fullfile(resdir,'band_sweep.txt'),'w');
hdr = sprintf('%8s %9s %10s %12s %14s\n', 'band[deg]','success','t_hit[s]','RMSmiss[m]','slew re-entries');
fprintf(hdr); fprintf(fid, hdr);
for b = 1:numel(bands_deg)
    ok = 0; th = nan(numel(cases),1); rm = nan(numel(cases),1); re = zeros(numel(cases),1);
    for c = 1:numel(cases)
        L = simulate_case(cases{c}, struct('acq_band', deg2rad(bands_deg(b))));
        m = L.metrics; ok = ok + m.success; th(c) = m.t_first_hit; rm(c) = m.rms_miss_acq_buf;
        e = abs(L.thref - L.th) > deg2rad(bands_deg(b));         % true while in slew mode
        k = find(~e, 1);                                         % first hand-off
        if ~isempty(k), re(c) = sum(diff(e(k:end)) == 1); end    % later re-entries into slew
    end
    line = sprintf('%8.1f %6d/%-2d %10.2f %12.2f %14d\n', bands_deg(b), ok, numel(cases), mean(th,'omitnan'), mean(rm,'omitnan'), sum(re));
    fprintf(line); fprintf(fid, line);
end
fclose(fid);
disp('run_band_sweep done.');
