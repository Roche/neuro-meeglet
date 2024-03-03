% Create unit test data for comparison between matlab and Python

addpath('./matlab')

fprintf('MATLAB version\n')

[stauts, data_path] = system('python -c "import mne; print(mne.datasets.testing.data_path(), end=str())"')


fname = strcat(data_path, '/MEG', '/sample', '/sample_audvis_trunc_raw.fif');

raw = fiff_setup_read_raw(fname);

picks_eeg = fiff_pick_types(raw.info, 0, 1);
picks_eeg(53) = [];

[data, times] = fiff_read_raw_segment(raw, raw.first_samp, raw.last_samp, picks_eeg);

data = data - mean(data, 1);  %average ref for comparability with Python example
data = data * 1e6

plot(data')

set(gcf, 'Units', 'inches');
screenPosition = get(gcf, 'Position');
set(gcf, 'Position', [screenPosition(1), screenPosition(2), 8, 4]);
ylabel('EEG [${\mu V}$]', 'Interpreter','latex')
xlabel('Time [ms]')
xlim([0, 6000])

cfg1.fsample = raw.info.sfreq;
cfg1.output = {'pow', 'csd', 'cov', 'coh', 'icoh', 'plv', 'pli', 'dwpli', 'r_plain', 'r_orth', 'gim'};
cfg2.density = 'oct';
cfg2 = cfg;
cfg2.density = 'Hz';
cfg1
cfg2


out1 = ro_freq_meeglet(data, cfg1);

out2 = ro_freq_meeglet(data, cfg2);

save('./nbs/api/data/mne_meeglet_testing_data.mat', 'out1', 'out2')