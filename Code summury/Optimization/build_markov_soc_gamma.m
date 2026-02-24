function gamma_M = build_markov_soc_gamma(params, gamma_drive)
% ——HGV——
def = struct( ...
    'S',21, ...
    's_min',0.20, ...          % 触发阈值（SoC<=s_min 需快充）
    'Eveh_kWh',450, ...        % HGV 电池 400–600 kWh
    'e_kWh_per_km',1.6, ...    % HGV 能耗 1.2–2.0 kWh/km
    'mean_day_km',280, ...     % HGV 日里程
    'dt_hours',1, ...
    'soc_init_mu',0.80, 'soc_init_sigma',0.05, ...
    's_refill',0.99 ...        % 回充目标 SoC
);
fn = fieldnames(def);
for k=1:numel(fn), if ~isfield(params,fn{k}), params.(fn{k}) = def.(fn{k}); end, end

H = numel(gamma_drive);   S = params.S;
soc = linspace(0,1,S).';
p   = normpdf(soc, params.soc_init_mu, params.soc_init_sigma);
p   = p / sum(p);

km_h   = params.mean_day_km * gamma_drive(:);
dSoC_h = (params.e_kWh_per_km .* km_h) / params.Eveh_kWh;

arr = zeros(H,1);
for h = 1:H
    d = dSoC_h(h);
    soc_after = max(0, soc - d);
    p_drive = interp1(soc, p, soc_after, 'linear', 0);
    p_drive = p_drive / max(sum(p_drive),eps);

    low = (soc <= params.s_min);
    arr(h) = sum(p_drive(low));              % 本小时触发快充
    mass = arr(h);

    p_keep = p_drive; p_keep(low)=0; p_keep = p_keep/max(sum(p_keep),1);
    [~,k0] = min(abs(soc - params.s_refill));
    add = zeros(S,1); add(k0)=mass;

    p = p_keep*(1-mass) + add;
    p = p / max(sum(p),eps);
end


S_arr = sum(arr);
if S_arr <= 1e-6 || nnz(arr > 1e-6) <= 3
    gamma_M = gamma_drive(:) / sum(gamma_drive);
else
    gamma_M = arr / S_arr;
end
end
