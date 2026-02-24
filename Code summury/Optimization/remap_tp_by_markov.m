% 用 NTS0501 + Markov(HGV参数)
% 写出 traffic_profiles_fast_new.csv
baseDir = 'C:\2023Research\Emobility\transportationmap\optimization3';
cd(baseDir); addpath(baseDir);
p = @(f) fullfile(baseDir,f);

% 1)  TP
TPold = readtable(p('traffic_profiles_fast.csv'), 'VariableNamingRule','preserve');
assert(all(ismember(["r_id","h","lambda_rh"], string(TPold.Properties.VariableNames))), ...
    'traffic_profiles_fast.csv 需包含列 r_id,h,lambda_rh');

% 2) NTS 
g  = read_nts0501_hour_shape(p('nts0501.xlsx'),[],[],[]);

% 3) HGV Markov 参数
params = struct();                   
gm = build_markov_soc_gamma(params, g);


% 4) 把日总分配到 24 小时
Rlist = unique(TPold.r_id);
rows_r = []; rows_h = []; rows_lam = [];
for rr = Rlist.'
    rows   = TPold(TPold.r_id==rr, :);
    totday = nansum(double(rows.lambda_rh)); if ~isfinite(totday), totday = 0; end
    lam_h  = totday * gm(:);                  % 24×1
    rows_r = [rows_r; repmat(rr,24,1)];
    rows_h = [rows_h; (1:24).'];
    rows_lam = [rows_lam; lam_h];
end

TPnew = table(rows_r, rows_h, rows_lam, 'VariableNames', {'r_id','h','lambda_rh'});
writetable(TPnew, p('traffic_profiles_fast_new.csv'));
fprintf('[Remap] DONE → traffic_profiles_fast_new.csv | rows=%d | unique r=%d\n', ...
        height(TPnew), numel(Rlist));


net_h = groupsummary(TPnew,"h","sum","lambda_rh");
disp(net_h(:,["h","sum_lambda_rh"]));

