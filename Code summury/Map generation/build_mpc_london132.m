%% 
% 输入: london_132kv.mat 
% 输出: mpc_london132.mat 

clear; clc;

%% 参数
% mat_path = 'london_132kv.mat';
mat_path = 'london_132kv_bigger2.mat';
BASE_KV  = 132;    % 基准电压
S_BASE   = 100;    % 系统基准 MVA
F_HZ     = 50;
MIN_LEN_KM = 0.05; 
DROP_SINGLE_ENDED = true;

%% 读取数据
S = load(mat_path);

subs = cell_pick(S.substations, {'id','name','lat','lon','voltages_kV'});

L1   = cell_pick(S.lines, {'id','name','from_sub','to_sub', ...
                           'voltage_kV','circuits','is_cable','length_km'});

Lx   = cell_pick(S.bridges, {'id','name','from_sub','to_sub', ...
                             'voltage_kV','circuits','is_cable','length_km'});   % ← 加了 circuits


L1 = clean_edges(L1, MIN_LEN_KM, DROP_SINGLE_ENDED);
Lx = clean_edges(Lx, MIN_LEN_KM, DROP_SINGLE_ENDED);

keep_ids = keep_components_with_132(L1, Lx, {subs.id});
keep_ids = cellstr(keep_ids);

subs = subs(ismember({subs.id}, keep_ids));
L1   = L1(ismember({L1.from_sub}, keep_ids) & ismember({L1.to_sub}, keep_ids));
Lx   = Lx(ismember({Lx.from_sub}, keep_ids) & ismember({Lx.to_sub}, keep_ids));

%% 建立 bus 索引
bus_id = cellfun(@char, {subs.id}, 'UniformOutput', false);
nbus   = numel(bus_id);
bus_map = containers.Map(bus_id, num2cell(1:nbus));


% R,X [ohm/km], C [uF/km]
param.overhead.v11  = [0.30 0.60 0.02];
param.overhead.v33  = [0.25 0.55 0.02];
param.overhead.v66  = [0.20 0.45 0.02];
param.overhead.v110 = [0.12 0.40 0.015];
param.overhead.v132 = [0.08 0.35 0.010];
param.overhead.v220 = [0.06 0.30 0.010];
param.overhead.v275 = [0.04 0.25 0.010];
param.overhead.v400 = [0.03 0.24 0.010];

param.cable.v11  = [0.15 0.30 0.20];
param.cable.v33  = [0.10 0.25 0.20];
param.cable.v66  = [0.08 0.22 0.25];
param.cable.v110 = [0.07 0.20 0.30];
param.cable.v132 = [0.06 0.27 0.35];
param.cable.v220 = [0.05 0.25 0.40];
param.cable.v275 = [0.04 0.22 0.45];
param.cable.v400 = [0.03 0.20 0.50];

Zbase = (BASE_KV^2)/S_BASE;
Ybase = 1/Zbase;

%% 组装 branch 
branches = [mark_core(L1,true); mark_core(Lx,false)];
branch = [];
kvs_ref = double([11 33 66 110 132 220 275 400]); 

for k = 1:numel(branches)
    b = branches(k);

    % 电压（kV）
    kv = round(asdouble(getf(b,'voltage_kV'), BASE_KV));
    if ~(isfinite(kv) && kv > 0), kv = BASE_KV; end

    typ = tern(is_true(getf(b,'is_cable')), 'cable', 'overhead');

    key = sprintf('v%d', kv);
    if ~isfield(param.(typ), key)
        [~,ix] = min(abs(kvs_ref - kv));
        kv = kvs_ref(ix);
        key = sprintf('v%d', kv);
    end
    pkm = param.(typ).(key);  % [R X C] per km

    % 线路长度（km）
    len_km = asdouble(getf(b,'length_km'), 0);

    R_ohm = pkm(1) * len_km;
    X_ohm = pkm(2) * len_km;
    C_uF  = pkm(3) * len_km;

    % 折算到 132kV 标幺
    B_S  = 2*pi*F_HZ * (C_uF*1e-6);  % siemens
    r    = R_ohm / Zbase;
    x    = X_ohm / Zbase;
    b_sh = B_S  / Ybase;

    % 端点母线号
    f = bus_map(char(getf(b,'from_sub')));
    t = bus_map(char(getf(b,'to_sub')));

    % MATPOWER branch 
    branch = [branch; [f t r x b_sh 999 999 999 0 0 1 -360 360]]; 
    branches(k).voltage_kV = kv;
    branches(k).length_km  = len_km;
    branches(k).is_cable   = is_true(getf(b,'is_cable'));
end

%% —— 给 branch 做映射（1..M）——
M = size(branch,1);
branch_map_bigger2_tbl = table( ...
    (1:M)', ...                         % branch_index
    branch(:,1), branch(:,2), ...       % fbus, tbus (已是 1..N 的母线编号)
    string({branches.from_sub}'), ...   % from_sub_osm
    string({branches.to_sub}'), ...     % to_sub_osm
    string({branches.id}'), ...         
    string({branches.name}'), ...       
    [branches.voltage_kV]', ...         % 电压等级 (kV)
    logical([branches.is_cable]'), ...  % 是否电缆
    logical([branches.is_core]'), ...   % 是否 132kV 主网
    [branches.length_km]', ...          % 长度 (km)
    'VariableNames', {'branch_index','fbus','tbus', ...
        'from_sub_osm','to_sub_osm','osm_line_id','name', ...
        'voltage_kV','is_cable','is_core','length_km'} ...
);
save branch_mapping_bigger2.mat branch_map_bigger2_tbl
disp(head(branch_map_bigger2_tbl,10))


%% bus 矩阵
bus = zeros(nbus,13);
for i=1:nbus
    bus(i,:) = [i 1 0 0 0 0 1 1.0 0 BASE_KV 1 1.1 0.9];
end
bus(1,2)=3; % slack

%% === generators ===
gen = [];
pf = 0.95;  

if isfield(S, 'generators') && ~isempty(S.generators)
    G = cell_pick(S.generators, {'at_sub','Pmax_MW'});
    % 聚合到母线
    cap_map = containers.Map('KeyType','double','ValueType','double');
    for i = 1:numel(G)
        sub_id = char(G(i).at_sub);
        if isempty(sub_id), continue; end
        if ~isKey(bus_map, sub_id), continue; end
        bi = bus_map(sub_id);
        Pmax = asdouble(G(i).Pmax_MW, 0);
        if ~(isfinite(Pmax) && Pmax > 0), continue; end
        if isKey(cap_map, bi), cap_map(bi) = cap_map(bi) + Pmax;
        else, cap_map(bi) = Pmax; end
    end

    if cap_map.Count > 0
        keys_bus = cell2mat(cap_map.keys);
        Pmaxs    = cell2mat(cap_map.values);     % MW
        Smaxs    = Pmaxs ./ max(pf, 1e-3);       % 近似 MVA 上限
        Qmaxs    = sqrt(max(Smaxs.^2 - Pmaxs.^2, 0));  % 无功上限
        
        gen = zeros(numel(keys_bus), 10);
        for k = 1:numel(keys_bus)
            bi = keys_bus(k);
            % [GEN_BUS PG QG QMAX QMIN VG MBASE GEN_STATUS PMAX PMIN]
            gen(k,:) = [bi, 0, 0, Qmaxs(k), -Qmaxs(k), 1.0, 100, 1, Pmaxs(k), 0];
        end
        bus(:,2) = 1;                        % PQ
        slack_bus = keys_bus(1);
        bus(slack_bus,2) = 3;
        others = setdiff(keys_bus, slack_bus);
        if ~isempty(others), bus(others,2) = 2; end
    end
end

if isempty(gen)
    bus(:,2) = 1; bus(1,2) = 3;
    gen = [1 0 0 999 -999 1.0 100 1 999 0];
end


%% 生成 mpc
mpc = struct();
mpc.baseMVA = S_BASE;
mpc.bus     = bus;
mpc.branch  = branch;
mpc.gen     = gen;


%% ===== graph =====
G = graph(mpc.branch(:,1), mpc.branch(:,2));
comp = conncomp(G);

slack_bus = find(mpc.bus(:,2)==3, 1);
if isempty(slack_bus)
    [~, keep_comp] = max(histcounts(comp, 1:max(comp)+1));
else
    keep_comp = comp(slack_bus);
end

keep_bus_mask = (comp == keep_comp);

% 过滤支路
keep_branch = keep_bus_mask(mpc.branch(:,1)) & keep_bus_mask(mpc.branch(:,2));
mpc.branch = mpc.branch(keep_branch,:);

% 压缩母线索引（1..N）
old_idx = find(keep_bus_mask);                 % 原 bus 索引（升序）
newN    = numel(old_idx);
old2new = zeros(size(mpc.bus,1),1);
old2new(old_idx) = 1:newN;

mpc.bus = mpc.bus(keep_bus_mask,:);
mpc.branch(:,1:2) = old2new(mpc.branch(:,1:2));

mpc.bus(:,1) = (1:newN)';

if ~isempty(mpc.gen)
    keep_gen = keep_bus_mask(mpc.gen(:,1));
    mpc.gen  = mpc.gen(keep_gen,:);
    if ~isempty(mpc.gen)
        mpc.gen(:,1) = old2new(mpc.gen(:,1));
    end
end

if ~any(mpc.bus(:,2)==3)
    mpc.bus(:,2) = 1;     % 全部先设 PQ
    mpc.bus(1,2) = 3;     % 第一个母线设 Slack
end


%% === 分配负荷 ===
kept_subs = subs(keep_bus_mask);               
bus_map_tbl = table( (1:size(mpc.bus,1))', ...
                     string({kept_subs.id})', string({kept_subs.name})', ...
                     'VariableNames',{'bus_index','osm_id','name'});

PD = zeros(size(mpc.bus,1),1);
try
    Ssub = load('london_132kv.mat','substations_table');   
    ST   = struct2table(Ssub.substations_table);
    J = outerjoin(bus_map_tbl, ST, 'LeftKeys','osm_id','RightKeys','id', ...
                  'Type','left','MergeKeys',true);

    if ismember('firm_capacity_MVA', J.Properties.VariableNames) && ...
       ismember('demand_headroom_MVA', J.Properties.VariableNames)
        PD = max(J.firm_capacity_MVA - J.demand_headroom_MVA, 0);
    elseif ismember('ecr_import_MW', J.Properties.VariableNames)
        PD = J.ecr_import_MW;
    elseif ismember('firm_capacity_MVA_est', J.Properties.VariableNames)
        PD = 0.60 * J.firm_capacity_MVA_est;
    else
        PD = 10*ones(size(PD)); 
    end
    PD(isnan(PD)) = 0;
catch ME
    warning('Load synthesis fell back to default: %s', ME.message);
    PD = 10*ones(size(PD));
end

target_total_MW = 1200;
s = sum(PD); if s>0, PD = PD * (target_total_MW/s); else, PD(:) = target_total_MW/numel(PD); end

pf = 0.96;
QD = PD * tan(acos(pf));

% 写入 mpc.bus
mpc.bus(:,3) = PD;
mpc.bus(:,4) = QD;

% 按 PMAX 比例分配计划出力
if ~isempty(mpc.gen)
    pmx = mpc.gen(:,9); pmx(pmx<=0) = 0;
    if sum(pmx)>0
        Ptot = sum(PD);
        share = pmx / sum(pmx);
        pg = min(pmx, Ptot*share);
        if size(mpc.gen,1)>=1
            slack_share = pg(1); others = pg; others(1)=0;
            mpc.gen(:,2) = others;                     
        end
    end
end






%% ===== 画图 =====

plot_mpc132_color(mpc);

% 保存 
save mpc_london132_bigger2.mat mpc

% 统计
fprintf('构建完成: bus=%d, branch=%d, gen=%d\n', size(mpc.bus,1), size(mpc.branch,1), size(mpc.gen,1));

%% 保存映射表
kept_subs = subs(keep_bus_mask);  
bus_map_bigger2_tbl = table( (1:size(mpc.bus,1))', ...
                     string({kept_subs.id})', string({kept_subs.name})', ...
                     'VariableNames',{'bus_index','osm_id','name'});
save bus_mapping_bigger2.mat bus_map_bigger2_tbl
disp(head(bus_map_tbl,45))



function A = cell_pick(C, fields)
if isempty(C), A = struct([]); return; end
if isstruct(C), C = num2cell(C); end
n = numel(C);
tmpl = cell2struct(cell(size(fields)), fields, 2);  
A    = repmat(tmpl, n, 1);
for i = 1:n
    s = C{i};
    if ~isstruct(s), s = struct(); end
    for f = 1:numel(fields)
        fn = fields{f};
        if isfield(s, fn)
            A(i).(fn) = s.(fn);
        else
            A(i).(fn) = default_for_field(fn);
        end
    end
end
end

function v = default_for_field(fn)
if ismember(fn, {'id','name','from_sub','to_sub','voltages_kV'})
    v = "";                      
elseif ismember(fn, {'is_cable'})
    v = false;                 
else
    v = 0;                       
end
end


function E = clean_edges(E, min_len_km, drop_single)
if isempty(E), return; end
keep = true(size(E));
for i = 1:numel(E)
    a = getf(E(i),'from_sub'); b = getf(E(i),'to_sub');
    if isempty(a) || isempty(b)
        keep(i) = ~drop_single && (~isempty(a) || ~isempty(b));
        continue
    end
    if strcmp(a,b), keep(i) = false; continue; end
    if isfield(E,'length_km') && ~isempty(E(i).length_km) && E(i).length_km < min_len_km
        keep(i) = false;
    end
    if ~isfield(E,'voltage_kV'), E(i).voltage_kV = 132; end
    if ~isfield(E,'is_core'),    E(i).is_core    = false; end
    if ~isfield(E,'is_cable'),   E(i).is_cable   = false; end
end
E = E(keep);
end


function x = asdouble(v, defaultVal)
if nargin < 2, defaultVal = NaN; end
if isnumeric(v)
    x = double(v);
elseif ischar(v) || (isstring(v) && isscalar(v))
    s = regexprep(char(v), '[^0-9.\-eE]', ''); 
    if isempty(s), x = defaultVal; else, x = str2double(s); if isnan(x), x = defaultVal; end, end
else
    x = defaultVal;
end
end

function tf = is_true(v)
if islogical(v)
    tf = v;
elseif isnumeric(v)
    tf = any(v ~= 0);
elseif ischar(v) || (isstring(v) && isscalar(v))
    s = lower(strtrim(char(v)));
    tf = any(strcmp(s, {'1','true','yes','y'}));
else
    tf = false;
end
end


function keep_ids = keep_components_with_132(L1, Lx, all_ids)

% 1) 边表
edges = {};
for i = 1:numel(L1)
    a = getf(L1(i),'from_sub'); b = getf(L1(i),'to_sub');
    if ~isempty(a) && ~isempty(b), edges(end+1,:) = {char(a), char(b)}; end 
for i = 1:numel(Lx)
    a = getf(Lx(i),'from_sub'); b = getf(Lx(i),'to_sub');
    if ~isempty(a) && ~isempty(b), edges(end+1,:) = {char(a), char(b)}; end 
end

if isempty(edges)
    if iscell(all_ids), keep_ids = all_ids(:)'; else, keep_ids = cellstr(all_ids); end
    return
end

% 2) 建图
Gtbl = cell2table(edges, 'VariableNames', {'s','t'});
uu = [Gtbl.s; Gtbl.t];
if ~iscellstr(uu)
    for k = 1:numel(uu)
        if ~ischar(uu{k}), uu{k} = char(uu{k}); end
    end
end
[uids,~,idx] = unique(uu, 'stable');
G = graph(idx(1:height(Gtbl)), idx(height(Gtbl)+1:end));
comp = conncomp(G);


edge_has132 = false(1, numedges(G));
for i = 1:numel(L1)
    a = char(getf(L1(i),'from_sub')); b = char(getf(L1(i),'to_sub'));
    ia = find(strcmp(uids, a), 1); ib = find(strcmp(uids, b), 1);
    if isempty(ia) || isempty(ib), continue; end
    e = findedge(G, ia, ib);
    if e > 0, edge_has132(e) = true; end
end


comp_has132 = false(1, max(comp));
[si, ti] = findedge(G);
for e = 1:numel(edge_has132)
    if edge_has132(e)
        comp_has132(comp(si(e))) = true;
        comp_has132(comp(ti(e))) = true;
    end
end
keep_nodes = uids(comp_has132(comp));
keep_ids = keep_nodes(:)';   % cellstr
end


function E = mark_core(E, iscore)
for i = 1:numel(E), E(i).is_core = iscore; end
end

function y = tern(cond, a, b)
if cond, y = a; else, y = b; end
end

function v = getf(s, f)
if isfield(s, f), v = s.(f); else, v = []; end
end




