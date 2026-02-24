
clc; close all;
baseDir = 'C:\2023Research\Emobility\transportationmap\optimization3';
p = @(f) fullfile(baseDir, f);

%% ==== 文件名 ====
fnCand   = 'candidates_new.csv';
fnRmid   = 'mapping_edge_midpoints.csv';
fnPaths  = 'paths_or_corridors_new.csv';
fnTP     = 'traffic_profiles_fast_new.csv';
fn03     = '03_network-nodes.csv';
fn04     = '04_network-edges.csv';
fnSeq    = 'edge_sequences.csv';   

%% ==== 参数 ====
Eveh_kWh = 450;
E_session_kWh = 200;               % HGV 单次快充能量
% 1) alpha_EV
alpha_EV = 0.077;   

% 2) p_stop(h)：hgv_stop_prob_by_hour.csv（24×1）
Tstop = readtable(p('hgv_stop_prob_by_hour.csv'));
p_stop_h = zeros(1,24);  p_stop_h(Tstop.h+1) = double(Tstop.prob(:));
p_stop_h = p_stop_h / sum(p_stop_h);   


% 3) p_public
p_public = 0.25;

% 4) 分时 phi(h)
phi_h = alpha_EV * p_public * p_stop_h;     % 24×1，和 = alpha_EV*p_public
PHI_day = mean(phi_h);                      % 作为“日均等效 φ”用于总量对齐

R_BATT_KM   = 300;   % 电池续航
R_DETOUR_KM = 5;     % 绕行半径

OBJ_Soft  = 0.1;

% edge_sequences
USE_COVERAGE_WINDOWS = isfile(p(fnSeq));
WINDOW_KM            = R_BATT_KM;  


%% ==== 读表 ====
C  = readtable(p(fnCand), 'VariableNamingRule','preserve', 'TextType','string');
TP = readtable(p(fnTP),   'VariableNamingRule','preserve', 'TextType','string');
Rmid = readtable(p(fnRmid), 'VariableNamingRule','preserve', 'TextType','string');
Nodes = readtable(p(fn03), 'VariableNamingRule','preserve');
Edges = readtable(p(fn04), 'VariableNamingRule','preserve');
C.cand_id = string(C.cand_id);
C.lon       = double(C.lon);
C.lat       = double(C.lat);
nmax        = double(C.nmax);                      % 枪数上限
mu          = double(C.Prated_kW) ./ E_session_kWh; % (veh/h)/plug

PRUNE_SMALL = true;         % 开关
y_min_default = 10;         % 全局最小枪数

if PRUNE_SMALL
    keepI = (nmax >= y_min_default);
    C     = C(keepI, :);        % 裁剪候选站表
    nmax  = nmax(keepI);        % 同步裁剪
    mu    = mu(keepI);          % 同步裁剪
end



% 1) 规范列名
names = string(Rmid.Properties.VariableNames);
names = regexprep(names, '^\xEF\xBB\xBF', '');  
names = strtrim(names);
names = replace(names, char(160), ' ');         
Rmid.Properties.VariableNames = cellstr(names);

% 2) edge_id
cands = ["edge_id","r_id","id","Network_Edge_ID","Network_EdgeId"];
hit = cands( find(ismember(cands, string(Rmid.Properties.VariableNames)), 1, 'first') );

eid_raw = Rmid.(hit);
if isstring(eid_raw) || iscellstr(eid_raw)
    eid = double(str2double(string(eid_raw)));
elseif isnumeric(eid_raw)
    eid = double(eid_raw);
end
Rmid.edge_id = eid;   

% 3) 03/04 
needLL = ~ismember("lon", Rmid.Properties.VariableNames) || ~ismember("lat", Rmid.Properties.VariableNames) ...
         || all(ismissing(Rmid.lon)) || all(ismissing(Rmid.lat));
if needLL
    Nodes = readtable(p(fn03), 'VariableNamingRule','preserve');
    Edges = readtable(p(fn04), 'VariableNamingRule','preserve');
    keyLon = containers.Map(double(Nodes.Network_Node_ID), double(Nodes.Network_Node_X));
    keyLat = containers.Map(double(Nodes.Network_Node_ID), double(Nodes.Network_Node_Y));
    u = double(Edges.Network_Node_A_ID); v = double(Edges.Network_Node_B_ID);
    ok = arrayfun(@(a,b) isKey(keyLon,a)&&isKey(keyLon,b)&&isKey(keyLat,a)&&isKey(keyLat,b), u, v);
    mid_lon = nan(height(Edges),1); mid_lat = mid_lon;
    mid_lon(ok) = 0.5*(arrayfun(@(a) keyLon(a), u(ok)) + arrayfun(@(b) keyLon(b), v(ok)));
    mid_lat(ok) = 0.5*(arrayfun(@(a) keyLat(a), u(ok)) + arrayfun(@(b) keyLat(b), v(ok)));
    Tmid = table(double(Edges.Network_Edge_ID), mid_lon, mid_lat, ...
                 'VariableNames', {'edge_id','lon','lat'});
    Rmid = outerjoin(Rmid, Tmid, 'Keys','edge_id', 'MergeKeys',true, 'Type','left');
    if ismember('lon_Tmid', Rmid.Properties.VariableNames)
        m = ismissing(Rmid.lon);  Rmid.lon(m) = Rmid.lon_Tmid(m);
        m = ismissing(Rmid.lat);  Rmid.lat(m) = Rmid.lat_Tmid(m);
        Rmid.lon_Tmid = []; Rmid.lat_Tmid = [];
    end
end

disp(string(Rmid.Properties.VariableNames).');
disp(Rmid(1:min(3,height(Rmid)), {'edge_id','lon','lat'}));


% paths
PathTbl = table();
if isfile(p(fnPaths))
    PathTbl = readtable(p(fnPaths), 'VariableNamingRule','preserve', 'TextType','string');
    if ~all(ismember(["r_id","cand_ids"], PathTbl.Properties.VariableNames))
        warning('paths_or_corridors.csv'); PathTbl = table();
    end
end

% 类型规整
C.cand_id = string(C.cand_id);  C.bus_id = string(C.bus_id);
for k=["r_id","h","lambda_rh"]; assert(ismember(k, TP.Properties.VariableNames)); end
TP.r_id = double(TP.r_id); TP.h = double(TP.h); TP.lambda_rh = double(TP.lambda_rh);
ridCol = "";  for k=["edge_id","r_id","id"], if ismember(k, Rmid.Properties.VariableNames), ridCol=k; break; end; end

if ~ismember("lon", Rmid.Properties.VariableNames) || ~ismember("lat", Rmid.Properties.VariableNames) ...
   || all(ismissing(Rmid.lon)) || all(ismissing(Rmid.lat))
    keyLon = containers.Map(double(Nodes.Network_Node_ID), double(Nodes.Network_Node_X));
    keyLat = containers.Map(double(Nodes.Network_Node_ID), double(Nodes.Network_Node_Y));
    u = double(Edges.Network_Node_A_ID); v = double(Edges.Network_Node_B_ID);
    ok = arrayfun(@(a,b) isKey(keyLon,a)&&isKey(keyLon,b)&&isKey(keyLat,a)&&isKey(keyLat,b), u, v);
    mid_lon = nan(height(Edges),1); mid_lat = mid_lon;
    mid_lon(ok) = 0.5*(arrayfun(@(a) keyLon(a), u(ok)) + arrayfun(@(b) keyLon(b), v(ok)));
    mid_lat(ok) = 0.5*(arrayfun(@(a) keyLat(a), u(ok)) + arrayfun(@(b) keyLat(b), v(ok)));
    Tmid = table(double(Edges.Network_Edge_ID), mid_lon, mid_lat, ...
                 'VariableNames', {'edge_id','lon','lat'});
    Rmid = outerjoin(Rmid, Tmid, 'LeftKeys', ridCol, 'RightKeys', 'edge_id', ...
                     'MergeKeys', false, 'Type', 'left');
    if ~ismember("lon", Rmid.Properties.VariableNames), Rmid.lon = Rmid.mid_lon; end
    if ~ismember("lat", Rmid.Properties.VariableNames), Rmid.lat = Rmid.mid_lat; end
end

R_use = unique(TP.r_id);
Rmid  = Rmid(ismember(double(Rmid.(ridCol)), R_use), :);
r_ids = double(Rmid.(ridCol));  NR = numel(r_ids);

hav = @(la1,lo1,la2,lo2) 2*6371000*asin( sqrt( sin(deg2rad((la2-la1)/2)).^2 + ...
                          cosd(la1).*cosd(la2).*sin(deg2rad((lo2-lo1)/2)).^2 ) );

%% ==== 构造 Ω====
NI   = height(C);
clon = C.lon;                     % double 向量
clat = C.lat;

cid2row = containers.Map(cellstr(string(C.cand_id)), 1:NI);
useaths = ~isempty(PathTbl);
if usePaths
    if ~isa(PathTbl.cand_ids,'string'), PathTbl.cand_ids = string(PathTbl.cand_ids); end
    PathTbl.r_id_num = double(PathTbl.r_id);
end

Rdetour_m = R_DETOUR_KM * 1000;   
Omega  = cell(NR,1);
n_empty=0; n_nopath=0; n_pathsempty=0;

for ir = 1:NR
    rid  = r_ids(ir);
    lonr = double(Rmid.lon(ir)); latr = double(Rmid.lat(ir));

    % 绕行半径可达
    if any(isnan([lonr,latr]))
        ok_detour = false(NI,1);
    else
        d_m = hav(latr,lonr,clat,clon);
        ok_detour = (d_m <= Rdetour_m);
    end
    ok = ok_detour;
    if usePaths
        rows = PathTbl(PathTbl.r_id_num == rid, :);
        if isempty(rows)
            n_nopath = n_nopath + 1;
        else
            sraw = rows.cand_ids(1);
            if ismissing(sraw) || strlength(sraw)==0
                n_pathsempty = n_pathsempty + 1;
            else
                ids_str = split(sraw, ",");
                ids_str = strtrim(ids_str);
                ids_str = ids_str(~ismissing(ids_str));
                ids_str = ids_str(ids_str~="" & ids_str~="0" & lower(ids_str)~="nan");
                if isempty(ids_str)
                    n_pathsempty = n_pathsempty + 1;
                else
                    ok_paths = false(NI,1);
                    for t = 1:numel(ids_str)
                        cs = char(ids_str(t));
                        if isKey(cid2row, cs), ok_paths(cid2row(cs)) = true; end
                    end
                    ok = ok & ok_paths;   
                end
            end
        end
    end

    Omega{ir} = find(ok);
    if isempty(Omega{ir}), n_empty = n_empty + 1; end
end

fprintf('[Ω] detour=%g km | Ω(r)空=%d/%d | paths缺=%d | paths空=%d\n', ...
        R_DETOUR_KM, n_empty, NR, n_nopath, n_pathsempty);
    




%% ==== λ_{r,h}====
Hset = unique(double(TP.h));  NH = numel(Hset);

TP_base = TP;                         % base veh/h
TP_base.r_id = double(TP_base.r_id);
TP_base.h    = double(TP_base.h);
TP_base.lambda_rh = double(TP_base.lambda_rh);

F_edge_day = double(Edges.Traffic_flow_trucks_2030) / 365;   % veh/day
eid04      = double(Edges.Network_Edge_ID);

[~,ia,ib] = intersect(eid04, r_ids); 
F_edge_use = F_edge_day(ia);

Fr_od = zeros(numel(r_ids),1);       % （veh/day）
for k = 1:height(TP_base)
    ir = find(r_ids == TP_base.r_id(k), 1);
    if ~isempty(ir)
        Fr_od(ir) = Fr_od(ir) + TP_base.lambda_rh(k); % 累加 veh/h
    end
end
Fr_od_day = Fr_od * 24 / NH;         % 变成 veh/day
sum_od_day   = nansum(Fr_od_day(ib)); % 与 04 交集的日总量
sum_edge_day = nansum(F_edge_use);
alpha_scale  = sum_edge_day / max(sum_od_day, eps);

fprintf('[CAL] α(pre-PHI)=%.6f | ΣOD(day)=%.1f | ΣEDGE(day)=%.1f\n', ...
        alpha_scale, sum_od_day, sum_edge_day);


%构造 lamRH
lamRH  = zeros(NR, NH);
haveQ  = false(NR, NH);              
Lmask  = false(NR, NH);               % 用于一致性校验

h2idx = containers.Map('KeyType','double','ValueType','double');
for t = 1:NH, h2idx(Hset(t)) = t; end

for k = 1:height(TP_base)
    ir = find(r_ids == TP_base.r_id(k), 1);
    if isempty(ir), continue; end
    if ~isKey(h2idx, TP_base.h(k)), continue; end
    hh = h2idx(TP_base.h(k));
    lamRH(ir, hh) = lamRH(ir, hh) + (alpha_scale * phi_h(hh)) * TP_base.lambda_rh(k); % 分时phi
    Lmask(ir, hh) = true;
end




BAN_SMALL = true;           
isLarge   = (nmax >= y_min_default);

if BAN_SMALL
    for ir = 1:NR
        if ~isempty(Omega{ir})
            Omega{ir} = Omega{ir}( isLarge(Omega{ir}) );  
        end
    end
end

% ===== 兜底Ω=====
K_NEAR = 3; R_FALLBACK_KM = 15; Rfb_m = R_FALLBACK_KM*1000;
hasDemandR = any(lamRH > 1e-12, 2);

if BAN_SMALL
    idxLarge = find(isLarge);
else
    idxLarge = (1:NI).';    
end

for ir = 1:NR
    if isempty(Omega{ir}) && hasDemandR(ir)
        lonr = double(Rmid.lon(ir)); latr = double(Rmid.lat(ir));
        if ~(isnan(lonr) || isnan(latr)) && ~isempty(idxLarge)
            d_m = hav(latr, lonr, clat(idxLarge), clon(idxLarge));
            [ds,ord] = sort(d_m,'ascend');
            pickLarge = idxLarge(ord(ds <= Rfb_m));
            if isempty(pickLarge)
                pickLarge = idxLarge(ord(1:min(K_NEAR, numel(ord))));
            else
                pickLarge = pickLarge(1:min(K_NEAR, numel(pickLarge)));
            end
            Omega{ir} = unique([Omega{ir}; pickLarge(:)]);
        end
    end
end

n_empty2 = sum(cellfun(@isempty, Omega));
fprintf('[Ω-fix] after fallback | Ω(r)空=%d/%d (K=%d, R=%.1f km)\n', n_empty2, NR, K_NEAR, R_FALLBACK_KM);



for ir = 1:NR
    for hh = 1:NH
        if lamRH(ir,hh) > 1e-12 && ~isempty(Omega{ir})
            haveQ(ir,hh) = true;
        end
    end
end


Hset = unique(double(TP.h));  NH = numel(Hset);
fprintf('[DBG] TP hours unique=%d; min=%g, max=%g\n', NH, min(Hset), max(Hset));


in_use = ismember(TP_base.r_id, r_ids) & ismember(TP_base.h, Hset);
TP_used = TP_base(in_use, :);


hour_sum = zeros(NH,1); 
for k = 1:height(TP_used)
    hh = h2idx(TP_used.h(k));                   % 1..NH
    hour_sum(hh) = hour_sum(hh) + TP_used.lambda_rh(k);   % veh/h
end
TP_used_tot = sum(hour_sum);      % 总基数

%EV 日总
sum_phi_h = sum(phi_h);                          % = alpha_EV * p_public
expected_ev = alpha_scale * sum( hour_sum(:) .* phi_h(:) );  % 逐小时乘后相加

%实际 lamRH 日总
dem_per_h   = sum(lamRH, 1);                     % 1×NH（veh/h）
tot_day_ev  = sum(dem_per_h);                    % veh/day


Q     = cell(NR, NH);
Qidx  = cell(NR, NH);
nBad  = 0;

for ir = 1:NR
    idxI = Omega{ir};
    for hh = 1:NH
        lam = lamRH(ir,hh);        
        if lam <= 1e-9, continue; end
        if isempty(idxI)
            nBad = nBad + 1;
            continue;
        end
        Q{ir,hh}    = sdpvar(numel(idxI),1,'full');
        Qidx{ir,hh} = idxI(:);
       
    end
end
fprintf('[DIAG] empty Ω 且 λ>0 的格子数: %d\n', nBad);








%% =====================================================================
%  power flow iteration

USE_LOSS_PROXY = true;         
loss_pf_ev     = 1.0;          
rho_loss       = 1.0;          
maxIter_loss   = 6;           
tol_loss_rel   = 2e-2;          % 收敛阈值

mpc0 = [];
PTDF = [];
k_loss = [];
slack_bus_i = [];
cand2busrow = [];   % cand_id(i) 
Mbus = [];          

if USE_LOSS_PROXY
    try
        % --- 读取 mpc---
        if exist(p('build_mpc_london132.m'), 'file')
            run(p('build_mpc_london132.m'));   
            
            Aih = sdpvar(NI, NH, 'full');
            Aih(:,:) = 0;
        end
        if exist('mpc','var')
            mpc0 = mpc;
        elseif exist(p('mpc_london132_bigger2.mat'),'file')
            tmp = load(p('mpc_london132_bigger2.mat'));
            if isfield(tmp,'mpc'), mpc0 = tmp.mpc; end
        end
  
        try
            define_constants;
            BUS_I = 1; BUS_TYPE = 2; REF = 3;
            BR_R = 3;
        catch
            BUS_I = 1; BUS_TYPE = 2; REF = 3;
            BR_R = 3;
        end

        nb = size(mpc0.bus,1);
        nl = size(mpc0.branch,1);

        % --- slack bus ---
        ref_rows = find(mpc0.bus(:,BUS_TYPE) == REF);
        if isempty(ref_rows), ref_rows = 1; end
        slack_bus_i = ref_rows(1);

        % --- 候选点 cand  ---
          bus_i = double(C.bus_id);
        cand2busrow = zeros(NI,1);
        for i = 1:NI
            r = find(mpc0.bus(:,BUS_I) == bus_i(i), 1);
            if ~isempty(r)
                cand2busrow(i) = r;
            elseif bus_i(i) >= 1 && bus_i(i) <= nb
                cand2busrow(i) = bus_i(i);
            else
                cand2busrow(i) = 0; % 无法映射
            end
        end

        Mbus = sparse(cand2busrow(cand2busrow>0), find(cand2busrow>0), 1, nb, NI);

        % --- PTDF ---
     
        try
            PTDF = makePTDF(mpc0);
        catch
            PTDF = makePTDF(mpc0.baseMVA, mpc0.bus, mpc0.branch, slack_bus_i);
        end

        % rho_loss 
        k_loss = max(0, mpc0.branch(:,BR_R));
  
    catch ME
         USE_LOSS_PROXY = false;
    end
end

%% ==== YALMIP 变量====
yalmip('clear');
x   = binvar(NI,1,'full');         % 建站
y   = intvar(NI,1,'full');         % 枪数（0..nmax）


%% ==== 约束 ====
Cons = [];

%% ==== Charger 约束 ====
% 参数
y_min_default = 2;          % 全局下限
y_max_global  = 80;         % 单站最大
BAN_SMALL     = true;       % 禁止 nmax<y_min_default 的小站被选

nmax = double(C.nmax);    
small = (nmax < y_min_default);


if BAN_SMALL
    Cons = [Cons, x(small) == 0];
end

y_min_vec = min(y_min_default * ones(NI,1), nmax);

for i = 1:NI
    Cons = [Cons, y_min_vec(i)*x(i) <= y(i) <= min(nmax(i), y_max_global)*x(i)];
end


% 需求守恒
for ir = 1:NR
    for hh = 1:NH
        if haveQ(ir,hh)
            Cons = [Cons, sum(Q{ir,hh}) == lamRH(ir,hh), Q{ir,hh} >= 0];
        end
    end
end

% 站容量
for i = 1:NI
    for hh = 1:NH
        expr = 0;
        for ir = 1:NR
            if haveQ(ir,hh)
                loc = find(Qidx{ir,hh} == i);
                if ~isempty(loc), expr = expr + Q{ir,hh}(loc); end
            end
        end
        Cons = [Cons, expr <= mu(i) * y(i)];
    end
end


% a_{i,h} <= mu(i) * nmax(i) 上限约束
for i = 1:NI
  for hh = 1:NH
    aih = 0;
    for ir = 1:NR
      if haveQ(ir,hh)
        loc = find(Qidx{ir,hh}==i);
        if ~isempty(loc), aih = aih + Q{ir,hh}(loc); end
      end
    end
    Cons = [Cons, aih <= mu(i) * nmax(i)];
  end
end


%% ==== 排队约束====
beta_wait = 1;     % 越大等待越短
delta     = 1e-3;    % 位移量（避免 s=0 斜率无穷）
N_KNOTS   = 8;       % 切线数量
EPS       = 0;       % 现在允许从 0 开始

% 
gamma_slack = 1e5;
slacks = sdpvar(0,1);  

for i = 1:NI
    smax_i = max(1, double(nmax(i)));              % s 的上界不超过枪数上限
    % 保有一条切线在 0 处
    s_knots = unique([0, exp(linspace(log(delta), log(smax_i+delta), N_KNOTS)) - delta]);
    for hh = 1:NH
        % a_{i,h} = Σ_r q_{r,h,i}
        aih = 0;
        for ir = 1:NR
            if haveQ(ir,hh)
                loc = find(Qidx{ir,hh} == i);
                if ~isempty(loc), aih = aih + Q{ir,hh}(loc); end
            end
        end

        % [LOSS] 
        Aih(i,hh) = aih;

        % 连续变量：s_{i,h}
        s_ih = sdpvar(1,1);
        t_ih = sdpvar(1,1);
        Cons = [Cons, s_ih >= 0, t_ih >= 0];

        % s >= a/mu(i)
        if mu(i) > 0
            Cons = [Cons, mu(i) * s_ih >= aih];
        else
            %无服务率时不允许有到达
            Cons = [Cons, aih == 0, s_ih == 0, t_ih == 0];
            continue
        end

        % 切线t >= m_k*s + b_k
        for kk = 1:numel(s_knots)
            sk = s_knots(kk);
            mk = 1/(2*sqrt(sk + delta));                        % g'(sk)
            bk = (sqrt(sk + delta) - sqrt(delta)) - mk*sk;      % 截距
            Cons = [Cons, t_ih >= mk*s_ih + bk];
        end
        
        % ... 在 for i/h 循环里
        u_ih = sdpvar(1,1);  Cons = [Cons, u_ih >= 0];
        Cons = [Cons, y(i) >= s_ih + beta_wait*t_ih - u_ih];
        slacks = [slacks; u_ih];
    end
end


Obj=0;
%% ==== 续航滑窗覆盖约束====
if USE_COVERAGE_WINDOWS
    fprintf('edge_sequences.csv ');
    Seq = readtable(p(fnSeq), 'VariableNamingRule','preserve','TextType','string'); % route_id, edge_ids
    win_km = WINDOW_KM;   % 窗口长度 = 续航（km）

    % 04 的长度列（单位：km）
    eid04   = double(Edges.Network_Edge_ID);
    dist_km = double(Edges.Distance);
    eid2row = containers.Map(num2cell(eid04), 1:numel(eid04));

    % 软松弛设置
    COV_SOFT  = true;     
    gamma_cov = 1e4;      

    nEmptyWin   = 0;      
    nNoLargeWin = 0;      
    nSkipWin    = 0;     

    for rr = 1:height(Seq)
        ids = double(str2double(split(Seq.edge_ids(rr), ',')));
        ids = ids(~isnan(ids));
        if isempty(ids), continue; end

        % 取距离与 r_ids 交集
        d = zeros(numel(ids),1);
        rid_mask = false(numel(ids),1);
        for k = 1:numel(ids)
            if isKey(eid2row, ids(k)), d(k) = dist_km(eid2row(ids(k))); else, d(k)=0; end
            rid_mask(k) = ismember(ids(k), r_ids);
        end
        ids2 = ids(rid_mask);
        d2   = d(rid_mask);
        if isempty(ids2), continue; end

        % 累计里程
        s = [0; cumsum(d2)];   % km

        for a = 1:numel(ids2)
            L = s(a);
            R = L + win_km;                       
            b = find(s <= R + 1e-9, 1, 'last');   
            bEdge = min(b-1, numel(ids2));
            if bEdge < a, continue; end
            idxC = [];
            for kk = a:bEdge
                rid_k = ids2(kk);
                ir_k  = find(r_ids == rid_k, 1);
                if ~isempty(ir_k) && ~isempty(Omega{ir_k})
                    idxC = union(idxC, Omega{ir_k});
                end
            end
            if isempty(idxC)
                nEmptyWin = nEmptyWin + 1;
                continue;
            end

            % 仅用大站优先满足覆盖
            idxC_good = idxC(nmax(idxC) >= y_min_default);

            if ~isempty(idxC_good)
                Cons = [Cons, sum(x(idxC_good)) >= 1];
            else
                % 没有大站，加松弛
                z = sdpvar(1,1);  Cons = [Cons, z >= 0, sum(x(idxC)) + z >= 1];
                Obj = Obj + gamma_cov * z;
            end
        end
    end
end



%% ==== 目标 ====
fi  = double(C.fi); 
ci  = double(C.ci); 
ct  = double(C.c_trench);
dBus= double(C.dist_bus_m);

Obj = Obj + fi.'*x + (ct.*dBus).'*x + ci.'*y;
Obj = Obj + gamma_slack * sum(slacks);

% 距离软惩罚
if OBJ_Soft > 0
    for ir = 1:NR
        lonr = double(Rmid.lon(ir)); latr = double(Rmid.lat(ir));
        if isnan(lonr) || isnan(latr), continue; end
        dkm = hav(latr,lonr,clat,clon)/1000;
        for hh = 1:NH
            if haveQ(ir,hh)
                Obj = Obj + OBJ_Soft * (dkm(Qidx{ir,hh}).' * Q{ir,hh});
            end
        end
    end
end

%% ==== 可行性 ====
cap_per_h  = sum((C.Prated_kW./E_session_kWh) .* C.nmax);  % 全网小时容量

dem_per_h = sum(lamRH,1);
fprintf('[CHECK] peak Dem=%.1f, max Cap=%.1f (veh/h)\n', max(dem_per_h), cap_per_h);
badH = find(dem_per_h > cap_per_h + 1e-9);
if ~isempty(badH)
    fprintf('[CHECK] infeasible hours: '); fprintf('%d ', Hset(badH)); fprintf('\n');
end
tot_day = sum(sum(lamRH,2));
fprintf('[CHK] System total (veh/hour aggregated)=%.0f → ≈ %.0f veh/day\n', sum(dem_per_h), sum(dem_per_h));

dem_per_h = sum(lamRH,1);              % 1×24
tot_day_ev = sum(dem_per_h);           % ≈ 18640.9
peak_dem   = max(dem_per_h);           % ≈ 1650
cap_per_h  = sum((C.Prated_kW./E_session_kWh) .* C.nmax);  % 全网小时容量
fprintf('[CHK] EV tot/day=%.1f | peakDem=%.1f | cap=%.1f (veh/h)\n', ...
        tot_day_ev, peak_dem, cap_per_h);


%% ==== 求解（外层迭代总网损最小OPF）====
ops = sdpsettings('solver','cplex','verbose',1);
ops.cplex.mip.tolerances.mipgap = 1e-3;
ops.cplex.timelimit             = 1800;
ops.cplex.display               = 2;

% -------- 网损--------
LossProxyTotal = 0;
if USE_LOSS_PROXY
    % Pi_MW(i,h) = E_session_kWh * a_{i,h}  (kW)  =>  /1000 => MW
    Pi_MW = (E_session_kWh/1000) * Aih;  % NI x NH 
    nb = size(mpc0.bus,1);

    for hh = 1:NH
        % 有功注入变化
        dp_full = - Mbus * Pi_MW(:,hh);

        % slack 
        dp_full(slack_bus_i) = -sum(dp_full);

        % 线路潮流
        df = PTDF * dp_full;

        % rho_loss
        LossProxyTotal = LossProxyTotal + sum( k_loss .* (df.^2) );
    end
end

Obj_base = Obj;     % 目标（建站/挖沟/等待/覆盖等）
diag = [];
best = struct('iter',0,'diag',[],'Obj',inf,'LossTrue',inf,'LossProxy',inf,'xv',[],'yv',[]);

loss_hist = table('Size',[0 6], ...
    'VariableTypes',{'double','double','double','double','double','double'}, ...
    'VariableNames',{'iter','rho_loss','ObjVal','LossProxyVal','LossTrueMW','isSuccess'});

for iter = 1:maxIter_loss
    Obj_iter = Obj_base + rho_loss * LossProxyTotal;

    fprintf('\n[ITER %d] rho_loss=%.4g | Solving MIQP...\n', iter, rho_loss);
    diag = optimize(Cons, Obj_iter, ops);

    if diag.problem ~= 0
        warning('[ITER %d] optimize fail：%s', iter, diag.info);
        loss_hist = [loss_hist; {iter, rho_loss, NaN, NaN, NaN, 0}];
        break
    end

    ObjVal = value(Obj_iter);
    LossProxyVal = value(LossProxyTotal);

    % ---- OPF 计算总网损----
    [LossTrueMW, isSuccess] = eval_total_loss_opf_from_solution(mpc0, Q, Qidx, cand2busrow, ...
        E_session_kWh, NH, loss_pf_ev);
    fprintf('[ITER %d] Obj=%.3f | LossProxy=%.3f | LossTrue=%.3f MW\n', ...
        iter, ObjVal, LossProxyVal, LossTrueMW);

    loss_hist = [loss_hist; {iter, rho_loss, ObjVal, LossProxyVal, LossTrueMW, double(isSuccess)}];

    % 保存
    if isSuccess && (LossTrueMW < best.LossTrue - 1e-6 || (abs(LossTrueMW-best.LossTrue)<1e-6 && ObjVal < best.Obj))
        best.iter = iter;
        best.diag = diag;
        best.Obj  = ObjVal;
        best.LossTrue  = LossTrueMW;
        best.LossProxy = LossProxyVal;
        best.xv = round(value(x));
        best.yv = round(value(y));
    end

    % ---- 终止判据 ----
    if iter >= 2
        prev = loss_hist.LossTrueMW(end-1);
        if isfinite(prev) && prev > 0
            rel_impr = abs(prev - LossTrueMW) / prev;
            if rel_impr < tol_loss_rel
                fprintf('[STOP]', rel_impr, tol_loss_rel);
                break
            end
        end
    end

    % ---- 权重校准----
    if USE_LOSS_PROXY && isfinite(LossProxyVal) && LossProxyVal > 1e-9 && isfinite(LossTrueMW) && LossTrueMW > 0
        ratio = LossTrueMW / LossProxyVal;
        ratio = min(max(ratio, 0.1), 10);   % 防止剧烈震荡
        rho_loss = rho_loss * ratio;
    else
        rho_loss = rho_loss * 2;
    end
end

if best.iter > 0
    fprintf('\n[FINAL] Use best iter=%d | LossTrue=%.3f MW | Obj=%.3f\n', best.iter, best.LossTrue, best.Obj);
    xv = best.xv;  yv = best.yv;
end

% 保存
try
    writetable(loss_hist, p('loss_iter_history.csv'));
catch
end


%% ==== 导出 ====
if diag.problem==0 || diag.problem==3 || diag.problem==4
    SolSite = table(C.cand_id, C.bus_id, C.type, C.lon, C.lat, xv, yv, ...
        'VariableNames', {'cand_id','bus_id','type','lon','lat','x','y'});
    writetable(SolSite, p('solution_sites_minimal_batt_new_final.csv'));
    disp('solution_sites_minimal_batt.csv');

    rows_r = []; rows_i = []; rows_h = []; rows_q = [];
    for ir = 1:NR
        rid = r_ids(ir);
        for hh = 1:NH
            if ~haveQ(ir,hh), continue; end
            qv = value(Q{ir,hh});
            keep = qv > 1e-6;
            if any(keep)
                rows_r = [rows_r; repmat(rid, nnz(keep),1)];
                rows_i = [rows_i; C.cand_id(Qidx{ir,hh}(keep))];
                rows_h = [rows_h; repmat(Hset(hh), nnz(keep),1)];
                rows_q = [rows_q; qv(keep)];
            end
        end
    end
    SolAssign = table(rows_r, rows_i, rows_h, rows_q, ...
        'VariableNames', {'r_id','cand_id','h','veh_per_h'});
    writetable(SolAssign, p('solution_assignment_batt_new_final.csv'));
    disp('solution_assignment_batt.csv');

    if ~isempty(SolAssign)
        Agg = groupsummary(SolAssign, {'cand_id','h'}, 'sum', 'veh_per_h');
        Agg.P_kW = Agg.sum_veh_per_h * E_session_kWh;
        Agg.sum_veh_per_h = [];
        writetable(Agg, p('station_hourly_load_kW_new_final.csv'));
        disp('station_hourly_load_kW.csv');
     end
else
    warning('无可行解');
end

try
    clear mpc;  
    run(p('build_mpc_london132.m'));
    mpc0 = mpc;   

 
    busN = size(mpc0.bus,1);
    BUS_I = 1; PD = 3; QD = 4;

    busNum = mpc0.bus(:, BUS_I);
    busNum2row = containers.Map('KeyType','double','ValueType','double');
    for rr = 1:busN
        busNum2row(busNum(rr)) = rr;
    end

    cand2row = containers.Map('KeyType','char','ValueType','double');
    for i = 1:height(C)
        cid = char(string(C.cand_id(i)));
        bi  = str2double(string(C.bus_id(i)));

        if isfinite(bi)
            if isKey(busNum2row, bi)
                cand2row(cid) = busNum2row(bi);          % 按母线编号映射
            elseif bi >= 1 && bi <= busN
                cand2row(cid) = bi;                      % 按行号映射
            else
            end
        end
    end

      Agg.bus_row = zeros(height(Agg),1);
    for r = 1:height(Agg)
        cid = char(string(Agg.cand_id(r)));
        if isKey(cand2row, cid)
            Agg.bus_row(r) = cand2row(cid);
        else
            Agg.bus_row(r) = 0; 
        end
    end
    Agg = Agg(Agg.bus_row>0, :);
    Agg.P_MW = Agg.P_kW / 1000;   % kW -> MW

    % 按 (bus_row, h) 汇总
    Agg2 = groupsummary(Agg, {'bus_row','h'}, 'sum', 'P_MW');
    Agg2.Pev_MW = Agg2.sum_P_MW;
    Agg2.sum_P_MW = [];

    %% ---- OPF  ----
    mpopt = mpoption('verbose', 0, 'out.all', 0);
    pf_ev = 1.00;
    if pf_ev < 1
        tanphi = tan(acos(pf_ev));
    else
        tanphi = 0;
    end

    rowsH = unique(Agg2.h);
    summary = table('Size',[NH 7], 'VariableTypes', ...
        {'double','double','double','double','double','double','double'}, ...
        'VariableNames', {'h','Pev_MW','Loss_MW','nV_viol','Vmin','Vmax','success'});

    bus_out = [];
    br_out  = [];

    for hh = 1:NH
        mpc_t = mpc0;

        % 该小时的 EV 负荷叠加到 bus PD/QD
        sel = (Agg2.h == hh);
        pev_bus = zeros(busN,1);
        pev_bus(Agg2.bus_row(sel)) = Agg2.Pev_MW(sel);

        mpc_t.bus(:,PD) = mpc_t.bus(:,PD) + pev_bus;
        mpc_t.bus(:,QD) = mpc_t.bus(:,QD) + pev_bus * tanphi;

        % OPF
        res = runopf(mpc_t, mpopt);

        summary.h(hh)      = hh;
        summary.Pev_MW(hh) = sum(pev_bus);

        if isfield(res,'success') && res.success
            % 网损（MW）
            Loss_MW = sum(res.gen(:,2)) - sum(res.bus(:,PD));

            % 电压越限统计
            VM = res.bus(:,8);
            Vmin = min(VM); Vmax = max(VM);
            nV  = nnz(VM < 0.95 | VM > 1.05);

            summary.Loss_MW(hh)  = Loss_MW;
            summary.nV_viol(hh)  = nV;
            summary.Vmin(hh)     = Vmin;
            summary.Vmax(hh)     = Vmax;
            summary.success(hh)  = 1;

            % 记录 bus 结果
            Tbus = table();
            Tbus.h      = repmat(hh, busN,1);
            Tbus.bus_i  = res.bus(:,BUS_I);
            Tbus.PD_MW  = res.bus(:,PD);
            Tbus.QD_MVAr= res.bus(:,QD);
            Tbus.VM_pu  = res.bus(:,8);
            Tbus.VA_deg = res.bus(:,9);
            bus_out = [bus_out; Tbus]; 

            % 记录 branch 结果
            Pf = res.branch(:,14);
            Pt = res.branch(:,16);
            loss_br = Pf + Pt;

            Tbr = table();
            Tbr.h = repmat(hh, size(res.branch,1), 1);
            Tbr.fbus = res.branch(:,1);
            Tbr.tbus = res.branch(:,2);
            Tbr.Pf_MW = Pf;
            Tbr.Pt_MW = Pt;
            Tbr.Loss_MW = loss_br;
            br_out = [br_out; Tbr];

        else
            summary.Loss_MW(hh)  = NaN;
            summary.nV_viol(hh)  = NaN;
            summary.Vmin(hh)     = NaN;
            summary.Vmax(hh)     = NaN;
            summary.success(hh)  = 0;
        end
    end

    %% ---- 结果 ----
    writetable(summary, p('opf_hourly_summary.csv'));
    if ~isempty(bus_out), writetable(bus_out, p('opf_bus_hourly.csv')); end
    if ~isempty(br_out),  writetable(br_out,  p('opf_branch_hourly.csv')); end


end
