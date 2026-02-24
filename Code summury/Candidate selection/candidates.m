%% ===================== candidates.m =====================
% 候选点为高速路网节点；按流量筛 + 接入距离 + 最小间距 贪心瘦身
% 输出：candidates.csv  

clc; clear; close all;
baseDir = 'C:\2023Research\Emobility\transportationmap\optimization2';
p = @(f) fullfile(baseDir, f);

%% --- 读数据 ---
Nodes = readtable(p('03_network-nodes.csv'),       'VariableNamingRule','preserve'); % Network_Node_ID, Network_Node_X, Network_Node_Y
Edges = readtable(p('04_network-edges.csv'),       'VariableNamingRule','preserve'); % Network_Edge_ID, Network_Node_A_ID, Network_Node_B_ID, Distance, Traffic_flow_trucks_2030
Subs  = readtable(p('london_132kv_substations_bigger2_main.csv'),'VariableNamingRule','preserve'); % id, lon, lat
BusE  = readtable(p('bus_energy_annual.csv'),      'VariableNamingRule','preserve'); % id|bus_id, MWh_day_2030
MapNd = readtable(p('mapping_nodes_to_bus.csv'),   'VariableNamingRule','preserve'); % id, lon, lat, bus_id

%% --- 参数 ---
FI = 300e3; CI = 20e3; C_TRENCH = 500; PRATED = 350;    % 成本/额定功率
DMAX_M = 12000;     
MIN_SPACING_M = 3000;  
TRUCKS_QMIN = 1e5;             
TOP_PER_BUS = 8;             

%% --- 经纬度映射---
Nodes.NodeID = double(Nodes.Network_Node_ID);
keyLon = containers.Map(Nodes.NodeID, double(Nodes.Network_Node_X));
keyLat = containers.Map(Nodes.NodeID, double(Nodes.Network_Node_Y));

%% --- 1) 选高速的边 ---
% 
cands = {'Traffic_flow_trucks_2030','Traffic_flow_trucks_2030_','Traffic_flow_trucks_2030 '};
flowCol = '';
for k = 1:numel(cands)
    if ismember(cands{k}, Edges.Properties.VariableNames), flowCol = cands{k}; break; end
end

if ~isempty(flowCol)
    Ekeep = Edges(Edges.(flowCol) >= TRUCKS_QMIN, :);
end

nid = unique([double(Ekeep.Network_Node_A_ID); double(Ekeep.Network_Node_B_ID)]);


hasXY = arrayfun(@(k) isKey(keyLon,k) && isKey(keyLat,k), nid);
nid   = nid(hasXY);
lon   = arrayfun(@(k) keyLon(k), nid);
lat   = arrayfun(@(k) keyLat(k), nid);

Cand0 = table(string(nid), lon(:), lat(:), 'VariableNames',{'cand_id','lon','lat'});
% —— 按伦敦范围过滤候选——
% 伦敦BBOX
LON = [-0.58, 0.42];   % west, east
LAT = [51.28, 51.75];  % south, north
inLDN = Cand0.lon >= LON(1) & Cand0.lon <= LON(2) & ...
        Cand0.lat >= LAT(1) & Cand0.lat <= LAT(2);
Cand0 = Cand0(inLDN, :);


%% --- 2) mapping_nodes_to_bus ---
MapNd.id     = string(MapNd.id);
MapNd.bus_id = string(MapNd.bus_id);

Cand0 = outerjoin(Cand0, MapNd(:,{'id','bus_id'}), ...
                  'LeftKeys','cand_id', 'RightKeys','id', ...
                  'MergeKeys',true, 'Type','left');
if ismember('bus_id_MapNd', Cand0.Properties.VariableNames)
    Cand0.bus_id = Cand0.bus_id_MapNd; Cand0.bus_id_MapNd = [];
end


Subs.bus_id = string(Subs.id);
hav = @(la1,lo1,la2,lo2) 2*6371000*asin( sqrt( sin(deg2rad((la2-la1)/2)).^2 + ...
                          cosd(la1).*cosd(la2).*sin(deg2rad((lo2-lo1)/2)).^2 ) );
need_bus = (Cand0.bus_id=="" | ismissing(Cand0.bus_id));
if any(need_bus)
    for t = find(need_bus).'
        [~, idx] = min(hav(Cand0.lat(t), Cand0.lon(t), Subs.lat, Subs.lon));
        Cand0.bus_id(t) = Subs.bus_id(idx);
    end
end

% 合并母线坐标算接入距离
SubsSlim = Subs(:,{'bus_id','lon','lat'}); SubsSlim.Properties.VariableNames = {'bus_id','bus_lon','bus_lat'};
Cand0 = outerjoin(Cand0, SubsSlim, 'Keys','bus_id', 'MergeKeys',true, 'Type','left');


Cand0.dist_bus_m = hav(Cand0.lat, Cand0.lon, Cand0.bus_lat, Cand0.bus_lon);
fprintf('[LOG]  dist_bus_m: min=%.0f p50=%.0f p90=%.0f max=%.0f (m)\n', ...
    min(Cand0.dist_bus_m,[],'omitnan'), prctile(Cand0.dist_bus_m,50), ...
    prctile(Cand0.dist_bus_m,90), max(Cand0.dist_bus_m,[],'omitnan'));

% 接入距离门槛
Cand1 = Cand0(Cand0.dist_bus_m <= DMAX_M & ~isnan(Cand0.dist_bus_m), :);

if ~ismember('cand_id', Cand1.Properties.VariableNames)
    alt = Cand1.Properties.VariableNames(startsWith(Cand1.Properties.VariableNames, 'cand_id'));
    if ~isempty(alt)
        Cand1.cand_id = Cand1.(alt{1});
    end
end
Cand1.cand_id = string(Cand1.cand_id);

%% --- 3) 重要度排序 + 最小间距贪心 + 总量上限 ---
nid2w = containers.Map('KeyType','char','ValueType','double');
if ~isempty(flowCol)
    for r = 1:height(Ekeep)
        na = string(Ekeep.Network_Node_A_ID(r));
        nb = string(Ekeep.Network_Node_B_ID(r));
        w  = double(Ekeep.(flowCol)(r));
        if ~isKey(nid2w,na), nid2w(na)=0; end
        if ~isKey(nid2w,nb), nid2w(nb)=0; end
        nid2w(na) = nid2w(na) + w;
        nid2w(nb) = nid2w(nb) + w;
    end
else
    for k = 1:height(Cand1), nid2w(Cand1.cand_id(k)) = 1; end
end


if ~ismember('cand_id', Cand1.Properties.VariableNames)
    alt = Cand1.Properties.VariableNames(startsWith(Cand1.Properties.VariableNames,'cand_id'));
    if ~isempty(alt), Cand1.cand_id = Cand1.(alt{1}); end
end
Cand1.cand_id = string(Cand1.cand_id);


nid2w = containers.Map('KeyType','char','ValueType','double');
if ~isempty(flowCol)
    for r = 1:height(Ekeep)
        na = string(Ekeep.Network_Node_A_ID(r));
        nb = string(Ekeep.Network_Node_B_ID(r));
        w  = double(Ekeep.(flowCol)(r));
        if ~isKey(nid2w,na), nid2w(na)=0; end
        if ~isKey(nid2w,nb), nid2w(nb)=0; end
        nid2w(na) = nid2w(na) + w;
        nid2w(nb) = nid2w(nb) + w;
    end
else
    for k = 1:height(Cand1), nid2w(Cand1.cand_id(k)) = 1; end
end

% 按母线分组
[grp, busKeys] = findgroups(Cand1.bus_id);
Cand2 = table();

for g = 1:numel(busKeys)
    T = Cand1(grp == g, :);      % 该母线下所有候选

    w = zeros(height(T),1);
    for t = 1:height(T)
        k = T.cand_id(t);
        if isKey(nid2w, k), w(t) = nid2w(k); else, w(t) = 0; end
    end
    [~,ord] = sort(w,'descend'); T = T(ord,:);

    % 最小间距贪心
    keep = false(height(T),1);
    pickedLat = []; pickedLon = [];
    for t = 1:height(T)
        if isempty(pickedLat)
            dmin = inf;
        else
            H = 2*6371000*asin( sqrt( sin(deg2rad((pickedLat - T.lat(t))/2)).^2 + ...
                cosd(T.lat(t)).*cosd(pickedLat).*sin(deg2rad((pickedLon - T.lon(t))/2)).^2 ) );
            dmin = min(H);
        end
        if dmin >= MIN_SPACING_M
            keep(t) = true;
            pickedLat(end+1) = T.lat(t);  
            pickedLon(end+1) = T.lon(t);  
        end
        if sum(keep) >= TOP_PER_BUS, break; end
    end

    Cand2 = [Cand2; T(keep,:)];  %#ok<AGROW>
end


% 去重
[~, ia] = unique(Cand2.cand_id);
dup = setdiff(1:height(Cand2), ia);


leftBus = setdiff(unique(Cand1.bus_id), unique(Cand2.bus_id));
if ~isempty(leftBus)
    addRows = table();
    for b = leftBus.'
        T = Cand1(Cand1.bus_id==b, :);
        [~,i0] = min(T.dist_bus_m);
        addRows = [addRows; T(i0,:)];
    end
    Cand2 = [Cand2; addRows];

end

%% --- 4) 成本/功率/能量---
n = height(Cand2);
Cand2.type       = repmat("highway_node", n, 1);
Cand2.fi         = repmat(FI,       n, 1);
Cand2.ci         = repmat(CI,       n, 1);
Cand2.c_trench   = repmat(C_TRENCH, n, 1);
Cand2.Prated_kW  = repmat(PRATED,   n, 1);
Cand2.trucks_2030    = zeros(n,1);
Cand2.truck_km_2030  = zeros(n,1);
Cand2.MWh_2030       = zeros(n,1);


BusE.Properties.VariableNames = strrep(BusE.Properties.VariableNames, 'id', 'bus_id');
BusE.bus_id = string(BusE.bus_id);
Cand2.bus_id = string(Cand2.bus_id);

Cand2 = outerjoin(Cand2, BusE(:,{'bus_id','MWh_day_2030'}), ...
                  'Keys','bus_id','MergeKeys',true,'Type','left');

Cand2.MWh_day_2030(isnan(Cand2.MWh_day_2030)) = 0;
H_eff = 8; U = 0.6;
Ppeak_kW = 1e3 * (Cand2.MWh_day_2030./H_eff)./U;
[g,~] = findgroups(Cand2.bus_id);
grp_counts = accumarray(g,1);
share = 1 ./ grp_counts(g);
Cand2.nmax = max(2, ceil((Ppeak_kW .* share) ./ Cand2.Prated_kW));





%% --- 5) 导出 ---
bad = isnan(Cand2.lon) | isnan(Cand2.lat) | Cand2.bus_id=="" | ismissing(Cand2.bus_id);
out = Cand2(:,{'cand_id','type','lon','lat','bus_id','dist_bus_m','fi','ci','c_trench','Prated_kW','nmax', ...
               'trucks_2030','truck_km_2030','MWh_2030'});
writetable(out, p('candidates.csv'));
fprintf('[LOG] candidates.csv 已生成：%s | 行数=%d\n', p('candidates.csv'), height(out));
