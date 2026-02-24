% ========= paths_or_corridors.csv 生成=========
clc;
baseDir = 'C:\2023Research\Emobility\transportationmap\optimization3';
p = @(f) fullfile(baseDir, f);

% 读表
C   = readtable(p('candidates_new.csv'),             'VariableNamingRule','preserve','TextType','string');   % cand_id,bus_id,lon,lat,dist_bus_m
Map = readtable(p('mapping_edge_midpoints.csv'), 'VariableNamingRule','preserve','TextType','string');   % edge_id,bus_id,(可无 lon/lat)

% 统一类型
if ~ismember('bus_id', Map.Properties.VariableNames) && ismember('id', Map.Properties.VariableNames)
    Map.bus_id = Map.id;
end
Map.bus_id = string(Map.bus_id);
C.bus_id   = string(C.bus_id);
if ~isa(C.cand_id,'string'), C.cand_id = string(C.cand_id); end

haveLL = ismember('lon', Map.Properties.VariableNames) && ismember('lat', Map.Properties.VariableNames);
if ~haveLL || all(ismissing(Map.lon)) || all(ismissing(Map.lat))
    Nodes = readtable(p('03_network-nodes.csv'), 'VariableNamingRule','preserve');
    Edges = readtable(p('04_network-edges.csv'), 'VariableNamingRule','preserve');

    keyLon = containers.Map(double(Nodes.Network_Node_ID), double(Nodes.Network_Node_X));
    keyLat = containers.Map(double(Nodes.Network_Node_ID), double(Nodes.Network_Node_Y));
    u = double(Edges.Network_Node_A_ID);
    v = double(Edges.Network_Node_B_ID);
    ok = arrayfun(@(a,b) isKey(keyLon,a)&&isKey(keyLon,b)&&isKey(keyLat,a)&&isKey(keyLat,b), u, v);

    mid_lon = nan(height(Edges),1); mid_lat = mid_lon;
    mid_lon(ok) = 0.5*(arrayfun(@(a) keyLon(a), u(ok)) + arrayfun(@(b) keyLon(b), v(ok)));
    mid_lat(ok) = 0.5*(arrayfun(@(a) keyLat(a), u(ok)) + arrayfun(@(b) keyLat(b), v(ok)));

    Tmid = table(string(Edges.Network_Edge_ID), mid_lon, mid_lat, ...
                 'VariableNames', {'edge_id','lon','lat'});

    Map = outerjoin(Map, Tmid, 'Keys','edge_id', 'MergeKeys',true, 'Type','left');
end

[~, ia] = unique(Map.edge_id, 'stable');
Map = Map(ia, :);

% ---------- 参数 ----------
R_WIN_M = 5000;   % 滑窗半径
K       = 8;      % 每条 r 最多候选数

% Haversine
hav = @(la1,lo1,la2,lo2) 2*6371000*asin( sqrt( sin(deg2rad((la2-la1)/2)).^2 + ...
                          cosd(la1).*cosd(la2).*sin(deg2rad((lo2-lo1)/2)).^2 ) );

% ---------- 为每个母线准备candidate set ----------
[gb, busKeys] = findgroups(C.bus_id);
poolByBus = containers.Map('KeyType','char','ValueType','any');
for k = 1:numel(busKeys)
    Tb = C(gb==k, {'cand_id','bus_id','lon','lat','dist_bus_m'});
    bad = (Tb.cand_id=="" | Tb.cand_id=="0" | Tb.cand_id=="NaN" | ismissing(Tb.cand_id));
    Tb(bad,:) = [];
    [~,ord] = sort(Tb.dist_bus_m, 'ascend', 'MissingPlacement','last');
    Tb = Tb(ord,:);
    [~,iu] = unique(Tb.cand_id, 'stable');
    Tb = Tb(iu,:);
    poolByBus(char(busKeys(k))) = Tb;
end

% ---------- 逐 r 生成候选（滑窗 + 距离最近回退） ----------
rows = cell(height(Map), 2);
nFallback = 0; nEmpty = 0;

for i = 1:height(Map)
    rows{i,1} = Map.edge_id(i);

    b = char(Map.bus_id(i));
    if ~isKey(poolByBus, b)
        K = 5; R_FB_M = 15000;  
        d = hav(Map.lat(i), Map.lon(i), C.lat, C.lon);
        [ds,ord] = sort(d,'ascend');
        pick = ord(ds <= R_FB_M);
        if isempty(pick), pick = ord(1:min(K, numel(ord))); else, pick = pick(1:min(K,numel(pick))); end
        rows{i,2} = strjoin(string(C.cand_id(pick)).', ',');
        nEmpty = nEmpty + 1;   
        continue
    end
    Tb = poolByBus(b);

    if any(ismissing([Map.lat(i), Map.lon(i)]))
        pick = 1:min(K, height(Tb));
        rows{i,2} = strjoin(string(Tb.cand_id(pick)).', ',');
        nFallback = nFallback + 1;
        continue
    end

    d = hav(Map.lat(i), Map.lon(i), Tb.lat, Tb.lon);   
    inwin = d <= R_WIN_M;

    if any(inwin)
        [~, o2] = sort(d(inwin), 'ascend');
        idx = find(inwin);
        pick = idx(o2(1:min(K, numel(o2))));
    else
        [~, oAll] = sort(d, 'ascend');
        pick = oAll(1:min(K, numel(oAll)));
        nFallback = nFallback + 1;
    end

    rows{i,2} = strjoin(string(Tb.cand_id(pick)).', ',');
end

P = cell2table(rows, 'VariableNames', {'r_id','cand_ids'});

% 1) 统计空行
emptyMask = (P.cand_ids=="" | ismissing(P.cand_ids));
fprintf('[FIX] empty rows before fallback = %d\n', nnz(emptyMask));

if nnz(emptyMask) > 0
    % 全局候选
    Cgood = C(~ismissing(C.cand_id) & C.cand_id~="" & C.cand_id~="0", ...
              {'cand_id','bus_id','lon','lat','dist_bus_m'});
    [~,iu] = unique(Cgood.cand_id,'stable'); Cgood = Cgood(iu,:);

    % r 的 lon/lat
    R = Map(:,{'edge_id','bus_id','lon','lat'}); R.Properties.VariableNames = {'r_id','bus_id','lon','lat'};
    P = outerjoin(P,R,'Keys','r_id','MergeKeys',true,'Type','left');  

    for i = find(emptyMask).'
        if ~(~ismissing(P.lon(i)) && ~ismissing(P.lat(i)))
            Tb = Cgood(string(Cgood.bus_id)==string(P.bus_id(i)),:);
            if ~isempty(Tb)
                P.lon(i) = mean(Tb.lon,'omitnan'); 
                P.lat(i) = mean(Tb.lat,'omitnan');
            end
        end
        if ~(~ismissing(P.lon(i)) && ~ismissing(P.lat(i)))
            continue
        end

        % 计算与全局候选的距离，取最近 K 个
        d = hav(P.lat(i), P.lon(i), Cgood.lat, Cgood.lon);
        [~,o] = sort(d, 'ascend');
        pick = o(1:min(K, numel(o)));
        P.cand_ids(i) = strjoin(string(Cgood.cand_id(pick)).', ',');
    end

    P.lon = []; P.lat = []; P.bus_id = [];
end

% 二次统计
emptyMask = (P.cand_ids=="" | ismissing(P.cand_ids));
fprintf('[FIX] empty rows after  fallback = %d\n', nnz(emptyMask));

writetable(P, p('paths_or_corridors_new_check.csv'), 'QuoteStrings', true);
fprintf('[CHK] rows=%d | nFallback=%d | empty=%d | R_WIN=%.0f m | K=%d\n', ...
        height(P), nFallback, nEmpty, R_WIN_M, K);

