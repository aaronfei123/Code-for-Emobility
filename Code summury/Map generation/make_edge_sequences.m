function make_edge_sequences()
%% ===================== make_edge_sequences.m =====================
% 从 01_Trucktrafficflow.csv 生成 edge_sequences.csv
% 列：route_id, edge_ids（逗号分隔, 按顺序）, len_km（序列总里程）
% 仅供“续航滑窗覆盖约束”使用；不参与任何需求加总

clc;
baseDir = 'C:\2023Research\Emobility\transportationmap\optimization2';
p = @(f) fullfile(baseDir, f);

%% ---- 加速参数 ----
CHUNK = 5e4;                  % 每块读取
PROGRESS_EVERY = 5e4;
USE_WAITBAR = true;

%% ---- 白名单边（限定建模范围） ----
Map = readtable(p('mapping_edge_midpoints.csv'), ...
    'VariableNamingRule','preserve','TextType','string');
v = string(Map.Properties.VariableNames); v = regexprep(v,'^\xEF\xBB\xBF','');
Map.Properties.VariableNames = cellstr(v);

eidCol = "";
for c = ["edge_id","id","r_id"]
    if ismember(c,Map.Properties.VariableNames), eidCol=c; break; end
end
assert(eidCol~="", 'mapping_edge_midpoints 缺 edge_id/id/r_id');

allowEdges = unique(double(Map.(eidCol)));
allowEdges = allowEdges(~isnan(allowEdges));
allowEdges = sort(allowEdges(:));     % 向量+排序，配合 ismember 二分

% anchors（粗筛，加速用）
RAND_ANCHORS = min(300, numel(allowEdges));
anchors = string(allowEdges(randperm(numel(allowEdges), RAND_ANCHORS)));
anchors = "," + anchors + ",";

%% ---- 04：距离映射（用于 len_km） ----
E = readtable(p('04_network-edges.csv'), 'VariableNamingRule','preserve');
eid04  = double(E.Network_Edge_ID);
lenMap = containers.Map(num2cell(eid04), num2cell(double(E.Distance)/1000)); % km

%% ---- 读 OD（分块） ----
OD_FILE = p('01_Trucktrafficflow.csv');
opts = detectImportOptions(OD_FILE,'TextType','string','PreserveVariableNames',true);
names = string(opts.VariableNames);
canon = @(s) lower(regexprep(string(s),'[\s_]+',''));
pick  = @(cands) names( find(ismember(canon(names), canon(string(cands))), 1,'first') );
pathCol = pick({'Edge_path_E_road','edge_path_e_road','edge_path','path_e_road','path'});
assert(~isempty(pathCol), 'OD 文件缺少路径列。');

ds = tabularTextDatastore(OD_FILE, 'TextType','string');
ds.SelectedVariableNames = pathCol;
ds.ReadSize = CHUNK;

seqs = strings(0,1);     % 原始序列（去重前）
block=0; totalRows=0; kept=0; skip=0;
if USE_WAITBAR, hwb = waitbar(0,'Parsing paths ...'); end
tblock = tic;

while hasdata(ds)
    T = read(ds); block=block+1;
    nrows = height(T); totalRows = totalRows + nrows;

    S = T.(pathCol);

    % anchors 粗筛
    mask_keep = false(nrows,1);
    for a = 1:numel(anchors)
        mask_keep = mask_keep | contains(S, anchors(a), 'IgnoreCase', false);
    end
    idxKeep = find(mask_keep);
    if isempty(idxKeep), skip = skip + nrows; 
        % 进度日志
        if mod(totalRows, PROGRESS_EVERY)==0
            fprintf('[PATH] rows=%d | kept=%d | skip=%d | t=%.1fs\n', totalRows, kept, skip, toc(tblock));
            if USE_WAITBAR
                waitbar(min(0.99, totalRows/1e7), hwb, sprintf('Parsed rows: %d', totalRows));
            end
            tblock = tic;
        end
        continue; 
    end

    for ii = 1:numel(idxKeep)
        i = idxKeep(ii);
        s = S(i);
        if ismissing(s) || strlength(s)==0, skip=skip+1; continue; end
        c = char(s); c(c=='[' | c==']' | c==' ') = [];
        ids = sscanf(c, '%f,');
        if isempty(ids), skip=skip+1; continue; end

        % 精确过滤 + 相邻去重
        mask = ismember(ids, allowEdges);
        if ~any(mask), skip=skip+1; continue; end
        ids = ids(mask);
        ids = ids([true; diff(ids)~=0]);

        seqs(end+1,1) = strjoin(string(ids), ','); 
        kept = kept + 1;
    end

    % 进度日志
    if mod(totalRows, PROGRESS_EVERY)==0
        fprintf('[PATH] rows=%d | kept=%d | skip=%d | t=%.1fs\n', totalRows, kept, skip, toc(tblock));
        if USE_WAITBAR
            waitbar(min(0.99, totalRows/1e7), hwb, sprintf('Parsed rows: %d', totalRows));
        end
        tblock = tic;
    end
end
if USE_WAITBAR, close(hwb); end

fprintf('[OD] 完成解析：kept=%d | skip=%d | totalRows=%d | raw unique=%d\n', ...
        kept, skip, totalRows, numel(unique(seqs)));

%% ---- 合并完全相同的序列（去重） ----
[uniqSeq, ~] = unique(seqs, 'stable');    % 不需要权重

%% ---- 计算 len_km ----
len_each = zeros(numel(uniqSeq),1);
for k = 1:numel(uniqSeq)
    ids = double(str2double(strsplit(uniqSeq(k), ',')));
    L = 0; 
    for e = ids(:).'
        if isKey(lenMap,e), L = L + lenMap(e); end
    end
    len_each(k) = L;
end

%% ---- 导出 edge_sequences.csv ----
route_id = (1:numel(uniqSeq)).';
SeqTbl = table(route_id, uniqSeq, len_each, 'VariableNames', {'route_id','edge_ids','len_km'});
writetable(SeqTbl, p('edge_sequences.csv'));

fprintf('[OUT] edge_sequences.csv 写出 rows=%d | 平均长度=%.1f km | 最长=%.1f km\n', ...
        height(SeqTbl), mean(len_each,'omitnan'), max(len_each,[],'omitnan'));
end
