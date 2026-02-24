function build_hgv_stop_prob_csv()
% 生成 24×1 的 hgv_stop_prob_by_hour.csv 
baseDir = 'C:\2023Research\Emobility\transportationmap\optimization3';
p = @(f) fullfile(baseDir,f);


USE_EVENT_LOG = false;             
EVENT_FILE    = p('hgv_stop_events.csv');  

USE_FLOW_SHARE = true;            
FLOW_SHARE_FILE = p('hgv_hourly_flow_share.csv'); 

OUT_FILE = p('hgv_stop_prob_by_hour.csv');

if USE_EVENT_LOG && isfile(EVENT_FILE)
    T = readtable(EVENT_FILE);
    assert(ismember('ts', T.Properties.VariableNames), 'include ts');
    ts = T.ts;
    if ~isdatetime(ts), ts = datetime(ts,'InputFormat','yyyy-MM-dd HH:mm:ss','TimeZone','local'); end
    h = hour(ts);
    cnt = accumarray(h+1, 1, [24,1], @sum, 0);   
    prob = cnt / max(sum(cnt),1);
else
  
    S = readtable(FLOW_SHARE_FILE);
    assert(ismember('h',S.Properties.VariableNames) && ismember('share',S.Properties.VariableNames), ...
        'FLOW_SHARE_FILE 需包含列 h, share');
    share = zeros(24,1); share(S.h+1) = double(S.share);
    share = share / sum(share + eps);          

    % 驾时规则
    Hday_drive = 9.0;        % 日驾驶时长（常见 9h；有时 10h）
    block_drive= 4.5;        % 每段连续驾驶最长 4.5h
    rest_min   = 0.75;       % 45分钟休息=0.75h
    nBlocks    = ceil(Hday_drive / block_drive); % 一天通常 2 段 4.5h
    long_rest  = 10;         % 夜间长休（小时，影响跨午夜模拟）

    % 蒙特卡洛参数
    Nveh = 50000;            % 模拟车辆天数
    stop_hist = zeros(24,1);

    % 服务/装卸停靠概率
    use_service_stops = true;
    p_service_per_h   = 0.10;    % 每小时有 10% 概率因装卸/排队临时停靠一次
    service_dur_h     = 0.5;     % 平均 30 分钟

    for n = 1:Nveh
        h0 = randsample(0:23, 1, true, share);
        t  = h0 + rand();         
        % 驾驶-休息序列
        drive_left = Hday_drive;
        while drive_left > 0
            d = min(block_drive, drive_left);
            % 这一段驾驶过程中，按小时步进，期间可能发生服务停靠
            t_end = t + d;
            if use_service_stops
                % 每小时掷一次硬币
                hh = floor(t):floor(t_end);
                for hh1 = hh
                    if rand() < p_service_per_h
                        stop_hist(mod(hh1,24)+1) = stop_hist(mod(hh1,24)+1) + 1;
                        t = min(t + service_dur_h, t_end);  % 插入停靠时间
                    end
                end
            end
            t = t_end;
            drive_left = drive_left - d;
            if drive_left <= 1e-6, break; end
            % 安排45分钟休息
            stop_hist(mod(floor(t),24)+1) = stop_hist(mod(floor(t),24)+1) + 1;
            t = t + rest_min;
        end
        % 夜间长休（不计入可快充停靠，通常回家）
        t = t + long_rest;
    end
    prob = stop_hist / max(sum(stop_hist),1);
end

% 输出CSV
h = (0:23).';
Tout = table(h, prob, 'VariableNames', {'h','prob'});
writetable(Tout, OUT_FILE);
fprintf('已生成 %s\n', OUT_FILE);
end
