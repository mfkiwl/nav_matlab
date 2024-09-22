classdef navigation_plots
    methods(Static)
        function axis_comparison(varargin)
            % 绘制数据对比图
            % 用法:
            % navigation_plots.axis_comparison('XLabel', '时间(s)', ...
            %                                  'YLabel', {{'East (m)', 'North (m)', 'Up (m)'}, {'Std East (m)', 'Std North (m)', 'Std Up (m)'}}, ...
            %                                  'Title', '位置对比', ...
            %                                  'Time', {data.imu.tow, data.imu.tow, data.imu.tow, data.imu.tow}, ...
            %                                  'Values', {log.pos, data.dev.pos_enu, log.P(:,7:9), data.dev.kf_p_pos}, ...
            %                                  'Labels', {'MATLAB', 'DEV', 'MATLAB_STD', 'DEV_STD'}, ...
            %                                  'Subplots', [1, 1, 2, 2], ...
            %                                  'LineStyles', {'-', '--', ':', '-.'});

            p = inputParser;
            addParameter(p, 'XLabel', '时间(s)', @ischar);
            addParameter(p, 'YLabel', {{'Y1', 'Y2', 'Y3'}}, @iscell);
            addParameter(p, 'Title', '', @ischar);
            addParameter(p, 'Time', {}, @iscell);
            addParameter(p, 'Values', {}, @iscell);
            addParameter(p, 'Labels', {}, @iscell);
            addParameter(p, 'Subplots', [], @isnumeric);
            addParameter(p, 'FontSize', 12, @isnumeric);
            addParameter(p, 'LineStyles', {}, @iscell);
            addParameter(p, 'Marker', 'none', @ischar);

            parse(p, varargin{:});

            xlabel_text = p.Results.XLabel;
            ylabel_text = p.Results.YLabel;
            title_text = p.Results.Title;
            time_data = p.Results.Time;
            values_data = p.Results.Values;
            labels = p.Results.Labels;
            subplots = p.Results.Subplots;
            font_size = p.Results.FontSize;
            line_styles = p.Results.LineStyles;
            marker = p.Results.Marker;

            % 输入验证
            assert(~isempty(time_data) && ~isempty(values_data), 'Time and Values data must be provided');
            assert(all(cellfun(@(x) length(x) == length(time_data{1}), time_data)), 'All time data must have the same length');
            assert(all(cellfun(@(x) size(x,1) == length(time_data{1}), values_data)), 'All value data must have the same number of rows as time data');
            assert(length(subplots) == length(values_data), 'Subplots must have the same length as the number of data sets');

            if isempty(labels)
                labels = arrayfun(@(x) sprintf('Data %d', x), 1:length(values_data), 'UniformOutput', false);
            end

            % 如果未提供LineStyles，则为所有数据集使用默认线型 '-'
            if isempty(line_styles)
                line_styles = repmat({'-'}, 1, length(values_data));
            elseif length(line_styles) < length(values_data)
                % 如果提供的LineStyles不足，用默认线型 '-' 补齐
                line_styles = [line_styles, repmat({'-'}, 1, length(values_data) - length(line_styles))];
            end

            num_subplots = max(subplots);

            figure('Name', title_text);

            colors = get(gca, 'ColorOrder');
            num_datasets = length(values_data);

            % 确保颜色足够
            if size(colors, 1) < num_datasets
                colors = repmat(colors, ceil(num_datasets / size(colors, 1)), 1);
            end
            colors = colors(1:num_datasets, :);

            for i = 1:3
                for sp = 1:num_subplots
                    ax = subplot(num_subplots, 3, i + (sp-1)*3);
                    hold(ax, 'on');
                    set(ax, 'UserData', {});  % 初始化UserData以存储图例条目

                    % 使用 for 循环替代 cellfun
                    for j = 1:num_datasets
                        navigation_plots.plot_data(ax, time_data{j}, values_data{j}(:,i), ...
                            labels{j}, subplots(j), colors(j,:), sp, line_styles{j}, marker);
                    end

                    title(ax, ylabel_text{sp}{i});
                    xlabel(ax, xlabel_text);
                    ylabel(ax, ylabel_text{sp}{i});
                    legend_entries = get(ax, 'UserData');
                    legend(ax, legend_entries{:});
                    grid(ax, 'on');
                end
            end

            sgtitle(title_text);
            set(gcf, 'Units', 'normalized', 'Position', [0.025, 0.05, 0.95, 0.85]);

            % 设置字体大小
            set(findall(gcf,'-property','FontSize'),'FontSize', font_size)
        end

        % 辅助函数用于绘制单个数据集
        function plot_data(ax, time, values, label, subplot, color, target_subplot, line_style, marker)
            if subplot == target_subplot
                line(ax, time, values, 'Color', color, 'Marker', marker, 'LineStyle', line_style, 'linewidth', 1.5);
                legend_entries = get(ax, 'UserData');
                legend_entries{end+1} = label;
                set(ax, 'UserData', legend_entries);
            end
        end

function trajectory_2d_plot(enu_data, labels)
    % 绘制二维轨迹
    % 参数:
    % enu_data: 包含多个轨迹数据的cell数组
    % labels: 图例标签

    figure('Name', '2D轨迹');
    
    % 定义不同的颜色
    colors = {'b', 'r', 'g', 'm', 'c', 'y', 'k', [0.5 0.5 0.5]};
    
    for i = 1:length(enu_data)
        % 使用循环索引选择标记和颜色，如果超出范围就循环使用
        color = colors{mod(i-1, length(colors)) + 1};
        
        % 绘制线条
        plot(enu_data{i}(:,1), enu_data{i}(:,2), ...
            '.-', ...
             'Color', color, ...
             'LineWidth', 0.5);
        hold on;
    end
     
    legend(labels, 'Location', 'best');
    axis equal;
    xlabel('East (m)');
    ylabel('North (m)');
    title('2D Trajectory');
    grid on;
    
    set(gcf, 'Units', 'normalized', 'Position', [0.025, 0.05, 0.95, 0.85]);
end


    end
end
