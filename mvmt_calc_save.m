function fig = mvmt_calc_save(mvmts, saveto, varargin)

softwr = 'FSL';
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'software'}
                softwr = varargin{i+1};
        end
    end
end

% save plot
fig = create_figure('mvmt', 2, 1);
subplot(3,1,1);
plot(mvmts(:,4:6));
legend('x', 'y', 'z');

subplot(3,1,2);
plot(mvmts(:,1:3));
legend('pitch', 'roll', 'yaw');

subplot(3,1,3);
FDs = calc_FD(mvmts, softwr);
plot(FDs); hold on;
text(1, 0.2, sprintf('mean FD: %.2f', mean(FDs)));
text(1, 0.22, sprintf('max FD: %.2f', max(FDs)));
ylim([-0.01 0.3]); legend('FD');

sz = get(0, 'screensize'); % Wani added two lines to make this visible (but it depends on the size of the monitor)
set(fig, 'Position', [sz(3)*.02 sz(4)*.05 sz(3) *.45 sz(4)*.85]);

saveas(fig, saveto);
close all;

end