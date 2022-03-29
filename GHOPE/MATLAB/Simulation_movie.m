load('/home/dankas/TEST_FOR_BUILDS/Debug/output/state_data.dat');
load('/home/dankas/TEST_FOR_BUILDS/Debug/output/debug/execution_times.dat');

AU = 1.4960e11;
state_data(:,3:5) = state_data(:,3:5)./AU;

exec_time_m = mean(execution_times(:,2))
exec_time_a = execution_times(end,3)/execution_times(end,1)

IDS = unique(state_data(:,2));
t_v = sort(unique(state_data(:,1)));

N = length(t_v);

figure(1);
for i=1:N

plot3(state_data(state_data(:,1) == t_v(i),3),state_data(state_data(:,1) == t_v(i),4),state_data(state_data(:,1) == t_v(i),5),'ob')
axis([-2 2 -2 2 -2 2]);
drawnow

pause(0.1)

end