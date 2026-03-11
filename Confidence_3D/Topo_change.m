% 1. 读入外部数据data.txt
data = load('topo.dat');

% 2. 将数据体内第三列全部置为0
data(:, 3) = 0;

% 3. 以topo.txt文件形式输出
dlmwrite('topo_0.dat', data, 'delimiter', '\t', 'precision', 6);
