function outputInfo(obj)
fprintf('iter:%d, \tdiff. rel. MSE:%2.2f\n',obj.iter,obj.optimParams.relMSEdiff(end));
if obj.iter>3
    fprintf('\t\tMSE mov. avg.:%2.2f\n',obj.optimParams.totalCost(end));
end
fprintf('\n');
end