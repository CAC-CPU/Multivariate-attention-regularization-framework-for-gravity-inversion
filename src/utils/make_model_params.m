function model = make_model_params(params)
nx = params.xGridsNum;
nz = params.zGridsNum;
fileName = params.fileName;
model = zeros(nz,nx);
if params.makeModelFlag == 1
    writematrix(model,'model.xlsx');
    system('model.xlsx');
end
model = readmatrix('model.xlsx');
saveDirectory = 'data/';
filePath = fullfile(saveDirectory,fileName);
save(filePath,'model','-ascii');
end