function ndata = normalizeSignal(data, dim)
dataMean = mean(data, dim);
ShapeDataMean = ones(1,length(size(data)));
ShapeDataMean(dim) = size(data, dim);
ndata = data - repmat(dataMean, ShapeDataMean);
dataStd = std(ndata, 0, dim);
ndata = ndata./repmat(dataStd, ShapeDataMean);