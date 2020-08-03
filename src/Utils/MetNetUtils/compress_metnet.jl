compress_model(model::MetNet) = MetNet(model; reshape = false, S = sparse(model.S), b = sparsevec(model.b));

uncompress_model(model::MetNet) = MetNet(model; reshape = false, S = Matrix(model.S), b = Vector(model.b));