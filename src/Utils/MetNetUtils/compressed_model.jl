compressed_model(model::MetNet; sparsity_th = 0.66) = 
    MetNet(compressed_copy(struct_to_dict(model); 
        sparsity_th = sparsity_th))
uncompressed_model(model::MetNet) = 
    MetNet(uncompressed_copy(struct_to_dict(model)))