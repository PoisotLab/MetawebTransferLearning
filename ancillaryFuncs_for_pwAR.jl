function leaf_traits_reconstruction(traits,tree)
    ancestral_rec = ancestralStateReconstruction(traits, tree);
    predint_rec = predint(ancestral_rec, level = 0.9)
    lower = predint_rec[:,1]
    upper = predint_rec[:,2]
    mean_trait = mean(predint_rec,dims = 2)[:,1]
    return [lower,upper,mean_trait]
end