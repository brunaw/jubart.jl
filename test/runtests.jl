using jubart
using Test


@testset "jubart.jl" begin
    # Testing the friedman_simulation_function
    @test isa(sim_friedman_simple(10,5,0.01),Dict)
    @test isa(sim_friedman_simple(10,5,0.01)["x"],Array{Float64,2})
    @test isa(sim_friedman_simple(10,5,0.01)["y"],Array{Float64,1})
    @test isa(sim_friedman_simple(10,5,0.01)["true_mean"],Array{Float64,1})

    # Testing the Julia model
    data_test=sim_friedman_simple(100,5,0.001)
    min_node_size_test=10
    number_tree=5
    model_test=juBart(data_test["x"],data_test["y"],number_tree,min_node_size_test)

    # Testing just the output
    @test isa(model_test,Dict)

    # Test for check the min node size
    node_size_count=0 # Counter of min_node_size

    for mh_tree in model_test["trees"]
        for m in 1:number_tree
            tree_nodes=mh_tree[m].nodes
            for run_tree_nodes in tree_nodes
                if(length(run_tree_nodes.observations_index)<min_node_size_test)
                    node_size_bol+=1
                end
            end
        end
    end

    if(node_size_count==0) node_size_bol=true else node_size_bol=false end

    # Veryfing the node boolean
    @test node_size_bol

    # Veryfing NA in test model
    @test !any(isnan,model_test["predictions"])


    # Testing the predict function
    prediction_test_model_all=predict_juBart(model_test,rand(100,10),"all") # Creating the all matrix
    prediction_test_model_mean=predict_juBart(model_test,rand(100,10)) # Creating the prediction vector

    # Verfying size
    @test isa(prediction_test_model_all,Array{Union{Nothing,Float64},2})
    @test isa(prediction_test_model_mean,Array{Float64,2})

    # Veryfing NA
    @test !any(isnan,prediction_test_model_all)
    @test !any(isnan,prediction_test_model_mean)


end
