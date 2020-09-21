using Statistics
using Distributions
using ProgressMeter

# Defining the abstract node type
abstract type node end

# Defining structures
mutable struct tree_node <: node
        node_label               :: Int #Node Label
        observations_index       :: Array{Int64,1} # Index of Observations
        depth                    :: Int64 # Depth of the Tree
        node_var                 :: Int64 # Index of Splitted var
        node_var_split           :: Float64 # Split var (first we will considerate continuous case)
        terminal                 :: Bool  # If is terminal or not
        parent_node              :: Int64 # Index of parent node
        mu                       :: Float64 # Mu value of distribution of observations that compose that node
        left_node                :: Bool
        right_node               :: Bool
end


mutable struct root <: node
        x                        :: AbstractMatrix # Index of Observations
        depth                    :: Int64 # Depth of the Tree
        mu                       :: Float64 # Mu value of distribution of observations that compose that node
end



mutable struct tree
        root    :: node
        nodes   :: Array{node,1}
end

#Changing print messages

function Base.show(io::IO,r::root)
    print(io,"Root node with ",size(r.x)[1]," observations.\n",
    " μ: ",r.mu,"\n")
end

function Base.show(io::IO,t::tree)
    print(io,"Tree with ",length(t.nodes)," nodes")
end

# Think if need to print more observations
function Base.show(io::IO,n::tree_node)
    print(io,"Node Label: ",n.node_label,"\n",
            " Observations: ",n.observations_index,"\n",
            " Depth: ",n.depth,"\n",
            " Node Var: ",n.node_var,"\n",
            " Node Var Split: ",n.node_var_split,"\n",
            " Parent Node: ",n.parent_node,"\n",
            " Is terminal: ",n.terminal,"\n",
            " μ: ",n.mu,"\n"
    )
end

#==== BART FUNCTIONS ====#

# Grow for root
function grow(node_to_be_spliten::root,min_node_size)

    # Creating a mecanism to reject trees that produce nodes with n<min_node_size
    bad_tree=true
    count_bad_tree=0

    # Creating the left and right node
    left_node_index=Int64[]
    right_node_index=Int64[]
    p_split = Int64
    p_split_value= Float64

    while bad_tree

        #Cleaning the nodes
        left_node_index=[]
        right_node_index=[]


        p_split=rand(1:size(node_to_be_spliten.x)[2],1)[1] # Choose what variable will
                                                                              # be used to split the root node
        p_split_value=rand(node_to_be_spliten.x[:,p_split],1)[1] # Selecting the splitting rule for that var

        # Splitting the node
        for i in 1:size(node_to_be_spliten.x)[1]
            if node_to_be_spliten.x[i,p_split]<p_split_value
                push!(left_node_index,i)
            else
                push!(right_node_index,i)
            end
        end

        #Verying the length of the new terminal nodes
        if (length(left_node_index)<min_node_size) || (length(right_node_index)<min_node_size)
            count_bad_tree+=1
        else
            bad_tree=false
        end

        # To not be trapped in a loop
        if count_bad_tree==2
            return tree(node_to_be_spliten,tree_node[])
        end



    end

    # Defining the children
    left_child_node=tree_node(1,left_node_index,1,p_split,p_split_value,true,0,0,true,false)
    right_child_node=tree_node(2,right_node_index,1,p_split,p_split_value,true,0,0,false,true)

    return tree(node_to_be_spliten,node[left_child_node,right_child_node])

end

# Grow for tree
function grow(tree_to_be_grown::tree,min_node_size)

    # Creating the new tree
    new_tree=deepcopy(tree_to_be_grown)

    terminal_nodes=node[]
    terminal_node_label=Int
    # Sweep all tree nodes
    for i in new_tree.nodes
        if i.terminal #Verify if it is a terminal node
            push!(terminal_nodes,i) # Add the terminal node
        end
    end

    terminal_to_be_spliten_index=rand(1:length(terminal_nodes),1)

    terminal_node_to_be_spliten=terminal_nodes[terminal_to_be_spliten_index][1]


    # Creating a mecanism to reject trees that produce nodes with n<min_node_size
    bad_tree=true
    count_bad_tree=0

    # Creating the left and right node
    left_node_index=Int64[]
    right_node_index=Int64[]
    p_split = Int64
    p_split_value= Float64

    while bad_tree

        #Cleaning the nodes
        left_node_index=[]
        right_node_index=[]


        p_split=rand(1:size(new_tree.root.x)[2],1)[1] # Choose what variable will
                                                                              # be used to split the root node
        # Selecting the splitting rule for that var
        p_split_value=rand(new_tree.root.x[terminal_node_to_be_spliten.observations_index,p_split],1)[1]

        # Splitting the node
        for i in terminal_node_to_be_spliten.observations_index
            if new_tree.root.x[i,p_split]<p_split_value
                push!(left_node_index,i)
            else
                push!(right_node_index,i)
            end
        end

        #Verying the length of the new terminal nodes
        if (length(left_node_index)<min_node_size) || (length(right_node_index)<min_node_size)
            count_bad_tree+=1
        else
            bad_tree=false
        end

        # To not be trapped in a loop
        if count_bad_tree==2
            return tree_to_be_grown
        end



    end

    # Defining the children
    left_child_node=tree_node(length(tree_to_be_grown.nodes)+1,
                                left_node_index,terminal_node_to_be_spliten.depth+1,p_split,
                                p_split_value,
                                true,terminal_node_to_be_spliten.node_label,0,true,false)
    right_child_node=tree_node(length(tree_to_be_grown.nodes)+2,
                                right_node_index,terminal_node_to_be_spliten.depth+1,p_split,
                                p_split_value,
                                true,terminal_node_to_be_spliten.node_label,0,false,true)


    # Updating the splited node
     for i in new_tree.nodes
        if(i.node_label==terminal_node_to_be_spliten.node_label)
            i.terminal=false
            break
        end
     end

    # Adding the new nodes
    push!(new_tree.nodes,left_child_node,right_child_node)

    return new_tree

end

# Prune

function prune(tree_to_be_pruned::tree)

     # Loading the new tree
    new_tree=deepcopy(tree_to_be_pruned)

    # Case of a root node only
    if(length(new_tree.nodes)==0)
        return new_tree
    end

    # In case of a shallow tree (just two terminal nodes)
    if(length(new_tree.nodes)==2)
        return tree(new_tree.root,tree_node[])
    end

    parent_nodes_indexes=Int[]
    parent_to_be_pruned_index=Int
    child_node_pruned=Int
    children_index=Int[0,0]

    # Sweep all tree nodes
    for i in 1:length(new_tree.nodes)
        if !new_tree.nodes[i].terminal #Verify if it is not a terminal node
            push!(parent_nodes_indexes,i) # Add the parent_node
        end
    end

    bad_to_prune=true
    # Veryfing if its children are both terminal nodes

    while(bad_to_prune)
        terminal_counter=0
        children_index=[0,0]

        parent_to_be_pruned_index=rand(parent_nodes_indexes,1)[1]


        children=get_children(new_tree,new_tree.nodes[parent_to_be_pruned_index])

        # Here I verify if the parent node is the chosen one, and if the node is terminal
        if(length(children)==2)
            children_index[1]=children[1]
            children_index[2]=children[2]
            bad_to_prune=false
        end

    end

    # Transforming back as terminal
    new_tree.nodes[parent_to_be_pruned_index].terminal=true

    # Removing the children

    deleteat!(new_tree.nodes,children_index)

    return new_tree

end

# Change function

function change(tree_to_be_changed::tree,min_node_size)
    new_tree=deepcopy(tree_to_be_changed)

    # If the tree is just composed by its terminal node

    if(length(tree_to_be_changed.nodes)==0)
        return new_tree
    end

    parent_nodes_indexes=Int[]
    parent_to_be_pruned_index=Int
    child_node_pruned=Int
    children_index=Int[]
    # Sweep all tree nodes
    for i in 1:length(new_tree.nodes)
        if !new_tree.nodes[i].terminal #Verify if it is not a terminal node
            push!(parent_nodes_indexes,i) # Add the parent_node
        end
    end

    # Creating a mecanism to reject trees that produce nodes with n<min_node_size
    bad_tree=true
    count_bad_tree=0
    new_tree_updated_list=Any[]

    while(bad_tree)

        # Selecting the node to be changed
        node_to_be_changed_label=rand(push!(parent_nodes_indexes,0),1)[1]

        #print("Node changed: ",node_to_be_changed_label,'\n')

        if(node_to_be_changed_label==0)
            current_node= new_tree.root
        else
            current_node=new_tree.nodes[node_to_be_changed_label]
        end


        # Choose what variable will be used to split the root node
        p_split=rand(1:size(new_tree.root.x)[2],1)[1]

        # Selecting the splitting rule for that var
        if(node_to_be_changed_label==0)
            p_split_value=rand(new_tree.root.x[:,p_split],1)[1]
        else
            p_split_value=rand(new_tree.root.x[current_node.observations_index,p_split],1)[1]
        end

        # Get the chidren from node that will be changed
        children_from_changed_node=get_children(new_tree, current_node)

        # Change the node and split rules
        new_tree.nodes[  children_from_changed_node[1]  ].node_var=p_split # Changing the node var
        new_tree.nodes[  children_from_changed_node[2]  ].node_var=p_split # Changing the node var

        new_tree.nodes[  children_from_changed_node[1]  ].node_var_split =p_split_value  # Changing the node split
        new_tree.nodes[  children_from_changed_node[2]  ].node_var_split =p_split_value #Changing the node split


        # Update the tree with the new splitting rule
        new_tree_updated_list=update_tree_function(new_tree,children_from_changed_node,min_node_size)

        # Verify if its a bad tree
        if(new_tree_updated_list[2])
            count_bad_tree+=1
            new_tree=deepcopy(tree_to_be_changed)
        else
            # Setting bad tree as false therefore go out of the loop
            bad_tree=new_tree_updated_list[2]
        end

        # Limit of tries
        if(count_bad_tree==2)
                return new_tree
        end

    end

    return new_tree_updated_list[1]


end

# Swap

function swap(tree_to_be_swaped::tree,min_node_size)

    new_tree=deepcopy(tree_to_be_swaped)

    parent_nodes_indexes=Int[]
    parent_to_be_pruned_index=Int
    child_node_pruned=Int
    children_index=Int[]

    # Count non-terminal
    for i in 1:length(new_tree.nodes)
        if (!new_tree.nodes[i].terminal && new_tree.nodes[i].parent_node!=0) #Verify if is not terminal child node, where the parent is not the root
            if(!new_tree.nodes[(new_tree.nodes[i].parent_node)].terminal)#Verify if its father is also a terminal node
                push!(parent_nodes_indexes,i) # Add the parent_node
            end
        end
    end

    parent_nodes_indexes=unique(parent_nodes_indexes)


    # If there is no terminal nodes enough
    if(length(parent_nodes_indexes)==0)
        return new_tree
    end

    #print(parent_nodes_indexes)

    # Creating a mecanism to reject trees that produce nodes with n<min_node_size
    bad_tree=true
    count_bad_tree=0
    new_tree_updated_list=Any[]

    while(bad_tree)


        # Selecting the node to be changed
        node_to_be_swaped_label=rand(parent_nodes_indexes,1)[1]

        #print("Node changed: ",node_to_be_swaped_label,'\n')

        # Getting the both internal nodes
        current_node_child=new_tree.nodes[node_to_be_swaped_label]
        current_node_father=new_tree.nodes[current_node_child.parent_node]

        #print(current_node_child,'\n')
        #print(current_node_father,'\n')

        # Get the chidren from node that will be swaped
        children_from_swaped_node_child=get_children(new_tree, current_node_child)[[1,2]]
        children_from_swaped_node_father=get_children(new_tree, current_node_father)

        #print(children_from_swaped_node_child,children_from_swaped_node_father)
        # Swapping the nodes

        new_tree.nodes[ children_from_swaped_node_father[1]].node_var, new_tree.nodes[children_from_swaped_node_child[1]].node_var= new_tree.nodes[children_from_swaped_node_child[1]].node_var,new_tree.nodes[ children_from_swaped_node_father[1]].node_var
        new_tree.nodes[ children_from_swaped_node_father[1]].node_var_split, new_tree.nodes[children_from_swaped_node_child[1]].node_var_split= new_tree.nodes[children_from_swaped_node_child[1]].node_var_split,new_tree.nodes[ children_from_swaped_node_father[1]].node_var_split

        new_tree.nodes[ children_from_swaped_node_father[2]].node_var, new_tree.nodes[children_from_swaped_node_child[2]].node_var= new_tree.nodes[children_from_swaped_node_child[2]].node_var,new_tree.nodes[ children_from_swaped_node_father[2]].node_var
        new_tree.nodes[ children_from_swaped_node_father[2]].node_var_split, new_tree.nodes[children_from_swaped_node_child[2]].node_var_split= new_tree.nodes[children_from_swaped_node_child[2]].node_var_split,new_tree.nodes[ children_from_swaped_node_father[2]].node_var_split

        # Update the tree with the new splitting rule
        new_tree_updated_list=update_tree_function(new_tree,children_from_swaped_node_father,min_node_size)

        #print(new_tree_updated_list[1].nodes)
        # Verify if its a bad tree
        if(new_tree_updated_list[2])
            count_bad_tree+=1
            new_tree=deepcopy(tree_to_be_swaped)
        else
            # Setting bad tree as false therefore go out of the loop
            bad_tree=new_tree_updated_list[2]
        end

        # Limit of tries
        if(count_bad_tree==2)
            #print("BAAAD TREEE")
                return tree_to_be_swaped
        end

    end

    return new_tree_updated_list[1]


end

# Auxiliar functions

# Get the children from a tree node
function get_children(update_tree::tree,parent_node::tree_node,all_children_node_copy=[])

    # Count aux
    count_aux=0
    all_children_node=[]
    for i in 1:length(update_tree.nodes)
    #all_chidren_node=node[]
        if(update_tree.nodes[i].parent_node==parent_node.node_label)
            push!(all_children_node,i)
            #print("Node number ",i,"\n")
            count_aux+=1
        end
    end


    all_children_node_copy=deepcopy(all_children_node)

        if(!update_tree.nodes[all_children_node[end]].terminal)
            push!(all_children_node_copy,get_children(update_tree,update_tree.nodes[all_children_node[end]],all_children_node_copy))
        end

        if(!update_tree.nodes[all_children_node[end-1]].terminal)
            push!(all_children_node_copy,get_children(update_tree,update_tree.nodes[all_children_node[end-1]],all_children_node_copy))
        end

    return collect(Iterators.flatten(all_children_node_copy))

end

# Get children from a root node

function get_children(update_tree::tree,parent_node::root,all_children_node_copy=[])

    # Count aux
    count_aux=0
    all_children_node=[]

    for i in 1:length(update_tree.nodes)
        if(update_tree.nodes[i].parent_node==0)
            push!(all_children_node,i)
            count_aux+=1
        end
    end

    all_children_node_copy=deepcopy(all_children_node)

        if(!update_tree.nodes[all_children_node[end]].terminal)
            push!(all_children_node_copy,get_children(update_tree,update_tree.nodes[all_children_node[end]],all_children_node_copy))
        end

        if(!update_tree.nodes[all_children_node[end-1]].terminal)
            push!(all_children_node_copy,get_children(update_tree,update_tree.nodes[all_children_node[end-1]],all_children_node_copy))
        end

    return collect(Iterators.flatten(all_children_node_copy))

end

# Update a tree

function update_tree_function(tree_to_be_updated::tree,list_of_nodes::Array{Int64,1},min_node_size)


    # Creating a copy of tree that will be updated
    new_tree=deepcopy(tree_to_be_updated)

    # Creating a boolean to indicate if will produce a 'bad_tree' (i.e: few number of observations in terminal node)
    is_a_bad_tree=false


    for i in 1:length(list_of_nodes)

        #print(" Updating node: ",list_of_nodes[i],"\n")
        # Creating the aux of each index
        left_index_aux=Int64[]
        right_index_aux=Int64[]

        if(new_tree.nodes[list_of_nodes[i]].left_node)
            # Consider the case of the parent node is the root node
            if(new_tree.nodes[list_of_nodes[i]].parent_node==0)
                for k in 1:size(new_tree.root.x)[1]
                    if(new_tree.root.x[k,new_tree.nodes[list_of_nodes[i]].node_var]<new_tree.nodes[list_of_nodes[i]].node_var_split)
                        push!(left_index_aux,k)
                    end
                end

                # Left index
                new_tree.nodes[list_of_nodes[i]].observations_index=left_index_aux


                # Validating the tree
                if(length(left_index_aux)<min_node_size)
                    is_a_bad_tree=true
                end

            else
            # Consider the other nodes
            for k in new_tree.nodes[new_tree.nodes[list_of_nodes[i]].parent_node].observations_index
                    if(new_tree.root.x[k,new_tree.nodes[list_of_nodes[i]].node_var]<new_tree.nodes[list_of_nodes[i]].node_var_split)
                        push!(left_index_aux,k)
                    end
            end

            # Left-index
            new_tree.nodes[list_of_nodes[i]].observations_index=left_index_aux

            #print("Left index",left_index_aux,"\n")

            # Validating the tree
            if(length(left_index_aux)<min_node_size)
                is_a_bad_tree=true
            end

            end
        else
            # Consider the case of the parent node is the root node
            if(new_tree.nodes[list_of_nodes[i]].parent_node==0)
                for k in 1:size(new_tree.root.x)[1]
                    if(new_tree.root.x[k,new_tree.nodes[list_of_nodes[i]].node_var]>=new_tree.nodes[list_of_nodes[i]].node_var_split)
                        push!(right_index_aux,k)
                    end
                end

                # Right index
                new_tree.nodes[list_of_nodes[i]].observations_index=right_index_aux

                #print("Right index",left_index_aux,"\n")

                # Validating the tree
                if(length(right_index_aux)<min_node_size)
                    is_a_bad_tree=true
                end
            else
            # Consider the other nodes
                for k in new_tree.nodes[new_tree.nodes[list_of_nodes[i]].parent_node].observations_index
                        if(new_tree.root.x[k,new_tree.nodes[list_of_nodes[i]].node_var]>=new_tree.nodes[list_of_nodes[i]].node_var_split)
                            push!(right_index_aux,k)
                        end
                end

                new_tree.nodes[list_of_nodes[i]].observations_index=right_index_aux

                #print("Right index",left_index_aux,"\n")

                if(length(right_index_aux)<min_node_size)
                        is_a_bad_tree=true
                end

            end
        end

    end
    return [new_tree,is_a_bad_tree]
end

# Moving a tree
function move_tree(tree_to_be_moved::tree,move,min_node_size)

    new_tree=deepcopy(tree_to_be_moved)

    if(move=="grow")

        # In case of a tree with just the root node
        if(length(tree_to_be_moved.nodes)==0)
            new_tree=grow(new_tree.root,min_node_size)
        else
            new_tree=grow(new_tree,min_node_size)
        end

    elseif(move=="prune")
        new_tree=prune(tree_to_be_moved)

    elseif(move=="change")
        new_tree=change(tree_to_be_moved,min_node_size)

    else
        new_tree=swap(tree_to_be_moved,min_node_size)

    end

    return [new_tree,move]

end

# Calculating full conditional
function tree_full_conditional(conditional_tree::tree,residuals::Array{Float64,1},tau::Number,tau_mu::Float64)

    # Consider the case of just the single root node
    if(length(conditional_tree.nodes)==0)
       #Seleting the terminal nodes
           sum_residuals_node_squared=sum(residuals.^2)
           sum_residuals_node=sum(residuals)
           terminal_nodes_size=size(conditional_tree.root.x)[1]

           # See the Math RBART .pdf Equation 1
           log_posterior=0.5*length(residuals)*log(tau) + 0.5*(sum(log.(tau_mu ./(tau_mu.+terminal_nodes_size.*tau))))-
                         0.5*tau*sum(sum_residuals_node_squared)+0.5*(tau^2)*sum((sum_residuals_node.^2) ./(tau_mu.+terminal_nodes_size.*tau))

           return log_posterior

    else

        # Seleting the terminal nodes
        terminal_nodes=tree_node[]
        sum_residuals_node_squared=Float64[]
        sum_residuals_node=Float64[]
        terminal_nodes_size=Int[]

        for i in conditional_tree.nodes
            if(i.terminal)
                push!(terminal_nodes,i)
                push!(sum_residuals_node_squared,sum(residuals[i.observations_index].^2))
                push!(sum_residuals_node,sum(residuals[i.observations_index]))
                push!(terminal_nodes_size,length(i.observations_index))
            end
        end

        # Obtatining the number of terminal nodes
        n_terminal_node=length(terminal_nodes)

        # See the Math RBART .pdf Equation 1
        log_posterior=0.5*length(residuals)*log(tau) + 0.5*(sum(log.(tau_mu ./(tau_mu.+terminal_nodes_size.*tau))))-
                      0.5*tau*sum(sum_residuals_node_squared)+0.5*(tau^2)*sum((sum_residuals_node.^2) ./(tau_mu.+terminal_nodes_size.*tau))

        return log_posterior

    end

end


# Updating the mu values

function update_mu(tree_object::tree,residuals,tau,tau_mu)

    # Creating the dummy tree
    new_tree=deepcopy(tree_object)

    # Case of just have the root node
    if(length(new_tree.nodes)==0)
        mean_values=(tau*sum(residuals))/(size(new_tree.root.x)[1]*tau +tau_mu)
        sd_values=sqrt(1/(size(new_tree.root.x)[1]*tau +tau_mu))
        mu_value=randn(1)[1]*sd_values+mean_values
        new_tree.root.mu=mu_value
        return new_tree

    else

         # Seleting the terminal nodes
        terminal_nodes_index=Int[]
        sum_residuals_node=Float64[]
        terminal_nodes_size=Int[]

        for i in 1:length(new_tree.nodes)
            if(new_tree.nodes[i].terminal)
                push!(terminal_nodes_index,i)
                push!(sum_residuals_node,sum(residuals[new_tree.nodes[i].observations_index]))
                push!(terminal_nodes_size,length(new_tree.nodes[i].observations_index))
            end
        end

        # Obtatining the number of terminal nodes
        n_terminal_node=length(terminal_nodes_index)


        # See the Math RBART .pdf Equation 1
        mean_values=(tau.*sum_residuals_node) ./(terminal_nodes_size*tau .+tau_mu)
        sd_values=sqrt.(1 ./(terminal_nodes_size .*tau .+tau_mu))
        mu_array=randn(n_terminal_node) .*sd_values .+mean_values

        #print(mu_array,"\n")

        for i in 1:length(terminal_nodes_index)
            new_tree.nodes[terminal_nodes_index[i]].mu=mu_array[i]
        end

        return new_tree

    end
end

# Updating the tau value

function update_tau(S,nu,lambda,n)

    # Updating tau value
    tau= rand(Gamma((n+nu)/2,  # Shape parameter
                    ((S+nu*lambda)/2)^(-1) ), # Rate Parameter (th's why it is the inverse)
         1)

    return tau
end

# Update tree prior

function tree_prior(tree_object::tree,alpha::Float64,beta::Float64)

    # Initilaizing the log_tree_prior
    log_tree_prior=0

    # Probability over the root node
    if(length(tree_object.nodes)==0)
            log_tree_prior+=log(1-alpha)

        # Prob. of being not-terminal
    else
            log_tree_prior+=log(alpha)

    end

    # Considering the other nodes

    for i in tree_object.nodes

        # Prob. of being terminal
        if(i.terminal)
            log_tree_prior+=log(1-alpha*(1+i.depth)^(-beta))

        # Prob. of being not-terminal
        else
            log_tree_prior+=log(alpha)-beta*log(1+i.depth)

        end
    end

    return log_tree_prior
end

# Get prediction
function get_prediction(multiple_trees_object::Array,x,single_tree=true)

    # Verify if it is a single tree
    if(length(multiple_trees_object)==1)
        single_tree=true
    end

    # Count the numbr of trees
    n_trees=length(multiple_trees_object)


    # Create an empty vector of the predictions
    predictions=repeat([0.0],size(x)[1])

    # Select just one tree
    if(single_tree)

            if(length(multiple_trees_object[1].nodes)==0)
                predictions.=multiple_trees_object[1].root.mu
            else
                # Getting the terminal nodes index
                for node in multiple_trees_object[1].nodes
                    if(node.terminal)
                            predictions[node.observations_index].=node.mu
                    end
                end
            end
    else
        # Recursion over the other trees
        frac_trees=deepcopy(multiple_trees_object)
        frac_trees=deleteat!(frac_trees,1)
        predictions=get_prediction([multiple_trees_object[1]],x,true)+ #Getting the values of the tree which was blanked
                    get_prediction(frac_trees,x, length(multiple_trees_object)==1)

    end

    return predictions
end


# Main juBART function

function juBart(x::Array, # Matrix of observations
               y::Array{Float64,1}, # Target variable
               number_trees::Int64, # Number of trees
               node_min_size::Int64=5, # Minimum size of observations in each terminal nodes
               alpha::Float64=0.95, # Alpha paramteter from tree prior
               beta::Float64=2.0, # Beta parameter from treee prior
               tau_mu::Float64=1.0, # Tau from mu prior
               nu::Float64=3.0, # Parameter from tau prior
               tau::Float64=0.1, # Initial value to tau
               lambda::Float64=0.1, # Lambda value
               tau_zero::Float64=0.1, # Tau initial value

               n_iter::Int64=1000,
               burn_in::Int64=250,
               thin::Int64=1)

    # Scaling hte y
    y_mean=mean(y)
    y_sd=std(y)
    y_scale=(y.-y_mean) ./y_sd

    sigma=0.0
    # Initializing the algorithm
    n=length(y)
    current_trees=Array{Any}(nothing,number_trees)

    # Think in a error message from number_tree=0

    # Storage containers
    store_size=Int((n_iter-burn_in)/thin)
    tree_store=Array{Any}(nothing,store_size)
    sigma_store=Array{Union{Nothing,Float64},1}(nothing,store_size)
    y_hat_store=Array{Union{Nothing,Float64}}(nothing,store_size,n)
    log_like_store=Array{Union{Nothing,Float64},1}(nothing,store_size)
    full_cond_store=Array{Union{Nothing,Float64},2}(nothing,store_size,number_trees)

    # Creating the stumps
    current_trees=fill(tree(root(x,0,0),tree_node[]),number_trees )


    # Getting the predictions
    predictions=get_prediction(current_trees,x,number_trees==1)
    #print(predictions,"\n")

    # Iterating over MH runs
    @showprogress 1 "MH Runs..." for i in 1:n_iter

        # Creating the storage
        if((i > burn_in) && (i%thin==0))
            curr=convert(Int,round((i-burn_in)/thin,digits=1))
            tree_store[curr]=deepcopy(current_trees)
            sigma_store[curr]=sigma
            y_hat_store[curr,:]=predictions
            #log_like_store[curr]=log_like
        end

        # Iterating over the trees
        for j in 1:number_trees

            # Calculating the residuals for each
            if(number_trees>1)
                partial_trees=deepcopy(current_trees)
                partial_trees=deleteat!(partial_trees,j)

                current_partial_residuals=y_scale-
                                          get_prediction(partial_trees,x, number_trees==1)
            else
                current_partial_residuals=y_scale
            end

            move=rand(["grow","prune","change","swap"],1)[1]
            if i < max(floor(0.1*burn_in), 10) move="grow" end

            #print(move,"\n")
            # Getting the new tree
            new_trees=deepcopy(current_trees)
            new_trees[j]=move_tree(new_trees[j],move,node_min_size)[1]

            #print(new_trees[j].nodes,"\n")
            #print(current_trees[j].root,"\n")

            #print(current_trees,"\n")
            #print(new_trees,"\n")



            #print(move_tree(new_trees[j],move,node_min_size)[1],"\n")

            # Calculating the likelihood of each set
            l_new=tree_full_conditional(new_trees[j],
                                                current_partial_residuals,tau,tau_mu)+
                      tree_prior(new_trees[j], alpha, beta)


            l_old=tree_full_conditional(current_trees[j],
                                            current_partial_residuals,tau,tau_mu)+
                  tree_prior(current_trees[j], alpha, beta)

                  #print("L_old: ",l_old,"\n")
                  #print("L_new: ",l_new,"\n")
            #Probability of accept the new proposed tree
            acceptance=exp(l_new-l_old)
            #print(acceptance,"\n")
            #If Storage or not based on thin and burn parameters

            if((i > burn_in) && ((i %thin) == 0) )
                    full_cond_store[curr, j] = l_old
            end

            if(rand(1)[1]<acceptance)
                    #print("ACCEPT! \n")
                    # Make changes if accept
                    #print(new_trees,"\n")
                    #print(current_trees,"\n")
                    current_trees = deepcopy(new_trees)
                    #print(current_trees,"\n")
                    #print(new_trees,"\n")
            end # End of accept if statement

            #To update the mu values
            current_trees[j]=update_mu(current_trees[j],
                                           current_partial_residuals,
                                           tau,tau_mu)


        end #End of iterations over the trees

     # Get the predicitions from updated trees
     predictions=get_prediction(current_trees,x,number_trees==1)
     S=sum((y_scale-predictions).^2)

    # Update tau
    tau=update_tau(S,nu,lambda,n)[1]

    sigma=1/sqrt(tau)
    #print(sigma)

    # NEED TO THINK HOW TO CALCULATE IT USING JULIA
    #log_like=loglikelihood(Normal(predictions,sigma),y_scale)

    end # End the MH iterations

    return Dict("trees"=>tree_store,
            "sigma"=>sigma_store,
            "predictions"=>y_hat_store,
            "full_conditionals"=>full_cond_store,
            "y"=>y,"X"=>x,"n_iter"=>n_iter,"burn_in"=>burn_in,
            "thin"=>thin,"store_size"=>store_size,"number_trees"=>number_trees)
end


# Function to update a tree with a new X
function predict_tree_model(tree_object::tree,x_new::Array{Float64,2})

    # Creating the new tree
    new_tree=deepcopy(tree_object)

    # Updating the root node
    new_tree.root.x=x_new

    # Getting the index of nodes to be updated
    nodes_index=collect(1:length(new_tree.nodes))

    # Getting the new tree
    new_tree=update_tree_function(new_tree,nodes_index,0)[1]

    return new_tree
end


function predict_jubart(jubart_model::Dict,x_predict::Array{Float64,2},type::String="mean")

    new_jubart_model=deepcopy(jubart_model)
    n_iter=size(new_jubart_model["predictions"])[1]
    n_xpredict=size(x_predict)[1]
    y_hat_matrix=Array{Union{Nothing,Float64}}(nothing,n_iter,n_xpredict)

    @showprogress 1 "Predicting..." for i in 1:n_iter
    #for i in 1:n_iter

        for k in 1:jubart_model["number_trees"]
             new_jubart_model["trees"][i][k]=predict_tree_model(new_jubart_model["trees"][i][k],x_predict)
        end

        y_hat_matrix[i,:]=get_prediction(new_jubart_model["trees"][i],x_predict,jubart_model["number_trees"]==1)

    end

    if(type=="mean")
        return mean(y_hat_matrix,dims=1)
    elseif(type=="median")
        return median(y_hat_matrix,dims=1)
    else
        return y_hat_matrix
    end
end


function sim_friedman_simple(n,p=0,scale_err=1)
    # Simulate some data using a simple version of friedman
    # y= 10sin(πx1*x2)+20(x3=0.5)2+10*4+ 5x5 +ε

    X=rand(n,5+p)
    pars=[10,20,10,5]

    mean=pars[1] .* sin.(pi .*X[:,1] .*X[:,2]) .+ pars[2] .*(X[:,3] .-0.5) .^2 .+
         pars[3] .*X[:,4] .+ pars[4] .*X[:,5]

    y=randn(n).*scale_err .+mean

    return Dict("x"=>X,"y"=>y,"true_mean"=>mean,"true_scale"=>scale_err)
end
