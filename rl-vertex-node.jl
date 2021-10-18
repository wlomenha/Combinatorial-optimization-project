using Gurobi, JuMP, LightGraphs, DelimitedFiles, GraphPlot, Random

path = "D:\\GitHub - Projects\\vertex-packing\\instancia.txt" #instancia

path = "D:\\GitHub - Projects\\vertex-packing\\instancia2.txt" #instancia 2

path = "D:\\GitHub - Projects\\vertex-packing\\instancia3.txt" #instancia 3 com 32 vértices

path = "D:\\GitHub - Projects\\vertex-packing\\instancia50v03.txt"

path = "D:\\GitHub - Projects\\vertex-packing\\instancia50v05.txt"

function plotGraph(graph)
    nodelabel = 1:nv(graph)
    p = gplot(graph, nodelabel = nodelabel)
    display(p)
end

function clean_matrix(a,n)
    for i in 1:n 
        for j in 1:n
            if i > j
                a[i,j] = 0
            end
        end
    end
    return a
end

function leitura_arquivo(path)
    n = readdlm(path)[1,1] 
    m = readdlm(path)[1,2]

    a = zeros(n,n)
    for i in 2:n+1
        for  j in 1:n
            a[i-1,j] = readdlm(path)[i,j]
        end
    end
    b = clean_matrix(a,n)
    return n,m,b
end 

function generate_graph(a,n)

    G = LightGraphs.SimpleGraph(n)

    for i in 1:n
        for j in 1:n
            if a[i,j] == 1
                add_edge!(G,i,j)
            end
        end
    end
    return G
end

function obj_value(x,u,n)
    soma = 0
    for i in 1:n 
        soma += x[i]
    end 

    for i in 1:n 
        for j in i:n 
            if a[i,j] == 1
                soma += u[i,j]*(1-x[i]-x[j])
            end 
        end 
    end 
    return soma
end

function subproblema(n,u,a)
    m = Model(Gurobi.Optimizer)
    @variable(m, x[i in 1:n], Bin)
    @objective(m, Max, sum(x[i] for i in 1:n) + sum(u[i,j]*(1-x[i] - x[j]) for i in 1:n, j in i:n))
    optimize!(m)

    valor_otimo = objective_value(m)
    #x_val = .value(x)
    x_val = []
    for i in 1:n
        x_val = value(x[i])
    end
    println(x_val)
    # x = zeros(n)
    # soma = 0
    # for i in 1:n 
    #     for j in i:n 
    #         if a[i,j] == 1
    #             soma += u[i,j]
    #         end 
    #     end 

    #     if 1-soma > 0 
    #         x[i] = 1
    #         #for j in i:n 
    #             #if a[i,j] == 1 
    #                 #x[j] = 0
    #             #end 
    #         #end
    #     #end 
    #     else 
    #         x[i] = 0
    #       #  for j in i:n 
    #        #     if a[i,j] == 1 
    #         #        x[j] = 0
    #          #   end 
    #         #end 
    #     end  
    #     soma = 0
    # end 

    # valor_otimo = obj_value(x, u, n)
    return valor_otimo, x_val
end

function viable_solution(a,n,x)
    model = Model(Gurobi.Optimizer)

    @variable(model, x[i in 1:n], Bin)
    @objective(model, Max, sum(x[i] for i in 1:n))
    @constraint(model, rc[i in 1:n, j in i:n; a[i,j] == 1], x[i] + x[j] <= 1)

    optimize!(model)
    y = value.(x)
    return objective_value(model), y  
 
end

function guloso(n,a)
    
    G = generate_graph(a,n)
    #plotGraph(G)
    packing = []
    v = zeros(n)

    while true

    #percorrer a lista de graus, pegar o com menor grau diferente de 0

        min_degree = 9999
        aux_index = 0
        index = 0

        for degree in degree(G)
            aux_index +=1
            if v[aux_index] == 0 && degree < min_degree
                min_degree = degree
                index = aux_index 
            end
        end
        v[index] = 1
        println("o vértice escolhido para entrar no packing foi: ", index)
    #adicionar ao empacotamento para

        append!(packing,index)

    #remover todas as arestas de u e dos vizinhos de um

        auxiliar = all_neighbors(G,index)
        vizinhos = copy(auxiliar)
        for h in vizinhos #impossibilitando de adicionar os vizinhos no packing
            v[h] = 1
        end     
        append!(vizinhos,index)
        

        for i in vizinhos
            aux = neighbors(G, i)
            aux_vizinho = copy(aux)
            for j in aux_vizinho
                rem_edge!(G,i,j)
            end
        end

    #repetir os processos anteriores até que todos os graus sejam 0
        println(degree(G))
        count_zeros = 0
        for k in degree(G)
            if k == 0
                count_zeros += 1
            end
            if count_zeros == length(degree(G))
                return packing
            end
        end  
    end
end


function guloso_pseudo_random(n,a)
    
    G = generate_graph(a,n)
    #plotGraph(G)
    packing = []
    v = zeros(n)

    choice = rand(1:n)
    v[choice] = 1
    append!(packing,choice)

    auxiliar = all_neighbors(G,choice)
    vizinhos = copy(auxiliar)

    for h in vizinhos #impossibilitando de adicionar os vizinhos no packing
        v[h] = 1
    end     

    append!(vizinhos,choice)
    
    for i in vizinhos
        aux = neighbors(G, i)
        aux_vizinho = copy(aux)
        for j in aux_vizinho
            rem_edge!(G,i,j)
        end
    end


    while true

    #percorrer a lista de graus, pegar o com menor grau diferente de 0

        min_degree = 9999
        aux_index = 0
        index = 0

        for degree in degree(G)
            aux_index +=1
            if v[aux_index] == 0 && degree < min_degree
                min_degree = degree
                index = aux_index 
            end
        end
        v[index] = 1
    #adicionar ao empacotamento para

        append!(packing,index)

    #remover todas as arestas de u e dos vizinhos de um

        auxiliar = all_neighbors(G,index)
        vizinhos = copy(auxiliar)
        for h in vizinhos #impossibilitando de adicionar os vizinhos no packing
            v[h] = 1
        end     
        append!(vizinhos,index)
        

        for i in vizinhos
            aux = neighbors(G, i)
            aux_vizinho = copy(aux)
            for j in aux_vizinho
                rem_edge!(G,i,j)
            end
        end

    #repetir os processos anteriores até que todos os graus sejam 0
        count_zeros = 0
        for k in degree(G)
            if k == 0
                count_zeros += 1
            end
            if count_zeros == length(degree(G))
                return size(packing)[1]
            end
        end  
    end
end


function guloso_random(n,a)
    
    G = generate_graph(a,n)
    #plotGraph(G)
    packing = []
    v = zeros(n)

    while true

    #percorrer a lista de graus, pegar o com menor grau diferente de 0
        choice = rand(1:n)
        if v[choice] == 1
            while v[choice] == 1
                choice = rand(1:n)
            end
        end
        
        v[choice] = 1
        append!(packing,choice)

    #remover todas as arestas de u e dos vizinhos de um

        auxiliar = all_neighbors(G,choice)
        vizinhos = copy(auxiliar)
        for h in vizinhos #impossibilitando de adicionar os vizinhos no packing
            v[h] = 1
        end     
        append!(vizinhos,choice)
        

        for i in vizinhos
            aux = neighbors(G, i)
            aux_vizinho = copy(aux)
            for j in aux_vizinho
                rem_edge!(G,i,j)
            end
        end

    #repetir os processos anteriores até que todos os graus sejam 0
        println(degree(G))
        count_zeros = 0
        for k in degree(G)
            if k == 0
                count_zeros += 1
            end
            if count_zeros == n
                return size(packing)[1]
            end
        end  
    end
end


n,m,a = leitura_arquivo(path)

maxIter = 1000
p_i = 2
pi_min = 0.0001
eps = 0.1


u = zeros(n,n)

best_lim_inf = 0
best_lim_sup = 9999999
improve = 0
opt = 0

G = generate_graph(a,n)


for k in 1:maxIter
    z, x_sub = subproblema(n, u, a)
    lb = guloso_pseudo_random(n,a)

    if lb > best_lim_inf
        best_lim_inf = lb
    end

    println(z)

    if z < best_lim_sup
        best_lim_sup = z
        improve = 0
    else
        improve += 1
    end

    if best_lim_sup - best_lim_inf < 1
        opt = best_lim_sup
        println("saiu aqui!")
        #x_opt = ...
        #y_opt = ...
        break
    end

    if improve >= maxIter/20
        p_i = p_i/2
        improve = 0
        if p_i < pi_min
            break
        end
    end

    s = zeros(n,n)
    soma_s = 0

    for i in 1:n
        for j in i:n
            s[i,j] = 1 - (x_sub[i] + x_sub[j]) 
            soma_s += s[i,j]^2
        end
    end

    T = p_i*((z - lb)/soma_s)

    if T < pi_min
        break
    end

    for i in 1:n
        for j in i:n
            u[i,j] = max(0, u[i,j] + 1.2*T*s[i,j])
        end
    end
end

println(best_lim_inf)
println(best_lim_sup) 


n,m,a = leitura_arquivo(path)

guloso(n,a) 

guloso_pseudo_random(n,a)

guloso_random(n,a)



