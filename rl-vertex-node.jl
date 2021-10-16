using Gurobi, JuMP, LightGraphs, DelimitedFiles, GraphPlot

path = "C:\\Users\\Wanderson\\OneDrive\\Área de Trabalho\\vertex-packing\\instancia.txt" #instancia

path = "C:\\Users\\Wanderson\\OneDrive\\Área de Trabalho\\vertex-packing\\instancia2.txt" #instancia 2


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
    x = zeros(n)
    soma = 0
    for i in 1:n 
        for j in i:n 
            if a[i,j] == 1
                soma += u[i,j]
            end 
        end 

        if 1-soma > 0 
            x[i] = 1
            for j in i:n 
                if a[i,j] == 1 
                    x[j] = 0
                end 
            end 
        else 
            x[i] = 0
            for j in i:n 
                if a[i,j] == 1 
                    x[j] = 1
                end 
            end 
        end  
        soma = 0
    end 

    valor_otimo = obj_value(x, u, n)
    return valor_otimo, x
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



maxIter = 1000
p_i = 2
pi_min = 0.0001
eps = 0.1
n,m,a = leitura_arquivo(path)

u = zeros(n,n)

best_lim_inf = 0
best_lim_sup = 9999999
improve = 0
opt = 0

G = generate_graph(a,n)

# best_lim_inf = 3

for k in 1:maxIter
    z, x_sub = subproblema(n, u, a) 
    lb, x_down = viable_solution(a,n,x_sub)
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
        for j in 1:n
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
            u[i,j] = max(0, u[i,j] - 1.1*T*s[i,j])
        end
    end
end

println(best_lim_inf)
println(best_lim_sup) 



function guloso(n,a)

    G = generate_graph(a,n)
    packing = []
    vistos = [9999]

    while true

    #percorrer a lista de graus, pegar o com menor grau diferente de 0

        min_degree = 9999
        aux_index = 0
        index = 0

        for degree in degree(G)
            aux_index +=1
            for i in vistos
                if aux_index != i && degree < min_degree
                    println("entrei nos que tem indice diferente")
                    println("indice", aux_index)
                    min_degree = degree
                    index = aux_index 
                end
            end
        end

    #adicionar ao empacotamento para

        append!(packing,index)

    #remover todas as arestas de u e dos vizinhos de um

        auxiliar = all_neighbors(G,index)
        vizinhos = copy(auxiliar)
        append!(vizinhos,index)

        for l in vizinhos
            append!(vistos,l)
            println("to salvando os já vistos")
        end

        println(vistos)

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
                return packing
            end
        end  
    end
end




n,m,a = leitura_arquivo(path)
clean_matrix(a,n)

guloso(n,a)

println(degree(G))