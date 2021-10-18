using Gurobi, JuMP

model = Model(Gurobi.Optimizer)

@variable(model, x[i in 1:n], Bin)
@objective(model, Max, sum(x[i] for i in 1:n))
@constraint(model, clique[i in 1:n, j in 1:n;i != j && a[i,j] == 1], x[i] + x[j] <= 1)

optimize!(model)

objective_value(model)
