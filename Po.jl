using JuMP, GLPK;  # Need to say it whenever we use JuMP

using GLPKMathProgInterface; # Loading the GLPK module for using its solver

mutable struct Instance

    T::Int64
    cc1::Float64;
    cc2::Float64;
    CF::Float64;
    d::Array{Int64};
    cec1::Float64;
    cec2 ::Float64;
    cep1::Float64;
    cep2::Float64;
    Instance() = new();


end


function ReadInstance(file_name::AbstractString)
    
    instance = Instance();    

    open(file_name, "r") do file

        instance.T = parse(Int64, readline(file));
        instance.cc1 = parse(Float64, readline(file));
        instance.cc2 = parse(Float64, readline(file));
        instance.CF = parse(Float64, readline(file));
        instance.d = parse.(Int64,split(readline(file)));
        instance.cec1 = parse(Float64, readline(file));
        instance.cec2 = parse(Float64, readline(file));
        instance.cep1 = parse(Float64, readline(file));
        instance.cep2 = parse(Float64, readline(file));
    end
    return instance;
end

#MODEL CONSTRUCTION
#--------------------
function SolveMasterProductionPlan(instance::Instance)
    
    myModel = Model(with_optimizer(GLPK.Optimizer) )


 
    #VARIABLES
    #---------
    @variable(myModel, Xp1[i=1:instance.T] >= 0,Int)
    @variable(myModel, Xp2[i=1:instance.T] >= 0,Int)

    @variable(myModel, Yc1[i=1:instance.T] >= 0,Int)
    @variable(myModel, Yc2[i=1:instance.T] >= 0,Int)

    @variable(myModel, Ec1[i=1:instance.T] >= 0,Int)
    @variable(myModel, Ec2[i=1:instance.T] >= 0,Int)
    @variable(myModel, Ep1[i=1:instance.T] >= 0,Int)
    @variable(myModel, Ep2[i=1:instance.T] >= 0,Int)
    @variable(myModel, P[i=1:instance.T] ,Bin)
    
    #OBJECTIVE

    @objective(myModel, Min, sum( instance.cc1*Yc1[j] 
                            + instance.cc2*Yc2[j] 
                            + Ec1[j] * instance.cec1 
                            + Ec2[j] * instance.cec2 
                            + Ep1[j] * instance.cep1 
                            + Ep2[j] * instance.cep2  
                            + instance.CF * P[j] for j=1:instance.T) ) 

    #CONSTRAINTS
    #-----------

        @constraint(myModel, tempo_max[t in 1:instance.T], Xp1[t] + Xp2[t] <= 800) 

        @constraint(myModel, demanda_p1[t in 1:instance.T], sum(Xp1[j] - instance.d[j] for j = 1:t) >= 0) 

        @constraint(myModel, demanda_c1c2[t in 1:instance.T], (Yc1[t] + Yc2[t]) <= P[t]* sum(instance.d * 9) ) 
        
                
        @constraint(myModel, [t in 2:instance.T], 2*Xp1[1] <= Xp2[1] ) 

        @constraint(myModel, [t in 2:instance.T], 3*Xp1[1] + Xp2[1] <= Yc1[1]) 
            
        @constraint(myModel, [t in 2:instance.T], 2*Xp2[1] <= Yc2[1]) 
                
        
        @constraint(myModel, [t in 2:instance.T], 2*Xp1[t] <= Xp2[t] + Ep2[t-1]) 

        @constraint(myModel, [t in 2:instance.T], 3*Xp1[t] + Xp2[t] <= Yc1[t] + Ec1[t-1]) 
        
        @constraint(myModel, [t in 2:instance.T], 2*Xp2[t] <= Yc2[t] + Ec2[t-1]) 
            
        @constraint(myModel, [t in 1:instance.T], Ep1[t] == sum( Xp1[j] - instance.d[j] for j = 1:t))
        
        @constraint(myModel, [t in 1:instance.T], Ep2[t] == sum( Xp2[j] - Xp1[j] * 2 for j = 1:t) )
        
        @constraint(myModel, [t in 1:instance.T], Ec1[t] == sum( Yc1[j] - (Xp1[j]*3 + Xp2[j]) for j = 1:t))
        
        @constraint(myModel, [t in 1:instance.T], Ec2[t] == sum( Yc2[j] - (Xp2[j] * 2) for j = 1:t))


    #THE MODEL IN A HUMAN-READABLE FORMAT
    #------------------------------------
    println("The optimization problem to be solved is:")
    print(myModel) # Shows the model constructed in a human-readable form

    #SOLVE IT AND DISPLAY THE RESULTS
    #--------------------------------
    optimize!(myModel) # solves the model  

    println("Objective value: ", objective_value(myModel)) # getObjectiveValue(model_name) gives the optimum objective value
    
    println("\nCompra por semana:\n")

    for i =1:instance.T
        println("Y1[",i,"] = ", value(Yc1[i])," ","Y2[",i,"] = ", value(Yc2[i]),)
    end

    println("\nFabricação por semana:\n")

    for i =1:instance.T
        println("X1[",i,"] = ", value(Xp1[i])," ","X2[",i,"] = ", value(Xp2[i]),)
    end

    println("\nEstoque C1:\n")

    for i =1:instance.T
        println("Ec1[",i,"] = ", value(Ec1[i]),)
    end

    println("\nEstoque C2:\n")

    for i =1:instance.T
        println("Ec1[",i,"] = ", value(Ec2[i]),)
    end

    println("\nEstoque P1:\n")

    for i =1:instance.T
        println("Ep1[",i,"] = ", value(Ep1[i]),)
    end

    println("\nEstoque P2:\n")

    for i =1:instance.T
        println("Ep2[",i,"] = ", value(Ep2[i]),)
    end
    println();
    for i =1:instance.T
        println("P[",i,"] = ", value(P[i]),)
    end
end


function SolverMasterProductionPlan(file_name::AbstractString)

    instance = ReadInstance(file_name);

    SolveMasterProductionPlan(instance);

end

SolverMasterProductionPlan("test01.txt");



